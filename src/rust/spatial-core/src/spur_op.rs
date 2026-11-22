//! Whole-operation SPUR diagnostics, half-life confidence intervals, and spatial transforms.

use crate::calibration::{correlation_floor, solve_c0_from_avc};
use crate::distance::{distance_matrix, lower_triangle, Metric};
use crate::threads::{is_serial, run_with_threads};
use crate::{basis, checked_len, linalg, CoreError};
use rayon::prelude::*;

mod halflife;
mod transform;

pub use halflife::{run_spur_halflife, SpurHalfLifeOutput, SpurHalfLifeRequest};
use transform::apply_lbmgls;
pub use transform::{run_spur_transform, SpurTransformKind, SpurTransformRequest};

const SUM_CHUNK_SIZE: usize = 32_768;
const TILE_COLS: usize = 1024;
const OMEGA_GRID_BLK: usize = 256;
const MAX_PAIR_CONTRIB_DOUBLES: usize = 16_000_000;
const SZ_LEVELS: [f64; 3] = [0.01, 0.05, 0.10];
const RHO_BM: f64 = 0.999;
const RHO_I0_ALT: f64 = 0.001;
const RHO_GRID_MIN: f64 = 0.0001;
const RHO_GRID_MAX: f64 = 0.03;
const RHO_GRID_LEN: usize = 30;
const MAX_CBAR_BRACKET_STEPS: usize = 500;
const MAX_HA_BRACKET_STEPS: usize = 200;

struct SpurGeometry {
    n: usize,
    max_dist: f64,
    norm_distmat: Vec<f64>,
    lower_tri: Vec<f64>,
    correlation_floor: f64,
    cbar_cache: Vec<(f64, f64)>,
    mean_exp_cache: Vec<(u64, f64)>,
}

impl SpurGeometry {
    fn from_coords(coords: &[f64], n: usize, d: usize, metric: Metric) -> Result<Self, CoreError> {
        if n < 2 || d == 0 || coords.len() != checked_len(&[n, d])? {
            return Err(CoreError::invalid_input(
                "SPUR requires at least 2 observations and 1 coordinate column",
            ));
        }
        let mut norm_distmat = distance_matrix(coords, n, d, metric)?;
        let max_dist = norm_distmat.iter().copied().fold(0.0_f64, f64::max);
        if !max_dist.is_finite() || max_dist <= 0.0 {
            return Err(CoreError::invalid_input(
                "all coordinates are identical; normalization requires distinct locations",
            ));
        }
        for v in &mut norm_distmat {
            *v /= max_dist;
        }
        let lower_tri = lower_triangle(&norm_distmat, n);
        let correlation_floor = correlation_floor(&lower_tri)?;
        Ok(Self {
            n,
            max_dist,
            norm_distmat,
            lower_tri,
            correlation_floor,
            cbar_cache: Vec::new(),
            mean_exp_cache: Vec::new(),
        })
    }

    fn validate_rhobar(&self, rhobar: f64) -> Result<(), CoreError> {
        if !rhobar.is_finite() || rhobar <= self.correlation_floor || rhobar >= 1.0 {
            return Err(CoreError::calibration(
                "target correlation must exceed the co-located-pair fraction and be below one",
            ));
        }
        Ok(())
    }

    /// Cached SPUR decay calibration, preserving replication-scale rounding for a valid bracket.
    fn get_cbar(&mut self, rhobar: f64) -> Result<f64, CoreError> {
        self.validate_rhobar(rhobar)?;
        if let Some(&(_, cached)) = self
            .cbar_cache
            .iter()
            .find(|&&(r, _)| r.to_bits() == rhobar.to_bits())
        {
            return Ok(cached);
        }

        let cbar = solve_cbar_cached(
            &self.lower_tri,
            rhobar,
            !is_serial(),
            &mut self.mean_exp_cache,
        )?;
        self.cbar_cache.push((rhobar, cbar));
        Ok(cbar)
    }

    /// Solve `get_cbar` for a batch of `rhobar` values, reusing cache and running in parallel.
    fn get_cbar_batch(&mut self, rhos: &[f64]) -> Result<Vec<f64>, CoreError> {
        for &r in rhos {
            self.validate_rhobar(r)?;
        }
        let mut missing: Vec<f64> = Vec::new();
        for &r in rhos {
            let bits = r.to_bits();
            let in_cache = self.cbar_cache.iter().any(|&(cr, _)| cr.to_bits() == bits);
            let in_missing = missing.iter().any(|&mr| mr.to_bits() == bits);
            if !in_cache && !in_missing {
                missing.push(r);
            }
        }

        if !missing.is_empty() {
            let solved: Vec<(f64, f64)> = if !is_serial() && missing.len() >= 4 {
                let seed_cache = self.mean_exp_cache.clone();
                missing
                    .par_iter()
                    .map(|&r| {
                        let mut local_cache = seed_cache.clone();
                        solve_cbar_cached(&self.lower_tri, r, false, &mut local_cache)
                            .map(|c| (r, c))
                    })
                    .collect::<Result<Vec<_>, _>>()?
            } else {
                let par_chunks = !is_serial();
                let mut solved_vec = Vec::with_capacity(missing.len());
                for &r in &missing {
                    let c = solve_cbar_cached(
                        &self.lower_tri,
                        r,
                        par_chunks,
                        &mut self.mean_exp_cache,
                    )?;
                    solved_vec.push((r, c));
                }
                solved_vec
            };
            self.cbar_cache.extend(solved);
        }

        let mut out = Vec::with_capacity(rhos.len());
        for &r in rhos {
            let bits = r.to_bits();
            let &(_, val) = self
                .cbar_cache
                .iter()
                .find(|&&(cr, _)| cr.to_bits() == bits)
                .expect("cbar missing after batch solve");
            out.push(val);
        }
        Ok(out)
    }

    fn sigma_lbm_dm(&self) -> Vec<f64> {
        let n = self.n;
        let mut mat = vec![0.0_f64; n * n];
        for i in 0..n {
            let d_i0 = self.norm_distmat[i * n];
            let src = &self.norm_distmat[i * n..(i + 1) * n];
            let dst = &mut mat[i * n..(i + 1) * n];
            for j in 0..n {
                dst[j] = 0.5 * (d_i0 + self.norm_distmat[j] - src[j]);
            }
        }
        double_demean(&mut mat, n);
        mat
    }

    fn exp_kernel(&self, c: f64) -> Vec<f64> {
        let n = self.n;
        let mut l_exp = vec![0.0_f64; self.lower_tri.len()];
        linalg::exp_scaled(&self.lower_tri, -c, &mut l_exp);
        let mut kern = vec![0.0_f64; n * n];
        let mut idx = 0usize;
        for i in 0..n {
            kern[i * n + i] = 1.0;
            for j in 0..i {
                let v = l_exp[idx];
                idx += 1;
                kern[i * n + j] = v;
                kern[j * n + i] = v;
            }
        }
        kern
    }

    fn projected_omega_from_mr(&self, c: f64, mr: &[f64], q: usize) -> Vec<f64> {
        let n = self.n;
        let sig = self.exp_kernel(c);
        let mut sig_mr = vec![0.0_f64; n * q];
        linalg::multiply(&sig, mr, &mut sig_mr, n, n, q, false);
        let mut om = vec![0.0_f64; q * q];
        linalg::multiply(mr, &sig_mr, &mut om, q, n, q, true);
        symmetrize(&mut om, q);
        om
    }

    fn projected_omega_grid(&self, c_vals: &[f64], mr: &[f64], q: usize) -> Vec<Vec<f64>> {
        let nc = c_vals.len();
        if nc == 0 {
            return Vec::new();
        }
        let n = self.n;
        let qq = q * q;
        let blk = OMEGA_GRID_BLK;

        let mut pairs = Vec::new();
        for i0 in (0..n).step_by(blk) {
            for j0 in (i0..n).step_by(blk) {
                pairs.push((i0, j0));
            }
        }

        let eval_pair_into = |i0: usize,
                              j0: usize,
                              c_slice: &[f64],
                              d_blk: &mut Vec<f64>,
                              s_blk: &mut Vec<f64>,
                              sr_blk: &mut Vec<f64>,
                              delta: &mut Vec<f64>,
                              dst_flat: &mut [f64]| {
            let bi = blk.min(n - i0);
            let bj = blk.min(n - j0);
            let len = bi * bj;
            s_blk.resize(len, 0.0);
            sr_blk.resize(bi * q, 0.0);
            delta.resize(qq, 0.0);
            let r_i = &mr[i0 * q..(i0 + bi) * q];
            let r_j = &mr[j0 * q..(j0 + bj) * q];

            if i0 == j0 {
                let tri_len = bi * (bi + 1) / 2;
                d_blk.resize(tri_len * 2, 0.0);
                {
                    let mut k = 0usize;
                    for r in 0..bi {
                        let row = &self.norm_distmat[(i0 + r) * n + i0..(i0 + r) * n + i0 + bi];
                        for col in r..bi {
                            d_blk[k] = row[col];
                            k += 1;
                        }
                    }
                }
                for (m, &c) in c_slice.iter().enumerate() {
                    let (d_tri, s_tri) = d_blk.split_at_mut(tri_len);
                    linalg::exp_scaled(d_tri, -c, s_tri);
                    let mut k = 0usize;
                    for r in 0..bi {
                        s_blk[r * bi + r] = s_tri[k];
                        k += 1;
                        for col in (r + 1)..bi {
                            let v = s_tri[k];
                            k += 1;
                            s_blk[r * bi + col] = v;
                            s_blk[col * bi + r] = v;
                        }
                    }
                    linalg::multiply(&s_blk[..len], r_j, &mut sr_blk[..bi * q], bi, bj, q, false);
                    linalg::multiply(r_i, &sr_blk[..bi * q], &mut delta[..qq], q, bi, q, true);
                    let om_m = &mut dst_flat[m * qq..(m + 1) * qq];
                    for k_idx in 0..qq {
                        om_m[k_idx] += delta[k_idx];
                    }
                }
            } else {
                d_blk.resize(len, 0.0);
                for r in 0..bi {
                    d_blk[r * bj..(r + 1) * bj].copy_from_slice(
                        &self.norm_distmat[(i0 + r) * n + j0..(i0 + r) * n + j0 + bj],
                    );
                }
                for (m, &c) in c_slice.iter().enumerate() {
                    linalg::exp_scaled(&d_blk[..len], -c, &mut s_blk[..len]);
                    linalg::multiply(&s_blk[..len], r_j, &mut sr_blk[..bi * q], bi, bj, q, false);
                    linalg::multiply(r_i, &sr_blk[..bi * q], &mut delta[..qq], q, bi, q, true);
                    let om_m = &mut dst_flat[m * qq..(m + 1) * qq];
                    for a in 0..q {
                        for b in 0..q {
                            om_m[a * q + b] += delta[a * q + b] + delta[b * q + a];
                        }
                    }
                }
            }
        };

        let mut flat_out = vec![0.0_f64; nc * qq];
        if !is_serial() && pairs.len() >= 4 && pairs.len() * nc * qq <= MAX_PAIR_CONTRIB_DOUBLES {
            let stride = nc * qq;
            let mut pair_contribs = vec![0.0_f64; pairs.len() * stride];
            pair_contribs
                .par_chunks_mut(stride)
                .zip(pairs.par_iter())
                .for_each_init(
                    || (Vec::new(), Vec::new(), Vec::new(), Vec::new()),
                    |(d_blk, s_blk, sr_blk, delta), (dst, &(i0, j0))| {
                        eval_pair_into(i0, j0, c_vals, d_blk, s_blk, sr_blk, delta, dst);
                    },
                );
            for p_idx in 0..pairs.len() {
                let src = &pair_contribs[p_idx * stride..(p_idx + 1) * stride];
                for k in 0..stride {
                    flat_out[k] += src[k];
                }
            }
        } else if !is_serial() && nc >= 4 {
            let nthreads = rayon::current_num_threads().max(1);
            let chunk_len = nc.div_ceil(nthreads).max(1);
            flat_out
                .par_chunks_mut(chunk_len * qq)
                .zip(c_vals.par_chunks(chunk_len))
                .for_each_init(
                    || (Vec::new(), Vec::new(), Vec::new(), Vec::new()),
                    |(d_blk, s_blk, sr_blk, delta), (dst, c_chunk)| {
                        for &(i0, j0) in &pairs {
                            eval_pair_into(i0, j0, c_chunk, d_blk, s_blk, sr_blk, delta, dst);
                        }
                    },
                );
        } else {
            let mut d_blk = Vec::new();
            let mut s_blk = Vec::new();
            let mut sr_blk = Vec::new();
            let mut delta = Vec::new();
            for &(i0, j0) in &pairs {
                eval_pair_into(
                    i0,
                    j0,
                    c_vals,
                    &mut d_blk,
                    &mut s_blk,
                    &mut sr_blk,
                    &mut delta,
                    &mut flat_out,
                );
            }
        }

        flat_out
            .chunks_exact(qq)
            .map(|slice| {
                let mut om = slice.to_vec();
                symmetrize(&mut om, q);
                om
            })
            .collect()
    }

    fn residual_sigma_bm(&self, c_bm: f64, qx: &[f64], k_x: usize) -> Vec<f64> {
        let n = self.n;
        // Residualization removes the constant matrix. Start with the
        // exponential deviations so that the near-unit-root signal is never
        // computed by subtracting order-one projections.
        let mut sig = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..i {
                let value = (-c_bm * self.norm_distmat[i * n + j]).exp_m1();
                sig[i * n + j] = value;
                sig[j * n + i] = value;
            }
        }
        // Since M = I - Q_x Q_x' and sig is symmetric:
        // M * sig * M = sig - Q_x * (Q_x' * sig) - (sig * Q_x) * Q_x' + Q_x * (Q_x' * sig * Q_x) * Q_x'.
        let mut qt_sig = vec![0.0_f64; k_x * n];
        linalg::multiply(qx, &sig, &mut qt_sig, k_x, n, n, true);
        let mut qt_sig_q = vec![0.0_f64; k_x * k_x];
        linalg::multiply(&qt_sig, qx, &mut qt_sig_q, k_x, n, k_x, false);

        // b_adj = Q_x' * sig - (Q_x' * sig * Q_x) * Q_x'  (k_x x n)
        let mut adj = vec![0.0_f64; k_x * n];
        // Build Q_x' directly via multiply(qt_sig_q, Q_x')? Use b[l, i] * q[j, l]:
        // Store (qt_sig_q * Q_x') in adj:
        let mut qx_t = vec![0.0_f64; k_x * n];
        for i in 0..n {
            for l in 0..k_x {
                qx_t[l * n + i] = qx[i * k_x + l];
            }
        }
        linalg::multiply(&qt_sig_q, &qx_t, &mut adj, k_x, k_x, n, false);

        for i in 0..n {
            for j in 0..=i {
                let mut sub = 0.0_f64;
                for l in 0..k_x {
                    sub += qx[i * k_x + l] * qt_sig[l * n + j]
                        + qx[j * k_x + l] * (qt_sig[l * n + i] - adj[l * n + i]);
                }
                let val = sig[i * n + j] - sub;
                sig[i * n + j] = val;
                sig[j * n + i] = val;
            }
        }
        sig
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SpurTestRequest<'a> {
    pub coords: &'a [f64],
    pub n: usize,
    pub d: usize,
    pub metric: Metric,
    pub y: &'a [f64],
    pub n_y: usize,
    pub x_in: Option<&'a [f64]>,
    pub k_x: usize,
    pub emat: &'a [f64],
    pub q: usize,
    pub nrep: usize,
    pub nthreads: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct SpurI0Output {
    pub lr: Vec<f64>,
    pub pvalue: Vec<f64>,
    pub cvalue: Vec<f64>,
    pub ha_param: f64,
    pub pvalue_mat: Vec<f64>,
    pub cvalue_mat: Vec<f64>,
    pub rho_grid: Vec<f64>,
    pub tied: bool,
}

#[derive(Clone, Debug)]
pub struct SpurI1Output {
    pub lr: Vec<f64>,
    pub pvalue: Vec<f64>,
    pub cvalue: Vec<f64>,
    pub ha_param: f64,
    pub tied: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct SpurPipelineRequest<'a> {
    pub coords: &'a [f64],
    pub n: usize,
    pub d: usize,
    pub metric: Metric,
    pub y: &'a [f64],
    pub x_in: &'a [f64],
    pub k_x: usize,
    pub emat: &'a [f64],
    pub q: usize,
    pub nrep: usize,
    pub transform_vars: &'a [f64],
    pub n_transform_vars: usize,
    pub nthreads: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct SpurPipelineOutput {
    pub i0: SpurI0Output,
    pub i1: SpurI1Output,
    pub i0resid: SpurI0Output,
    pub i1resid: SpurI1Output,
    pub transformed_vars: Vec<f64>,
}

fn deterministic_mean_exp_mode(vd: &[f64], c: f64, allow_parallel: bool) -> f64 {
    if vd.is_empty() {
        return f64::NAN;
    }
    let chunk_sums: Vec<f64> = if allow_parallel && vd.len() >= 4 * SUM_CHUNK_SIZE {
        vd.par_chunks(SUM_CHUNK_SIZE)
            .map_init(Vec::new, |scratch, chunk| {
                scratch.resize(chunk.len(), 0.0);
                linalg::exp_scaled(chunk, -c, scratch);
                scratch.iter().sum::<f64>()
            })
            .collect()
    } else {
        let mut scratch = Vec::with_capacity(vd.len().min(SUM_CHUNK_SIZE));
        vd.chunks(SUM_CHUNK_SIZE)
            .map(|chunk| {
                scratch.resize(chunk.len(), 0.0);
                linalg::exp_scaled(chunk, -c, &mut scratch);
                scratch.iter().sum::<f64>()
            })
            .collect()
    };
    chunk_sums.into_iter().sum::<f64>() / vd.len() as f64
}

#[cfg(test)]
fn deterministic_mean_exp(vd: &[f64], c: f64) -> f64 {
    deterministic_mean_exp_mode(vd, c, !is_serial())
}

fn solve_cbar_cached(
    vd: &[f64],
    rhobar: f64,
    allow_parallel: bool,
    cache: &mut Vec<(u64, f64)>,
) -> Result<f64, CoreError> {
    let mut eval_c = |c: f64| -> f64 {
        let bits = c.to_bits();
        if let Some(&(_, v)) = cache.iter().find(|&&(cb, _)| cb == bits) {
            return v;
        }
        let v = deterministic_mean_exp_mode(vd, c, allow_parallel);
        cache.push((bits, v));
        v
    };

    let mut c0 = 10.0_f64;
    let mut c1 = 10.0_f64;

    let mut it = 0usize;
    loop {
        let v = eval_c(c0);
        if v > rhobar {
            break;
        }
        c1 = c0;
        c0 *= 0.5;
        it += 1;
        if it > MAX_CBAR_BRACKET_STEPS || c0 == 0.0 {
            return Err(CoreError::calibration("rhobar too large for get_cbar"));
        }
    }

    it = 0;
    loop {
        let v = eval_c(c1);
        if v < rhobar {
            break;
        }
        if c1 >= 10_240.0 {
            return if v == rhobar {
                Ok(c1)
            } else {
                solve_c0_from_avc(vd, rhobar)
            };
        }
        c0 = c1;
        c1 *= 2.0;
        it += 1;
        if it > MAX_CBAR_BRACKET_STEPS || !c1.is_finite() {
            return Err(CoreError::calibration("rhobar too small for get_cbar"));
        }
    }

    while (c1 - c0) > 0.001 {
        let cm = (c0 * c1).sqrt();
        let v = eval_c(cm);
        if v < rhobar {
            c1 = cm;
        } else {
            c0 = cm;
        }
    }

    Ok((c0 * c1).sqrt())
}

fn double_demean(mat: &mut [f64], n: usize) {
    let mut col_means = vec![0.0_f64; n];
    for row in mat.chunks_exact(n) {
        for j in 0..n {
            col_means[j] += row[j];
        }
    }
    for m in &mut col_means {
        *m /= n as f64;
    }
    for row in mat.chunks_exact_mut(n) {
        let mut row_sum = 0.0_f64;
        for j in 0..n {
            row[j] -= col_means[j];
            row_sum += row[j];
        }
        let row_mean = row_sum / n as f64;
        for v in row.iter_mut() {
            *v -= row_mean;
        }
    }
    symmetrize(mat, n);
}

fn symmetrize(mat: &mut [f64], n: usize) {
    for i in 0..n {
        for j in 0..i {
            let avg = 0.5 * (mat[i * n + j] + mat[j * n + i]);
            mat[i * n + j] = avg;
            mat[j * n + i] = avg;
        }
    }
}

fn center_columns(out: &mut [f64], n: usize, q: usize) {
    let mut col_means = vec![0.0_f64; q];
    for row in out.chunks_exact(q) {
        for j in 0..q {
            col_means[j] += row[j];
        }
    }
    for m in &mut col_means {
        *m /= n as f64;
    }
    for row in out.chunks_exact_mut(q) {
        for j in 0..q {
            row[j] -= col_means[j];
        }
    }
}

fn residual_project_matrix(src: &[f64], qx: &[f64], n: usize, k_x: usize, q: usize) -> Vec<f64> {
    let mut qt_src = vec![0.0_f64; k_x * q];
    linalg::multiply(qx, src, &mut qt_src, k_x, n, q, true);
    let mut q_qt_src = vec![0.0_f64; n * q];
    linalg::multiply(qx, &qt_src, &mut q_qt_src, n, k_x, q, false);
    src.iter().zip(&q_qt_src).map(|(&s, &p)| s - p).collect()
}

fn orthonormal_design_basis(x_in: &[f64], n: usize, k_x: usize) -> Result<Vec<f64>, CoreError> {
    if k_x == 0 || k_x >= n || x_in.len() != checked_len(&[n, k_x])? {
        return Err(CoreError::invalid_input(
            "invalid residual test design matrix dimensions",
        ));
    }
    if x_in.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input(
            "residual test design matrix must be finite",
        ));
    }
    let mut q_mat = vec![0.0_f64; n * k_x];
    // Scale each column to unit L2 norm first so rank detection and projection are
    // invariant to regressor units across arbitrarily small or large column scales.
    for j in 0..k_x {
        let max_abs = (0..n)
            .map(|i| x_in[i * k_x + j].abs())
            .fold(0.0_f64, f64::max);
        if !max_abs.is_finite() || max_abs == 0.0 {
            return Err(CoreError::invalid_input(
                "design matrix is rank-deficient; residual tests require full column rank",
            ));
        }
        let scaled_norm = (0..n)
            .map(|i| (x_in[i * k_x + j] / max_abs).powi(2))
            .sum::<f64>()
            .sqrt();
        if !scaled_norm.is_finite() || scaled_norm == 0.0 {
            return Err(CoreError::invalid_input(
                "design matrix is rank-deficient; residual tests require full column rank",
            ));
        }
        let inv = scaled_norm.recip();
        for i in 0..n {
            q_mat[i * k_x + j] = (x_in[i * k_x + j] / max_abs) * inv;
        }
    }
    for j in 0..k_x {
        for _ in 0..2 {
            for prev in 0..j {
                let dot: f64 = (0..n)
                    .map(|i| q_mat[i * k_x + prev] * q_mat[i * k_x + j])
                    .sum();
                for i in 0..n {
                    q_mat[i * k_x + j] -= dot * q_mat[i * k_x + prev];
                }
            }
        }
        let norm = (0..n)
            .map(|i| q_mat[i * k_x + j].powi(2))
            .sum::<f64>()
            .sqrt();
        if !norm.is_finite() || norm <= 1e-10 {
            return Err(CoreError::invalid_input(
                "design matrix is rank-deficient; residual tests require full column rank",
            ));
        }
        let inv = norm.recip();
        for i in 0..n {
            q_mat[i * k_x + j] *= inv;
        }
    }
    Ok(q_mat)
}

struct CholFactor {
    /// Upper triangular factor `U` with `U' * U = A` (`q x q`, row-major).
    upper: Vec<f64>,
    /// Lower triangular factor `L = U'` with `L * L' = A` (`q x q`, row-major).
    lower: Vec<f64>,
    /// `0.5 * ln(det(A)) = sum_i ln(U[i, i])`.
    log_sqrt_det: f64,
}

fn cholesky_factor(a: &[f64], q: usize) -> Result<CholFactor, CoreError> {
    let mut lower = vec![0.0_f64; q * q];
    let mut log_sqrt_det = 0.0_f64;
    for i in 0..q {
        for j in 0..=i {
            let mut sum = a[i * q + j];
            for k in 0..j {
                sum -= lower[i * q + k] * lower[j * q + k];
            }
            if i == j {
                if !sum.is_finite() || sum <= 0.0 {
                    return Err(CoreError::degenerate(
                        "SPUR covariance matrix is not positive definite",
                    ));
                }
                let diag = sum.sqrt();
                lower[i * q + j] = diag;
                log_sqrt_det += diag.ln();
            } else {
                lower[i * q + j] = sum / lower[j * q + j];
            }
        }
    }
    let mut upper = vec![0.0_f64; q * q];
    for i in 0..q {
        for j in 0..=i {
            upper[j * q + i] = lower[i * q + j];
        }
    }
    Ok(CholFactor {
        upper,
        lower,
        log_sqrt_det,
    })
}

fn invert_spd_from_lower(lower: &[f64], q: usize) -> Vec<f64> {
    // Invert lower triangular `L`, then `A^{-1} = L^{-T} L^{-1}`.
    let mut linv = vec![0.0_f64; q * q];
    for j in 0..q {
        linv[j * q + j] = lower[j * q + j].recip();
        for i in (j + 1)..q {
            let mut sum = 0.0_f64;
            for k in j..i {
                sum -= lower[i * q + k] * linv[k * q + j];
            }
            linv[i * q + j] = sum / lower[i * q + i];
        }
    }
    let mut inv = vec![0.0_f64; q * q];
    linalg::multiply(&linv, &linv, &mut inv, q, q, q, true);
    symmetrize(&mut inv, q);
    inv
}

struct OmegaPairFactor {
    ch_om: CholFactor,
    om_inv: Vec<f64>,
    ch_omi: CholFactor,
}

fn factor_omega(om: &[f64], q: usize) -> Result<OmegaPairFactor, CoreError> {
    let ch_om = cholesky_factor(om, q)?;
    let om_inv = invert_spd_from_lower(&ch_om.lower, q);
    let ch_omi = cholesky_factor(&om_inv, q)?;
    Ok(OmegaPairFactor {
        ch_om,
        om_inv,
        ch_omi,
    })
}

// Covariance scale is a nuisance: only trace-normalized shapes are informative.
fn covariance_shapes_differ(a: &[f64], b: &[f64], q: usize) -> Result<bool, CoreError> {
    let ta: f64 = (0..q).map(|i| a[i * q + i]).sum();
    let tb: f64 = (0..q).map(|i| b[i * q + i]).sum();
    if !ta.is_finite() || ta <= 0.0 || !tb.is_finite() || tb <= 0.0 {
        return Err(CoreError::degenerate(
            "SPUR covariance has non-positive trace",
        ));
    }
    Ok(a.iter()
        .zip(b)
        .any(|(&x, &y)| (x / ta - y / tb).abs() > 1e-12))
}

fn quantile_type7_inplace(vals: &mut [f64], prob: f64) -> f64 {
    let m = vals.len();
    if m == 0 {
        return f64::NAN;
    }
    if m == 1 {
        return vals[0];
    }
    let pos = prob.clamp(0.0, 1.0) * ((m - 1) as f64);
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let cmp = |a: &f64, b: &f64| a.total_cmp(b);
    let (_, lo_val, _) = vals.select_nth_unstable_by(lo, cmp);
    let lo_v = *lo_val;
    if lo == hi {
        return lo_v;
    }
    let (_, hi_val, _) = vals[lo + 1..].select_nth_unstable_by(hi - (lo + 1), cmp);
    let frac = pos - (lo as f64);
    lo_v + frac * (*hi_val - lo_v)
}

fn quantiles_type7(vals: &mut [f64], probs: &[f64]) -> Vec<f64> {
    let m = vals.len();
    if m == 0 {
        return vec![f64::NAN; probs.len()];
    }
    if m == 1 {
        return vec![vals[0]; probs.len()];
    }
    let cmp = |a: &f64, b: &f64| a.total_cmp(b);
    struct Req {
        lo: usize,
        hi: usize,
        frac: f64,
    }
    let mut reqs = Vec::with_capacity(probs.len());
    let mut needed = Vec::with_capacity(probs.len() * 2);
    for &prob in probs {
        let pos = prob.clamp(0.0, 1.0) * ((m - 1) as f64);
        let lo = pos.floor() as usize;
        let hi = pos.ceil() as usize;
        let frac = pos - (lo as f64);
        reqs.push(Req { lo, hi, frac });
        needed.push(lo);
        if hi != lo {
            needed.push(hi);
        }
    }
    needed.sort_unstable();
    needed.dedup();
    let mut start = 0usize;
    for &idx in &needed {
        if idx >= start {
            vals[start..].select_nth_unstable_by(idx - start, cmp);
            start = idx + 1;
        }
    }
    reqs.into_iter()
        .map(|r| vals[r.lo] + r.frac * (vals[r.hi] - vals[r.lo]))
        .collect()
}

struct PowQfBase<'a> {
    q: usize,
    nrep: usize,
    emat: &'a [f64],
    qe: Vec<f64>,
    ch_om0_lower: Vec<f64>,
    ch_om0i_upper: Vec<f64>,
}

impl<'a> PowQfBase<'a> {
    fn new(om0: &[f64], emat: &'a [f64], q: usize, nrep: usize) -> Result<Self, CoreError> {
        let base = factor_omega(om0, q)?;
        let mut qe = vec![0.0_f64; nrep];
        for r in 0..q {
            let row = &emat[r * nrep..(r + 1) * nrep];
            for k in 0..nrep {
                qe[k] += row[k] * row[k];
            }
        }
        Ok(Self {
            q,
            nrep,
            emat,
            qe,
            ch_om0_lower: base.ch_om.lower,
            ch_om0i_upper: base.ch_omi.upper,
        })
    }

    fn eval(&self, om1: &[f64]) -> Result<f64, CoreError> {
        let q = self.q;
        let nrep = self.nrep;
        let alt = factor_omega(om1, q)?;
        let mut ho = vec![0.0_f64; q * q];
        let mut ha = vec![0.0_f64; q * q];
        linalg::multiply(
            &alt.ch_omi.upper,
            &self.ch_om0_lower,
            &mut ho,
            q,
            q,
            q,
            false,
        );
        linalg::multiply(
            &self.ch_om0i_upper,
            &alt.ch_om.lower,
            &mut ha,
            q,
            q,
            q,
            false,
        );

        let mut qa_o = vec![0.0_f64; nrep];
        let mut qo_a = vec![0.0_f64; nrep];
        let mut ya_buf = vec![0.0_f64; q * TILE_COLS.min(nrep)];
        let mut yo_buf = vec![0.0_f64; q * TILE_COLS.min(nrep)];
        let mut lr_o = vec![0.0_f64; nrep];
        let mut tile_e = if nrep > TILE_COLS {
            vec![0.0_f64; q * TILE_COLS]
        } else {
            Vec::new()
        };

        for start in (0..nrep).step_by(TILE_COLS) {
            let cols = (nrep - start).min(TILE_COLS);
            let e_slice: &[f64] = if cols == nrep {
                &self.emat[..q * cols]
            } else {
                for r in 0..q {
                    tile_e[r * cols..(r + 1) * cols]
                        .copy_from_slice(&self.emat[r * nrep + start..r * nrep + start + cols]);
                }
                &tile_e[..q * cols]
            };
            linalg::multiply(&ho, e_slice, &mut ya_buf[..q * cols], q, q, cols, false);
            linalg::multiply(&ha, e_slice, &mut yo_buf[..q * cols], q, q, cols, false);
            for r in 0..q {
                let ra = &ya_buf[r * cols..(r + 1) * cols];
                let ro = &yo_buf[r * cols..(r + 1) * cols];
                for c in 0..cols {
                    qa_o[start + c] += ra[c] * ra[c];
                    qo_a[start + c] += ro[c] * ro[c];
                }
            }
        }

        for k in 0..nrep {
            lr_o[k] = self.qe[k] / qa_o[k];
        }
        let cv = quantile_type7_inplace(&mut lr_o, 0.95);
        let exceed = (0..nrep).filter(|&k| (qo_a[k] / self.qe[k]) > cv).count();
        Ok((exceed as f64) / (nrep as f64))
    }
}

fn select_leading_basis(sigdm: &[f64], n: usize, q: usize) -> Result<(Vec<f64>, bool), CoreError> {
    let eigen = basis::symmetric_eigenpairs(sigdm, n, q)?;
    let lead = eigen.values.first().copied().unwrap_or(0.0);
    let tail = eigen.values.get(q - 1).copied().unwrap_or(0.0);
    if !lead.is_finite() || lead <= 0.0 || !tail.is_finite() || tail <= 1e-12 * lead {
        return Err(CoreError::invalid_input(
            "SPUR covariance basis has insufficient positive rank for requested q",
        ));
    }
    let mut r_mat = vec![0.0_f64; n * q];
    for j in 0..q {
        let src = &eigen.vectors[j * n..(j + 1) * n];
        // Small solver roundoff can swap equal-magnitude entries. Use the
        // earliest entry within a relative tolerance of the maximum.
        let max_abs = src.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let tolerance = 1e-12 * max_abs;
        let max_idx = src
            .iter()
            .position(|v| max_abs - v.abs() <= tolerance)
            .expect("an eigenvector has a largest-magnitude entry");
        let sign = if src[max_idx] < 0.0 { -1.0 } else { 1.0 };
        for i in 0..n {
            r_mat[i * q + j] = src[i] * sign;
        }
    }
    Ok((r_mat, eigen.tied))
}

fn project_scaled_y_cols(
    y_scaled: &[f64],
    y_norms: &[f64],
    n: usize,
    n_y: usize,
    basis_r: &[f64],
    q: usize,
) -> Result<Vec<f64>, CoreError> {
    let mut p_mat = vec![0.0_f64; q * n_y];
    linalg::multiply(basis_r, y_scaled, &mut p_mat, q, n, n_y, true);
    for j in 0..n_y {
        let p_max = (0..q)
            .map(|r| p_mat[r * n_y + j].abs())
            .fold(0.0_f64, f64::max);
        if !p_max.is_finite() || p_max == 0.0 {
            return Err(CoreError::degenerate("projected SPUR vector has zero norm"));
        }
        let p_norm = p_max
            * (0..q)
                .map(|r| (p_mat[r * n_y + j] / p_max).powi(2))
                .sum::<f64>()
                .sqrt();
        if !p_norm.is_finite() || p_norm <= 1e-12 * y_norms[j] {
            return Err(CoreError::degenerate("projected SPUR vector has zero norm"));
        }
    }
    Ok(p_mat)
}

fn prepare_y_projections(
    y: &[f64],
    n: usize,
    n_y: usize,
    basis_r: &[f64],
    q: usize,
) -> Result<Vec<f64>, CoreError> {
    if n_y == 0 || y.len() != checked_len(&[n, n_y])? {
        return Err(CoreError::invalid_input("invalid response Y dimensions"));
    }
    if y.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input("response Y values must be finite"));
    }
    let mut y_min = vec![f64::INFINITY; n_y];
    let mut y_max = vec![f64::NEG_INFINITY; n_y];
    for row in y.chunks_exact(n_y) {
        for j in 0..n_y {
            let v = row[j];
            y_min[j] = y_min[j].min(v);
            y_max[j] = y_max[j].max(v);
        }
    }
    for j in 0..n_y {
        if y_max[j] - y_min[j] == 0.0 {
            return Err(CoreError::degenerate(
                "response Y has zero variation; SPUR inference requires non-constant Y",
            ));
        }
    }
    let mut y_dm = y.to_vec();
    center_columns(&mut y_dm, n, n_y);
    let mut y_dm_max = vec![0.0_f64; n_y];
    for row in y_dm.chunks_exact(n_y) {
        for j in 0..n_y {
            y_dm_max[j] = y_dm_max[j].max(row[j].abs());
        }
    }
    let mut y_dm_norm = vec![0.0_f64; n_y];
    for j in 0..n_y {
        let max_d = y_dm_max[j];
        if !max_d.is_finite() || max_d == 0.0 {
            return Err(CoreError::degenerate(
                "response Y has zero variation; SPUR inference requires non-constant Y",
            ));
        }
        for i in 0..n {
            y_dm[i * n_y + j] /= max_d;
        }
        y_dm_norm[j] = (0..n)
            .map(|i| y_dm[i * n_y + j].powi(2))
            .sum::<f64>()
            .sqrt();
    }
    project_scaled_y_cols(&y_dm, &y_dm_norm, n, n_y, basis_r, q)
}

#[allow(clippy::too_many_arguments)]
fn prepare_residual_y_projections(
    y: &[f64],
    n: usize,
    n_y: usize,
    basis_r: &[f64],
    qx: &[f64],
    k_x: usize,
    q: usize,
) -> Result<Vec<f64>, CoreError> {
    if n_y == 0 || y.len() != checked_len(&[n, n_y])? {
        return Err(CoreError::invalid_input("invalid response Y dimensions"));
    }
    if y.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input("response Y values must be finite"));
    }
    let mut y_res = vec![0.0_f64; n * n_y];
    let mut y_res_norm = vec![0.0_f64; n_y];
    let mut col = vec![0.0_f64; n];

    for j in 0..n_y {
        let mut col_min = f64::INFINITY;
        let mut col_max_val = f64::NEG_INFINITY;
        let mut raw_max = 0.0_f64;
        for i in 0..n {
            let v = y[i * n_y + j];
            col_min = col_min.min(v);
            col_max_val = col_max_val.max(v);
            raw_max = raw_max.max(v.abs());
        }
        if col_max_val == col_min {
            return Err(CoreError::degenerate(
                "response Y has zero variation; SPUR inference requires non-constant Y",
            ));
        }
        // The validated nuisance design spans constants. Subtract in the
        // original units to preserve small differences near a large offset.
        let anchor = y[j];
        for i in 0..n {
            col[i] = y[i * n_y + j] - anchor;
        }
        let mut scale = col.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if !scale.is_finite() {
            let scaled_anchor = anchor / raw_max;
            for i in 0..n {
                col[i] = y[i * n_y + j] / raw_max - scaled_anchor;
            }
            scale = col.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        }
        for v in &mut col {
            *v /= scale;
        }
        let col_norm = col.iter().map(|v| v * v).sum::<f64>().sqrt();
        let u_once = residual_project_matrix(&col, qx, n, k_x, 1);
        let u = residual_project_matrix(&u_once, qx, n, k_x, 1);
        let u_max = u.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if !u_max.is_finite() || u_max == 0.0 {
            return Err(CoreError::degenerate(
                "residualized response Y has zero variation; SPUR inference requires variation orthogonal to X",
            ));
        }
        let u_norm = u_max * u.iter().map(|v| (v / u_max).powi(2)).sum::<f64>().sqrt();
        if !u_norm.is_finite() || u_norm <= 1e-12 * col_norm {
            return Err(CoreError::degenerate(
                "residualized response Y has zero variation; SPUR inference requires variation orthogonal to X",
            ));
        }
        for i in 0..n {
            y_res[i * n_y + j] = u[i] / u_max;
        }
        y_res_norm[j] = u_norm / u_max;
    }

    project_scaled_y_cols(&y_res, &y_res_norm, n, n_y, basis_r, q)
}

fn validate_emat(emat: &[f64], n: usize, q: usize, nrep: usize) -> Result<(), CoreError> {
    if q < 2 || q >= n {
        return Err(CoreError::invalid_input(
            "SPUR inference requires 2 <= q < n; q = 1 has no covariance-shape information",
        ));
    }
    if nrep == 0 || emat.len() != checked_len(&[q, nrep])? {
        return Err(CoreError::invalid_input(
            "emat must have dimensions q x nrep with nrep >= 1",
        ));
    }
    if emat.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input("emat entries must be finite"));
    }
    Ok(())
}

fn calibrate_ha_i0(
    om_ho: &[f64],
    om_i0: &[f64],
    om_bm: &[f64],
    emat: &[f64],
    q: usize,
    nrep: usize,
) -> Result<f64, CoreError> {
    let base = PowQfBase::new(om_ho, emat, q, nrep)?;
    let mix = |g: f64| -> Vec<f64> { om_i0.iter().zip(om_bm).map(|(&a, &b)| a + g * b).collect() };

    let mut pow = 1.0_f64;
    let mut gtry = 1.0_f64;
    let mut g = gtry;
    let mut it = 0usize;
    while pow > 0.5 {
        g = gtry;
        pow = base.eval(&mix(g))?;
        gtry = g * 0.5;
        it += 1;
        if it >= MAX_HA_BRACKET_STEPS && pow > 0.5 {
            return Err(CoreError::calibration(
                "get_ha_param_i0 step-1 bracketing did not converge",
            ));
        }
    }
    let mut g1 = g;

    pow = 0.0;
    gtry = 30.0;
    it = 0;
    while pow < 0.5 {
        g = gtry;
        pow = base.eval(&mix(g))?;
        gtry = g * 2.0;
        it += 1;
        if it >= MAX_HA_BRACKET_STEPS && pow < 0.5 {
            return Err(CoreError::calibration(
                "get_ha_param_i0 step-2 bracketing did not converge",
            ));
        }
    }
    let mut g2 = g;

    let mut ii = 1usize;
    while (pow - 0.5).abs() > 0.01 {
        g = 0.5 * (g1 + g2);
        pow = base.eval(&mix(g))?;
        if pow > 0.5 {
            g2 = g;
        } else if pow < 0.5 {
            g1 = g;
        }
        ii += 1;
        if ii > 20 {
            break;
        }
    }
    Ok(g)
}

fn calibrate_ha_i1(
    geom: &mut SpurGeometry,
    om_ho: &[f64],
    mr: &[f64],
    emat: &[f64],
    q: usize,
    nrep: usize,
) -> Result<f64, CoreError> {
    let base = PowQfBase::new(om_ho, emat, q, nrep)?;
    let mut pow = 1.0_f64;
    let mut ctry = geom.get_cbar(0.95)?;
    let mut c = ctry;
    let mut it = 0usize;
    while pow > 0.5 {
        c = ctry;
        let om_c = geom.projected_omega_from_mr(c, mr, q);
        pow = base.eval(&om_c)?;
        ctry = c * 0.5;
        it += 1;
        if it >= MAX_HA_BRACKET_STEPS && pow > 0.5 {
            return Err(CoreError::calibration(
                "get_ha_param_i1 step-1 bracketing did not converge",
            ));
        }
    }
    let mut c1 = c;

    pow = 0.0;
    ctry = geom.get_cbar(0.01)?;
    it = 0;
    while pow < 0.5 {
        c = ctry;
        let om_c = geom.projected_omega_from_mr(c, mr, q);
        pow = base.eval(&om_c)?;
        ctry = 2.0 * c;
        it += 1;
        if it >= MAX_HA_BRACKET_STEPS && pow < 0.5 {
            return Err(CoreError::calibration(
                "get_ha_param_i1 step-2 bracketing did not converge",
            ));
        }
    }
    let mut c2 = c;

    let mut ii = 1usize;
    while (pow - 0.5).abs() > 0.01 {
        c = 0.5 * (c1 + c2);
        let om_c = geom.projected_omega_from_mr(c, mr, q);
        pow = base.eval(&om_c)?;
        if pow > 0.5 {
            c2 = c;
        } else if pow < 0.5 {
            c1 = c;
        }
        ii += 1;
        if ii > 20 {
            break;
        }
    }
    Ok(c)
}

fn eval_sim_lr_ratio(
    ch_om_ho_lower: &[f64],
    ch_omi_ho_upper: &[f64],
    ch_omi_ha_upper: &[f64],
    emat: &[f64],
    q: usize,
    nrep: usize,
) -> Vec<f64> {
    let mut a_ho = vec![0.0_f64; q * q];
    let mut a_ha = vec![0.0_f64; q * q];
    linalg::multiply(ch_omi_ho_upper, ch_om_ho_lower, &mut a_ho, q, q, q, false);
    linalg::multiply(ch_omi_ha_upper, ch_om_ho_lower, &mut a_ha, q, q, q, false);

    let mut lr_ho = vec![0.0_f64; nrep];
    let mut tile_e = if nrep > TILE_COLS {
        vec![0.0_f64; q * TILE_COLS]
    } else {
        Vec::new()
    };
    let mut buf_ho = vec![0.0_f64; q * TILE_COLS.min(nrep)];
    let mut buf_ha = vec![0.0_f64; q * TILE_COLS.min(nrep)];
    let mut q_ho = vec![0.0_f64; TILE_COLS.min(nrep)];
    let mut q_ha = vec![0.0_f64; TILE_COLS.min(nrep)];

    for start in (0..nrep).step_by(TILE_COLS) {
        let cols = (nrep - start).min(TILE_COLS);
        let e_slice: &[f64] = if cols == nrep {
            &emat[..q * cols]
        } else {
            for r in 0..q {
                tile_e[r * cols..(r + 1) * cols]
                    .copy_from_slice(&emat[r * nrep + start..r * nrep + start + cols]);
            }
            &tile_e[..q * cols]
        };
        linalg::multiply(&a_ho, e_slice, &mut buf_ho[..q * cols], q, q, cols, false);
        linalg::multiply(&a_ha, e_slice, &mut buf_ha[..q * cols], q, q, cols, false);
        q_ho[..cols].fill(0.0);
        q_ha[..cols].fill(0.0);
        for r in 0..q {
            let rho = &buf_ho[r * cols..(r + 1) * cols];
            let rha = &buf_ha[r * cols..(r + 1) * cols];
            for c in 0..cols {
                q_ho[c] += rho[c] * rho[c];
                q_ha[c] += rha[c] * rha[c];
            }
        }
        for c in 0..cols {
            lr_ho[start + c] = q_ho[c] / q_ha[c];
        }
    }
    lr_ho
}

fn eval_observed_lr(
    om_ho_inv: &[f64],
    om_ha_inv: &[f64],
    ch_omi_ho_upper: &[f64],
    ch_omi_ha_upper: &[f64],
    p_mat: &[f64],
    q: usize,
    n_y: usize,
    use_chol_norm: bool,
) -> Result<Vec<f64>, CoreError> {
    let mut lr = vec![0.0_f64; n_y];
    if use_chol_norm {
        let mut y_ho = vec![0.0_f64; q * n_y];
        let mut y_ha = vec![0.0_f64; q * n_y];
        linalg::multiply(ch_omi_ho_upper, p_mat, &mut y_ho, q, q, n_y, false);
        linalg::multiply(ch_omi_ha_upper, p_mat, &mut y_ha, q, q, n_y, false);
        for j in 0..n_y {
            let num: f64 = (0..q).map(|r| y_ho[r * n_y + j].powi(2)).sum();
            let den: f64 = (0..q).map(|r| y_ha[r * n_y + j].powi(2)).sum();
            if !num.is_finite() || !den.is_finite() || den <= 0.0 {
                return Err(CoreError::degenerate("SPUR LR denominator is non-positive"));
            }
            lr[j] = num / den;
        }
    } else {
        let mut m_ho = vec![0.0_f64; q * n_y];
        let mut m_ha = vec![0.0_f64; q * n_y];
        linalg::multiply(om_ho_inv, p_mat, &mut m_ho, q, q, n_y, false);
        linalg::multiply(om_ha_inv, p_mat, &mut m_ha, q, q, n_y, false);
        for j in 0..n_y {
            let num: f64 = (0..q).map(|r| p_mat[r * n_y + j] * m_ho[r * n_y + j]).sum();
            let den: f64 = (0..q).map(|r| p_mat[r * n_y + j] * m_ha[r * n_y + j]).sum();
            if !num.is_finite() || !den.is_finite() || den <= 0.0 {
                return Err(CoreError::degenerate("SPUR LR denominator is non-positive"));
            }
            lr[j] = num / den;
        }
    }
    Ok(lr)
}

fn build_rho_grid() -> Vec<f64> {
    (0..RHO_GRID_LEN)
        .map(|i| {
            RHO_GRID_MIN + (RHO_GRID_MAX - RHO_GRID_MIN) * (i as f64) / ((RHO_GRID_LEN - 1) as f64)
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn run_spur_i0_from_basis(
    geom: &mut SpurGeometry,
    tied: bool,
    mr: &[f64],
    om_bm: &[f64],
    p_mat: &[f64],
    n_y: usize,
    emat: &[f64],
    q: usize,
    nrep: usize,
) -> Result<SpurI0Output, CoreError> {
    let c_rho = geom.get_cbar(RHO_I0_ALT)?;
    let om_rho = geom.projected_omega_from_mr(c_rho, mr, q);
    if !covariance_shapes_differ(&om_rho, om_bm, q)? {
        return Err(CoreError::degenerate(
            "SPUR experiment is uninformative: null and alternative covariance shapes are proportional",
        ));
    }

    let ha_param = calibrate_ha_i0(&om_rho, &om_rho, om_bm, emat, q, nrep)?;
    let om_ha: Vec<f64> = om_rho
        .iter()
        .zip(om_bm)
        .map(|(&a, &b)| a + ha_param * b)
        .collect();

    let ho_f = factor_omega(&om_rho, q)?;
    let ha_f = factor_omega(&om_ha, q)?;
    let lr = eval_observed_lr(
        &ho_f.om_inv,
        &ha_f.om_inv,
        &ho_f.ch_omi.upper,
        &ha_f.ch_omi.upper,
        p_mat,
        q,
        n_y,
        true,
    )?;

    let rho_grid = build_rho_grid();
    let pos_rhos: Vec<f64> = rho_grid.iter().copied().filter(|&r| r > 0.0).collect();
    let pos_cs = geom.get_cbar_batch(&pos_rhos)?;
    let pos_oms = geom.projected_omega_grid(&pos_cs, mr, q);
    let mut om_grid: Vec<Option<Vec<f64>>> = Vec::with_capacity(RHO_GRID_LEN);
    let mut pos_iter = pos_oms.into_iter();
    for &rho_i in &rho_grid {
        if rho_i > 0.0 {
            om_grid.push(Some(pos_iter.next().expect("missing pos_oms entry")));
        } else {
            om_grid.push(None);
        }
    }

    let probs = [1.0 - SZ_LEVELS[0], 1.0 - SZ_LEVELS[1], 1.0 - SZ_LEVELS[2]];
    let eval_grid_point = |om_opt: Option<&Vec<f64>>| -> Result<(Vec<f64>, Vec<f64>), CoreError> {
        let ch_om_lower = if let Some(om_i) = om_opt {
            cholesky_factor(om_i, q)?.lower
        } else {
            let mut eye = vec![0.0_f64; q * q];
            for d_i in 0..q {
                eye[d_i * q + d_i] = 1.0;
            }
            eye
        };
        let mut lr_ho = eval_sim_lr_ratio(
            &ch_om_lower,
            &ho_f.ch_omi.upper,
            &ha_f.ch_omi.upper,
            emat,
            q,
            nrep,
        );
        let pvals: Vec<f64> = lr
            .iter()
            .map(|&obs_lr| (lr_ho.iter().filter(|&&v| v > obs_lr).count() as f64) / (nrep as f64))
            .collect();
        let cvs = quantiles_type7(&mut lr_ho, &probs);
        Ok((cvs, pvals))
    };

    let grid_results: Vec<(Vec<f64>, Vec<f64>)> = if is_serial() {
        om_grid
            .iter()
            .map(|om_opt| eval_grid_point(om_opt.as_ref()))
            .collect::<Result<Vec<_>, _>>()?
    } else {
        om_grid
            .par_iter()
            .map(|om_opt| eval_grid_point(om_opt.as_ref()))
            .collect::<Result<Vec<_>, _>>()?
    };

    let mut cvalue_mat = vec![0.0_f64; RHO_GRID_LEN * 3];
    let mut pvalue_mat = vec![0.0_f64; RHO_GRID_LEN * n_y];
    let mut cvalue = vec![f64::NEG_INFINITY; 3];
    let mut pvalue = vec![f64::NEG_INFINITY; n_y];

    for (ir, (cvs, pvals)) in grid_results.into_iter().enumerate() {
        for k in 0..3 {
            cvalue_mat[ir * 3 + k] = cvs[k];
            cvalue[k] = cvalue[k].max(cvs[k]);
        }
        for j in 0..n_y {
            pvalue_mat[ir * n_y + j] = pvals[j];
            pvalue[j] = pvalue[j].max(pvals[j]);
        }
    }

    Ok(SpurI0Output {
        lr,
        pvalue,
        cvalue,
        ha_param,
        pvalue_mat,
        cvalue_mat,
        rho_grid,
        tied,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_spur_i1_from_basis(
    geom: &mut SpurGeometry,
    tied: bool,
    mr: &[f64],
    om_ho: &[f64],
    p_mat: &[f64],
    n_y: usize,
    emat: &[f64],
    q: usize,
    nrep: usize,
) -> Result<SpurI1Output, CoreError> {
    let c_lo = geom.get_cbar(0.01)?;
    let c_mid = geom.get_cbar(0.5)?;
    let c_hi = geom.get_cbar(0.95)?;
    let om_lo = geom.projected_omega_from_mr(c_lo, mr, q);
    let om_mid = geom.projected_omega_from_mr(c_mid, mr, q);
    let om_hi = geom.projected_omega_from_mr(c_hi, mr, q);
    if !covariance_shapes_differ(&om_lo, &om_mid, q)?
        && !covariance_shapes_differ(&om_lo, &om_hi, q)?
    {
        return Err(CoreError::degenerate(
            "SPUR experiment is uninformative: covariance shapes are proportional",
        ));
    }
    let ha_param = calibrate_ha_i1(geom, om_ho, mr, emat, q, nrep)?;
    let om_ha = geom.projected_omega_from_mr(ha_param, mr, q);

    let ho_f = factor_omega(om_ho, q)?;
    let ha_f = factor_omega(&om_ha, q)?;
    let mut lr_ho = eval_sim_lr_ratio(
        &ho_f.ch_om.lower,
        &ho_f.ch_omi.upper,
        &ha_f.ch_omi.upper,
        emat,
        q,
        nrep,
    );

    let lr = eval_observed_lr(
        &ho_f.om_inv,
        &ha_f.om_inv,
        &ho_f.ch_omi.upper,
        &ha_f.ch_omi.upper,
        p_mat,
        q,
        n_y,
        false,
    )?;
    let pvalue: Vec<f64> = lr
        .iter()
        .map(|&obs_lr| (lr_ho.iter().filter(|&&v| v > obs_lr).count() as f64) / (nrep as f64))
        .collect();
    let probs = [1.0 - SZ_LEVELS[0], 1.0 - SZ_LEVELS[1], 1.0 - SZ_LEVELS[2]];
    let cvalue = quantiles_type7(&mut lr_ho, &probs);

    Ok(SpurI1Output {
        lr,
        pvalue,
        cvalue,
        ha_param,
        tied,
    })
}

fn prepare_lbm_basis(
    sigdm_bm: &[f64],
    n: usize,
    q: usize,
) -> Result<(Vec<f64>, bool, Vec<f64>, Vec<f64>), CoreError> {
    let (basis_r, tied) = select_leading_basis(sigdm_bm, n, q)?;
    let mut mr = basis_r.clone();
    center_columns(&mut mr, n, q);
    let mut sig_r = vec![0.0_f64; n * q];
    linalg::multiply(sigdm_bm, &basis_r, &mut sig_r, n, n, q, false);
    let mut om_bm = vec![0.0_f64; q * q];
    linalg::multiply(&basis_r, &sig_r, &mut om_bm, q, n, q, true);
    symmetrize(&mut om_bm, q);
    Ok((basis_r, tied, mr, om_bm))
}

struct ResidBasis {
    basis_r: Vec<f64>,
    tied: bool,
    mr: Vec<f64>,
    om_bm: Vec<f64>,
    qx: Vec<f64>,
}

fn prepare_resid_design(
    x_in: &[f64],
    n: usize,
    k_x: usize,
    q: usize,
) -> Result<Vec<f64>, CoreError> {
    let qx = orthonormal_design_basis(x_in, n, k_x)?;
    if q > n - k_x {
        return Err(CoreError::invalid_input(
            "q cannot exceed residual degrees of freedom (n - rank(X))",
        ));
    }
    let unit_constant = vec![(n as f64).sqrt().recip(); n];
    let unit_once = residual_project_matrix(&unit_constant, &qx, n, k_x, 1);
    let unit_residual = residual_project_matrix(&unit_once, &qx, n, k_x, 1);
    if unit_residual.iter().map(|v| v * v).sum::<f64>().sqrt() > 64.0 * f64::EPSILON {
        return Err(CoreError::invalid_input(
            "SPUR residual tests require a design that spans a constant",
        ));
    }
    Ok(qx)
}

fn prepare_resid_basis(
    geom: &mut SpurGeometry,
    qx: Vec<f64>,
    k_x: usize,
    q: usize,
) -> Result<ResidBasis, CoreError> {
    let n = geom.n;
    let c_bm = geom.get_cbar(RHO_BM)?;
    let sigdm_bm = geom.residual_sigma_bm(c_bm, &qx, k_x);
    let (basis_r, tied) = select_leading_basis(&sigdm_bm, n, q)?;
    let mr = residual_project_matrix(&basis_r, &qx, n, k_x, q);
    let mut sig_r = vec![0.0_f64; n * q];
    linalg::multiply(&sigdm_bm, &basis_r, &mut sig_r, n, n, q, false);
    let mut om_bm = vec![0.0_f64; q * q];
    linalg::multiply(&basis_r, &sig_r, &mut om_bm, q, n, q, true);
    symmetrize(&mut om_bm, q);
    Ok(ResidBasis {
        basis_r,
        tied,
        mr,
        om_bm,
        qx,
    })
}

pub fn run_spur_i0(req: &SpurTestRequest<'_>) -> Result<SpurI0Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    run_with_threads(req.nthreads, || {
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_r, tied, mr, om_bm) = prepare_lbm_basis(&geom.sigma_lbm_dm(), req.n, req.q)?;
        let p_mat = prepare_y_projections(req.y, req.n, req.n_y, &basis_r, req.q)?;
        run_spur_i0_from_basis(
            &mut geom, tied, &mr, &om_bm, &p_mat, req.n_y, req.emat, req.q, req.nrep,
        )
    })?
}

pub fn run_spur_i1(req: &SpurTestRequest<'_>) -> Result<SpurI1Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    run_with_threads(req.nthreads, || {
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_r, tied, mr, om_bm) = prepare_lbm_basis(&geom.sigma_lbm_dm(), req.n, req.q)?;
        let p_mat = prepare_y_projections(req.y, req.n, req.n_y, &basis_r, req.q)?;
        run_spur_i1_from_basis(
            &mut geom, tied, &mr, &om_bm, &p_mat, req.n_y, req.emat, req.q, req.nrep,
        )
    })?
}

pub fn run_spur_i0_resid(req: &SpurTestRequest<'_>) -> Result<SpurI0Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    let x_in = req
        .x_in
        .ok_or_else(|| CoreError::invalid_input("residual test requires X_in design matrix"))?;
    run_with_threads(req.nthreads, || {
        let qx = prepare_resid_design(x_in, req.n, req.k_x, req.q)?;
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let resid = prepare_resid_basis(&mut geom, qx, req.k_x, req.q)?;
        let p_mat = prepare_residual_y_projections(
            req.y,
            req.n,
            req.n_y,
            &resid.basis_r,
            &resid.qx,
            req.k_x,
            req.q,
        )?;
        run_spur_i0_from_basis(
            &mut geom,
            resid.tied,
            &resid.mr,
            &resid.om_bm,
            &p_mat,
            req.n_y,
            req.emat,
            req.q,
            req.nrep,
        )
    })?
}

pub fn run_spur_i1_resid(req: &SpurTestRequest<'_>) -> Result<SpurI1Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    let x_in = req
        .x_in
        .ok_or_else(|| CoreError::invalid_input("residual test requires X_in design matrix"))?;
    run_with_threads(req.nthreads, || {
        let qx = prepare_resid_design(x_in, req.n, req.k_x, req.q)?;
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let resid = prepare_resid_basis(&mut geom, qx, req.k_x, req.q)?;
        let p_mat = prepare_residual_y_projections(
            req.y,
            req.n,
            req.n_y,
            &resid.basis_r,
            &resid.qx,
            req.k_x,
            req.q,
        )?;
        run_spur_i1_from_basis(
            &mut geom,
            resid.tied,
            &resid.mr,
            &resid.om_bm,
            &p_mat,
            req.n_y,
            req.emat,
            req.q,
            req.nrep,
        )
    })?
}

impl SpurGeometry {
    fn kernel_dm(&self, c: f64) -> Vec<f64> {
        let mut kern = self.exp_kernel(c);
        double_demean(&mut kern, self.n);
        kern
    }
}

pub fn run_spur_pipeline(req: &SpurPipelineRequest<'_>) -> Result<SpurPipelineOutput, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    if req.y.len() != req.n || req.y.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input(
            "invalid y vector for SPUR pipeline",
        ));
    }
    if req.n_transform_vars == 0
        || req.transform_vars.len() != checked_len(&[req.n, req.n_transform_vars])?
        || req.transform_vars.iter().any(|v| !v.is_finite())
    {
        return Err(CoreError::invalid_input(
            "invalid transform variables for SPUR pipeline",
        ));
    }

    run_with_threads(req.nthreads, || {
        let qx = prepare_resid_design(req.x_in, req.n, req.k_x, req.q)?;
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let sigdm_lbm = geom.sigma_lbm_dm();
        let (basis_lbm, tied_lbm, mr_lbm, om_lbm) = prepare_lbm_basis(&sigdm_lbm, req.n, req.q)?;
        let p_lbm = prepare_y_projections(req.y, req.n, 1, &basis_lbm, req.q)?;
        let i0 = run_spur_i0_from_basis(
            &mut geom, tied_lbm, &mr_lbm, &om_lbm, &p_lbm, 1, req.emat, req.q, req.nrep,
        )?;
        let i1 = run_spur_i1_from_basis(
            &mut geom, tied_lbm, &mr_lbm, &om_lbm, &p_lbm, 1, req.emat, req.q, req.nrep,
        )?;
        let transformed_vars =
            apply_lbmgls(sigdm_lbm, req.n, req.transform_vars, req.n_transform_vars)?;

        let resid = prepare_resid_basis(&mut geom, qx, req.k_x, req.q)?;
        let p_res = prepare_residual_y_projections(
            req.y,
            req.n,
            1,
            &resid.basis_r,
            &resid.qx,
            req.k_x,
            req.q,
        )?;
        let i0resid = run_spur_i0_from_basis(
            &mut geom,
            resid.tied,
            &resid.mr,
            &resid.om_bm,
            &p_res,
            1,
            req.emat,
            req.q,
            req.nrep,
        )?;
        let i1resid = run_spur_i1_from_basis(
            &mut geom,
            resid.tied,
            &resid.mr,
            &resid.om_bm,
            &p_res,
            1,
            req.emat,
            req.q,
            req.nrep,
        )?;

        Ok(SpurPipelineOutput {
            i0,
            i1,
            i0resid,
            i1resid,
            transformed_vars,
        })
    })?
}

#[cfg(test)]
mod tests;
