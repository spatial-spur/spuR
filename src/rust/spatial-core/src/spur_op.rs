//! Whole-operation SPUR diagnostics, half-life confidence intervals, and spatial transforms.

use crate::distance::{distance_matrix, lower_triangle, Distances, Metric};
use crate::threads::{is_serial, run_with_threads};
use crate::{basis, checked_len, linalg, CoreError};
use rayon::prelude::*;

const SUM_CHUNK_SIZE: usize = 32_768;
const TILE_COLS: usize = 1024;
const MAX_HALFLIFE_TILE_BYTES: usize = 4 * 1024 * 1024;
const OMEGA_GRID_BLK: usize = 256;
const MAX_PAIR_CONTRIB_DOUBLES: usize = 16_000_000;
const SZ_LEVELS: [f64; 3] = [0.01, 0.05, 0.10];
const RHO_BM: f64 = 0.999;
const RHO_I0_ALT: f64 = 0.001;
const RHO_GRID_MIN: f64 = 0.0001;
const RHO_GRID_MAX: f64 = 0.03;
const RHO_GRID_LEN: usize = 30;
const HL_HO_SMALL_LEN: usize = 100;
const HL_HO_LARGE_LEN: usize = 30;
const HL_HO_LEN: usize = HL_HO_SMALL_LEN + HL_HO_LARGE_LEN + 1;
const HL_HA_LEN: usize = 50;
const MAX_CBAR_BRACKET_STEPS: usize = 500;
const MAX_HA_BRACKET_STEPS: usize = 200;
const LBM_GLS_EIG_TOL: f64 = 1e-10;

#[derive(Clone, Debug)]
pub struct SpurGeometry {
    pub n: usize,
    pub raw_distmat: Vec<f64>,
    pub max_dist: f64,
    pub norm_distmat: Vec<f64>,
    lower_tri: Vec<f64>,
    cbar_cache: Vec<(f64, f64)>,
    mean_exp_cache: Vec<(u64, f64)>,
}

impl SpurGeometry {
    pub fn from_coords(
        coords: &[f64],
        n: usize,
        d: usize,
        metric: Metric,
    ) -> Result<Self, CoreError> {
        if n < 2 || d == 0 || coords.len() != checked_len(&[n, d])? {
            return Err(CoreError::invalid_input(
                "SPUR requires at least 2 observations and 1 coordinate column",
            ));
        }
        let raw_distmat = distance_matrix(coords, n, d, metric)?;
        let max_dist = raw_distmat.iter().copied().fold(0.0_f64, f64::max);
        if !max_dist.is_finite() || max_dist <= 0.0 {
            return Err(CoreError::invalid_input(
                "all coordinates are identical; normalization requires distinct locations",
            ));
        }
        let norm_distmat: Vec<f64> = raw_distmat.iter().map(|&v| v / max_dist).collect();
        let lower_tri = lower_triangle(&norm_distmat, n);
        Ok(Self {
            n,
            raw_distmat,
            max_dist,
            norm_distmat,
            lower_tri,
            cbar_cache: Vec::new(),
            mean_exp_cache: Vec::new(),
        })
    }

    pub fn from_normalized_distmat(norm_distmat: &[f64], n: usize) -> Result<Self, CoreError> {
        if n < 2 || norm_distmat.len() != checked_len(&[n, n])? {
            return Err(CoreError::invalid_input(
                "invalid normalized distance matrix dimensions",
            ));
        }
        if norm_distmat.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(CoreError::invalid_input(
                "normalized distance matrix must be finite and nonnegative",
            ));
        }
        let norm_vec = norm_distmat.to_vec();
        let lower_tri = lower_triangle(&norm_vec, n);
        Ok(Self {
            n,
            raw_distmat: norm_vec.clone(),
            max_dist: 1.0,
            norm_distmat: norm_vec,
            lower_tri,
            cbar_cache: Vec::new(),
            mean_exp_cache: Vec::new(),
        })
    }

    /// Cached SPUR decay calibration `get_cbar(rhobar, distmat)` matching pinned bisection.
    pub fn get_cbar(&mut self, rhobar: f64) -> Result<f64, CoreError> {
        if !rhobar.is_finite() || rhobar <= 0.0 || rhobar >= 1.0 {
            return Err(CoreError::calibration(
                "rhobar must lie strictly between 0 and 1",
            ));
        }
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
    pub fn get_cbar_batch(&mut self, rhos: &[f64]) -> Result<Vec<f64>, CoreError> {
        for &r in rhos {
            if !r.is_finite() || r <= 0.0 || r >= 1.0 {
                return Err(CoreError::calibration(
                    "rhobar must lie strictly between 0 and 1",
                ));
            }
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
        let diag_zero = (0..n).all(|i| self.norm_distmat[i * n + i] == 0.0);
        if diag_zero {
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
        } else {
            let mut kern = vec![0.0_f64; self.norm_distmat.len()];
            linalg::exp_scaled(&self.norm_distmat, -c, &mut kern);
            kern
        }
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
        let mut sig = self.exp_kernel(c_bm);
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
pub struct SpurHalfLifeRequest<'a> {
    pub coords: &'a [f64],
    pub n: usize,
    pub d: usize,
    pub metric: Metric,
    pub y: &'a [f64],
    pub emat: &'a [f64],
    pub q: usize,
    pub nrep: usize,
    pub level: f64,
    pub nthreads: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct SpurHalfLifeOutput {
    pub ci_lower: f64,
    pub ci_upper: f64,
    pub max_dist: f64,
    pub pv_mat: Vec<f64>,
    pub hl_grid_ho: Vec<f64>,
    pub tied: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpurTransformKind {
    Nn { exclude_zero_distances: bool },
    Iso { radius: f64, strict_lt: bool },
    Cluster,
    LbmGls,
}

#[derive(Clone, Copy, Debug)]
pub struct SpurTransformRequest<'a> {
    pub coords: Option<&'a [f64]>,
    pub n: usize,
    pub d: usize,
    pub metric: Metric,
    pub kind: SpurTransformKind,
    pub cluster_id: Option<&'a [usize]>,
    pub y: Option<&'a [f64]>,
    pub n_y: usize,
    pub return_matrix: bool,
    pub demean_lbmgls: bool,
    pub nthreads: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct SpurTransformOutput {
    pub matrix: Option<Vec<f64>>,
    pub transformed: Option<Vec<f64>>,
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
    pub demean_lbmgls: bool,
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
        if v < rhobar || c1 > 10_000.0 {
            break;
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

fn col_demean_matrix(src: &[f64], n: usize, q: usize) -> Vec<f64> {
    let mut out = src.to_vec();
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
    out
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

fn select_leading_basis(
    sigdm: &mut [f64],
    n: usize,
    q: usize,
) -> Result<(Vec<f64>, bool), CoreError> {
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
        let tolerance = 64.0 * f64::EPSILON * max_abs;
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
    let mut y_dm = col_demean_matrix(y, n, n_y);
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
    let mut p_mat = vec![0.0_f64; q * n_y];
    linalg::multiply(basis_r, &y_dm, &mut p_mat, q, n, n_y, true);
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
        if !p_norm.is_finite() || p_norm <= 1e-12 * y_dm_norm[j] {
            return Err(CoreError::degenerate("projected SPUR vector has zero norm"));
        }
    }
    Ok(p_mat)
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

fn run_spur_i0_from_basis(
    geom: &mut SpurGeometry,
    basis_r: &[f64],
    tied: bool,
    mr: &[f64],
    om_bm: &[f64],
    y: &[f64],
    n_y: usize,
    emat: &[f64],
    q: usize,
    nrep: usize,
) -> Result<SpurI0Output, CoreError> {
    let n = geom.n;
    let p_mat = prepare_y_projections(y, n, n_y, basis_r, q)?;
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
        &p_mat,
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

fn run_spur_i1_from_basis(
    geom: &mut SpurGeometry,
    basis_r: &[f64],
    tied: bool,
    mr: &[f64],
    om_ho: &[f64],
    y: &[f64],
    n_y: usize,
    emat: &[f64],
    q: usize,
    nrep: usize,
) -> Result<SpurI1Output, CoreError> {
    let n = geom.n;
    let p_mat = prepare_y_projections(y, n, n_y, basis_r, q)?;
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
        &p_mat,
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
    geom: &SpurGeometry,
    q: usize,
) -> Result<(Vec<f64>, bool, Vec<f64>, Vec<f64>), CoreError> {
    let n = geom.n;
    let sigdm_bm = geom.sigma_lbm_dm();
    let mut sigdm_work = sigdm_bm.clone();
    let (basis_r, tied) = select_leading_basis(&mut sigdm_work, n, q)?;
    let mr = col_demean_matrix(&basis_r, n, q);
    let mut sig_r = vec![0.0_f64; n * q];
    linalg::multiply(&sigdm_bm, &basis_r, &mut sig_r, n, n, q, false);
    let mut om_bm = vec![0.0_f64; q * q];
    linalg::multiply(&basis_r, &sig_r, &mut om_bm, q, n, q, true);
    symmetrize(&mut om_bm, q);
    Ok((basis_r, tied, mr, om_bm))
}

fn prepare_resid_basis(
    geom: &mut SpurGeometry,
    x_in: &[f64],
    k_x: usize,
    q: usize,
) -> Result<(Vec<f64>, bool, Vec<f64>, Vec<f64>), CoreError> {
    let n = geom.n;
    let qx = orthonormal_design_basis(x_in, n, k_x)?;
    if q > n - k_x {
        return Err(CoreError::invalid_input(
            "q cannot exceed residual degrees of freedom (n - rank(X))",
        ));
    }
    let c_bm = geom.get_cbar(RHO_BM)?;
    let sigdm_bm = geom.residual_sigma_bm(c_bm, &qx, k_x);
    let mut sigdm_work = sigdm_bm.clone();
    let (basis_r, tied) = select_leading_basis(&mut sigdm_work, n, q)?;
    let mr = residual_project_matrix(&basis_r, &qx, n, k_x, q);
    let mut sig_r = vec![0.0_f64; n * q];
    linalg::multiply(&sigdm_bm, &basis_r, &mut sig_r, n, n, q, false);
    let mut om_bm = vec![0.0_f64; q * q];
    linalg::multiply(&basis_r, &sig_r, &mut om_bm, q, n, q, true);
    symmetrize(&mut om_bm, q);
    Ok((basis_r, tied, mr, om_bm))
}

pub fn run_spur_i0(req: &SpurTestRequest<'_>) -> Result<SpurI0Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    run_with_threads(req.nthreads, || {
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_r, tied, mr, om_bm) = prepare_lbm_basis(&geom, req.q)?;
        run_spur_i0_from_basis(
            &mut geom, &basis_r, tied, &mr, &om_bm, req.y, req.n_y, req.emat, req.q, req.nrep,
        )
    })?
}

pub fn run_spur_i1(req: &SpurTestRequest<'_>) -> Result<SpurI1Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    run_with_threads(req.nthreads, || {
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_r, tied, mr, om_bm) = prepare_lbm_basis(&geom, req.q)?;
        run_spur_i1_from_basis(
            &mut geom, &basis_r, tied, &mr, &om_bm, req.y, req.n_y, req.emat, req.q, req.nrep,
        )
    })?
}

pub fn run_spur_i0_resid(req: &SpurTestRequest<'_>) -> Result<SpurI0Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    let x_in = req
        .x_in
        .ok_or_else(|| CoreError::invalid_input("residual test requires X_in design matrix"))?;
    run_with_threads(req.nthreads, || {
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_r, tied, mr, om_bm) = prepare_resid_basis(&mut geom, x_in, req.k_x, req.q)?;
        run_spur_i0_from_basis(
            &mut geom, &basis_r, tied, &mr, &om_bm, req.y, req.n_y, req.emat, req.q, req.nrep,
        )
    })?
}

pub fn run_spur_i1_resid(req: &SpurTestRequest<'_>) -> Result<SpurI1Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    let x_in = req
        .x_in
        .ok_or_else(|| CoreError::invalid_input("residual test requires X_in design matrix"))?;
    run_with_threads(req.nthreads, || {
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_r, tied, mr, om_bm) = prepare_resid_basis(&mut geom, x_in, req.k_x, req.q)?;
        run_spur_i1_from_basis(
            &mut geom, &basis_r, tied, &mr, &om_bm, req.y, req.n_y, req.emat, req.q, req.nrep,
        )
    })?
}

#[inline(always)]
fn fast_neg_half_pow(x: f64, q: usize) -> f64 {
    if q == 8 {
        let x2 = x * x;
        return (x2 * x2).recip();
    }
    if q == 15 {
        let x2 = x * x;
        let x4 = x2 * x2;
        let x7 = x4 * x2 * x;
        return (x7 * x.sqrt()).recip();
    }
    let mut m = q >> 1;
    let mut acc = 1.0_f64;
    let mut base = x;
    while m > 0 {
        if (m & 1) != 0 {
            acc *= base;
        }
        base *= base;
        m >>= 1;
    }
    if (q & 1) != 0 {
        acc *= x.sqrt();
    }
    acc.recip()
}

fn halflife_tile_size(q: usize, p_pairs: usize, nrep: usize) -> (usize, bool) {
    // Packed pair products dominate scratch for large q. If even one packed
    // draw would exceed the budget, use the factor-based path for every draw.
    let base_per_draw = 2 * q + HL_HA_LEN + 3;
    let budget = MAX_HALFLIFE_TILE_BYTES / size_of::<f64>();
    let packed = base_per_draw + p_pairs <= budget;
    let per_draw = base_per_draw + if packed { p_pairs } else { 0 };
    ((budget / per_draw).clamp(1, TILE_COLS).min(nrep), packed)
}

fn quadratic_form_from_upper(upper: &[f64], draws: &[f64], cols: usize, c: usize) -> f64 {
    let q = draws.len() / cols;
    upper
        .chunks_exact(q)
        .enumerate()
        .map(|(i, row)| {
            (i..q)
                .map(|j| row[j] * draws[j * cols + c])
                .sum::<f64>()
                .powi(2)
        })
        .sum()
}

fn log_sum_exp(terms: &[f64]) -> f64 {
    let max = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let scaled_sum = terms.iter().map(|&term| (term - max).exp()).sum::<f64>();
    max + scaled_sum.ln()
}

fn accepts_halflife_null(exceedances: usize, nrep: usize, level: f64) -> bool {
    // Compare the complement from the integer count so an exact grid boundary
    // is not accepted due to rounding 1 - level.
    (nrep - exceedances) as f64 / (nrep as f64) < level
}

pub fn run_spur_halflife(req: &SpurHalfLifeRequest<'_>) -> Result<SpurHalfLifeOutput, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    if !req.level.is_finite() || req.level <= 0.0 || req.level >= 1.0 {
        return Err(CoreError::invalid_input(
            "confidence level must lie strictly in (0, 1)",
        ));
    }
    if req.y.len() != req.n || req.y.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input(
            "response y must have length n with finite values",
        ));
    }
    let (y_min, y_max) = req
        .y
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(mn, mx), &v| {
            (mn.min(v), mx.max(v))
        });
    if y_max - y_min == 0.0 {
        return Err(CoreError::degenerate(
            "response y is constant; half-life inference requires variation in y",
        ));
    }

    run_with_threads(req.nthreads, || {
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let q = req.q;
        let n = req.n;
        let nrep = req.nrep;

        let c_bm = geom.get_cbar(RHO_BM)?;
        let mut sigdm_bm = geom.kernel_dm(c_bm);
        let (basis_r, tied) = select_leading_basis(&mut sigdm_bm, n, q)?;
        let mr = col_demean_matrix(&basis_r, n, q);
        let p_vec = prepare_y_projections(req.y, n, 1, &basis_r, q)?;

        let mut hl_grid_ho = Vec::with_capacity(HL_HO_LEN);
        for i in 0..HL_HO_SMALL_LEN {
            hl_grid_ho.push(0.001 + (1.0 - 0.001) * (i as f64) / ((HL_HO_SMALL_LEN - 1) as f64));
        }
        for i in 0..HL_HO_LARGE_LEN {
            hl_grid_ho.push(1.01 + (3.0 - 1.01) * (i as f64) / ((HL_HO_LARGE_LEN - 1) as f64));
        }
        hl_grid_ho.push(100.0);

        let hl_grid_ha: Vec<f64> = (0..HL_HA_LEN)
            .map(|i| 0.001 + (1.0 - 0.001) * (i as f64) / ((HL_HA_LEN - 1) as f64))
            .collect();

        let neg_ln_half = -(0.5_f64).ln();
        let c_grid_ho: Vec<f64> = hl_grid_ho.iter().map(|&h| neg_ln_half / h).collect();
        let c_grid_ha: Vec<f64> = hl_grid_ha.iter().map(|&h| neg_ln_half / h).collect();

        let mut c_all = Vec::with_capacity(HL_HO_LEN + HL_HA_LEN);
        c_all.extend_from_slice(&c_grid_ho);
        c_all.extend_from_slice(&c_grid_ha);
        let om_all = geom.projected_omega_grid(&c_all, &mr, q);

        // If every tested covariance has the same shape, the angular likelihood
        // is constant. Strict simulation exceedances then have no meaningful
        // interpretation, even when q > 1 (e.g. an equilateral triangle).
        let reference_shape = &om_all[0];
        let mut informative = false;
        for candidate in &om_all[1..] {
            if covariance_shapes_differ(candidate, reference_shape, q)? {
                informative = true;
                break;
            }
        }
        if !informative {
            return Err(CoreError::degenerate(
                "half-life experiment is uninformative: retained covariance shapes are proportional",
            ));
        }

        let (om_ho_slice, om_ha_slice) = om_all.split_at(HL_HO_LEN);
        let ho_factors: Vec<OmegaPairFactor> = if is_serial() {
            om_ho_slice
                .iter()
                .map(|om| factor_omega(om, q))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            om_ho_slice
                .par_iter()
                .map(|om| factor_omega(om, q))
                .collect::<Result<Vec<_>, _>>()?
        };
        let ha_factors: Vec<OmegaPairFactor> = if is_serial() {
            om_ha_slice
                .iter()
                .map(|om| factor_omega(om, q))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            om_ha_slice
                .par_iter()
                .map(|om| factor_omega(om, q))
                .collect::<Result<Vec<_>, _>>()?
        };

        let q_f = q as f64;
        let ha_x_log: Vec<f64> = ha_factors
            .iter()
            .map(|ha| {
                let mut xc = vec![0.0_f64; q];
                linalg::multiply(&ha.ch_omi.upper, &p_vec, &mut xc, q, q, 1, false);
                let qnorm: f64 = xc.iter().map(|&v| v * v).sum();
                if !qnorm.is_finite() || qnorm <= 0.0 {
                    return Err(CoreError::degenerate(
                        "projected SPUR vector has zero quadratic form",
                    ));
                }
                Ok(ha.ch_omi.log_sqrt_det - 0.5 * q_f * qnorm.ln())
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut emat_sq = vec![0.0_f64; nrep];
        for r in 0..q {
            let row = &req.emat[r * nrep..(r + 1) * nrep];
            for c in 0..nrep {
                emat_sq[c] += row[c] * row[c];
            }
        }
        if emat_sq.iter().any(|&v| !v.is_finite() || v <= 0.0) {
            return Err(CoreError::invalid_input(
                "half-life simulation draws must have finite positive squared norms",
            ));
        }

        let p_pairs = q * (q + 1) / 2;
        let (max_tile, packed) = halflife_tile_size(q, p_pairs, nrep);
        let mut w_ha = Vec::new();
        let mut cancellation_scale = Vec::new();
        if packed {
            w_ha.resize(HL_HA_LEN * p_pairs, 0.0);
            for (j, ha) in ha_factors.iter().enumerate() {
                let mut kp = 0usize;
                for b in 0..q {
                    for a in 0..=b {
                        let val = ha.om_inv[a * q + b];
                        w_ha[j * p_pairs + kp] = if a == b { val } else { 2.0 * val };
                        kp += 1;
                    }
                }
                // Sum of absolute terms is bounded by ||A||_F * ||z||².
                // Flag positive forms too when cancellation can dominate.
                let frob = ha.om_inv.iter().map(|v| v * v).sum::<f64>().sqrt();
                cancellation_scale.push(8.0 * f64::EPSILON * (p_pairs as f64) * frob);
            }
        }

        let eval_null_hl = |ho: &OmegaPairFactor| -> Result<usize, CoreError> {
            let mut xc_obs = vec![0.0_f64; q];
            linalg::multiply(&ho.ch_omi.upper, &p_vec, &mut xc_obs, q, q, 1, false);
            let qnorm_obs: f64 = xc_obs.iter().map(|&v| v * v).sum();
            if !qnorm_obs.is_finite() || qnorm_obs <= 0.0 {
                return Err(CoreError::degenerate(
                    "projected SPUR vector has zero quadratic form",
                ));
            }
            let log_ho_obs = ho.ch_omi.log_sqrt_det - 0.5 * q_f * qnorm_obs.ln();
            let log_obs_sum = log_sum_exp(&ha_x_log) - log_ho_obs;
            let obs_sum = log_obs_sum.exp();

            let det_ratios: Vec<f64> = ha_factors
                .iter()
                .map(|ha| (ha.ch_omi.log_sqrt_det - ho.ch_omi.log_sqrt_det).exp())
                .collect();
            let use_fast_pow = packed
                && q <= 60
                && obs_sum.is_finite()
                && obs_sum >= f64::MIN_POSITIVE
                && det_ratios
                    .iter()
                    .all(|&v| v.is_finite() && v > 1e-150 && v < 1e150);

            let mut exceed = 0usize;
            let mut tile_e = vec![0.0_f64; q * max_tile];
            let mut e_tile = vec![0.0_f64; q * max_tile];
            let mut u_tile = vec![0.0_f64; if packed { p_pairs * max_tile } else { 0 }];
            let mut q_tile = vec![0.0_f64; if packed { HL_HA_LEN * max_tile } else { 0 }];
            let mut inv_q_ho = vec![0.0_f64; if use_fast_pow { max_tile } else { 0 }];
            let mut sum_ha_over_ho = vec![0.0_f64; if use_fast_pow { max_tile } else { 0 }];
            let mut col_sq = vec![0.0_f64; if packed { max_tile } else { 0 }];

            for start in (0..nrep).step_by(max_tile) {
                let cols = (nrep - start).min(max_tile);
                for r in 0..q {
                    tile_e[r * cols..(r + 1) * cols]
                        .copy_from_slice(&req.emat[r * nrep + start..r * nrep + start + cols]);
                }

                // e_tile = L_ho * tile_e has covariance om_ho and exact null quadratic form ||tile_e||^2.
                linalg::multiply(
                    &ho.ch_om.lower,
                    &tile_e[..q * cols],
                    &mut e_tile[..q * cols],
                    q,
                    q,
                    cols,
                    false,
                );

                if packed {
                    col_sq[..cols].fill(0.0);
                    let mut kp = 0usize;
                    for b in 0..q {
                        let row_b = &e_tile[b * cols..(b + 1) * cols];
                        for a in 0..=b {
                            let row_a = &e_tile[a * cols..(a + 1) * cols];
                            let dst = &mut u_tile[kp * cols..(kp + 1) * cols];
                            for c in 0..cols {
                                dst[c] = row_a[c] * row_b[c];
                            }
                            if a == b {
                                for (sq, &term) in col_sq[..cols].iter_mut().zip(dst.iter()) {
                                    *sq += term;
                                }
                            }
                            kp += 1;
                        }
                    }
                    linalg::multiply(
                        &w_ha,
                        &u_tile[..p_pairs * cols],
                        &mut q_tile[..HL_HA_LEN * cols],
                        HL_HA_LEN,
                        p_pairs,
                        cols,
                        false,
                    );
                    for j in 0..HL_HA_LEN {
                        let row = &mut q_tile[j * cols..(j + 1) * cols];
                        let thresh = cancellation_scale[j];
                        // The normal case needs only the reduction. A second
                        // pass repairs suspect forms without branching on each
                        // simulation draw in the common path.
                        let needs_repair = row
                            .iter()
                            .zip(&col_sq[..cols])
                            .fold(false, |bad, (&form, &norm)| {
                                bad | !form.is_finite() | (form <= thresh * norm)
                            });
                        if !needs_repair {
                            continue;
                        }
                        for c in 0..cols {
                            let mut form = row[c];
                            if !form.is_finite() || form <= thresh * col_sq[c] {
                                form = quadratic_form_from_upper(
                                    &ha_factors[j].ch_omi.upper,
                                    &e_tile[..q * cols],
                                    cols,
                                    c,
                                );
                                if !form.is_finite() || form <= 0.0 {
                                    return Err(CoreError::degenerate(
                                        "half-life simulation draw has an invalid covariance quadratic form",
                                    ));
                                }
                                row[c] = form;
                            }
                        }
                    }
                }

                let log_draw_sum = |c: usize| -> Result<f64, CoreError> {
                    let log_ho = ho.ch_omi.log_sqrt_det - 0.5 * q_f * emat_sq[start + c].ln();
                    let mut log_ha = [0.0; HL_HA_LEN];
                    if packed {
                        for (j, ha) in ha_factors.iter().enumerate() {
                            let q_ha = q_tile[j * cols + c];
                            log_ha[j] = ha.ch_omi.log_sqrt_det - 0.5 * q_f * q_ha.ln();
                        }
                    } else {
                        for (j, ha) in ha_factors.iter().enumerate() {
                            let q_ha = quadratic_form_from_upper(
                                &ha.ch_omi.upper,
                                &e_tile[..q * cols],
                                cols,
                                c,
                            );
                            if !q_ha.is_finite() || q_ha <= 0.0 {
                                return Err(CoreError::degenerate(
                                    "half-life simulation draw has an invalid covariance quadratic form",
                                ));
                            }
                            log_ha[j] = ha.ch_omi.log_sqrt_det - 0.5 * q_f * q_ha.ln();
                        }
                    }
                    Ok(log_sum_exp(&log_ha) - log_ho)
                };

                if use_fast_pow {
                    for c in 0..cols {
                        inv_q_ho[c] = emat_sq[start + c].recip();
                    }
                    sum_ha_over_ho[..cols].fill(0.0);
                    for j in 0..HL_HA_LEN {
                        let coeff = det_ratios[j];
                        let row = &q_tile[j * cols..(j + 1) * cols];
                        for c in 0..cols {
                            // Each form and original draw norm was validated.
                            // A non-normal power marks this draw for log-space
                            // recomputation, even if other terms are normal.
                            let ratio = row[c] * inv_q_ho[c];
                            let pow = fast_neg_half_pow(ratio, q);
                            let term = if pow >= f64::MIN_POSITIVE {
                                coeff * pow
                            } else {
                                f64::NAN
                            };
                            sum_ha_over_ho[c] += term;
                        }
                    }
                    for c in 0..cols {
                        let exceeds = if sum_ha_over_ho[c].is_finite()
                            && sum_ha_over_ho[c] >= f64::MIN_POSITIVE
                            && inv_q_ho[c].is_finite()
                        {
                            sum_ha_over_ho[c] > obs_sum
                        } else {
                            log_draw_sum(c)? > log_obs_sum
                        };
                        if exceeds {
                            exceed += 1;
                        }
                    }
                } else {
                    for c in 0..cols {
                        if log_draw_sum(c)? > log_obs_sum {
                            exceed += 1;
                        }
                    }
                }
            }
            Ok(exceed)
        };

        let exceedances = if is_serial() {
            ho_factors
                .iter()
                .map(eval_null_hl)
                .collect::<Result<Vec<_>, _>>()?
        } else {
            ho_factors
                .par_iter()
                .map(eval_null_hl)
                .collect::<Result<Vec<_>, _>>()?
        };

        let mut ci_lower = f64::INFINITY;
        let mut ci_upper = f64::NEG_INFINITY;
        let mut any_accepted = false;
        for (&count, &hl) in exceedances.iter().zip(&hl_grid_ho) {
            if accepts_halflife_null(count, nrep, req.level) {
                any_accepted = true;
                ci_lower = ci_lower.min(hl);
                ci_upper = ci_upper.max(hl);
            }
        }
        if !any_accepted {
            ci_lower = f64::NAN;
            ci_upper = f64::NAN;
        }
        let pv_mat = exceedances
            .into_iter()
            .map(|count| count as f64 / nrep as f64)
            .collect();

        Ok(SpurHalfLifeOutput {
            ci_lower,
            ci_upper,
            max_dist: geom.max_dist,
            pv_mat,
            hl_grid_ho,
            tied,
        })
    })?
}

impl SpurGeometry {
    fn kernel_dm(&self, c: f64) -> Vec<f64> {
        let mut kern = self.exp_kernel(c);
        double_demean(&mut kern, self.n);
        kern
    }
}

fn transform_by_neighbors(
    req: &SpurTransformRequest<'_>,
    select: impl Fn(usize, &mut Vec<usize>) -> Result<(), CoreError> + Sync,
) -> Result<SpurTransformOutput, CoreError> {
    let n = req.n;
    let (width, y) = if req.return_matrix {
        (n, None)
    } else {
        (req.y.map_or(0, |_| req.n_y), req.y)
    };
    let mut values = vec![0.0; checked_len(&[n, width])?];
    let write_row = |neighbors: &mut Vec<usize>, (i, row): (usize, &mut [f64])| {
        neighbors.clear();
        select(i, neighbors)?;
        if neighbors.is_empty() || row.is_empty() {
            return Ok(());
        }
        let weight = (neighbors.len() as f64).recip();
        if let Some(y) = y {
            row.copy_from_slice(&y[i * width..(i + 1) * width]);
            for &j in neighbors.iter() {
                for (out, &value) in row.iter_mut().zip(&y[j * width..(j + 1) * width]) {
                    *out -= weight * value;
                }
            }
        } else {
            row[i] = 1.0;
            for &j in neighbors.iter() {
                row[j] = -weight;
            }
        }
        Ok(())
    };
    if width == 0 {
        let mut neighbors = Vec::new();
        for i in 0..n {
            write_row(&mut neighbors, (i, &mut []))?;
        }
    } else if is_serial() || n < 64 {
        let mut neighbors = Vec::new();
        for row in values.chunks_exact_mut(width).enumerate() {
            write_row(&mut neighbors, row)?;
        }
    } else {
        values
            .par_chunks_mut(width)
            .enumerate()
            .try_for_each_init(Vec::new, write_row)?;
    }

    if req.return_matrix {
        let transformed = req.y.map(|y| {
            let mut out = vec![0.0; y.len()];
            linalg::multiply(&values, y, &mut out, n, n, req.n_y, false);
            out
        });
        Ok(SpurTransformOutput {
            matrix: Some(values),
            transformed,
        })
    } else {
        Ok(SpurTransformOutput {
            matrix: None,
            transformed: y.map(|_| values),
        })
    }
}

fn neighbor_transform(
    req: &SpurTransformRequest<'_>,
    coords: &[f64],
) -> Result<SpurTransformOutput, CoreError> {
    let n = req.n;
    let distances = Distances::new(coords, n, req.d, req.metric)?;
    // Cache pairwise distances when the caller already needs a dense matrix.
    // A transformed-only request keeps only the current row's neighbor indices.
    let dense = req
        .return_matrix
        .then(|| distance_matrix(coords, n, req.d, req.metric))
        .transpose()?;
    let distance = |i: usize, j: usize| {
        let value = match &dense {
            Some(matrix) => matrix[i * n + j],
            None => distances.pair(i.max(j), i.min(j)),
        };
        if !value.is_finite() {
            return Err(CoreError::invalid_input(
                "pairwise distance exceeds the floating-point range",
            ));
        }
        Ok(value)
    };
    match req.kind {
        SpurTransformKind::Nn {
            exclude_zero_distances,
        } => {
            // Preserve the previous normalization before equality checks. For
            // all-colocated data, every off-diagonal distance remains zero.
            let max_in_row =
                |i: usize| (0..i).try_fold(0.0_f64, |max, j| Ok(max.max(distance(i, j)?)));
            let diameter = if let Some(matrix) = &dense {
                matrix.iter().copied().fold(0.0, f64::max)
            } else if is_serial() || n < 64 {
                (0..n).try_fold(0.0_f64, |max, i| Ok(max.max(max_in_row(i)?)))?
            } else {
                (0..n)
                    .into_par_iter()
                    .map(max_in_row)
                    .try_reduce(|| 0.0_f64, |a, b| Ok(a.max(b)))?
            };
            if exclude_zero_distances && diameter == 0.0 {
                return Err(CoreError::invalid_input(
                    "all coordinates are identical; normalization requires distinct locations",
                ));
            }
            let scale = if diameter > 0.0 { diameter } else { 1.0 };
            transform_by_neighbors(req, |i, neighbors| {
                let mut nearest = f64::INFINITY;
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let mut d = distance(i, j)? / scale;
                    if exclude_zero_distances && d == 0.0 {
                        d = 1e10;
                    }
                    if d < nearest {
                        nearest = d;
                        neighbors.clear();
                    }
                    if d == nearest {
                        neighbors.push(j);
                    }
                }
                Ok(())
            })
        }
        SpurTransformKind::Iso { radius, strict_lt } => {
            transform_by_neighbors(req, |i, neighbors| {
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let d = distance(i, j)?;
                    if d > 0.0 && if strict_lt { d < radius } else { d <= radius } {
                        neighbors.push(j);
                    }
                }
                Ok(())
            })
        }
        _ => unreachable!("neighbor transform requires NN or ISO"),
    }
}

pub fn run_spur_transform(
    req: &SpurTransformRequest<'_>,
) -> Result<SpurTransformOutput, CoreError> {
    let n = req.n;
    if n < 2 {
        return Err(CoreError::invalid_input(
            "spatial transform requires at least 2 observations",
        ));
    }
    if let Some(y) = req.y {
        if req.n_y == 0 || y.len() != checked_len(&[n, req.n_y])? {
            return Err(CoreError::invalid_input(
                "invalid Y dimensions for transform",
            ));
        }
        if y.iter().any(|v| !v.is_finite()) {
            return Err(CoreError::invalid_input("Y values must be finite"));
        }
    }

    run_with_threads(req.nthreads, || match req.kind {
        SpurTransformKind::Cluster => {
            let cid = req.cluster_id.ok_or_else(|| {
                CoreError::invalid_input("cluster transformation requires cluster_id")
            })?;
            if cid.len() != n {
                return Err(CoreError::invalid_input("cluster_id length mismatch"));
            }
            let max_c = cid.iter().copied().max().unwrap_or(0) + 1;
            let mut counts = vec![0usize; max_c];
            for &g in cid {
                counts[g] += 1;
            }
            let transformed = req.y.map(|y| {
                let ny = req.n_y;
                let mut sums = vec![0.0_f64; max_c * ny];
                for (i, &g) in cid.iter().enumerate() {
                    for col in 0..ny {
                        sums[g * ny + col] += y[i * ny + col];
                    }
                }
                let mut out = vec![0.0_f64; n * ny];
                for (i, &g) in cid.iter().enumerate() {
                    let inv_cnt = (counts[g] as f64).recip();
                    for col in 0..ny {
                        out[i * ny + col] = y[i * ny + col] - sums[g * ny + col] * inv_cnt;
                    }
                }
                out
            });
            let matrix = if req.return_matrix {
                let mut mat = vec![0.0_f64; n * n];
                for i in 0..n {
                    let gi = cid[i];
                    let w = (counts[gi] as f64).recip();
                    for j in 0..n {
                        if cid[j] == gi {
                            mat[i * n + j] = -w;
                        }
                    }
                    mat[i * n + i] += 1.0;
                }
                Some(mat)
            } else {
                None
            };
            Ok(SpurTransformOutput {
                matrix,
                transformed,
            })
        }
        SpurTransformKind::Nn { .. } => {
            let coords = req
                .coords
                .ok_or_else(|| CoreError::invalid_input("nn transform requires coordinates"))?;
            neighbor_transform(req, coords)
        }
        SpurTransformKind::Iso { radius, .. } => {
            if !radius.is_finite() || radius <= 0.0 {
                return Err(CoreError::invalid_input(
                    "iso transform radius must be positive and finite",
                ));
            }
            let coords = req
                .coords
                .ok_or_else(|| CoreError::invalid_input("iso transform requires coordinates"))?;
            neighbor_transform(req, coords)
        }
        SpurTransformKind::LbmGls => {
            let coords = req
                .coords
                .ok_or_else(|| CoreError::invalid_input("lbmgls transform requires coordinates"))?;
            let geom = SpurGeometry::from_coords(coords, n, req.d, req.metric)?;
            apply_lbmgls_with_geom(&geom, req.y, req.n_y, req.return_matrix, req.demean_lbmgls)
        }
    })?
}

fn apply_lbmgls_with_geom(
    geom: &SpurGeometry,
    y: Option<&[f64]>,
    n_y: usize,
    return_matrix: bool,
    demean_lbmgls: bool,
) -> Result<SpurTransformOutput, CoreError> {
    let n = geom.n;
    let mut sigdm = geom.sigma_lbm_dm();
    let pairs = basis::dense_eigenpairs(&mut sigdm, n, n)?;
    let keep: Vec<usize> = pairs
        .values
        .iter()
        .enumerate()
        .filter(|(_, &v)| v > LBM_GLS_EIG_TOL)
        .map(|(idx, _)| idx)
        .collect();
    let rank = keep.len();
    if rank == 0 {
        return Err(CoreError::degenerate(
            "LBM-GLS covariance has no eigenvalues above 1e-10",
        ));
    }

    // Build V_keep (n x rank) and V_scaled = V_keep * diag(eval^{-1/2}) (n x rank).
    let mut v_keep = vec![0.0_f64; n * rank];
    let mut v_scaled = vec![0.0_f64; n * rank];
    for (col_idx, &eig_idx) in keep.iter().enumerate() {
        let inv_sqrt = pairs.values[eig_idx].sqrt().recip();
        let vec_slice = &pairs.vectors[eig_idx * n..(eig_idx + 1) * n];
        for i in 0..n {
            let v = vec_slice[i];
            v_keep[i * rank + col_idx] = v;
            v_scaled[i * rank + col_idx] = v * inv_sqrt;
        }
    }

    let matrix = if return_matrix {
        let mut v_keep_t = vec![0.0_f64; rank * n];
        for i in 0..n {
            for r in 0..rank {
                v_keep_t[r * n + i] = v_keep[i * rank + r];
            }
        }
        let mut mat = vec![0.0_f64; n * n];
        linalg::multiply(&v_scaled, &v_keep_t, &mut mat, n, rank, n, false);
        symmetrize(&mut mat, n);
        Some(mat)
    } else {
        None
    };

    let transformed = y.map(|y_mat| {
        let mut vt_y = vec![0.0_f64; rank * n_y];
        linalg::multiply(&v_keep, y_mat, &mut vt_y, rank, n, n_y, true);
        let mut out = vec![0.0_f64; n * n_y];
        linalg::multiply(&v_scaled, &vt_y, &mut out, n, rank, n_y, false);
        if demean_lbmgls {
            col_demean_matrix(&out, n, n_y)
        } else {
            out
        }
    });

    Ok(SpurTransformOutput {
        matrix,
        transformed,
    })
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
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_lbm, tied_lbm, mr_lbm, om_lbm) = prepare_lbm_basis(&geom, req.q)?;
        let i0 = run_spur_i0_from_basis(
            &mut geom, &basis_lbm, tied_lbm, &mr_lbm, &om_lbm, req.y, 1, req.emat, req.q, req.nrep,
        )?;
        let i1 = run_spur_i1_from_basis(
            &mut geom, &basis_lbm, tied_lbm, &mr_lbm, &om_lbm, req.y, 1, req.emat, req.q, req.nrep,
        )?;

        let (basis_res, tied_res, mr_res, om_res) =
            prepare_resid_basis(&mut geom, req.x_in, req.k_x, req.q)?;
        let i0resid = run_spur_i0_from_basis(
            &mut geom, &basis_res, tied_res, &mr_res, &om_res, req.y, 1, req.emat, req.q, req.nrep,
        )?;
        let i1resid = run_spur_i1_from_basis(
            &mut geom, &basis_res, tied_res, &mr_res, &om_res, req.y, 1, req.emat, req.q, req.nrep,
        )?;

        let tr = apply_lbmgls_with_geom(
            &geom,
            Some(req.transform_vars),
            req.n_transform_vars,
            false,
            req.demean_lbmgls,
        )?;
        Ok(SpurPipelineOutput {
            i0,
            i1,
            i0resid,
            i1resid,
            transformed_vars: tr.transformed.expect("transformed vars requested"),
        })
    })?
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ErrorKind;

    fn load_bench_200_spur() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let text = include_str!("../tests/fixtures/bench_200.csv");
        let mut y = Vec::new();
        let mut x_in = Vec::new();
        let mut coords = Vec::new();
        for line in text.lines().skip(1) {
            if line.trim().is_empty() {
                continue;
            }
            let cols: Vec<f64> = line
                .split(',')
                .map(|s| s.trim().parse::<f64>().unwrap())
                .collect();
            y.push(cols[1]);
            x_in.extend_from_slice(&[1.0, cols[2], cols[3]]);
            // (s_1, s_2)
            coords.extend_from_slice(&[cols[6], cols[7]]);
        }
        (y, x_in, coords)
    }

    #[test]
    fn streamed_neighbor_transforms_match_matrix_paths() {
        let coords = [0.0, 0.0, 1.0, 2.0, 6.0];
        let y = [1.0, 3.0, 5.0, 7.0, 9.0];
        let nn = SpurTransformKind::Nn {
            exclude_zero_distances: false,
        };
        let nn_excluded = SpurTransformKind::Nn {
            exclude_zero_distances: true,
        };
        let iso = SpurTransformKind::Iso {
            radius: 1.0,
            strict_lt: false,
        };
        let iso_strict = SpurTransformKind::Iso {
            radius: 1.0,
            strict_lt: true,
        };
        let base = SpurTransformRequest {
            coords: Some(&coords),
            n: coords.len(),
            d: 1,
            metric: Metric::Euclidean,
            kind: nn,
            cluster_id: None,
            y: Some(&y),
            n_y: 1,
            return_matrix: false,
            demean_lbmgls: false,
            nthreads: Some(1),
        };
        let compare = |base: SpurTransformRequest<'_>| {
            let dense = run_spur_transform(&SpurTransformRequest {
                return_matrix: true,
                ..base
            })
            .unwrap();
            let streamed = run_spur_transform(&base).unwrap();
            assert!(dense.matrix.is_some());
            assert!(streamed.matrix.is_none());
            let values = streamed.transformed.unwrap();
            for (&actual, expected) in values.iter().zip(dense.transformed.unwrap()) {
                assert!((actual - expected).abs() < 1e-12);
            }
            values
        };
        for (kind, expected) in [
            (nn, [-2.0, 2.0, 4.0 / 3.0, 2.0, 2.0]),
            (nn_excluded, [-4.0, -2.0, 4.0 / 3.0, 2.0, 2.0]),
            (iso, [-4.0, -2.0, 4.0 / 3.0, 2.0, 0.0]),
            (iso_strict, [0.0; 5]),
        ] {
            let values = compare(SpurTransformRequest { kind, ..base });
            for (actual, expected) in values.into_iter().zip(expected) {
                assert!((actual - expected).abs() < 1e-12);
            }
        }
        assert_eq!(
            compare(SpurTransformRequest {
                coords: Some(&[5.0; 3]),
                n: 3,
                y: Some(&[1.0, 3.0, 5.0]),
                ..base
            }),
            vec![-3.0, 0.0, 3.0]
        );
        for return_matrix in [false, true] {
            let err = run_spur_transform(&SpurTransformRequest {
                coords: Some(&[5.0; 3]),
                n: 3,
                kind: nn_excluded,
                y: Some(&[1.0, 3.0, 5.0]),
                return_matrix,
                ..base
            })
            .unwrap_err();
            assert_eq!(err.kind(), ErrorKind::InvalidInput);
        }
        let coords_many: Vec<f64> = (0..80).flat_map(|i| [i as f64, (i % 5) as f64]).collect();
        let y_many: Vec<f64> = (0..160).map(|i| (i as f64).sin()).collect();
        for kind in [nn, iso] {
            compare(SpurTransformRequest {
                coords: Some(&coords_many),
                n: 80,
                d: 2,
                metric: Metric::GreatCircleFractionLatLon,
                kind,
                y: Some(&y_many),
                n_y: 2,
                nthreads: Some(4),
                ..base
            });
        }
    }

    #[test]
    fn spur_basis_orientation_uses_the_first_nearly_maximal_entry() {
        for (offset, first_is_positive) in [(8.0, true), (1024.0, false)] {
            let a = -1.0;
            let b = 1.0 + offset * f64::EPSILON;
            let mut covariance = [a * a, a * b, a * b, b * b];
            let (basis, tied) = select_leading_basis(&mut covariance, 2, 1).unwrap();
            assert!(!tied);
            assert_eq!(basis[0].is_sign_positive(), first_is_positive);
            assert_eq!(basis[1].is_sign_negative(), first_is_positive);
        }
    }

    #[test]
    fn halflife_worker_scratch_and_likelihood_terms_are_stable() {
        for q in [4, 64, 500, 1200] {
            let pairs = q * (q + 1) / 2;
            let (cols, packed) = halflife_tile_size(q, pairs, 10_000);
            assert!(cols <= TILE_COLS);
            assert_eq!(packed, q <= 500);
            let tile_doubles = cols * (2 * q + HL_HA_LEN + 3 + if packed { pairs } else { 0 });
            assert!(tile_doubles * size_of::<f64>() <= MAX_HALFLIFE_TILE_BYTES);
        }
        let upper = [1.0, -1.0, 0.0, 1e-7];
        let draws = [1.0, 2.0, 1.0, 3.0];
        let form = quadratic_form_from_upper(&upper, &draws, 2, 0);
        assert!((form / 1e-14 - 1.0).abs() < 4.0 * f64::EPSILON);
        let log_total = log_sum_exp(&[1000.0, 1000.0]);
        assert!((log_total - 1000.0 - 2.0_f64.ln()).abs() < 1e-13);
    }

    #[test]
    fn halflife_confidence_set_excludes_exact_monte_carlo_boundaries() {
        for (level, boundary) in [(0.90, 10), (0.95, 5)] {
            assert!(!accepts_halflife_null(boundary - 1, 100, level));
            assert!(!accepts_halflife_null(boundary, 100, level));
            assert!(accepts_halflife_null(boundary + 1, 100, level));
        }
    }

    fn make_deterministic_emat(q: usize, nrep: usize) -> Vec<f64> {
        // Box-Muller from deterministic LCG for self-contained Rust tests.
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut next_u = || -> f64 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let bits = (state >> 11) as f64 * (1.0 / (1_u64 << 53) as f64);
            bits.clamp(1e-12, 1.0 - 1e-12)
        };
        let mut out = vec![0.0_f64; q * nrep];
        for chunk in out.chunks_mut(2) {
            let u1 = next_u();
            let u2 = next_u();
            let r = (-2.0 * u1.ln()).sqrt();
            let th = 2.0 * std::f64::consts::PI * u2;
            chunk[0] = r * th.cos();
            if chunk.len() > 1 {
                chunk[1] = r * th.sin();
            }
        }
        out
    }

    #[test]
    fn spur_cbar_and_diagnostics_are_deterministic_and_thread_invariant() {
        let (y, x_in, coords) = load_bench_200_spur();
        let n = y.len();
        let q = 8usize;
        let nrep = 1000usize;
        let emat = make_deterministic_emat(q, nrep);

        let mut geom = SpurGeometry::from_coords(&coords, n, 2, Metric::Euclidean).unwrap();
        let c099 = geom.get_cbar(0.999).unwrap();
        let c099_cached = geom.get_cbar(0.999).unwrap();
        assert_eq!(c099.to_bits(), c099_cached.to_bits());
        assert!((c099 - 0.002_662_372_394_202_289_4).abs() < 1e-12);
        let mean_rho = deterministic_mean_exp(&geom.lower_tri, c099);
        assert!((mean_rho - 0.999).abs() < 1e-3, "mean_rho={mean_rho}");

        // Pipeline vs individual tests match to machine precision and are thread invariant.
        let pipe_1 = run_spur_pipeline(&SpurPipelineRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            x_in: &x_in,
            k_x: 3,
            emat: &emat,
            q,
            nrep,
            transform_vars: &x_in,
            n_transform_vars: 3,
            demean_lbmgls: false,
            nthreads: Some(1),
        })
        .unwrap();
        let pipe_4 = run_spur_pipeline(&SpurPipelineRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            x_in: &x_in,
            k_x: 3,
            emat: &emat,
            q,
            nrep,
            transform_vars: &x_in,
            n_transform_vars: 3,
            demean_lbmgls: false,
            nthreads: Some(4),
        })
        .unwrap();
        let i0_ind = run_spur_i0(&SpurTestRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            n_y: 1,
            x_in: None,
            k_x: 0,
            emat: &emat,
            q,
            nrep,
            nthreads: Some(2),
        })
        .unwrap();
        let i0r_ind = run_spur_i0_resid(&SpurTestRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            n_y: 1,
            x_in: Some(&x_in),
            k_x: 3,
            emat: &emat,
            q,
            nrep,
            nthreads: Some(2),
        })
        .unwrap();

        // Exact parity against canonical-sign spur-python with frozen deterministic emat:
        assert!((pipe_1.i0.ha_param - 1.90625).abs() < 1e-12);
        assert!((pipe_1.i0.lr[0] - 14.931_257_343_240_508).abs() < 1e-11);
        assert!((pipe_1.i0.pvalue[0] - 0.0).abs() < 1e-12);
        assert!((pipe_1.i1.ha_param - 27.469_699_236_099_583).abs() < 1e-11);
        assert!((pipe_1.i1.lr[0] - 1.622_386_080_187_783_6).abs() < 1e-11);
        assert!((pipe_1.i1.pvalue[0] - 0.807).abs() < 1e-12);
        assert!((pipe_1.i1.cvalue[0] - 7.928_843_817_974_141).abs() < 1e-10);
        assert!((pipe_1.i1.cvalue[1] - 6.426_857_918_629_026).abs() < 1e-10);
        assert!((pipe_1.i1.cvalue[2] - 5.480_387_716_657_392).abs() < 1e-10);
        assert!((pipe_1.i0resid.ha_param - 600.375).abs() < 1e-12);
        assert!((pipe_1.i0resid.lr[0] - 8.492_145_848_011_605).abs() < 1e-9);
        assert!((pipe_1.i0resid.pvalue[0] - 0.121).abs() < 1e-12);
        assert!((pipe_1.i1resid.ha_param - 875.033_754_868_347_6).abs() < 1e-10);
        assert!((pipe_1.i1resid.lr[0] - 223.478_786_080_475_8).abs() < 1e-8);
        assert!((pipe_1.i1resid.pvalue[0] - 0.24).abs() < 1e-12);
        assert!((pipe_1.i1resid.cvalue[0] - 456.320_954_571_764_9).abs() < 1e-9);
        assert!((pipe_1.i1resid.cvalue[1] - 359.524_332_941_563_05).abs() < 1e-9);
        assert!((pipe_1.i1resid.cvalue[2] - 307.072_910_917_401_26).abs() < 1e-9);

        assert!((pipe_1.i0.lr[0] - i0_ind.lr[0]).abs() < 1e-12);
        assert!((pipe_1.i0.pvalue[0] - i0_ind.pvalue[0]).abs() < 1e-12);
        assert!((pipe_1.i0.ha_param - i0_ind.ha_param).abs() < 1e-12);
        assert!((pipe_1.i0resid.lr[0] - i0r_ind.lr[0]).abs() < 1e-12);
        assert!((pipe_1.i0resid.pvalue[0] - i0r_ind.pvalue[0]).abs() < 1e-12);

        // Check I(0) grid maxima invariant
        assert_eq!(pipe_1.i0.rho_grid.len(), 30);
        let max_p = pipe_1.i0.pvalue_mat.iter().copied().fold(0.0_f64, f64::max);
        assert!((pipe_1.i0.pvalue[0] - max_p).abs() < 1e-15);

        // Thread count invariance
        assert!((pipe_1.i0.lr[0] - pipe_4.i0.lr[0]).abs() < 1e-13);
        assert!((pipe_1.i1.lr[0] - pipe_4.i1.lr[0]).abs() < 1e-13);
        assert!((pipe_1.i0resid.lr[0] - pipe_4.i0resid.lr[0]).abs() < 1e-13);
        assert!((pipe_1.i1resid.lr[0] - pipe_4.i1resid.lr[0]).abs() < 1e-13);
        for (&a, &b) in pipe_1.transformed_vars.iter().zip(&pipe_4.transformed_vars) {
            assert!((a - b).abs() < 1e-13);
        }
    }

    #[test]
    fn spur_halflife_and_lbmgls_invariance_and_colocation() {
        let (y, _, coords) = load_bench_200_spur();
        let n = y.len();
        let q = 8usize;
        let nrep = 300usize;
        let emat = make_deterministic_emat(q, nrep);

        let hl_1 = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            emat: &emat,
            q,
            nrep,
            level: 0.95,
            nthreads: Some(1),
        })
        .unwrap();
        let hl_2 = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            emat: &emat,
            q,
            nrep,
            level: 0.95,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(hl_1.pv_mat.len(), 131);
        assert_eq!(hl_1.hl_grid_ho.len(), 131);
        assert!(
            (hl_1.ci_lower - 0.152_363_636_363_636_37).abs() < 1e-12,
            "got ci_lower={}, ci_upper={}",
            hl_1.ci_lower,
            hl_1.ci_upper
        );
        assert!((hl_1.ci_upper - 100.0).abs() < 1e-12);
        for (&a, &b) in hl_1.pv_mat.iter().zip(&hl_2.pv_mat) {
            assert!((a - b).abs() < 1e-14);
        }

        // Affine invariance: y -> -2.5 * y + 17 preserves angular likelihood.
        let y_aff: Vec<f64> = y.iter().map(|&v| -2.5 * v + 17.0).collect();
        let hl_aff = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y_aff,
            emat: &emat,
            q,
            nrep,
            level: 0.95,
            nthreads: Some(2),
        })
        .unwrap();
        for (&a, &b) in hl_1.pv_mat.iter().zip(&hl_aff.pv_mat) {
            assert!((a - b).abs() < 1e-12, "p_orig={a}, p_aff={b}");
        }

        // Co-located LBM-GLS rank < n - 1 and zero constant vector action:
        let coords_col = vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let mut y_col = vec![1.0, 2.0, -1.0, 0.5, 3.0];
        let ones = vec![1.0_f64; 5];
        y_col.extend(ones);
        // Pass 5 x 2 row-major matrix [y_col, ones]:
        let mut y_mat = vec![0.0_f64; 10];
        for i in 0..5 {
            y_mat[i * 2] = [1.0, 2.0, -1.0, 0.5, 3.0][i];
            y_mat[i * 2 + 1] = 1.0;
        }
        let lbm = run_spur_transform(&SpurTransformRequest {
            coords: Some(&coords_col),
            n: 5,
            d: 2,
            metric: Metric::Euclidean,
            kind: SpurTransformKind::LbmGls,
            cluster_id: None,
            y: Some(&y_mat),
            n_y: 2,
            return_matrix: true,
            demean_lbmgls: false,
            nthreads: Some(1),
        })
        .unwrap();
        let h_mat = lbm.matrix.unwrap();
        // Verify H is symmetric and removes constants, without adding 11'/n.
        for i in 0..5 {
            assert!((lbm.transformed.as_ref().unwrap()[i * 2 + 1]).abs() < 1e-12);
            // Co-located rows 0 and 1 must have equal rows in H:
            for j in 0..5 {
                assert!((h_mat[i * 5 + j] - h_mat[j * 5 + i]).abs() < 1e-12);
            }
        }
        for j in 0..5 {
            assert!((h_mat[j] - h_mat[5 + j]).abs() < 1e-12);
        }

        // Iso isolated node stays zero row:
        let iso = run_spur_transform(&SpurTransformRequest {
            coords: Some(&coords_col),
            n: 5,
            d: 2,
            metric: Metric::Euclidean,
            kind: SpurTransformKind::Iso {
                radius: 0.5,
                strict_lt: true,
            },
            cluster_id: None,
            y: None,
            n_y: 0,
            return_matrix: true,
            demean_lbmgls: false,
            nthreads: None,
        })
        .unwrap();
        let iso_m = iso.matrix.unwrap();
        // Points (1,0), (0,1), (1,1) are >= 1.0 apart, so rows 2, 3, 4 have zero neighbors -> 0 rows.
        for r in 2..5 {
            assert!(iso_m[r * 5..(r + 1) * 5].iter().all(|&v| v == 0.0));
        }
    }

    #[test]
    fn spur_validates_residual_and_halflife_dimensions_and_scales_invariantly() {
        let n = 12usize;
        let k_x = 3usize;
        let coords: Vec<f64> = (0..n)
            .flat_map(|i| {
                let fi = i as f64;
                [(fi * 0.7).sin() + 0.1 * fi, (fi * 1.1).cos()]
            })
            .collect();
        let y: Vec<f64> = (0..n).map(|i| ((i as f64) * 0.53).sin() + 0.2).collect();
        let mut x_in = Vec::with_capacity(n * k_x);
        for i in 0..n {
            let fi = i as f64;
            x_in.extend_from_slice(&[1.0, (fi * 0.37).cos(), (fi * 0.89).sin() + 0.1 * fi]);
        }

        // Finding 3a: q > n - k_x is rejected cleanly (q=10, 11 for n=12, k_x=3 -> dof=9).
        for q_bad in [n - k_x + 1, n - 1] {
            let emat_bad = make_deterministic_emat(q_bad, 50);
            let err0 = run_spur_i0_resid(&SpurTestRequest {
                coords: &coords,
                n,
                d: 2,
                metric: Metric::Euclidean,
                y: &y,
                n_y: 1,
                x_in: Some(&x_in),
                k_x,
                emat: &emat_bad,
                q: q_bad,
                nrep: 50,
                nthreads: Some(1),
            })
            .unwrap_err();
            assert_eq!(err0.kind(), ErrorKind::InvalidInput);
            let err1 = run_spur_i1_resid(&SpurTestRequest {
                coords: &coords,
                n,
                d: 2,
                metric: Metric::Euclidean,
                y: &y,
                n_y: 1,
                x_in: Some(&x_in),
                k_x,
                emat: &emat_bad,
                q: q_bad,
                nrep: 50,
                nthreads: Some(1),
            })
            .unwrap_err();
            assert_eq!(err1.kind(), ErrorKind::InvalidInput);
        }

        // Legal boundary q == n - k_x succeeds cleanly.
        let q_legal = n - k_x;
        let emat_legal = make_deterministic_emat(q_legal, 60);
        let ok_res = run_spur_i0_resid(&SpurTestRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            n_y: 1,
            x_in: Some(&x_in),
            k_x,
            emat: &emat_legal,
            q: q_legal,
            nrep: 60,
            nthreads: Some(1),
        })
        .unwrap();
        assert!(ok_res.lr[0].is_finite() && ok_res.pvalue[0].is_finite());

        // Co-located points reducing covariance rank below q are rejected explicitly.
        let mut coords_col = Vec::with_capacity(n * 2);
        for i in 0..n {
            let site = (i % 6) as f64;
            coords_col.extend_from_slice(&[site, (site * 1.3).sin()]);
        }
        let emat_rank = make_deterministic_emat(7, 40);
        let err_rank = run_spur_i0_resid(&SpurTestRequest {
            coords: &coords_col,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            n_y: 1,
            x_in: Some(&x_in),
            k_x,
            emat: &emat_rank,
            q: 7,
            nrep: 40,
            nthreads: Some(1),
        })
        .unwrap_err();
        assert_eq!(err_rank.kind(), ErrorKind::InvalidInput);

        // Finding 3b: Half-life rejects q=1, constant y, and zero projection, while accepting n=3, q=2.
        let emat_1 = make_deterministic_emat(1, 50);
        let err_q1 = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            emat: &emat_1,
            q: 1,
            nrep: 50,
            level: 0.95,
            nthreads: Some(1),
        })
        .unwrap_err();
        assert_eq!(err_q1.kind(), ErrorKind::InvalidInput);
        let q1_test = SpurTestRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            n_y: 1,
            x_in: Some(&x_in),
            k_x,
            emat: &emat_1,
            q: 1,
            nrep: 50,
            nthreads: Some(1),
        };
        assert_eq!(
            run_spur_i0(&q1_test).unwrap_err().kind(),
            ErrorKind::InvalidInput
        );
        assert_eq!(
            run_spur_i1(&q1_test).unwrap_err().kind(),
            ErrorKind::InvalidInput
        );
        assert_eq!(
            run_spur_i0_resid(&q1_test).unwrap_err().kind(),
            ErrorKind::InvalidInput
        );
        assert_eq!(
            run_spur_i1_resid(&q1_test).unwrap_err().kind(),
            ErrorKind::InvalidInput
        );

        let emat_2 = make_deterministic_emat(2, 50);
        // All centered covariance shapes on an equilateral triangle are scalar
        // multiples of the identity, despite retaining two positive directions.
        let equilateral = [0.0, 0.0, 1.0, 0.0, 0.5, 3.0_f64.sqrt() / 2.0];
        let err_proportional = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &equilateral,
            n: 3,
            d: 2,
            metric: Metric::Euclidean,
            y: &[0.0, 1.0, -0.5],
            emat: &emat_2,
            q: 2,
            nrep: 50,
            level: 0.95,
            nthreads: Some(1),
        })
        .unwrap_err();
        assert_eq!(err_proportional.kind(), ErrorKind::DegenerateVariance);
        let equilateral_y = [0.0, 1.0, -0.5];
        let intercept = [1.0, 1.0, 1.0];
        let equilateral_test = SpurTestRequest {
            coords: &equilateral,
            n: 3,
            d: 2,
            metric: Metric::Euclidean,
            y: &equilateral_y,
            n_y: 1,
            x_in: Some(&intercept),
            k_x: 1,
            emat: &emat_2,
            q: 2,
            nrep: 50,
            nthreads: Some(1),
        };
        assert_eq!(
            run_spur_i0(&equilateral_test).unwrap_err().kind(),
            ErrorKind::DegenerateVariance
        );
        assert_eq!(
            run_spur_i1(&equilateral_test).unwrap_err().kind(),
            ErrorKind::DegenerateVariance
        );
        assert_eq!(
            run_spur_i0_resid(&equilateral_test).unwrap_err().kind(),
            ErrorKind::DegenerateVariance
        );
        assert_eq!(
            run_spur_i1_resid(&equilateral_test).unwrap_err().kind(),
            ErrorKind::DegenerateVariance
        );
        let const_y = vec![5.0; n];
        let err_const = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &const_y,
            emat: &emat_2,
            q: 2,
            nrep: 50,
            level: 0.95,
            nthreads: Some(1),
        })
        .unwrap_err();
        assert_eq!(err_const.kind(), ErrorKind::DegenerateVariance);

        // 3-site legal half-life check (n=3, q=2).
        let coords_3 = [0.0, 0.0, 1.0, 0.0, 0.3, 0.9];
        let y_3 = [0.0, 1.0, -0.5];
        let hl_3 = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &coords_3,
            n: 3,
            d: 2,
            metric: Metric::Euclidean,
            y: &y_3,
            emat: &emat_2,
            q: 2,
            nrep: 50,
            level: 0.95,
            nthreads: Some(1),
        })
        .unwrap();
        assert!(hl_3.ci_lower.is_finite() && hl_3.ci_upper.is_finite());

        // Finding 5: Column scaling of X_in (1e-12 and 1e12) leaves projection and inference invariant,
        // while true rank deficiency is still rejected.
        let (_, _, coords_200) = load_bench_200_spur();
        let n200 = 100usize;
        let coords_sub = &coords_200[..n200 * 2];
        let y_sub: Vec<f64> = (0..n200).map(|i| ((i as f64) * 0.31).sin()).collect();
        let mut x_orig = Vec::with_capacity(n200 * 2);
        let mut x_tiny = Vec::with_capacity(n200 * 2);
        let mut x_huge = Vec::with_capacity(n200 * 2);
        for i in 0..n200 {
            let v = ((i as f64) * 0.77).cos() + 0.1 * ((i as f64) * 0.13).sin();
            x_orig.extend_from_slice(&[1.0, v]);
            x_tiny.extend_from_slice(&[1.0, v * 1e-12]);
            x_huge.extend_from_slice(&[1.0, v * 1e12]);
        }
        let emat_s = make_deterministic_emat(8, 200);
        let run_scaled = |xm: &[f64]| {
            run_spur_i0_resid(&SpurTestRequest {
                coords: coords_sub,
                n: n200,
                d: 2,
                metric: Metric::Euclidean,
                y: &y_sub,
                n_y: 1,
                x_in: Some(xm),
                k_x: 2,
                emat: &emat_s,
                q: 8,
                nrep: 200,
                nthreads: Some(1),
            })
            .unwrap()
        };
        let r_orig = run_scaled(&x_orig);
        let r_tiny = run_scaled(&x_tiny);
        let r_huge = run_scaled(&x_huge);
        assert!((r_orig.lr[0] - r_tiny.lr[0]).abs() < 1e-11);
        assert!((r_orig.lr[0] - r_huge.lr[0]).abs() < 1e-11);
        assert!((r_orig.ha_param - r_tiny.ha_param).abs() < 1e-12);
        assert!((r_orig.pvalue[0] - r_tiny.pvalue[0]).abs() < 1e-12);

        // Collinear column still rejected despite scaling.
        let mut x_dep = Vec::with_capacity(n200 * 3);
        for i in 0..n200 {
            let v = x_orig[i * 2 + 1];
            x_dep.extend_from_slice(&[1.0, v * 1e-12, v * 3e-12]);
        }
        let err_dep = orthonormal_design_basis(&x_dep, n200, 3).unwrap_err();
        assert_eq!(err_dep.kind(), ErrorKind::InvalidInput);

        // Division by the column maximum remains valid even when its reciprocal
        // would overflow, and response location does not determine degeneracy.
        let tiny = [1e-320, 0.0, 0.0, 1e-320, -1e-320, -1e-320];
        assert!(orthonormal_design_basis(&tiny, 3, 2).is_ok());
        let basis = [1.0 / 2.0_f64.sqrt(), -1.0 / 2.0_f64.sqrt(), 0.0];
        let p = prepare_y_projections(&[1e13, 1e13 + 1.0, 1e13 - 1.0], 3, 1, &basis, 1).unwrap();
        assert!(p[0].is_finite() && p[0] != 0.0);
    }

    #[test]
    fn optimized_omega_grid_quantiles_halflife_and_tiny_coords_match_reference() {
        // 1. projected_omega_grid parity for n < 256 and n > 256 (serial & parallel)
        let (_, _, coords_200) = load_bench_200_spur();
        let geom_200 = SpurGeometry::from_coords(&coords_200, 200, 2, Metric::Euclidean).unwrap();
        let (basis_200, _, mr_200, _) = prepare_lbm_basis(&geom_200, 8).unwrap();
        let c_vals = [0.25, 1.5, 8.0, 42.0, 120.0];
        let grid_serial = run_with_threads(Some(1), || {
            geom_200.projected_omega_grid(&c_vals, &mr_200, 8)
        })
        .unwrap();
        let grid_par = run_with_threads(Some(4), || {
            geom_200.projected_omega_grid(&c_vals, &mr_200, 8)
        })
        .unwrap();
        for (idx, &c) in c_vals.iter().enumerate() {
            let ref_om = geom_200.projected_omega_from_mr(c, &mr_200, 8);
            for k in 0..64 {
                assert!(
                    (grid_serial[idx][k] - ref_om[k]).abs() < 1e-12,
                    "n=200 serial diff at c={c}: {}",
                    (grid_serial[idx][k] - ref_om[k]).abs()
                );
                assert_eq!(grid_serial[idx][k].to_bits(), grid_par[idx][k].to_bits());
            }
        }
        let _ = basis_200;

        let n_big = 513usize;
        let coords_big: Vec<f64> = (0..n_big)
            .flat_map(|i| {
                let f = i as f64;
                [(f * 0.37).sin() + 0.05 * f, (f * 0.83).cos()]
            })
            .collect();
        let geom_big = SpurGeometry::from_coords(&coords_big, n_big, 2, Metric::Euclidean).unwrap();
        let (_, _, mr_big, _) = prepare_lbm_basis(&geom_big, 6).unwrap();
        let big_serial = run_with_threads(Some(1), || {
            geom_big.projected_omega_grid(&c_vals, &mr_big, 6)
        })
        .unwrap();
        let big_par = run_with_threads(Some(4), || {
            geom_big.projected_omega_grid(&c_vals, &mr_big, 6)
        })
        .unwrap();
        for (idx, &c) in c_vals.iter().enumerate() {
            let ref_om = geom_big.projected_omega_from_mr(c, &mr_big, 6);
            for k in 0..36 {
                assert!((big_serial[idx][k] - ref_om[k]).abs() < 1e-12);
                assert_eq!(big_serial[idx][k].to_bits(), big_par[idx][k].to_bits());
            }
        }

        // 2. quantiles_type7 linear selection matches naive sort bitwise
        let raw_vals: Vec<f64> = (0..503)
            .map(|i| ((i as f64) * 1.618_033_988_7).sin() * 100.0)
            .collect();
        let probs = [0.0, 0.01, 0.25, 0.5, 0.9, 0.95, 0.99, 1.0];
        let mut sorted = raw_vals.clone();
        sorted.sort_unstable_by(|a, b| a.total_cmp(b));
        let expected_q: Vec<f64> = probs
            .iter()
            .map(|&p| {
                let pos = p * ((sorted.len() - 1) as f64);
                let lo = pos.floor() as usize;
                let hi = pos.ceil() as usize;
                sorted[lo] + (pos - lo as f64) * (sorted[hi] - sorted[lo])
            })
            .collect();
        let mut work = raw_vals;
        let got_q = quantiles_type7(&mut work, &probs);
        for (&g, &e) in got_q.iter().zip(&expected_q) {
            assert_eq!(g.to_bits(), e.to_bits());
        }

        // 3. Compare both single- and multiple-tile, fast and log-scale
        // half-life paths to a factor-based calculation.
        let n_hl = 64usize;
        let coords_hl = &coords_200[..n_hl * 2];
        let y_hl: Vec<f64> = (0..n_hl).map(|i| ((i as f64) * 0.41).sin() + 0.1).collect();
        for (q_hl, nrep_hl) in [(4, TILE_COLS + 3), (8, 80), (15, 80), (60, 7), (61, 7)] {
            let emat_hl = make_deterministic_emat(q_hl, nrep_hl);
            let out_fast = run_spur_halflife(&SpurHalfLifeRequest {
                coords: coords_hl,
                n: n_hl,
                d: 2,
                metric: Metric::Euclidean,
                y: &y_hl,
                emat: &emat_hl,
                q: q_hl,
                nrep: nrep_hl,
                level: 0.95,
                nthreads: Some(1),
            })
            .unwrap();

            // Reference log-sum-exp evaluation without W_ha / fast_neg_half_pow
            let mut g = SpurGeometry::from_coords(coords_hl, n_hl, 2, Metric::Euclidean).unwrap();
            let c_bm = g.get_cbar(RHO_BM).unwrap();
            let mut sigdm_bm = g.kernel_dm(c_bm);
            let (basis_r, _) = select_leading_basis(&mut sigdm_bm, n_hl, q_hl).unwrap();
            let mr = col_demean_matrix(&basis_r, n_hl, q_hl);
            let p_vec = prepare_y_projections(&y_hl, n_hl, 1, &basis_r, q_hl).unwrap();
            let neg_ln_half = -(0.5_f64).ln();
            let ho_f: Vec<_> = out_fast
                .hl_grid_ho
                .iter()
                .map(|&h| {
                    factor_omega(&g.projected_omega_from_mr(neg_ln_half / h, &mr, q_hl), q_hl)
                        .unwrap()
                })
                .collect();
            let ha_f: Vec<_> = (0..HL_HA_LEN)
                .map(|i| {
                    let h = 0.001 + (1.0 - 0.001) * (i as f64) / ((HL_HA_LEN - 1) as f64);
                    factor_omega(&g.projected_omega_from_mr(neg_ln_half / h, &mr, q_hl), q_hl)
                        .unwrap()
                })
                .collect();
            let q_f = q_hl as f64;
            for (idx, ho) in ho_f.iter().enumerate() {
                let mut xc_obs = vec![0.0_f64; q_hl];
                linalg::multiply(&ho.ch_omi.upper, &p_vec, &mut xc_obs, q_hl, q_hl, 1, false);
                let log_ho_obs = ho.ch_omi.log_sqrt_det
                    - 0.5 * q_f * xc_obs.iter().map(|v| v * v).sum::<f64>().ln();
                let lr_obs: f64 = ha_f
                    .iter()
                    .map(|ha| {
                        let mut xc = vec![0.0_f64; q_hl];
                        linalg::multiply(&ha.ch_omi.upper, &p_vec, &mut xc, q_hl, q_hl, 1, false);
                        let lha = ha.ch_omi.log_sqrt_det
                            - 0.5 * q_f * xc.iter().map(|v| v * v).sum::<f64>().ln();
                        (lha - log_ho_obs).exp()
                    })
                    .sum::<f64>()
                    / (HL_HA_LEN as f64);
                let mut exceed = 0usize;
                for c in 0..nrep_hl {
                    let e_col: Vec<f64> = (0..q_hl).map(|r| emat_hl[r * nrep_hl + c]).collect();
                    let mut x = vec![0.0_f64; q_hl];
                    linalg::multiply(&ho.ch_om.lower, &e_col, &mut x, q_hl, q_hl, 1, false);
                    let mut z0 = vec![0.0_f64; q_hl];
                    linalg::multiply(&ho.ch_omi.upper, &x, &mut z0, q_hl, q_hl, 1, false);
                    let log_ho = ho.ch_omi.log_sqrt_det
                        - 0.5 * q_f * z0.iter().map(|v| v * v).sum::<f64>().ln();
                    let avg: f64 = ha_f
                        .iter()
                        .map(|ha| {
                            let mut z1 = vec![0.0_f64; q_hl];
                            linalg::multiply(&ha.ch_omi.upper, &x, &mut z1, q_hl, q_hl, 1, false);
                            let lha = ha.ch_omi.log_sqrt_det
                                - 0.5 * q_f * z1.iter().map(|v| v * v).sum::<f64>().ln();
                            (lha - log_ho).exp()
                        })
                        .sum::<f64>()
                        / (HL_HA_LEN as f64);
                    if avg > lr_obs {
                        exceed += 1;
                    }
                }
                let ref_p = (exceed as f64) / (nrep_hl as f64);
                assert!(
                    (out_fast.pv_mat[idx] - ref_p).abs() < 1e-12,
                    "q={q_hl} idx={idx}: fast={} ref={ref_p}",
                    out_fast.pv_mat[idx]
                );
            }
        }

        // 4. Tiny positive coordinate scale invariance (coords * 1e-12)
        let coords_tiny: Vec<f64> = coords_hl.iter().map(|&v| v * 1e-12).collect();
        let emat_tiny = make_deterministic_emat(6, 100);
        let hl_orig = run_spur_halflife(&SpurHalfLifeRequest {
            coords: coords_hl,
            n: n_hl,
            d: 2,
            metric: Metric::Euclidean,
            y: &y_hl,
            emat: &emat_tiny,
            q: 6,
            nrep: 100,
            level: 0.95,
            nthreads: Some(1),
        })
        .unwrap();
        let hl_tiny = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &coords_tiny,
            n: n_hl,
            d: 2,
            metric: Metric::Euclidean,
            y: &y_hl,
            emat: &emat_tiny,
            q: 6,
            nrep: 100,
            level: 0.95,
            nthreads: Some(1),
        })
        .unwrap();
        assert!((hl_orig.ci_lower - hl_tiny.ci_lower).abs() < 1e-12);
        assert!((hl_orig.ci_upper - hl_tiny.ci_upper).abs() < 1e-12);
        assert!((hl_tiny.max_dist - hl_orig.max_dist * 1e-12).abs() < 1e-22);

        // 5. fast_neg_half_pow overflows ratio^(q/2) -> 0.0 at q=60, ratio=1e11; verify
        // underflow guard forces log-space fallback instead of zeroing a dominant term.
        let pow_underflow = fast_neg_half_pow(1e11, 60);
        assert_eq!(pow_underflow, 0.0);
        assert!(pow_underflow < f64::MIN_POSITIVE);
        let dominant_log = 140.0_f64 * 10.0_f64.ln() - 30.0 * (1e11_f64).ln();
        let secondary_log = 0.0_f64 - 30.0 * (1e7_f64).ln();
        let true_log_sum = log_sum_exp(&[dominant_log, secondary_log]);
        assert!(true_log_sum > -200.0 * 10.0_f64.ln());
        assert!(true_log_sum - secondary_log > 40.0);
    }
}
