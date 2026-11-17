//! Whole-operation SPUR diagnostics, half-life confidence intervals, and spatial transforms.

use crate::distance::{distance_matrix, lower_triangle, Metric};
use crate::threads::{is_serial, run_with_threads};
use crate::{basis, checked_len, linalg, CoreError};
use rayon::prelude::*;

const SUM_CHUNK_SIZE: usize = 32_768;
const TILE_COLS: usize = 4_096;
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
        if !max_dist.is_finite() || max_dist <= 1e-10 {
            return Err(CoreError::invalid_input(
                "all coordinates are identical (or nearly so); normalization requires distinct locations",
            ));
        }
        let inv_max = max_dist.recip();
        let norm_distmat: Vec<f64> = raw_distmat.iter().map(|&v| v * inv_max).collect();
        let lower_tri = lower_triangle(&norm_distmat, n);
        Ok(Self {
            n,
            raw_distmat,
            max_dist,
            norm_distmat,
            lower_tri,
            cbar_cache: Vec::new(),
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

        let vd = &self.lower_tri;
        let mut c0 = 10.0_f64;
        let mut c1 = 10.0_f64;

        let mut it = 0usize;
        loop {
            let v = deterministic_mean_exp(vd, c0);
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
            let v = deterministic_mean_exp(vd, c1);
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
            let v = deterministic_mean_exp(vd, cm);
            if v < rhobar {
                c1 = cm;
            } else {
                c0 = cm;
            }
        }

        let cbar = (c0 * c1).sqrt();
        self.cbar_cache.push((rhobar, cbar));
        Ok(cbar)
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
        let mut kern = vec![0.0_f64; self.norm_distmat.len()];
        linalg::exp_scaled(&self.norm_distmat, -c, &mut kern);
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

fn deterministic_mean_exp(vd: &[f64], c: f64) -> f64 {
    if vd.is_empty() {
        return f64::NAN;
    }
    let chunk_sums: Vec<f64> = if !is_serial() && vd.len() >= 4 * SUM_CHUNK_SIZE {
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
    src.iter()
        .zip(&q_qt_src)
        .map(|(&s, &p)| s - p)
        .collect()
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
        return Err(CoreError::degenerate("SPUR covariance has non-positive trace"));
    }
    Ok(a.iter().zip(b).any(|(&x, &y)| (x / ta - y / tb).abs() > 1e-12))
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
    vals.sort_unstable_by(|a, b| a.total_cmp(b));
    probs
        .iter()
        .map(|&prob| {
            if m == 1 {
                vals[0]
            } else {
                let pos = prob.clamp(0.0, 1.0) * ((m - 1) as f64);
                let lo = pos.floor() as usize;
                let hi = pos.ceil() as usize;
                let frac = pos - (lo as f64);
                vals[lo] + frac * (vals[hi] - vals[lo])
            }
        })
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
        linalg::multiply(&alt.ch_omi.upper, &self.ch_om0_lower, &mut ho, q, q, q, false);
        linalg::multiply(&self.ch_om0i_upper, &alt.ch_om.lower, &mut ha, q, q, q, false);

        let mut qa_o = vec![0.0_f64; nrep];
        let mut qo_a = vec![0.0_f64; nrep];
        let mut ya_buf = vec![0.0_f64; q * TILE_COLS.min(nrep)];
        let mut yo_buf = vec![0.0_f64; q * TILE_COLS.min(nrep)];
        let mut lr_o = vec![0.0_f64; nrep];
        let mut tile_e = vec![0.0_f64; q * TILE_COLS.min(nrep)];

        for start in (0..nrep).step_by(TILE_COLS) {
            let cols = (nrep - start).min(TILE_COLS);
            for r in 0..q {
                tile_e[r * cols..(r + 1) * cols]
                    .copy_from_slice(&self.emat[r * nrep + start..r * nrep + start + cols]);
            }
            linalg::multiply(&ho, &tile_e[..q * cols], &mut ya_buf[..q * cols], q, q, cols, false);
            linalg::multiply(&ha, &tile_e[..q * cols], &mut yo_buf[..q * cols], q, q, cols, false);
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
        let exceed = (0..nrep)
            .filter(|&k| (qo_a[k] / self.qe[k]) > cv)
            .count();
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
        // Canonicalize each eigenvector's sign across LAPACK/Jacobi/host backends
        // (largest-magnitude entry positive, earliest index breaks exact ties).
        let mut max_idx = 0usize;
        let mut max_abs = -1.0_f64;
        for (i, &val) in src.iter().enumerate() {
            let a = val.abs();
            if a > max_abs {
                max_abs = a;
                max_idx = i;
            }
        }
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
        y_dm_norm[j] = (0..n).map(|i| y_dm[i * n_y + j].powi(2)).sum::<f64>().sqrt();
    }
    let mut p_mat = vec![0.0_f64; q * n_y];
    linalg::multiply(basis_r, &y_dm, &mut p_mat, q, n, n_y, true);
    for j in 0..n_y {
        let p_max = (0..q).map(|r| p_mat[r * n_y + j].abs()).fold(0.0_f64, f64::max);
        if !p_max.is_finite() || p_max == 0.0 {
            return Err(CoreError::degenerate(
                "projected SPUR vector has zero norm",
            ));
        }
        let p_norm = p_max
            * (0..q)
                .map(|r| (p_mat[r * n_y + j] / p_max).powi(2))
                .sum::<f64>()
                .sqrt();
        if !p_norm.is_finite() || p_norm <= 1e-12 * y_dm_norm[j] {
            return Err(CoreError::degenerate(
                "projected SPUR vector has zero norm",
            ));
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
    let mix = |g: f64| -> Vec<f64> {
        om_i0
            .iter()
            .zip(om_bm)
            .map(|(&a, &b)| a + g * b)
            .collect()
    };

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
    let mut tile_e = vec![0.0_f64; q * TILE_COLS.min(nrep)];
    let mut buf_ho = vec![0.0_f64; q * TILE_COLS.min(nrep)];
    let mut buf_ha = vec![0.0_f64; q * TILE_COLS.min(nrep)];
    let mut q_ho = vec![0.0_f64; TILE_COLS.min(nrep)];
    let mut q_ha = vec![0.0_f64; TILE_COLS.min(nrep)];

    for start in (0..nrep).step_by(TILE_COLS) {
        let cols = (nrep - start).min(TILE_COLS);
        for r in 0..q {
            tile_e[r * cols..(r + 1) * cols]
                .copy_from_slice(&emat[r * nrep + start..r * nrep + start + cols]);
        }
        linalg::multiply(&a_ho, &tile_e[..q * cols], &mut buf_ho[..q * cols], q, q, cols, false);
        linalg::multiply(&a_ha, &tile_e[..q * cols], &mut buf_ha[..q * cols], q, q, cols, false);
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
            RHO_GRID_MIN
                + (RHO_GRID_MAX - RHO_GRID_MIN) * (i as f64) / ((RHO_GRID_LEN - 1) as f64)
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
    let mut c_grid = Vec::with_capacity(RHO_GRID_LEN);
    for &rho_i in &rho_grid {
        c_grid.push(if rho_i > 0.0 {
            Some(geom.get_cbar(rho_i)?)
        } else {
            None
        });
    }

    let geom_ref: &SpurGeometry = geom;
    let probs = [1.0 - SZ_LEVELS[0], 1.0 - SZ_LEVELS[1], 1.0 - SZ_LEVELS[2]];
    let eval_grid_point = |c_opt: Option<f64>| -> Result<(Vec<f64>, Vec<f64>), CoreError> {
        let ch_om_lower = if let Some(c_i) = c_opt {
            let om_i = geom_ref.projected_omega_from_mr(c_i, mr, q);
            cholesky_factor(&om_i, q)?.lower
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
        c_grid
            .iter()
            .map(|&c_opt| eval_grid_point(c_opt))
            .collect::<Result<Vec<_>, _>>()?
    } else {
        c_grid
            .par_iter()
            .map(|&c_opt| eval_grid_point(c_opt))
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
            &mut geom,
            &basis_r,
            tied,
            &mr,
            &om_bm,
            req.y,
            req.n_y,
            req.emat,
            req.q,
            req.nrep,
        )
    })?
}

pub fn run_spur_i1(req: &SpurTestRequest<'_>) -> Result<SpurI1Output, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    run_with_threads(req.nthreads, || {
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_r, tied, mr, om_bm) = prepare_lbm_basis(&geom, req.q)?;
        run_spur_i1_from_basis(
            &mut geom,
            &basis_r,
            tied,
            &mr,
            &om_bm,
            req.y,
            req.n_y,
            req.emat,
            req.q,
            req.nrep,
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
            &mut geom,
            &basis_r,
            tied,
            &mr,
            &om_bm,
            req.y,
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
        let mut geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let (basis_r, tied, mr, om_bm) = prepare_resid_basis(&mut geom, x_in, req.k_x, req.q)?;
        run_spur_i1_from_basis(
            &mut geom,
            &basis_r,
            tied,
            &mr,
            &om_bm,
            req.y,
            req.n_y,
            req.emat,
            req.q,
            req.nrep,
        )
    })?
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

        // If every tested covariance has the same shape, the angular likelihood
        // is constant. Strict simulation exceedances then have no meaningful
        // interpretation, even when q > 1 (e.g. an equilateral triangle).
        let reference_shape = geom.projected_omega_from_mr(c_grid_ho[0], &mr, q);
        let mut informative = false;
        for &c in c_grid_ho.iter().chain(&c_grid_ha) {
            let candidate = geom.projected_omega_from_mr(c, &mr, q);
            if covariance_shapes_differ(&candidate, &reference_shape, q)? {
                informative = true;
                break;
            }
        }
        if !informative {
            return Err(CoreError::degenerate(
                "half-life experiment is uninformative: retained covariance shapes are proportional",
            ));
        }

        let build_factor = |c: f64| -> Result<OmegaPairFactor, CoreError> {
            let om = geom.projected_omega_from_mr(c, &mr, q);
            factor_omega(&om, q)
        };

        let ho_factors: Vec<OmegaPairFactor> = if is_serial() {
            c_grid_ho
                .iter()
                .map(|&c| build_factor(c))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            c_grid_ho
                .par_iter()
                .map(|&c| build_factor(c))
                .collect::<Result<Vec<_>, _>>()?
        };
        let ha_factors: Vec<OmegaPairFactor> = if is_serial() {
            c_grid_ha
                .iter()
                .map(|&c| build_factor(c))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            c_grid_ha
                .par_iter()
                .map(|&c| build_factor(c))
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

        let eval_null_hl = |ho: &OmegaPairFactor| -> Result<f64, CoreError> {
            let mut xc_obs = vec![0.0_f64; q];
            linalg::multiply(&ho.ch_omi.upper, &p_vec, &mut xc_obs, q, q, 1, false);
            let qnorm_obs: f64 = xc_obs.iter().map(|&v| v * v).sum();
            if !qnorm_obs.is_finite() || qnorm_obs <= 0.0 {
                return Err(CoreError::degenerate(
                    "projected SPUR vector has zero quadratic form",
                ));
            }
            let log_ho_obs = ho.ch_omi.log_sqrt_det - 0.5 * q_f * qnorm_obs.ln();
            let lr_obs: f64 = ha_x_log
                .iter()
                .map(|&lha| (lha - log_ho_obs).exp())
                .sum::<f64>()
                / (HL_HA_LEN as f64);

            let mut a_ho = vec![0.0_f64; q * q];
            linalg::multiply(&ho.ch_omi.upper, &ho.ch_om.lower, &mut a_ho, q, q, q, false);
            let a_ha_list: Vec<Vec<f64>> = ha_factors
                .iter()
                .map(|ha| {
                    let mut m = vec![0.0_f64; q * q];
                    linalg::multiply(&ha.ch_omi.upper, &ho.ch_om.lower, &mut m, q, q, q, false);
                    m
                })
                .collect();

            let mut exceed = 0usize;
            let max_tile = TILE_COLS.min(nrep);
            let mut tile_e = vec![0.0_f64; q * max_tile];
            let mut buf = vec![0.0_f64; q * max_tile];
            let mut log_ho = vec![0.0_f64; max_tile];
            let mut sum_ha_over_ho = vec![0.0_f64; max_tile];

            for start in (0..nrep).step_by(TILE_COLS) {
                let cols = (nrep - start).min(TILE_COLS);
                for r in 0..q {
                    tile_e[r * cols..(r + 1) * cols]
                        .copy_from_slice(&req.emat[r * nrep + start..r * nrep + start + cols]);
                }
                linalg::multiply(&a_ho, &tile_e[..q * cols], &mut buf[..q * cols], q, q, cols, false);
                log_ho[..cols].fill(0.0);
                for r in 0..q {
                    let row = &buf[r * cols..(r + 1) * cols];
                    for c in 0..cols {
                        log_ho[c] += row[c] * row[c];
                    }
                }
                for c in 0..cols {
                    log_ho[c] = ho.ch_omi.log_sqrt_det - 0.5 * q_f * log_ho[c].ln();
                }

                sum_ha_over_ho[..cols].fill(0.0);
                for (j, a_ha) in a_ha_list.iter().enumerate() {
                    linalg::multiply(a_ha, &tile_e[..q * cols], &mut buf[..q * cols], q, q, cols, false);
                    let log_c_ha = ha_factors[j].ch_omi.log_sqrt_det;
                    for c in 0..cols {
                        let mut qsum = 0.0_f64;
                        for r in 0..q {
                            let v = buf[r * cols + c];
                            qsum += v * v;
                        }
                        let log_ha = log_c_ha - 0.5 * q_f * qsum.ln();
                        sum_ha_over_ho[c] += (log_ha - log_ho[c]).exp();
                    }
                }
                let inv_nha = 1.0 / (HL_HA_LEN as f64);
                for c in 0..cols {
                    if (sum_ha_over_ho[c] * inv_nha) > lr_obs {
                        exceed += 1;
                    }
                }
            }
            Ok((exceed as f64) / (nrep as f64))
        };

        let pv_mat: Vec<f64> = if is_serial() {
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

        let alpha = 1.0 - req.level;
        let mut ci_lower = f64::INFINITY;
        let mut ci_upper = f64::NEG_INFINITY;
        let mut any_accepted = false;
        for (&p, &hl) in pv_mat.iter().zip(&hl_grid_ho) {
            if p > alpha {
                any_accepted = true;
                ci_lower = ci_lower.min(hl);
                ci_upper = ci_upper.max(hl);
            }
        }
        if !any_accepted {
            ci_lower = f64::NAN;
            ci_upper = f64::NAN;
        }

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

pub fn run_spur_transform(req: &SpurTransformRequest<'_>) -> Result<SpurTransformOutput, CoreError> {
    let n = req.n;
    if n < 2 {
        return Err(CoreError::invalid_input(
            "spatial transform requires at least 2 observations",
        ));
    }
    if let Some(y) = req.y {
        if req.n_y == 0 || y.len() != checked_len(&[n, req.n_y])? {
            return Err(CoreError::invalid_input("invalid Y dimensions for transform"));
        }
        if y.iter().any(|v| !v.is_finite()) {
            return Err(CoreError::invalid_input("Y values must be finite"));
        }
    }

    run_with_threads(req.nthreads, || {
        match req.kind {
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
                Ok(SpurTransformOutput { matrix, transformed })
            }
            SpurTransformKind::Nn {
                exclude_zero_distances,
            } => {
                let coords = req.coords.ok_or_else(|| {
                    CoreError::invalid_input("nn transform requires coordinates")
                })?;
                let geom = SpurGeometry::from_coords(coords, n, req.d, req.metric)?;
                let mut mat = vec![0.0_f64; n * n];
                for i in 0..n {
                    let row_d = &geom.norm_distmat[i * n..(i + 1) * n];
                    let mut min_d = f64::INFINITY;
                    for j in 0..n {
                        if i == j {
                            continue;
                        }
                        let dval = if exclude_zero_distances && row_d[j] == 0.0 {
                            1e10
                        } else {
                            row_d[j]
                        };
                        min_d = min_d.min(dval);
                    }
                    let mut cnt = 0usize;
                    for j in 0..n {
                        if i == j {
                            continue;
                        }
                        let dval = if exclude_zero_distances && row_d[j] == 0.0 {
                            1e10
                        } else {
                            row_d[j]
                        };
                        if dval == min_d {
                            cnt += 1;
                        }
                    }
                    let w = if cnt > 0 { (cnt as f64).recip() } else { 0.0 };
                    for j in 0..n {
                        if i == j {
                            mat[i * n + j] = 1.0;
                        } else {
                            let dval = if exclude_zero_distances && row_d[j] == 0.0 {
                                1e10
                            } else {
                                row_d[j]
                            };
                            if dval == min_d {
                                mat[i * n + j] = -w;
                            }
                        }
                    }
                }
                let transformed = req.y.map(|y| {
                    let mut out = vec![0.0_f64; n * req.n_y];
                    linalg::multiply(&mat, y, &mut out, n, n, req.n_y, false);
                    out
                });
                Ok(SpurTransformOutput {
                    matrix: if req.return_matrix { Some(mat) } else { None },
                    transformed,
                })
            }
            SpurTransformKind::Iso { radius, strict_lt } => {
                if !radius.is_finite() || radius <= 0.0 {
                    return Err(CoreError::invalid_input(
                        "iso transform radius must be positive and finite",
                    ));
                }
                let coords = req.coords.ok_or_else(|| {
                    CoreError::invalid_input("iso transform requires coordinates")
                })?;
                let raw_distmat = distance_matrix(coords, n, req.d, req.metric)?;
                let mut mat = vec![0.0_f64; n * n];
                for i in 0..n {
                    let row_d = &raw_distmat[i * n..(i + 1) * n];
                    let in_rad = |j: usize| {
                        let d = row_d[j];
                        i != j && d > 0.0 && if strict_lt { d < radius } else { d <= radius }
                    };
                    let cnt = (0..n).filter(|&j| in_rad(j)).count();
                    if cnt == 0 {
                        // Isolated observations stay zero row (matches Stata/R/Python iso_matrix).
                        continue;
                    }
                    let w = (cnt as f64).recip();
                    for j in 0..n {
                        if i == j {
                            mat[i * n + j] = 1.0;
                        } else if in_rad(j) {
                            mat[i * n + j] = -w;
                        }
                    }
                }
                let transformed = req.y.map(|y| {
                    let mut out = vec![0.0_f64; n * req.n_y];
                    linalg::multiply(&mat, y, &mut out, n, n, req.n_y, false);
                    out
                });
                Ok(SpurTransformOutput {
                    matrix: if req.return_matrix { Some(mat) } else { None },
                    transformed,
                })
            }
            SpurTransformKind::LbmGls => {
                let coords = req.coords.ok_or_else(|| {
                    CoreError::invalid_input("lbmgls transform requires coordinates")
                })?;
                let geom = SpurGeometry::from_coords(coords, n, req.d, req.metric)?;
                apply_lbmgls_with_geom(
                    &geom,
                    req.y,
                    req.n_y,
                    req.return_matrix,
                    req.demean_lbmgls,
                )
            }
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

    Ok(SpurTransformOutput { matrix, transformed })
}

pub fn run_spur_pipeline(req: &SpurPipelineRequest<'_>) -> Result<SpurPipelineOutput, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep)?;
    if req.y.len() != req.n || req.y.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input("invalid y vector for SPUR pipeline"));
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
        let path = format!(
            "{}/../../../verification/data/bench_200.csv",
            env!("CARGO_MANIFEST_DIR")
        );
        let text = std::fs::read_to_string(path).expect("bench_200.csv");
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
        for (&a, &b) in pipe_1
            .transformed_vars
            .iter()
            .zip(&pipe_4.transformed_vars)
        {
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
        assert!((hl_1.ci_lower - 0.152_363_636_363_636_37).abs() < 1e-12, "got ci_lower={}, ci_upper={}", hl_1.ci_lower, hl_1.ci_upper);
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
            coords: &coords, n, d: 2, metric: Metric::Euclidean, y: &y, n_y: 1,
            x_in: Some(&x_in), k_x, emat: &emat_1, q: 1, nrep: 50, nthreads: Some(1),
        };
        assert_eq!(run_spur_i0(&q1_test).unwrap_err().kind(), ErrorKind::InvalidInput);
        assert_eq!(run_spur_i1(&q1_test).unwrap_err().kind(), ErrorKind::InvalidInput);
        assert_eq!(run_spur_i0_resid(&q1_test).unwrap_err().kind(), ErrorKind::InvalidInput);
        assert_eq!(run_spur_i1_resid(&q1_test).unwrap_err().kind(), ErrorKind::InvalidInput);

        let emat_2 = make_deterministic_emat(2, 50);
        // All centered covariance shapes on an equilateral triangle are scalar
        // multiples of the identity, despite retaining two positive directions.
        let equilateral = [0.0, 0.0, 1.0, 0.0, 0.5, 3.0_f64.sqrt() / 2.0];
        let err_proportional = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &equilateral, n: 3, d: 2, metric: Metric::Euclidean,
            y: &[0.0, 1.0, -0.5], emat: &emat_2, q: 2, nrep: 50,
            level: 0.95, nthreads: Some(1),
        }).unwrap_err();
        assert_eq!(err_proportional.kind(), ErrorKind::DegenerateVariance);
        let equilateral_y = [0.0, 1.0, -0.5];
        let intercept = [1.0, 1.0, 1.0];
        let equilateral_test = SpurTestRequest {
            coords: &equilateral, n: 3, d: 2, metric: Metric::Euclidean,
            y: &equilateral_y, n_y: 1, x_in: Some(&intercept), k_x: 1,
            emat: &emat_2, q: 2, nrep: 50, nthreads: Some(1),
        };
        assert_eq!(run_spur_i0(&equilateral_test).unwrap_err().kind(), ErrorKind::DegenerateVariance);
        assert_eq!(run_spur_i1(&equilateral_test).unwrap_err().kind(), ErrorKind::DegenerateVariance);
        assert_eq!(run_spur_i0_resid(&equilateral_test).unwrap_err().kind(), ErrorKind::DegenerateVariance);
        assert_eq!(run_spur_i1_resid(&equilateral_test).unwrap_err().kind(), ErrorKind::DegenerateVariance);
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
}
