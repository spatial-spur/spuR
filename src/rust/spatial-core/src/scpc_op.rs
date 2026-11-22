//! Whole-operation exact SCPC inference pipeline.

use crate::basis::{scpc_basis, ScpcBasis};
use crate::calibration::{correlation_floor, solve_c0_from_avc};
use crate::covariance;
use crate::distance::{distance_matrix, lower_triangle, Metric};
use crate::linalg;
use crate::nystrom::{
    is_psd_prefix_grid, landmark_basis, sampled_covariance_grid, LandmarkRequest, SampledPairs,
};
use crate::student_t::ln_gamma;
use crate::threads::{is_serial, run_with_threads};
use crate::{checked_len, CoreError};
use rayon::prelude::*;

mod probability;
use probability::ProbabilityWorkspace;

const MIN_AVC: f64 = 0.00001;
const CGRID_FAC: f64 = 1.2;
const DEFAULT_CV_TOL: f64 = 0.001;
const EXTRA_CV_LEVELS: [f64; 3] = [0.32, 0.10, 0.01];
const SCALE_REL_TOL: f64 = 1e-12;

/// Complete input bundle for one exact SCPC inference call.
///
/// All two-dimensional arrays are contiguous and row-major.
pub struct ExactScpcRequest<'a> {
    /// Coordinate matrix (`n x d`). For geodesic coordinates, select
    /// `Metric::GreatCircleFractionLonLat` or `Metric::GreatCircleFractionLatLon`.
    pub coords: &'a [f64],
    pub n: usize,
    pub d: usize,
    pub metric: Metric,
    pub avc: f64,
    pub uncond: bool,
    pub cvs: bool,
    /// Estimated regression coefficients (`k_use`).
    pub coef: &'a [f64],
    /// Coefficient influence contributions, scaled by the site count (`n x k_use`).
    pub wj: &'a [f64],
    /// Cluster membership indices in `0..n` for `m_obs` observations when clustered,
    /// or `None` when spatial sites equal observations (`m_obs == n`).
    pub cluster_id: Option<&'a [usize]>,
    /// Optional 1-way absorbed fixed-effect membership indices in `0..n_fe` for `m_obs`
    /// observations when `uncond == false`, or `None` when not demeaning by a 1-way FE.
    pub fe_id: Option<&'a [usize]>,
    /// Conditional regressor direction (`m_obs x k_use`): `xj` for both OLS and IV.
    /// Required when `uncond == false`.
    pub xj: Option<&'a [f64]>,
    /// Conditional scaling direction (`m_obs x k_use`): `sign(xj)` or within-cluster
    /// normalized `xj` for both OLS and IV. Required when `uncond == false`.
    pub xjs: Option<&'a [f64]>,
    /// Left factor `U` (`m_obs x proj_rank`) of conditional projection `A -> A - U (V' A)`
    /// (for OLS `Q_x`, for IV `V = P_Z X (X' P_Z X)^{-1}` corresponding to `M' = I - V U'`).
    /// Required when `uncond == false` and `proj_rank > 0`.
    pub proj_u: Option<&'a [f64]>,
    /// Right factor `V` (`m_obs x proj_rank`) of conditional projection `A -> A - U (V' A)`
    /// (for OLS `Q_x`, for IV `U = X` corresponding to `M' = I - V U'`).
    /// Required when `uncond == false` and `proj_rank > 0`.
    pub proj_v: Option<&'a [f64]>,
    pub proj_rank: usize,
    /// Critical value solver tolerance in ordinary t-units (released default `0.001`).
    pub cv_tolerance: f64,
    /// Optional thread budget (`Some(1)` runs sequentially on caller thread).
    pub nthreads: Option<usize>,
}

/// Complete input bundle for one approximate (Nyström) SCPC inference call.
///
/// Caller must pass `coords` already normalized via `nystrom::normalize_coords` (or host
/// equivalent), with `wj` and site-indexed / observation inputs aligned to that row permutation.
pub struct ApproxScpcRequest<'a> {
    pub coords: &'a [f64],
    pub n: usize,
    pub d: usize,
    pub metric: Metric,
    pub avc: f64,
    pub uncond: bool,
    pub cvs: bool,
    pub coef: &'a [f64],
    pub wj: &'a [f64],
    pub cluster_id: Option<&'a [usize]>,
    pub fe_id: Option<&'a [usize]>,
    pub xj: Option<&'a [f64]>,
    pub xjs: Option<&'a [f64]>,
    pub proj_u: Option<&'a [f64]>,
    pub proj_v: Option<&'a [f64]>,
    pub proj_rank: usize,
    pub cv_tolerance: f64,
    pub nthreads: Option<usize>,
    pub seed: u32,
    pub landmarks: usize,
    pub samples: usize,
    pub pair_count: usize,
}

/// Output of exact and approximate SCPC inference.
#[derive(Debug, Clone)]
pub struct ScpcOutput {
    /// Row-major `k_use x 6`: `[coef, se, tstat, pvalue, ci_lower, ci_upper]`.
    pub scpcstats: Vec<f64>,
    /// Optional row-major `k_use x 4`: two-sided critical values at `[32%, 10%, 5%, 1%]`.
    pub scpccvs: Option<Vec<f64>>,
    /// Selected spatial basis `n x (q + 1)` in row-major order.
    pub wfin: Vec<f64>,
    pub c0: f64,
    pub cmax: f64,
    pub cvfin: f64,
    pub q: usize,
    pub tied: bool,
}

fn map_indices<T: Send>(
    len: usize,
    f: impl Fn(usize) -> Result<T, CoreError> + Send + Sync,
) -> Result<Vec<T>, CoreError> {
    if is_serial() || len <= 1 {
        (0..len).map(f).collect()
    } else {
        (0..len).into_par_iter().map(f).collect()
    }
}

fn initial_qmax(avc: f64) -> usize {
    if avc >= 0.05 {
        10
    } else if avc >= 0.01 {
        20
    } else if avc >= 0.005 {
        60
    } else {
        120
    }
}

fn covariance_grids(
    dist: &[f64],
    n: usize,
    weights: &[&[f64]],
    widths: &[usize],
    scales: &[f64],
) -> Result<Vec<Vec<f64>>, CoreError> {
    if dist.len() != checked_len(&[n, n])? {
        return Err(CoreError::invalid_input(
            "distance matrix dimensions do not match",
        ));
    }
    covariance::project(n, weights, widths, scales, |i, j, ni, nj, tile| {
        for (row, values) in tile.chunks_exact_mut(nj).take(ni).enumerate() {
            let start = (i + row) * n + j;
            values.copy_from_slice(&dist[start..start + nj]);
        }
    })
}

fn select_final_q(
    oms: &[f64],
    width: usize,
    qmax: usize,
    tolerance: f64,
) -> Result<(usize, f64), CoreError> {
    let results = map_indices(qmax, |idx| {
        let q = idx + 1;
        let mut ws = ProbabilityWorkspace::new(q);
        let cv = ws.critical(oms, width, q, 0.05, tolerance)?;
        let qf = q as f64;
        let gamma_ratio = (ln_gamma(0.5 * (qf + 1.0)) - ln_gamma(0.5 * qf)).exp();
        let length = cv * gamma_ratio / qf.sqrt();
        Ok((q, cv, length))
    })?;
    let best = results
        .into_iter()
        .min_by(|a, b| a.2.total_cmp(&b.2).then_with(|| a.0.cmp(&b.0)))
        .ok_or(CoreError::invalid_input("empty qmax search range"))?;
    Ok((best.0, best.1))
}

fn scaled_norm(values: impl Iterator<Item = f64> + Clone) -> f64 {
    let scale = values.clone().map(f64::abs).fold(0.0_f64, f64::max);
    if scale == 0.0 || !scale.is_finite() {
        return scale;
    }
    scale
        * values
            .map(|value| (value / scale).powi(2))
            .sum::<f64>()
            .sqrt()
}

/// SCPC probabilities are homogeneous in all columns of a coefficient's W.
/// Equilibrate them together before forming the quadratic covariance grid.
fn normalize_covariance_basis(mut weights: Vec<f64>) -> Vec<f64> {
    let scale = weights
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
    if scale > 0.0 {
        for value in &mut weights {
            *value /= scale;
        }
    }
    weights
}

#[allow(clippy::too_many_arguments)]
fn build_conditional_basis(
    wfin: &[f64],
    n: usize,
    q: usize,
    cluster_id: Option<&[usize]>,
    fe_id: Option<&[usize]>,
    xj: &[f64],
    xjs: &[f64],
    k_use: usize,
    col_j: usize,
    proj_u: Option<&[f64]>,
    proj_v: Option<&[f64]>,
    proj_rank: usize,
) -> Result<Vec<f64>, CoreError> {
    let width = q + 1;
    let mut wx = vec![0.0; checked_len(&[n, width])?];
    let mut pre_max = 0.0_f64;

    if let Some(fe) = fe_id {
        let m_obs = fe.len();
        let n_fe = fe.iter().copied().max().map_or(0, |m| m + 1);
        let mut counts = vec![0usize; n_fe];
        for &f in fe {
            counts[f] += 1;
        }
        if let Some(cl_id) = cluster_id {
            let mut diag_g = vec![0.0_f64; n];
            for (i, &g) in cl_id.iter().enumerate() {
                diag_g[g] += xj[i * k_use + col_j] * xjs[i * k_use + col_j];
            }
            for g in 0..n {
                wx[g * width] = wfin[g * width] * diag_g[g];
            }
            for c in 0..q {
                let column = (0..n).map(|g| diag_g[g] * wfin[g * width + 1 + c]);
                pre_max = pre_max.max(scaled_norm(column));
            }
        } else {
            for i in 0..n {
                let xj_i = xj[i * k_use + col_j];
                let xjs_i = xjs[i * k_use + col_j];
                wx[i * width] = wfin[i * width] * xj_i * xjs_i;
            }
            for c in 0..q {
                let column = (0..n).map(|i| {
                    let idx = i * k_use + col_j;
                    wfin[i * width + 1 + c] * xj[idx] * xjs[idx]
                });
                pre_max = pre_max.max(scaled_norm(column));
            }
        }

        if q > 0 {
            let u_opt = if proj_rank > 0 {
                Some(proj_u.ok_or(CoreError::invalid_input(
                    "missing conditional projection U factor",
                ))?)
            } else {
                None
            };
            let v_opt = if proj_rank > 0 {
                Some(proj_v.ok_or(CoreError::invalid_input(
                    "missing conditional projection V factor",
                ))?)
            } else {
                None
            };
            let mut buf = vec![0.0_f64; m_obs];
            let mut fe_sum = vec![0.0_f64; n_fe];
            let mut proj_dot = vec![0.0_f64; proj_rank];

            for c in 0..q {
                if let Some(cl_id) = cluster_id {
                    for i in 0..m_obs {
                        buf[i] = wfin[cl_id[i] * width + 1 + c] * xj[i * k_use + col_j];
                    }
                } else {
                    for i in 0..m_obs {
                        buf[i] = wfin[i * width + 1 + c] * xj[i * k_use + col_j];
                    }
                }

                fe_sum.fill(0.0);
                for (i, &f) in fe.iter().enumerate() {
                    fe_sum[f] += buf[i];
                }
                for f in 0..n_fe {
                    fe_sum[f] /= counts[f] as f64;
                }
                for (i, &f) in fe.iter().enumerate() {
                    buf[i] -= fe_sum[f];
                }

                if let (Some(u), Some(v)) = (u_opt, v_opt) {
                    proj_dot.fill(0.0);
                    for i in 0..m_obs {
                        let v_row = &v[i * proj_rank..(i + 1) * proj_rank];
                        let val = buf[i];
                        for l in 0..proj_rank {
                            proj_dot[l] += v_row[l] * val;
                        }
                    }
                    for i in 0..m_obs {
                        let u_row = &u[i * proj_rank..(i + 1) * proj_rank];
                        let mut pval = 0.0_f64;
                        for l in 0..proj_rank {
                            pval += u_row[l] * proj_dot[l];
                        }
                        buf[i] -= pval;
                    }
                }

                if let Some(cl_id) = cluster_id {
                    for (i, &g) in cl_id.iter().enumerate() {
                        wx[g * width + 1 + c] += buf[i] * xjs[i * k_use + col_j];
                    }
                } else {
                    for i in 0..m_obs {
                        wx[i * width + 1 + c] = buf[i] * xjs[i * k_use + col_j];
                    }
                }
            }
        }
    } else if let Some(cl_id) = cluster_id {
        let mut diag_g = vec![0.0_f64; n];
        for (i, &g) in cl_id.iter().enumerate() {
            diag_g[g] += xj[i * k_use + col_j] * xjs[i * k_use + col_j];
        }
        for g in 0..n {
            wx[g * width] = wfin[g * width] * diag_g[g];
        }

        if q > 0 {
            let mut b_q = vec![0.0_f64; checked_len(&[n, q])?];
            for g in 0..n {
                b_q[g * q..(g + 1) * q].copy_from_slice(&wfin[g * width + 1..g * width + 1 + q]);
            }
            for c in 0..q {
                for g in 0..n {
                    wx[g * width + 1 + c] = diag_g[g] * b_q[g * q + c];
                }
                pre_max = pre_max.max(scaled_norm((0..n).map(|g| wx[g * width + 1 + c])));
            }

            if proj_rank > 0 {
                let u = proj_u.ok_or(CoreError::invalid_input(
                    "missing conditional projection U factor",
                ))?;
                let v = proj_v.ok_or(CoreError::invalid_input(
                    "missing conditional projection V factor",
                ))?;
                let mut u_cl = vec![0.0_f64; checked_len(&[n, proj_rank])?];
                let mut v_cl = vec![0.0_f64; checked_len(&[n, proj_rank])?];
                for (i, &g) in cl_id.iter().enumerate() {
                    let xj_i = xj[i * k_use + col_j];
                    let xjs_i = xjs[i * k_use + col_j];
                    let u_row = &u[i * proj_rank..(i + 1) * proj_rank];
                    let v_row = &v[i * proj_rank..(i + 1) * proj_rank];
                    let uc_row = &mut u_cl[g * proj_rank..(g + 1) * proj_rank];
                    let vc_row = &mut v_cl[g * proj_rank..(g + 1) * proj_rank];
                    for l in 0..proj_rank {
                        uc_row[l] += xjs_i * u_row[l];
                        vc_row[l] += xj_i * v_row[l];
                    }
                }
                let mut vt_bq = vec![0.0_f64; checked_len(&[proj_rank, q])?];
                let mut corr = vec![0.0_f64; checked_len(&[n, q])?];
                linalg::multiply(&v_cl, &b_q, &mut vt_bq, proj_rank, n, q, true);
                linalg::multiply(&u_cl, &vt_bq, &mut corr, n, proj_rank, q, false);
                for g in 0..n {
                    for c in 0..q {
                        wx[g * width + 1 + c] -= corr[g * q + c];
                    }
                }
            }
        }
    } else {
        for i in 0..n {
            let xj_i = xj[i * k_use + col_j];
            let xjs_i = xjs[i * k_use + col_j];
            wx[i * width] = wfin[i * width] * xj_i * xjs_i;
        }

        if q > 0 {
            // Denominator weights are S_x M_V D_x B_q.
            let mut rx = vec![0.0; checked_len(&[n, q])?];
            for c in 0..q {
                for i in 0..n {
                    let idx = i * k_use + col_j;
                    rx[i * q + c] = wfin[i * width + 1 + c] * xj[idx];
                }
                let column = (0..n).map(|i| rx[i * q + c] * xjs[i * k_use + col_j]);
                pre_max = pre_max.max(scaled_norm(column));
            }
            let rr = if proj_rank > 0 {
                let u = proj_u.ok_or(CoreError::invalid_input(
                    "missing conditional projection U factor",
                ))?;
                let v = proj_v.ok_or(CoreError::invalid_input(
                    "missing conditional projection V factor",
                ))?;
                let mut vt_rx = vec![0.0; checked_len(&[proj_rank, q])?];
                let mut u_vt_rx = vec![0.0; checked_len(&[n, q])?];
                linalg::multiply(v, &rx, &mut vt_rx, proj_rank, n, q, true);
                linalg::multiply(u, &vt_rx, &mut u_vt_rx, n, proj_rank, q, false);
                for (r, p) in rx.iter_mut().zip(u_vt_rx.iter()) {
                    *r -= *p;
                }
                rx
            } else {
                rx
            };
            for i in 0..n {
                let xjs_i = xjs[i * k_use + col_j];
                for c in 0..q {
                    wx[i * width + 1 + c] = rr[i * q + c] * xjs_i;
                }
            }
        }
    }

    if wx.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::degenerate(
            "conditional projection produced non-finite W values",
        ));
    }
    if width > 1 {
        let tol = SCALE_REL_TOL * pre_max;
        let all_collapsed = !pre_max.is_finite()
            || pre_max == 0.0
            || (1..width).all(|c| {
                let col_norm = scaled_norm((0..n).map(|i| wx[i * width + c]));
                !col_norm.is_finite() || col_norm <= tol
            });
        if all_collapsed {
            return Err(CoreError::degenerate(
                "conditional projection collapsed all non-constant W columns",
            ));
        }
    }
    Ok(wx)
}

#[allow(clippy::too_many_arguments)]
fn validate_common_inputs(
    coords: &[f64],
    n: usize,
    d: usize,
    avc: f64,
    uncond: bool,
    coef: &[f64],
    wj: &[f64],
    cluster_id: Option<&[usize]>,
    fe_id: Option<&[usize]>,
    xj: Option<&[f64]>,
    xjs: Option<&[f64]>,
    proj_u: Option<&[f64]>,
    proj_v: Option<&[f64]>,
    proj_rank: usize,
    cv_tolerance: f64,
) -> Result<f64, CoreError> {
    if n < 2 || d == 0 || coef.is_empty() {
        return Err(CoreError::invalid_input(
            "invalid SCPC sample or coefficient dimensions",
        ));
    }
    let k_use = coef.len();
    if coords.len() != checked_len(&[n, d])? || wj.len() != checked_len(&[n, k_use])? {
        return Err(CoreError::invalid_input(
            "SCPC input array dimensions do not match",
        ));
    }
    if coef.iter().any(|v| !v.is_finite()) || wj.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input("SCPC inputs must be finite"));
    }
    if !avc.is_finite() || avc <= 0.001 || avc >= 0.99 {
        return Err(CoreError::invalid_input(
            "avc must lie strictly in (0.001, 0.99)",
        ));
    }
    let tol = if cv_tolerance > 0.0 && cv_tolerance.is_finite() {
        cv_tolerance
    } else {
        DEFAULT_CV_TOL
    };
    let m_obs = if let Some(cl) = cluster_id {
        if cl.len() < n {
            return Err(CoreError::invalid_input(
                "cluster_id length must be at least the number of clusters n",
            ));
        }
        let mut counts = vec![0usize; n];
        for &g in cl {
            if g >= n {
                return Err(CoreError::invalid_input(
                    "cluster_id entry is out of bounds for n clusters",
                ));
            }
            counts[g] += 1;
        }
        if counts.contains(&0) {
            return Err(CoreError::invalid_input(
                "every cluster index in 0..n must appear in cluster_id",
            ));
        }
        cl.len()
    } else {
        n
    };
    if uncond {
        if fe_id.is_some() {
            return Err(CoreError::invalid_input(
                "unconditional SCPC does not accept absorbed fe_id",
            ));
        }
    } else {
        if let Some(fe) = fe_id {
            if fe.len() != m_obs {
                return Err(CoreError::invalid_input(
                    "fe_id length must match the number of model observations",
                ));
            }
            let max_fe = *fe.iter().max().ok_or(CoreError::invalid_input(
                "fe_id must contain at least one element",
            ))?;
            if max_fe >= m_obs {
                return Err(CoreError::invalid_input(
                    "fe_id entry is out of bounds for observation count",
                ));
            }
            let mut fe_counts = vec![0usize; max_fe + 1];
            for &f in fe {
                fe_counts[f] += 1;
            }
            if fe_counts.contains(&0) {
                return Err(CoreError::invalid_input(
                    "every fixed-effect index in 0..n_fe must appear in fe_id",
                ));
            }
        }
        let xj = xj.ok_or(CoreError::invalid_input(
            "conditional SCPC requires xj direction matrix",
        ))?;
        let xjs = xjs.ok_or(CoreError::invalid_input(
            "conditional SCPC requires xjs scaling matrix",
        ))?;
        if xj.len() != checked_len(&[m_obs, k_use])? || xjs.len() != checked_len(&[m_obs, k_use])? {
            return Err(CoreError::invalid_input(
                "SCPC conditional xj/xjs dimensions do not match",
            ));
        }
        if xj.iter().any(|v| !v.is_finite()) || xjs.iter().any(|v| !v.is_finite()) {
            return Err(CoreError::invalid_input(
                "SCPC conditional directions must be finite",
            ));
        }
        if proj_rank > 0 {
            let expected = checked_len(&[m_obs, proj_rank])?;
            let u = proj_u.ok_or(CoreError::invalid_input(
                "conditional SCPC requires projection U factor",
            ))?;
            let v = proj_v.ok_or(CoreError::invalid_input(
                "conditional SCPC requires projection V factor",
            ))?;
            if u.len() != expected || v.len() != expected {
                return Err(CoreError::invalid_input(
                    "conditional projection factor dimensions do not match",
                ));
            }
            if u.iter().any(|val| !val.is_finite()) || v.iter().any(|val| !val.is_finite()) {
                return Err(CoreError::invalid_input(
                    "conditional projection factors must be finite",
                ));
            }
        }
    }
    Ok(tol)
}

fn build_scale_grid(c0: f64, cmax: f64) -> Vec<f64> {
    let raw_nc = if c0 > 0.0 && cmax > c0 && c0.is_finite() && cmax.is_finite() {
        let ratio = cmax / c0;
        let log_ratio = if ratio.is_finite() {
            ratio.ln()
        } else {
            cmax.ln() - c0.ln()
        };
        (log_ratio / CGRID_FAC.ln()).ceil()
    } else {
        2.0
    };
    let nc = if raw_nc.is_finite() {
        raw_nc.clamp(2.0, 1_000_000.0) as usize
    } else {
        2
    };
    let mut scales = Vec::with_capacity(nc + 1);
    let mut c = c0;
    for _ in 0..nc {
        if !c.is_finite() {
            break;
        }
        scales.push(c);
        c *= CGRID_FAC;
    }
    if scales.is_empty() {
        scales.push(c0);
    }
    if cmax.is_finite() && cmax > *scales.last().expect("scale grid includes c0") {
        scales.push(cmax);
    }
    scales
}

fn trim_basis_columns(
    weights: &[f64],
    n: usize,
    full_width: usize,
    width: usize,
) -> Result<Vec<f64>, CoreError> {
    let mut wfin = vec![0.0; checked_len(&[n, width])?];
    for (dst, src) in wfin
        .chunks_exact_mut(width)
        .zip(weights.chunks_exact(full_width))
    {
        dst.copy_from_slice(&src[..width]);
    }
    Ok(wfin)
}

// SCPC weights are homogeneous within each coefficient. Work in bounded units
// before projections and restore the score column's units only on its SE.
fn normalized_columns(values: &[f64], width: usize) -> (Vec<f64>, Vec<f64>) {
    let mut scales = vec![0.0_f64; width];
    for row in values.chunks_exact(width) {
        for (scale, &value) in scales.iter_mut().zip(row) {
            *scale = scale.max(value.abs());
        }
    }
    let mut normalized = values.to_vec();
    for row in normalized.chunks_exact_mut(width) {
        for (value, &scale) in row.iter_mut().zip(&scales) {
            if scale > 0.0 {
                *value /= scale;
            }
        }
    }
    (normalized, scales)
}

fn compute_uncond_stats(
    wfin: &[f64],
    wj: &[f64],
    n: usize,
    q: usize,
    coef: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), CoreError> {
    let k_use = coef.len();
    let width = q + 1;
    let mut projected = vec![0.0; checked_len(&[width, k_use])?];
    let (mut wj, scales) = normalized_columns(wj, k_use);
    let mut means = vec![0.0; k_use];
    for row in wj.chunks_exact(k_use) {
        for (mean, &value) in means.iter_mut().zip(row) {
            *mean += value;
        }
    }
    for mean in &mut means {
        *mean /= n as f64;
    }
    for row in wj.chunks_exact_mut(k_use) {
        for (value, mean) in row.iter_mut().zip(&means) {
            *value -= mean;
        }
    }
    linalg::multiply(wfin, &wj, &mut projected, width, n, k_use, true);
    let q_root = (q as f64).sqrt();
    let se_denom = q_root * (n as f64).sqrt();

    let mut se_vec = Vec::with_capacity(k_use);
    let mut tau_vec = Vec::with_capacity(k_use);
    for j in 0..k_use {
        let wj_norm = scaled_norm((0..n).map(|i| wj[i * k_use + j]));
        let norm = scaled_norm((1..width).map(|row| projected[row * k_use + j]));
        let se = (norm / se_denom) * scales[j];
        if !se.is_finite() || se == 0.0 || norm <= SCALE_REL_TOL * wj_norm {
            return Err(CoreError::degenerate(
                "SCPC variance projection is degenerate",
            ));
        }
        tau_vec.push(coef[j] / se);
        se_vec.push(se);
    }
    Ok((se_vec, tau_vec))
}

#[allow(clippy::too_many_arguments)]
fn finalize_scpc_output(
    k_use: usize,
    q: usize,
    width: usize,
    oms_full: &[f64],
    full_width: usize,
    cvfin: f64,
    coef: &[f64],
    se_vec: &[f64],
    tau_vec: &[f64],
    omsx_list: Option<&[Vec<f64>]>,
    cvs: bool,
    cv_tolerance: f64,
    wfin: Vec<f64>,
    c0: f64,
    cmax: f64,
    tied: bool,
) -> Result<ScpcOutput, CoreError> {
    let q_root = (q as f64).sqrt();
    let p_u = map_indices(k_use, |j| {
        let mut ws = ProbabilityWorkspace::new(q);
        ws.maximum(oms_full, full_width, tau_vec[j].abs() / q_root, None)
            .map(|(p, _)| p)
    })?;

    let cvs_uncond: Option<[f64; 4]> = if cvs {
        let vals = map_indices(EXTRA_CV_LEVELS.len(), |idx| {
            let mut ws = ProbabilityWorkspace::new(q);
            ws.critical(oms_full, full_width, q, EXTRA_CV_LEVELS[idx], cv_tolerance)
        })?;
        Some([vals[0], vals[1], cvfin, vals[2]])
    } else {
        None
    };

    let mut scpcstats = vec![0.0; checked_len(&[k_use, 6])?];
    let mut scpccvs = if cvs {
        Some(vec![0.0; checked_len(&[k_use, 4])?])
    } else {
        None
    };

    if let Some(omsx_list) = omsx_list {
        let cond_results = map_indices(k_use, |j| {
            let omsx = &omsx_list[j];
            let mut ws = ProbabilityWorkspace::new(q);
            let (p_c, _) = ws.maximum(omsx, width, tau_vec[j].abs() / q_root, None)?;
            let cvx = ws.critical(omsx, width, q, 0.05, cv_tolerance)?;
            let cond_cvs = if cvs {
                let extra = map_indices(EXTRA_CV_LEVELS.len(), |idx| {
                    let mut inner = ProbabilityWorkspace::new(q);
                    inner.critical(omsx, width, q, EXTRA_CV_LEVELS[idx], cv_tolerance)
                })?;
                Some([extra[0], extra[1], cvx, extra[2]])
            } else {
                None
            };
            Ok((p_c, cvx, cond_cvs))
        })?;

        for j in 0..k_use {
            let (p_c, cvx, cond_cvs) = cond_results[j];
            let p_final = p_u[j].max(p_c);
            let cv = cvfin.max(cvx);
            let row = &mut scpcstats[j * 6..(j + 1) * 6];
            row[0] = coef[j];
            row[1] = se_vec[j];
            row[2] = tau_vec[j];
            row[3] = p_final;
            row[4] = coef[j] - cv * se_vec[j];
            row[5] = coef[j] + cv * se_vec[j];

            if let (Some(out_cvs), Some(uncond_cvs), Some(c_cvs)) =
                (scpccvs.as_mut(), cvs_uncond.as_ref(), cond_cvs.as_ref())
            {
                let row = &mut out_cvs[j * 4..(j + 1) * 4];
                for idx in 0..4 {
                    row[idx] = uncond_cvs[idx].max(c_cvs[idx]);
                }
            }
        }
    } else {
        for j in 0..k_use {
            let row = &mut scpcstats[j * 6..(j + 1) * 6];
            row[0] = coef[j];
            row[1] = se_vec[j];
            row[2] = tau_vec[j];
            row[3] = p_u[j];
            row[4] = coef[j] - cvfin * se_vec[j];
            row[5] = coef[j] + cvfin * se_vec[j];
        }
        if let (Some(out_cvs), Some(uncond_cvs)) = (scpccvs.as_mut(), cvs_uncond.as_ref()) {
            for row in out_cvs.chunks_exact_mut(4) {
                row.copy_from_slice(uncond_cvs);
            }
        }
    }

    Ok(ScpcOutput {
        scpcstats,
        scpccvs,
        wfin,
        c0,
        cmax,
        cvfin,
        q,
        tied,
    })
}

/// Run the full exact SCPC inference pipeline in one call.
pub fn run_exact_scpc(req: &ExactScpcRequest<'_>) -> Result<ScpcOutput, CoreError> {
    let cv_tolerance = validate_common_inputs(
        req.coords,
        req.n,
        req.d,
        req.avc,
        req.uncond,
        req.coef,
        req.wj,
        req.cluster_id,
        req.fe_id,
        req.xj,
        req.xjs,
        req.proj_u,
        req.proj_v,
        req.proj_rank,
        req.cv_tolerance,
    )?;
    let n = req.n;
    let k_use = req.coef.len();

    run_with_threads(req.nthreads, || {
        let distmat = distance_matrix(req.coords, n, req.d, req.metric)?;
        let (c0, cmax) = {
            let dist_vec = lower_triangle(&distmat, n);
            let floor = correlation_floor(&dist_vec)?;
            let min_avc_eff = floor + (1.0 - floor) * MIN_AVC;
            let c0 = solve_c0_from_avc(&dist_vec, req.avc)?;
            let cmax = solve_c0_from_avc(&dist_vec, min_avc_eff)?;
            (c0, cmax)
        };
        let scales = build_scale_grid(c0, cmax);

        let mut qmax = initial_qmax(req.avc);
        let (basis, oms_full, q, cvfin): (ScpcBasis, Vec<f64>, usize, f64) = loop {
            qmax = qmax.min(n - 1);
            let basis = scpc_basis(&distmat, n, c0, qmax)?;
            let mut grids = covariance_grids(&distmat, n, &[&basis.weights], &[qmax + 1], &scales)?;
            let oms = grids.remove(0);
            let (q_opt, cv_opt) = select_final_q(&oms, qmax + 1, qmax, cv_tolerance)?;
            if q_opt < qmax || qmax == n - 1 {
                break (basis, oms, q_opt, cv_opt);
            }
            qmax = ((qmax as f64) + (qmax as f64) * 0.5).round() as usize;
        };

        let full_width = qmax + 1;
        let width = q + 1;
        let wfin = trim_basis_columns(&basis.weights, n, full_width, width)?;
        let (se_vec, tau_vec) = compute_uncond_stats(&wfin, req.wj, n, q, req.coef)?;

        // Precompute conditional Oms before dropping the n x n distance matrix.
        let omsx_list = if req.uncond {
            drop(distmat);
            None
        } else {
            let (xj, _) = normalized_columns(req.xj.expect("validated conditional xj"), k_use);
            let xjs = req.xjs.expect("validated conditional xjs");
            let wx_list = map_indices(k_use, |j| {
                build_conditional_basis(
                    &wfin,
                    n,
                    q,
                    req.cluster_id,
                    req.fe_id,
                    &xj,
                    xjs,
                    k_use,
                    j,
                    req.proj_u,
                    req.proj_v,
                    req.proj_rank,
                )
                .map(normalize_covariance_basis)
            })?;
            let wx_refs: Vec<&[f64]> = wx_list.iter().map(Vec::as_slice).collect();
            let widths = vec![width; k_use];
            let grids = covariance_grids(&distmat, n, &wx_refs, &widths, &scales)?;
            drop(distmat);
            Some(grids)
        };

        finalize_scpc_output(
            k_use,
            q,
            width,
            &oms_full,
            full_width,
            cvfin,
            req.coef,
            &se_vec,
            &tau_vec,
            omsx_list.as_deref(),
            req.cvs,
            cv_tolerance,
            wfin,
            c0,
            cmax,
            basis.tied,
        )
    })?
}

/// Run the full approximate (Nyström) SCPC inference pipeline in one call.
pub fn run_approx_scpc(req: &ApproxScpcRequest<'_>) -> Result<ScpcOutput, CoreError> {
    let cv_tolerance = validate_common_inputs(
        req.coords,
        req.n,
        req.d,
        req.avc,
        req.uncond,
        req.coef,
        req.wj,
        req.cluster_id,
        req.fe_id,
        req.xj,
        req.xjs,
        req.proj_u,
        req.proj_v,
        req.proj_rank,
        req.cv_tolerance,
    )?;
    if req.landmarks < 2 || req.samples == 0 || req.pair_count == 0 {
        return Err(CoreError::invalid_input(
            "invalid Nystrom landmark or sampling parameters",
        ));
    }
    let n = req.n;
    let k_use = req.coef.len();
    let m = req.landmarks.min(n);
    let q_limit = m - 1;

    run_with_threads(req.nthreads, || {
        let mut qmax = initial_qmax(req.avc);
        let (lb, scales, oms_full, state_after_uncond, q, cvfin) = loop {
            qmax = qmax.min(q_limit);
            let lb = landmark_basis(
                req.coords,
                n,
                req.d,
                req.metric,
                &LandmarkRequest {
                    avc: req.avc,
                    min_avc: MIN_AVC,
                    components: qmax,
                    landmarks: m,
                    samples: req.samples,
                    seed: req.seed,
                },
            )?;
            let scales = build_scale_grid(lb.decay[0], lb.decay[1]);
            let full_width = qmax + 1;
            let (oms, state_after) = sampled_covariance_grid(
                req.coords,
                n,
                req.d,
                req.metric,
                &lb.weights,
                full_width,
                &scales,
                req.pair_count,
                lb.state,
            )?;
            if !is_psd_prefix_grid(&oms, full_width, qmax) {
                return Err(CoreError::invalid_sampled_covariance(
                    "sampled SCPC covariance grid is not positive semidefinite",
                ));
            }
            let (q_opt, cv_opt) = select_final_q(&oms, full_width, qmax, cv_tolerance)?;
            if q_opt < qmax || qmax >= q_limit {
                break (lb, scales, oms, state_after, q_opt, cv_opt);
            }
            qmax = ((qmax as f64) + (qmax as f64) * 0.5).round() as usize;
        };

        let [c0, cmax] = lb.decay;
        let full_width = qmax + 1;
        let width = q + 1;
        let wfin = trim_basis_columns(&lb.weights, n, full_width, width)?;
        let (se_vec, tau_vec) = compute_uncond_stats(&wfin, req.wj, n, q, req.coef)?;

        let omsx_list = if req.uncond {
            None
        } else {
            let (xj, _) = normalized_columns(req.xj.expect("validated conditional xj"), k_use);
            let xjs = req.xjs.expect("validated conditional xjs");
            // All reported coefficients use the same sample, starting immediately
            // after the unconditional grid. Their order or selection cannot change it.
            let pairs = SampledPairs::new(
                req.coords,
                n,
                req.d,
                req.metric,
                req.pair_count,
                state_after_uncond,
            )?;
            let grids = map_indices(k_use, |j| {
                let wx = build_conditional_basis(
                    &wfin,
                    n,
                    q,
                    req.cluster_id,
                    req.fe_id,
                    &xj,
                    xjs,
                    k_use,
                    j,
                    req.proj_u,
                    req.proj_v,
                    req.proj_rank,
                )?;
                let wx = normalize_covariance_basis(wx);
                let omsx = pairs.covariance_grid(&wx, width, &scales)?;
                if !is_psd_prefix_grid(&omsx, width, q) {
                    return Err(CoreError::invalid_sampled_covariance(
                        "sampled conditional SCPC covariance grid is not positive semidefinite",
                    ));
                }
                Ok(omsx)
            })?;
            Some(grids)
        };

        finalize_scpc_output(
            k_use,
            q,
            width,
            &oms_full,
            full_width,
            cvfin,
            req.coef,
            &se_vec,
            &tau_vec,
            omsx_list.as_deref(),
            req.cvs,
            cv_tolerance,
            wfin,
            c0,
            cmax,
            false,
        )
    })?
}

#[cfg(test)]
mod tests;
