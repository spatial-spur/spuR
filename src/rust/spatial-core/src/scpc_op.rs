//! Whole-operation exact SCPC inference pipeline.

use crate::basis::{scpc_basis, ScpcBasis};
use crate::calibration::{correlation_floor, solve_c0_from_avc};
use crate::covariance;
use crate::distance::{distance_matrix, lower_triangle, Metric};
use crate::linalg::{self, EigenWorkspace};
use crate::nystrom::{
    advance_pair_chain_state, is_psd_prefix_grid, landmark_basis, sampled_covariance_grid,
    LandmarkRequest,
};
use crate::quadrature::RULE;
use crate::student_t::{ln_gamma, student_t_ppf};
use crate::threads::{is_serial, run_with_threads};
use crate::{checked_len, CoreError};
use rayon::prelude::*;
use std::f64::consts::{FRAC_2_PI, FRAC_PI_2};
use std::sync::OnceLock;

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
    /// Coefficient influence matrix (`n x k_use`), column `j` is `n * (B^{-1}[j,.] @ S') + coef[j]`.
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

/// Complete output bundle returned by [`run_exact_scpc`].
#[derive(Debug, Clone)]
pub struct ExactScpcOutput {
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

fn log_one_plus_exp(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}

fn angular_product(cos2: f64, ratios: &[f64], logs: &[f64]) -> f64 {
    if cos2 <= 0.0 {
        return 0.0;
    }
    let product: f64 = ratios.iter().map(|value| cos2 / (cos2 + value)).product();
    if product >= f64::MIN_POSITIVE {
        return product.sqrt();
    }
    let log_cos = cos2.ln();
    (-0.5
        * logs
            .iter()
            .map(|ratio| log_one_plus_exp(ratio - log_cos))
            .sum::<f64>())
    .exp()
}

fn complement_integrand(sin2: f64, denom: f64, ratios: &[f64], scale: f64) -> f64 {
    let mut product = 1.0;
    let mut deficit = 0.0;
    for &ratio in ratios {
        let extra = (ratio / scale) * denom;
        let inv = 1.0 / (sin2 + extra);
        deficit += product * extra * inv;
        product *= sin2 * inv;
    }
    deficit / ((1.0 + product.sqrt()) * denom)
}

fn panel(f: &impl Fn(f64, f64) -> f64, left: f64, right: f64) -> f64 {
    const HALF: f64 = FRAC_PI_2 * 0.5;
    type Rule = [(f64, f64, f64); RULE.len()];
    static COMMON: OnceLock<[Rule; 3]> = OnceLock::new();
    let index = match (left, right) {
        (0.0, FRAC_PI_2) => Some(0),
        (0.0, HALF) => Some(1),
        (HALF, FRAC_PI_2) => Some(2),
        _ => None,
    };
    if let Some(index) = index {
        let rules = COMMON.get_or_init(|| {
            [(0.0, FRAC_PI_2), (0.0, HALF), (HALF, FRAC_PI_2)].map(|(lo, hi)| {
                RULE.map(|(node, weight)| {
                    let (sin, cos) = (lo + (hi - lo) * node).sin_cos();
                    (sin * sin, cos * cos, (hi - lo) * weight)
                })
            })
        });
        return rules[index]
            .iter()
            .map(|&(sin2, cos2, weight)| weight * f(sin2, cos2))
            .sum();
    }
    let width = right - left;
    width
        * RULE
            .iter()
            .map(|&(node, weight)| {
                let (sin, cos) = (left + width * node).sin_cos();
                weight * f(sin * sin, cos * cos)
            })
            .sum::<f64>()
}

fn adaptive(
    f: &impl Fn(f64, f64) -> f64,
    lo: f64,
    hi: f64,
    whole: f64,
    tol: f64,
    depth: u8,
) -> f64 {
    let mid = 0.5 * (lo + hi);
    let left = panel(f, lo, mid);
    let right = panel(f, mid, hi);
    if depth == 20 || (left + right - whole).abs() <= tol.max(1e-11 * whole.abs()) {
        return left + right;
    }
    adaptive(f, lo, mid, left, tol * 0.5, depth + 1)
        + adaptive(f, mid, hi, right, tol * 0.5, depth + 1)
}

fn integrate(f: &impl Fn(f64, f64) -> f64) -> f64 {
    let whole = panel(f, 0.0, FRAC_PI_2);
    let tolerance = if whole.abs() < 1e-8 {
        1e-10 * whole.abs()
    } else {
        1e-11
    };
    adaptive(f, 0.0, FRAC_PI_2, whole, tolerance, 0) * FRAC_2_PI
}

struct ProbabilityWorkspace {
    eigensolver: EigenWorkspace,
    ratios: Vec<f64>,
    logs: Vec<f64>,
}

impl ProbabilityWorkspace {
    fn new(q: usize) -> Self {
        Self {
            eigensolver: EigenWorkspace::new(q + 1),
            ratios: Vec::with_capacity(q),
            logs: Vec::with_capacity(q),
        }
    }

    /// Evaluate tail-stable adaptive angular/complement rejection probability.
    fn evaluate(&mut self, omega: &[f64], width: usize, cv: f64) -> Result<f64, CoreError> {
        if !cv.is_finite() {
            return Err(CoreError::invalid_input("SCPC critical value must be finite"));
        }
        if cv <= 0.0 {
            return Ok(1.0);
        }
        let (positive, negative) = if cv > 1.0 {
            (cv.recip().powi(2), -1.0)
        } else {
            (1.0, -cv.powi(2))
        };
        let eigenvalues = self
            .eigensolver
            .weighted_eigenvalues(omega, width, positive, negative)?;
        let max_positive = eigenvalues.iter().copied().fold(0.0, f64::max);
        if !max_positive.is_finite() || max_positive <= 0.0 {
            return Ok(0.0);
        }
        self.ratios.clear();
        self.logs.clear();
        let max_pos_ln = max_positive.ln();
        for &value in eigenvalues.iter().filter(|&&value| value < 0.0) {
            self.ratios.push(-value / max_positive);
            self.logs.push((-value).ln() - max_pos_ln);
        }
        if self.ratios.is_empty() {
            return Ok(1.0);
        }
        if self.ratios.len() == 1 {
            return Ok(FRAC_2_PI * (-0.5 * self.logs[0]).exp().atan());
        }

        let scale: f64 = self.ratios.iter().sum();
        if scale == 0.0 {
            return Ok(1.0);
        }
        let probability = if scale < 0.0625 {
            let deficit = integrate(&|sin2, cos2| {
                complement_integrand(sin2, cos2 + scale * sin2, &self.ratios, scale)
            });
            1.0 - scale.sqrt() * deficit
        } else {
            integrate(&|_, cos2| angular_product(cos2, &self.ratios, &self.logs))
        };
        Ok(probability.clamp(0.0, 1.0))
    }

    fn maximum(
        &mut self,
        omegas: &[f64],
        width: usize,
        cv: f64,
        known: Option<(f64, usize)>,
    ) -> Result<(f64, usize), CoreError> {
        let mut best = (-1.0, 0);
        for (i, omega) in omegas.chunks_exact(width * width).enumerate() {
            let p = if let Some((p, known_i)) = known.filter(|&(_, index)| index == i) {
                debug_assert_eq!(i, known_i);
                p
            } else {
                self.evaluate(omega, width, cv)?
            };
            if p > best.0 {
                best = (p, i);
            }
        }
        Ok(best)
    }

    /// Released `get_cv(Oms, q, level)` solver: steps out by `1/sqrt(q)` and
    /// bisects to `tol / sqrt(q)`, returning the upper bracket endpoint in `t` units.
    fn critical(
        &mut self,
        omegas: &[f64],
        width: usize,
        q: usize,
        level: f64,
        tolerance: f64,
    ) -> Result<f64, CoreError> {
        let root_q = (q as f64).sqrt();
        let step = 1.0 / root_q;
        let tol = tolerance * step;
        let matrix_size = width * width;
        let grid_count = omegas.len() / matrix_size;

        let mut rp = 1.0;
        let mut active = 0usize;
        let mut cv0 = student_t_ppf(1.0 - 0.5 * level, q)? / root_q;
        let mut cv1 = cv0;

        for _ in 0..=grid_count + 4 {
            if rp <= level {
                break;
            }
            let omega = &omegas[active * matrix_size..(active + 1) * matrix_size];
            cv1 = cv0;
            for iter in 0..1024 {
                if self.evaluate(omega, width, cv1)? > level {
                    cv0 = cv1;
                    cv1 += if iter < 32 { step } else { cv1 };
                } else {
                    break;
                }
            }
            if self.evaluate(omega, width, cv1)? > level {
                return Err(CoreError::calibration(
                    "SCPC critical value could not be bracketed",
                ));
            }

            let mut p_upper = self.evaluate(omega, width, cv1)?;
            while cv1 - cv0 > tol {
                let mid = 0.5 * (cv0 + cv1);
                if mid == cv0 || mid == cv1 {
                    break;
                }
                let p_mid = self.evaluate(omega, width, mid)?;
                if p_mid > level {
                    cv0 = mid;
                } else {
                    cv1 = mid;
                    p_upper = p_mid;
                }
            }

            let (max_p, max_i) = self.maximum(omegas, width, cv1, Some((p_upper, active)))?;
            if max_i == active {
                break;
            }
            active = max_i;
            rp = max_p;
            cv0 = cv1;
        }
        Ok(cv1 * root_q)
    }
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
                let mut col_sq = 0.0_f64;
                for g in 0..n {
                    let val = diag_g[g] * wfin[g * width + 1 + c];
                    col_sq += val * val;
                }
                pre_max = pre_max.max(col_sq.sqrt());
            }
        } else {
            for i in 0..n {
                let xj_i = xj[i * k_use + col_j];
                let xjs_i = xjs[i * k_use + col_j];
                wx[i * width] = wfin[i * width] * xj_i * xjs_i;
            }
            for c in 0..q {
                let mut col_sq = 0.0_f64;
                for i in 0..n {
                    let idx = i * k_use + col_j;
                    let b_ic = wfin[i * width + 1 + c];
                    let scale_i = xj[idx] * xjs[idx];
                    col_sq += (b_ic * scale_i).powi(2);
                }
                pre_max = pre_max.max(col_sq.sqrt());
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
                let mut col_sq = 0.0_f64;
                for g in 0..n {
                    let val = diag_g[g] * b_q[g * q + c];
                    wx[g * width + 1 + c] = val;
                    col_sq += val * val;
                }
                pre_max = pre_max.max(col_sq.sqrt());
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
                let mut col_sq = 0.0_f64;
                for i in 0..n {
                    let idx = i * k_use + col_j;
                    let b_ic = wfin[i * width + 1 + c];
                    let scale_i = xj[idx] * xjs[idx];
                    rx[i * q + c] = b_ic * xj[idx];
                    col_sq += (b_ic * scale_i).powi(2);
                }
                pre_max = pre_max.max(col_sq.sqrt());
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
                let col_norm = (0..n)
                    .map(|i| wx[i * width + c].powi(2))
                    .sum::<f64>()
                    .sqrt();
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
    let nc = ((cmax / c0).ln() / CGRID_FAC.ln()).ceil().max(2.0) as usize;
    let mut scales = Vec::with_capacity(nc);
    let mut c = c0;
    for _ in 0..nc {
        scales.push(c);
        c *= CGRID_FAC;
    }
    scales
}

fn trim_basis_columns(weights: &[f64], n: usize, full_width: usize, width: usize) -> Result<Vec<f64>, CoreError> {
    let mut wfin = vec![0.0; checked_len(&[n, width])?];
    for (dst, src) in wfin
        .chunks_exact_mut(width)
        .zip(weights.chunks_exact(full_width))
    {
        dst.copy_from_slice(&src[..width]);
    }
    Ok(wfin)
}

fn compute_uncond_stats(
    wfin: &[f64],
    wj: &[f64],
    n: usize,
    q: usize,
    k_use: usize,
) -> Result<(Vec<f64>, Vec<f64>), CoreError> {
    let width = q + 1;
    let mut projected = vec![0.0; checked_len(&[width, k_use])?];
    linalg::multiply(wfin, wj, &mut projected, width, n, k_use, true);
    let q_root = (q as f64).sqrt();
    let se_denom = q_root * (n as f64).sqrt();

    let mut se_vec = Vec::with_capacity(k_use);
    let mut tau_vec = Vec::with_capacity(k_use);
    for j in 0..k_use {
        let wj_norm = (0..n)
            .map(|i| wj[i * k_use + j].powi(2))
            .sum::<f64>()
            .sqrt();
        let norm = (1..width)
            .map(|row| projected[row * k_use + j].powi(2))
            .sum::<f64>()
            .sqrt();
        if !norm.is_finite()
            || !wj_norm.is_finite()
            || wj_norm == 0.0
            || norm <= SCALE_REL_TOL * wj_norm
        {
            return Err(CoreError::degenerate(
                "SCPC variance projection is degenerate",
            ));
        }
        tau_vec.push(q_root * projected[j] / norm);
        se_vec.push(norm / se_denom);
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
) -> Result<ExactScpcOutput, CoreError> {
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

            if let (Some(out_cvs), Some(uncond_cvs), Some(c_cvs)) = (
                scpccvs.as_mut(),
                cvs_uncond.as_ref(),
                cond_cvs.as_ref(),
            ) {
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

    Ok(ExactScpcOutput {
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
pub fn run_exact_scpc(req: &ExactScpcRequest<'_>) -> Result<ExactScpcOutput, CoreError> {
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
        let (se_vec, tau_vec) = compute_uncond_stats(&wfin, req.wj, n, q, k_use)?;

        // Precompute conditional Oms before dropping the n x n distance matrix.
        let omsx_list = if req.uncond {
            drop(distmat);
            None
        } else {
            let xj = req.xj.expect("validated conditional xj");
            let xjs = req.xjs.expect("validated conditional xjs");
            let wx_list = map_indices(k_use, |j| {
                build_conditional_basis(
                    &wfin,
                    n,
                    q,
                    req.cluster_id,
                    req.fe_id,
                    xj,
                    xjs,
                    k_use,
                    j,
                    req.proj_u,
                    req.proj_v,
                    req.proj_rank,
                )
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
pub fn run_approx_scpc(req: &ApproxScpcRequest<'_>) -> Result<ExactScpcOutput, CoreError> {
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
                return Err(CoreError::tied(
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
        let (se_vec, tau_vec) = compute_uncond_stats(&wfin, req.wj, n, q, k_use)?;

        let omsx_list = if req.uncond {
            None
        } else {
            let xj = req.xj.expect("validated conditional xj");
            let xjs = req.xjs.expect("validated conditional xjs");
            // Advance LCG in coefficient order before parallel conditional evaluation.
            let mut seeds = Vec::with_capacity(k_use);
            let mut state = state_after_uncond;
            for _ in 0..k_use {
                seeds.push(state);
                state = advance_pair_chain_state(req.pair_count, state);
            }
            let grids = map_indices(k_use, |j| {
                let wx = build_conditional_basis(
                    &wfin,
                    n,
                    q,
                    req.cluster_id,
                    req.fe_id,
                    xj,
                    xjs,
                    k_use,
                    j,
                    req.proj_u,
                    req.proj_v,
                    req.proj_rank,
                )?;
                let (omsx, _) = sampled_covariance_grid(
                    req.coords,
                    n,
                    req.d,
                    req.metric,
                    &wx,
                    width,
                    &scales,
                    req.pair_count,
                    seeds[j],
                )?;
                if !is_psd_prefix_grid(&omsx, width, q) {
                    return Err(CoreError::tied(
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
mod tests {
    use super::*;
    use crate::student_t::betainc;
    use crate::ErrorKind;

    fn solve_linear(a: &mut [f64], b: &mut [f64], p: usize) {
        for col in 0..p {
            let mut pivot = col;
            let mut max_val = a[col * p + col].abs();
            for row in (col + 1)..p {
                if a[row * p + col].abs() > max_val {
                    max_val = a[row * p + col].abs();
                    pivot = row;
                }
            }
            if pivot != col {
                for k in 0..p {
                    a.swap(col * p + k, pivot * p + k);
                }
                b.swap(col, pivot);
            }
            let diag = a[col * p + col];
            for k in col..p {
                a[col * p + k] /= diag;
            }
            b[col] /= diag;
            for row in 0..p {
                if row == col {
                    continue;
                }
                let factor = a[row * p + col];
                for k in col..p {
                    a[row * p + k] -= factor * a[col * p + k];
                }
                b[row] -= factor * b[col];
            }
        }
    }

    fn invert_symmetric(mat: &[f64], p: usize) -> Vec<f64> {
        let mut inv = vec![0.0; p * p];
        for col in 0..p {
            let mut a = mat.to_vec();
            let mut rhs = vec![0.0; p];
            rhs[col] = 1.0;
            solve_linear(&mut a, &mut rhs, p);
            for row in 0..p {
                inv[row * p + col] = rhs[row];
            }
        }
        inv
    }

    fn modified_gram_schmidt(x: &[f64], n: usize, p: usize) -> Vec<f64> {
        let mut q = vec![0.0; n * p];
        let mut rank = 0usize;
        for col in 0..p {
            let mut v: Vec<f64> = (0..n).map(|i| x[i * p + col]).collect();
            let orig_norm = v.iter().map(|val| val * val).sum::<f64>().sqrt();
            if orig_norm == 0.0 {
                continue;
            }
            for val in &mut v {
                *val /= orig_norm;
            }
            for _ in 0..2 {
                for r in 0..rank {
                    let dot: f64 = (0..n).map(|i| q[i * p + r] * v[i]).sum();
                    for i in 0..n {
                        v[i] -= dot * q[i * p + r];
                    }
                }
            }
            let norm = v.iter().map(|val| val * val).sum::<f64>().sqrt();
            if norm > 1e-12 {
                for i in 0..n {
                    q[i * p + rank] = v[i] / norm;
                }
                rank += 1;
            }
        }
        if rank < p {
            let mut compact = vec![0.0; n * rank];
            for i in 0..n {
                compact[i * rank..(i + 1) * rank].copy_from_slice(&q[i * p..i * p + rank]);
            }
            compact
        } else {
            q
        }
    }

    fn load_bench_csv(name: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let path = format!(
            "{}/../../../verification/data/{}",
            env!("CARGO_MANIFEST_DIR"),
            name
        );
        let text = std::fs::read_to_string(path).expect("fixture csv");
        let mut y = Vec::new();
        let mut x = Vec::new();
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
            x.extend_from_slice(&[1.0, cols[2], cols[3]]);
            coords.extend_from_slice(&[cols[6], cols[7]]);
        }
        (y, x, coords)
    }

    fn prepare_bench_ols(y: &[f64], x: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, usize) {
        let n = y.len();
        let p = 3usize;
        let mut xtx = vec![0.0; p * p];
        let mut xty = vec![0.0; p];
        for i in 0..n {
            for r in 0..p {
                xty[r] += x[i * p + r] * y[i];
                for c in 0..p {
                    xtx[r * p + c] += x[i * p + r] * x[i * p + c];
                }
            }
        }
        let inv_xtx = invert_symmetric(&xtx, p);
        let mut coef = xty.clone();
        let mut a = xtx.clone();
        solve_linear(&mut a, &mut coef, p);

        let mut wj = vec![0.0; n * p];
        let mut xj = vec![0.0; n * p];
        let mut xjs = vec![0.0; n * p];
        for i in 0..n {
            let yhat: f64 = (0..p).map(|k| x[i * p + k] * coef[k]).sum();
            let resid = y[i] - yhat;
            for j in 0..p {
                let x_dir: f64 =
                    (0..p).map(|k| (n as f64) * inv_xtx[j * p + k] * x[i * p + k]).sum();
                xj[i * p + j] = x_dir;
                xjs[i * p + j] = if x_dir > 0.0 {
                    1.0
                } else if x_dir < 0.0 {
                    -1.0
                } else {
                    0.0
                };
                wj[i * p + j] = x_dir * resid + coef[j];
            }
        }
        let q_mat = modified_gram_schmidt(x, n, p);
        let proj_rank = q_mat.len() / n;
        (coef, wj, xj, xjs, q_mat, proj_rank)
    }

    #[test]
    fn conditional_basis_weights_apply_regressor_before_projection() {
        let x = [-2.0, -1.0, 1.0, 2.0];
        let sign = [-1.0, -1.0, 1.0, 1.0];
        let norm = 10.0_f64.sqrt();
        let w: Vec<_> = x.iter().flat_map(|value| [0.5, -value / norm]).collect();
        let projection: Vec<_> = x.iter().flat_map(|value| [0.5, value / norm]).collect();
        let weights = build_conditional_basis(
            &w,
            4,
            1,
            None,
            None,
            &x,
            &sign,
            1,
            0,
            Some(&projection),
            Some(&projection),
            2,
        )
        .unwrap();
        let a = 1.5 / norm;
        let expected = [1.0, a, 0.5, -a, 0.5, a, 1.0, -a];
        for (got, expected) in weights.iter().zip(expected) {
            assert!((got - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn identity_covariance_probability_matches_student_t_and_near_zero_tail() {
        for q in [1usize, 2, 5, 10] {
            let width = q + 1;
            let mut identity = vec![0.0; width * width];
            for i in 0..width {
                identity[i * width + i] = 1.0;
            }
            let mut ws = ProbabilityWorkspace::new(q);
            for t in [1e-6, 0.1, 1.0, 2.5, 10.0] {
                let cv_norm = t / (q as f64).sqrt();
                let got = ws.evaluate(&identity, width, cv_norm).unwrap();
                let expected = if t * t < q as f64 {
                    1.0 - betainc(0.5, 0.5 * q as f64, (t * t) / (q as f64 + t * t))
                } else {
                    betainc(0.5 * q as f64, 0.5, (q as f64) / (q as f64 + t * t))
                };
                assert!(
                    (got - expected).abs() < 1e-10,
                    "q={q} t={t}: got={got}, expected={expected}"
                );
            }
        }

        // Near-zero test statistic reproducer from rust_core_handoff.md:
        // coords=[[0.0],[1.0]], coef=[1e-6], wj=[[1+1e-6],[-1+1e-6]], avc=0.03, uncond=true
        let coords = [0.0, 1.0];
        let coef = [1e-6];
        let wj = [1.0 + 1e-6, -1.0 + 1e-6];
        let out = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n: 2,
            d: 1,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: true,
            cvs: false,
            coef: &coef,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
            proj_rank: 0,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        assert_eq!(out.q, 1);
        let rho = 0.03_f64;
        let expected_p = 1.0 - FRAC_2_PI * (1e-6 * ((1.0 - rho) / (1.0 + rho)).sqrt()).atan();
        assert!(
            (out.scpcstats[3] - expected_p).abs() < 1e-12,
            "got p={}, expected={expected_p}",
            out.scpcstats[3]
        );
    }

    #[test]
    fn exact_scpc_matches_release_score_statistics_and_is_scale_invariant() {
        let (y, x, coords) = load_bench_csv("bench_200.csv");
        let n = y.len();
        let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);

        let out_euclid = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();

        assert_eq!(out_euclid.q, 19);
        assert!(
            (out_euclid.c0 - 3.500409719601247).abs() < 1e-3,
            "c0={}",
            out_euclid.c0
        );
        let dist_euclid = distance_matrix(&coords, n, 2, Metric::Euclidean).unwrap();
        let dvec_euclid = lower_triangle(&dist_euclid, n);
        assert!(
            (crate::calibration::avg_correlation(out_euclid.c0, &dvec_euclid) - 0.03).abs() < 1e-12
        );
        assert!(
            (out_euclid.cvfin - 2.1008365544083087).abs() < 1e-3,
            "cvfin={}",
            out_euclid.cvfin
        );
        let expected_stats = [
            [1.4834167407175964, 0.11454779987636342, 12.95019845268713],
            [1.5473280763431876, 0.09650902048915827, 16.032989128897157],
            [-0.6589528089605742, 0.10738682321940547, -6.13625386435211],
        ];
        for (row, expected) in out_euclid.scpcstats.chunks_exact(6).zip(expected_stats) {
            for (got, exp) in row.iter().zip(expected) {
                assert!((got - exp).abs() <= 1e-3 * exp.abs().max(1.0), "got={got}, exp={exp}");
            }
        }
        // Check 5% CV reuse consistency between scpcstats CI and scpccvs[, 2]
        let cvs = out_euclid.scpccvs.as_ref().unwrap();
        for j in 0..3 {
            let se = out_euclid.scpcstats[j * 6 + 1];
            let half_w = out_euclid.scpcstats[j * 6 + 5] - out_euclid.scpcstats[j * 6];
            assert!((half_w - cvs[j * 4 + 2] * se).abs() < 1e-12);
        }

        // Single-threaded serial execution on caller thread matches multi-threaded exactly
        let out_one = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: false,
            coef: &coef,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        assert_eq!(out_one.q, out_euclid.q);
        for (a, b) in out_one.scpcstats.iter().zip(out_euclid.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        // Positive outcome rescaling (`y * 1e-12`) preserves q, t, p and scales SE, CI proportionally
        let y_scale = 1e-12;
        let y_scaled: Vec<f64> = y.iter().map(|&v| v * y_scale).collect();
        let (coef_ys, wj_ys, _, _, _, _) = prepare_bench_ols(&y_scaled, &x);
        let out_ys = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef_ys,
            wj: &wj_ys,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        assert_eq!(out_ys.q, out_euclid.q);
        for j in 0..3 {
            assert!((out_ys.scpcstats[j * 6 + 1] - out_euclid.scpcstats[j * 6 + 1] * y_scale).abs() < 1e-22);
            assert!((out_ys.scpcstats[j * 6 + 2] - out_euclid.scpcstats[j * 6 + 2]).abs() < 1e-10);
            assert!((out_ys.scpcstats[j * 6 + 3] - out_euclid.scpcstats[j * 6 + 3]).abs() < 1e-10);
            assert!((out_ys.scpcstats[j * 6 + 4] - out_euclid.scpcstats[j * 6 + 4] * y_scale).abs() < 1e-22);
            assert!((out_ys.scpcstats[j * 6 + 5] - out_euclid.scpcstats[j * 6 + 5] * y_scale).abs() < 1e-22);
        }

        // Coordinate rescaling (`coords * 1e-6`) preserves q, se, t, p and scales c0 by 1e6
        let c_scale = 1e-6;
        let coords_scaled: Vec<f64> = coords.iter().map(|&v| v * c_scale).collect();
        let out_cs = run_exact_scpc(&ExactScpcRequest {
            coords: &coords_scaled,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: false,
            coef: &coef,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_cs.q, out_euclid.q);
        assert!((out_cs.c0 * c_scale - out_euclid.c0).abs() / out_euclid.c0 < 1e-10);
        for (a, b) in out_cs.scpcstats.iter().zip(out_euclid.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-10, "a={a}, b={b}");
        }

        // Exact zero wj produces structured DegenerateVariance error
        let wj_zero = vec![0.0; n * 3];
        let err = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: true,
            cvs: false,
            coef: &coef,
            wj: &wj_zero,
            cluster_id: None,
            fe_id: None,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
            proj_rank: 0,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap_err();
        assert_eq!(err.kind(), ErrorKind::DegenerateVariance);

        // Latlong conditional (`s_1 = lon, s_2 = lat`)
        let out_geo = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::GreatCircleFractionLonLat,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_geo.q, 18);
    }

    #[test]
    fn exact_scpc_matches_release_score_statistics_bench_1000() {
        let (y, x, coords) = load_bench_csv("bench_1000.csv");
        let n = y.len();
        let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);

        let out = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(4),
        })
        .unwrap();

        assert_eq!(out.q, 8);
        let dist_1000 = distance_matrix(&coords, n, 2, Metric::Euclidean).unwrap();
        let dvec_1000 = lower_triangle(&dist_1000, n);
        assert!((crate::calibration::avg_correlation(out.c0, &dvec_1000) - 0.03).abs() < 1e-12);
        assert!((out.c0 - 0.43896032977784594).abs() < 1e-3, "c0={}", out.c0);
        assert!(
            (out.cvfin - 2.365574447704166).abs() < 1e-3,
            "cvfin={}",
            out.cvfin
        );
        let expected_stats = [
            [1.5780823405600297, 0.13216765018777865, 11.940004519396016],
            [1.745699446088806, 0.0697326838936868, 25.03416401913148],
            [-0.9039531747843945, 0.037101498337623266, -24.3643306951765],
        ];
        for (row, expected) in out.scpcstats.chunks_exact(6).zip(expected_stats) {
            assert!((row[0] - expected[0]).abs() < 1e-12);
            for (got, exp) in row.iter().zip(expected) {
                assert!((got - exp).abs() <= 2e-3 * exp.abs().max(1.0), "got={got}, exp={exp}");
            }
        }
    }

    #[test]
    fn cluster_id_five_site_cutoffs_and_singleton_equivalence() {
        // Five-site clustered Gaussian cutoffs from test_conditional_ols.do:
        // s = [0,1,2,4,6], 2 obs/cluster, x_g = {x_0 - 0.25, x_0 + 0.25}, y = x + (s == 0)
        let s_sites = [0.0, 1.0, 2.0, 4.0, 6.0];
        let x0 = [-1.0, -4.0, 2.0, 2.0, 1.0];
        let mut s_obs = Vec::new();
        let mut x_obs = Vec::new();
        let mut y_obs = Vec::new();
        let mut cl_id = Vec::new();
        for (g, (&sg, &xg)) in s_sites.iter().zip(x0.iter()).enumerate() {
            for delta in [-0.25_f64, 0.25_f64] {
                let xval = xg + delta;
                s_obs.push(sg);
                x_obs.push(xval);
                y_obs.push(xval + if sg == 0.0 { 1.0 } else { 0.0 });
                cl_id.push(g);
            }
        }
        let g_sites = s_sites.len();
        let m_obs = y_obs.len();

        // Case 1: With intercept (y ~ x), cluster-aggregated scores, n_eff = G = 5
        let mut x_int = Vec::with_capacity(m_obs * 2);
        for &xv in &x_obs {
            x_int.extend_from_slice(&[1.0, xv]);
        }
        let mut xtx = [0.0; 4];
        let mut xty = [0.0; 2];
        for i in 0..m_obs {
            for r in 0..2 {
                xty[r] += x_int[i * 2 + r] * y_obs[i];
                for c in 0..2 {
                    xtx[r * 2 + c] += x_int[i * 2 + r] * x_int[i * 2 + c];
                }
            }
        }
        let bread = invert_symmetric(&xtx, 2);
        let beta: Vec<f64> = (0..2)
            .map(|r| (0..2).map(|c| bread[r * 2 + c] * xty[c]).sum())
            .collect();
        let mut s_cl = vec![0.0; g_sites * 2];
        for i in 0..m_obs {
            let resid = y_obs[i] - (beta[0] + beta[1] * x_obs[i]);
            let g = cl_id[i];
            s_cl[g * 2] += x_int[i * 2] * resid;
            s_cl[g * 2 + 1] += x_int[i * 2 + 1] * resid;
        }
        let mut wj_int = vec![0.0; g_sites];
        for g in 0..g_sites {
            let mut dot = 0.0;
            for c in 0..2 {
                dot += bread[2 + c] * s_cl[g * 2 + c];
            }
            wj_int[g] = (g_sites as f64) * dot + beta[1];
        }
        let mut xj_int = vec![0.0; m_obs];
        let mut sq_sum_int = vec![0.0; g_sites];
        for i in 0..m_obs {
            let val = (g_sites as f64) * (bread[2] + bread[3] * x_obs[i]);
            xj_int[i] = val;
            sq_sum_int[cl_id[i]] += val * val;
        }
        let xjs_int: Vec<f64> = (0..m_obs)
            .map(|i| {
                let d = sq_sum_int[cl_id[i]].sqrt();
                if d > 0.0 { xj_int[i] / d } else { 0.0 }
            })
            .collect();
        let q_int = modified_gram_schmidt(&x_int, m_obs, 2);
        let out_int = run_exact_scpc(&ExactScpcRequest {
            coords: &s_sites,
            n: g_sites,
            d: 1,
            metric: Metric::Euclidean,
            avc: 0.5,
            uncond: false,
            cvs: true,
            coef: &[beta[1]],
            wj: &wj_int,
            cluster_id: Some(&cl_id),
            fe_id: None,
            xj: Some(&xj_int),
            xjs: Some(&xjs_int),
            proj_u: Some(&q_int),
            proj_v: Some(&q_int),
            proj_rank: 2,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        let cv5_int = out_int.scpccvs.as_ref().unwrap()[2];
        assert!(
            (cv5_int - 17.380977).abs() < 0.002,
            "clustered int cv5={cv5_int}"
        );

        // Case 2: No intercept (y ~ x - 1)
        let denom_noint: f64 = x_obs.iter().map(|v| v * v).sum();
        let beta_noint: f64 = x_obs.iter().zip(y_obs.iter()).map(|(xv, yv)| xv * yv).sum::<f64>()
            / denom_noint;
        let mut s_cl_no = vec![0.0; g_sites];
        for i in 0..m_obs {
            let resid = y_obs[i] - beta_noint * x_obs[i];
            s_cl_no[cl_id[i]] += x_obs[i] * resid;
        }
        let wj_no: Vec<f64> = s_cl_no
            .iter()
            .map(|sg| (g_sites as f64) * (sg / denom_noint) + beta_noint)
            .collect();
        let xj_no: Vec<f64> = x_obs.iter().map(|xv| (g_sites as f64) * (xv / denom_noint)).collect();
        let mut sq_sum_no = vec![0.0; g_sites];
        for i in 0..m_obs {
            sq_sum_no[cl_id[i]] += xj_no[i] * xj_no[i];
        }
        let xjs_no: Vec<f64> = (0..m_obs)
            .map(|i| {
                let d = sq_sum_no[cl_id[i]].sqrt();
                if d > 0.0 { xj_no[i] / d } else { 0.0 }
            })
            .collect();
        let q_no = modified_gram_schmidt(&x_obs, m_obs, 1);
        let out_no = run_exact_scpc(&ExactScpcRequest {
            coords: &s_sites,
            n: g_sites,
            d: 1,
            metric: Metric::Euclidean,
            avc: 0.5,
            uncond: false,
            cvs: true,
            coef: &[beta_noint],
            wj: &wj_no,
            cluster_id: Some(&cl_id),
            fe_id: None,
            xj: Some(&xj_no),
            xjs: Some(&xjs_no),
            proj_u: Some(&q_no),
            proj_v: Some(&q_no),
            proj_rank: 1,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        let cv5_no = out_no.scpccvs.as_ref().unwrap()[2];
        assert!(
            (cv5_no - 9.902023).abs() < 0.002,
            "clustered noint cv5={cv5_no}"
        );

        // Case 3: Singleton clusters (cluster_id = 0..n-1) match unclustered cluster_id = None to < 1e-12
        let (y, x, coords) = load_bench_csv("bench_200.csv");
        let n = y.len();
        let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);
        let singleton: Vec<usize> = (0..n).collect();
        let out_none = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        let out_sing = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj,
            cluster_id: Some(&singleton),
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        for (a, b) in out_none.scpcstats.iter().zip(out_sing.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-12, "a={a}, b={b}");
        }
        for (a, b) in out_none
            .scpccvs
            .unwrap()
            .iter()
            .zip(out_sing.scpccvs.unwrap().iter())
        {
            assert!((a - b).abs() < 1e-12, "cv a={a}, b={b}");
        }
    }

    #[test]
    fn exact_scpc_iv_and_one_way_fe_match_explicit_dummies() {
        let (_, x_all, _) = load_bench_csv("bench_200.csv");
        let n_obs = 60usize;
        let g_sites = 30usize;
        let n_fe = 5usize;
        // Staggered clusters (pairs of obs) and fixed effects (crossing cluster boundaries)
        let cl_id: Vec<usize> = (0..n_obs).map(|i| i / 2).collect();
        let fe_id: Vec<usize> = (0..n_obs).map(|i| i % n_fe).collect();
        let mut s_cl = Vec::with_capacity(g_sites * 2);
        for g in 0..g_sites {
            s_cl.push((g as f64) * 0.13 + ((g % 5) as f64) * 0.21);
            s_cl.push(((g % 7) as f64) * 0.17 - (g as f64) * 0.04);
        }
        let s_obs: Vec<f64> = cl_id
            .iter()
            .enumerate()
            .flat_map(|(i, &g)| [s_cl[g * 2] + (i as f64) * 0.01, s_cl[g * 2 + 1]])
            .collect();

        // Build exogenous w, excluded instrument z, endogenous x1, and outcome y
        let mut w = Vec::with_capacity(n_obs);
        let mut z = Vec::with_capacity(n_obs);
        let mut x1 = Vec::with_capacity(n_obs);
        let mut y = Vec::with_capacity(n_obs);
        for i in 0..n_obs {
            let wi = x_all[i * 3 + 1];
            let zi = x_all[i * 3 + 2] + 0.3 * (i as f64).sin();
            let ui = 0.4 * ((i as f64) * 1.1).cos() + 0.15 * (fe_id[i] as f64);
            let xi = 0.85 * zi + 0.35 * wi + 0.55 * ui + 0.2 * ((fe_id[i] as f64) - 2.0);
            let yi = 1.2 * xi + 0.6 * wi + 0.9 * (fe_id[i] as f64) + ui;
            w.push(wi);
            z.push(zi);
            x1.push(xi);
            y.push(yi);
        }

        let demean = |col: &[f64]| -> Vec<f64> {
            let mut sums = vec![0.0; n_fe];
            let mut cnts = vec![0usize; n_fe];
            for (i, &val) in col.iter().enumerate() {
                sums[fe_id[i]] += val;
                cnts[fe_id[i]] += 1;
            }
            col.iter()
                .enumerate()
                .map(|(i, &val)| val - sums[fe_id[i]] / (cnts[fe_id[i]] as f64))
                .collect()
        };

        // --- Part 1: Unclustered & Clustered 1-way FE OLS (y ~ x1 + w | fe) ---
        let x_dm_c0 = demean(&x1);
        let x_dm_c1 = demean(&w);
        let y_dm = demean(&y);
        let mut x_dm = vec![0.0; n_obs * 2];
        let mut x_full_dum = vec![0.0; n_obs * (2 + n_fe)];
        for i in 0..n_obs {
            x_dm[i * 2] = x_dm_c0[i];
            x_dm[i * 2 + 1] = x_dm_c1[i];
            x_full_dum[i * (2 + n_fe)] = x1[i];
            x_full_dum[i * (2 + n_fe) + 1] = w[i];
            x_full_dum[i * (2 + n_fe) + 2 + fe_id[i]] = 1.0;
        }
        let solve_ols = |x_mat: &[f64], y_vec: &[f64], p_cols: usize| -> (Vec<f64>, Vec<f64>, Vec<f64>) {
            let mut xtx = vec![0.0; p_cols * p_cols];
            let mut xty = vec![0.0; p_cols];
            for i in 0..n_obs {
                for r in 0..p_cols {
                    xty[r] += x_mat[i * p_cols + r] * y_vec[i];
                    for c in 0..p_cols {
                        xtx[r * p_cols + c] += x_mat[i * p_cols + r] * x_mat[i * p_cols + c];
                    }
                }
            }
            let inv = invert_symmetric(&xtx, p_cols);
            let mut beta = xty.clone();
            let mut a = xtx;
            solve_linear(&mut a, &mut beta, p_cols);
            let bread: Vec<f64> = inv.iter().map(|v| v * (n_obs as f64)).collect();
            let mut res = vec![0.0; n_obs];
            for i in 0..n_obs {
                let fit: f64 = (0..p_cols).map(|c| x_mat[i * p_cols + c] * beta[c]).sum();
                res[i] = y_vec[i] - fit;
            }
            (beta, bread, res)
        };
        let (beta_dm, bread_dm, res_dm) = solve_ols(&x_dm, &y_dm, 2);
        let mut wj_ols_u = vec![0.0; n_obs * 2];
        let mut xj_ols = vec![0.0; n_obs * 2];
        for i in 0..n_obs {
            for j in 0..2 {
                let bdot: f64 = (0..2)
                    .map(|c| bread_dm[j * 2 + c] * x_dm[i * 2 + c])
                    .sum();
                wj_ols_u[i * 2 + j] = (n_obs as f64) * res_dm[i] * bdot + beta_dm[j];
                xj_ols[i * 2 + j] = (n_obs as f64) * bdot;
            }
        }
        let xjs_ols_u: Vec<f64> = xj_ols.iter().map(|&v| v.signum()).collect();
        let q_dm = modified_gram_schmidt(&x_dm, n_obs, 2);
        let q_dum = modified_gram_schmidt(&x_full_dum, n_obs, 2 + n_fe);

        let out_fe_ols = run_exact_scpc(&ExactScpcRequest {
            coords: &s_obs,
            n: n_obs,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.1,
            uncond: false,
            cvs: true,
            coef: &beta_dm,
            wj: &wj_ols_u,
            cluster_id: None,
            fe_id: Some(&fe_id),
            xj: Some(&xj_ols),
            xjs: Some(&xjs_ols_u),
            proj_u: Some(&q_dm),
            proj_v: Some(&q_dm),
            proj_rank: 2,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        let out_dum_ols = run_exact_scpc(&ExactScpcRequest {
            coords: &s_obs,
            n: n_obs,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.1,
            uncond: false,
            cvs: true,
            coef: &beta_dm,
            wj: &wj_ols_u,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj_ols),
            xjs: Some(&xjs_ols_u),
            proj_u: Some(&q_dum),
            proj_v: Some(&q_dum),
            proj_rank: 2 + n_fe,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        assert_eq!(out_fe_ols.q, out_dum_ols.q);
        for (a, b) in out_fe_ols.scpcstats.iter().zip(out_dum_ols.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-11, "uncl fe ols stats a={a}, b={b}");
        }
        for (a, b) in out_fe_ols.scpccvs.unwrap().iter().zip(out_dum_ols.scpccvs.unwrap().iter()) {
            assert!((a - b).abs() < 1e-11, "uncl fe ols cvs a={a}, b={b}");
        }

        // Clustered 1-way FE OLS
        let mut s_cl_ols = vec![0.0; g_sites * 2];
        for i in 0..n_obs {
            let g = cl_id[i];
            for c in 0..2 {
                s_cl_ols[g * 2 + c] += x_dm[i * 2 + c] * res_dm[i];
            }
        }
        let mut wj_ols_c = vec![0.0; g_sites * 2];
        for g in 0..g_sites {
            for j in 0..2 {
                let mut dot = 0.0;
                for c in 0..2 {
                    dot += (bread_dm[j * 2 + c] / (n_obs as f64)) * s_cl_ols[g * 2 + c];
                }
                wj_ols_c[g * 2 + j] = (g_sites as f64) * dot + beta_dm[j];
            }
        }
        let xj_ols_c: Vec<f64> = (0..n_obs)
            .flat_map(|i| {
                let x_dm = &x_dm;
                let bread_dm = &bread_dm;
                (0..2).map(move |j| {
                    let bdot: f64 = (0..2)
                        .map(|c| (bread_dm[j * 2 + c] / (n_obs as f64)) * x_dm[i * 2 + c])
                        .sum();
                    (g_sites as f64) * bdot
                })
            })
            .collect();
        let mut sq_ols_c = vec![0.0; g_sites * 2];
        for i in 0..n_obs {
            for j in 0..2 {
                sq_ols_c[cl_id[i] * 2 + j] += xj_ols_c[i * 2 + j].powi(2);
            }
        }
        let xjs_ols_c: Vec<f64> = (0..n_obs)
            .flat_map(|i| {
                let cl_id = &cl_id;
                let xj_ols_c = &xj_ols_c;
                let sq_ols_c = &sq_ols_c;
                (0..2).map(move |j| {
                    let d = sq_ols_c[cl_id[i] * 2 + j].sqrt();
                    if d > 0.0 { xj_ols_c[i * 2 + j] / d } else { 0.0 }
                })
            })
            .collect();

        let out_fe_ols_c = run_exact_scpc(&ExactScpcRequest {
            coords: &s_cl,
            n: g_sites,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.1,
            uncond: false,
            cvs: true,
            coef: &beta_dm,
            wj: &wj_ols_c,
            cluster_id: Some(&cl_id),
            fe_id: Some(&fe_id),
            xj: Some(&xj_ols_c),
            xjs: Some(&xjs_ols_c),
            proj_u: Some(&q_dm),
            proj_v: Some(&q_dm),
            proj_rank: 2,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();
        let out_dum_ols_c = run_exact_scpc(&ExactScpcRequest {
            coords: &s_cl,
            n: g_sites,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.1,
            uncond: false,
            cvs: true,
            coef: &beta_dm,
            wj: &wj_ols_c,
            cluster_id: Some(&cl_id),
            fe_id: None,
            xj: Some(&xj_ols_c),
            xjs: Some(&xjs_ols_c),
            proj_u: Some(&q_dum),
            proj_v: Some(&q_dum),
            proj_rank: 2 + n_fe,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();
        for (a, b) in out_fe_ols_c.scpcstats.iter().zip(out_dum_ols_c.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-11, "cl fe ols stats a={a}, b={b}");
        }

        // --- Part 2: Unclustered & Clustered 1-way FE IV (y ~ w | fe | x1 ~ z) ---
        let z_dm_c0 = demean(&z);
        let mut z_dm = vec![0.0; n_obs * 2];
        let mut z_full_dum = vec![0.0; n_obs * (2 + n_fe)];
        for i in 0..n_obs {
            z_dm[i * 2] = z_dm_c0[i];
            z_dm[i * 2 + 1] = x_dm_c1[i];
            z_full_dum[i * (2 + n_fe)] = z[i];
            z_full_dum[i * (2 + n_fe) + 1] = w[i];
            z_full_dum[i * (2 + n_fe) + 2 + fe_id[i]] = 1.0;
        }
        let build_iv_proj = |x_mat: &[f64], z_mat: &[f64], p_cols: usize| -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
            let qz = modified_gram_schmidt(z_mat, n_obs, p_cols);
            let rz = qz.len() / n_obs;
            let mut qzt_x = vec![0.0; rz * p_cols];
            linalg::multiply(&qz, x_mat, &mut qzt_x, rz, n_obs, p_cols, true);
            let mut pzx = vec![0.0; n_obs * p_cols];
            linalg::multiply(&qz, &qzt_x, &mut pzx, n_obs, rz, p_cols, false);
            let mut a = vec![0.0; p_cols * p_cols];
            linalg::multiply(x_mat, &pzx, &mut a, p_cols, n_obs, p_cols, true);
            let a_inv = invert_symmetric(&a, p_cols);
            let mut proj_v = vec![0.0; n_obs * p_cols];
            linalg::multiply(&pzx, &a_inv, &mut proj_v, n_obs, p_cols, p_cols, false);
            (pzx, a_inv, x_mat.to_vec(), proj_v)
        };
        let (pzx_dm, a_inv_dm, u_iv_dm, v_iv_dm) = build_iv_proj(&x_dm, &z_dm, 2);
        let (_, _, u_iv_dum, v_iv_dum) = build_iv_proj(&x_full_dum, &z_full_dum, 2 + n_fe);

        // Compute IV coefficients & influence functions
        let mut pzy = vec![0.0; 2];
        for i in 0..n_obs {
            pzy[0] += pzx_dm[i * 2] * y_dm[i];
            pzy[1] += pzx_dm[i * 2 + 1] * y_dm[i];
        }
        let mut beta_iv = vec![0.0; 2];
        for j in 0..2 {
            beta_iv[j] = a_inv_dm[j * 2] * pzy[0] + a_inv_dm[j * 2 + 1] * pzy[1];
        }
        let mut res_iv = vec![0.0; n_obs];
        for i in 0..n_obs {
            res_iv[i] = y_dm[i] - (x_dm[i * 2] * beta_iv[0] + x_dm[i * 2 + 1] * beta_iv[1]);
        }
        let mut wj_iv_c = vec![0.0; g_sites * 2];
        let mut s_cl_iv = vec![0.0; g_sites * 2];
        for i in 0..n_obs {
            let g = cl_id[i];
            for c in 0..2 {
                s_cl_iv[g * 2 + c] += pzx_dm[i * 2 + c] * res_iv[i];
            }
        }
        for g in 0..g_sites {
            for j in 0..2 {
                let dot: f64 = (0..2).map(|c| a_inv_dm[j * 2 + c] * s_cl_iv[g * 2 + c]).sum();
                wj_iv_c[g * 2 + j] = (g_sites as f64) * dot + beta_iv[j];
            }
        }
        let xj_iv_c: Vec<f64> = (0..n_obs)
            .flat_map(|i| {
                let pzx_dm = &pzx_dm;
                let a_inv_dm = &a_inv_dm;
                (0..2).map(move |j| {
                    let dot: f64 = (0..2)
                        .map(|c| a_inv_dm[j * 2 + c] * pzx_dm[i * 2 + c])
                        .sum();
                    (g_sites as f64) * dot
                })
            })
            .collect();
        let mut sq_iv_c = vec![0.0; g_sites * 2];
        for i in 0..n_obs {
            for j in 0..2 {
                sq_iv_c[cl_id[i] * 2 + j] += xj_iv_c[i * 2 + j].powi(2);
            }
        }
        let xjs_iv_c: Vec<f64> = (0..n_obs)
            .flat_map(|i| {
                let cl_id = &cl_id;
                let xj_iv_c = &xj_iv_c;
                let sq_iv_c = &sq_iv_c;
                (0..2).map(move |j| {
                    let d = sq_iv_c[cl_id[i] * 2 + j].sqrt();
                    if d > 0.0 { xj_iv_c[i * 2 + j] / d } else { 0.0 }
                })
            })
            .collect();

        // IV denominator weights are Xi_s' M' Xi R where M' = I - V U'.
        let out_fe_iv_c = run_exact_scpc(&ExactScpcRequest {
            coords: &s_cl,
            n: g_sites,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.1,
            uncond: false,
            cvs: true,
            coef: &beta_iv,
            wj: &wj_iv_c,
            cluster_id: Some(&cl_id),
            fe_id: Some(&fe_id),
            xj: Some(&xj_iv_c),
            xjs: Some(&xjs_iv_c),
            proj_u: Some(&v_iv_dm),
            proj_v: Some(&u_iv_dm),
            proj_rank: 2,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();
        let out_dum_iv_c = run_exact_scpc(&ExactScpcRequest {
            coords: &s_cl,
            n: g_sites,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.1,
            uncond: false,
            cvs: true,
            coef: &beta_iv,
            wj: &wj_iv_c,
            cluster_id: Some(&cl_id),
            fe_id: None,
            xj: Some(&xj_iv_c),
            xjs: Some(&xjs_iv_c),
            proj_u: Some(&v_iv_dum),
            proj_v: Some(&u_iv_dum),
            proj_rank: 2 + n_fe,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_fe_iv_c.q, out_dum_iv_c.q);
        for (a, b) in out_fe_iv_c.scpcstats.iter().zip(out_dum_iv_c.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-11, "cl fe iv stats a={a}, b={b}");
        }
        for (a, b) in out_fe_iv_c.scpccvs.unwrap().iter().zip(out_dum_iv_c.scpccvs.unwrap().iter()) {
            assert!((a - b).abs() < 1e-11, "cl fe iv cvs a={a}, b={b}");
        }
    }

    fn permute_rows(src: &[f64], perm: &[usize], cols: usize) -> Vec<f64> {
        let mut out = Vec::with_capacity(perm.len() * cols);
        for &src_idx in perm {
            out.extend_from_slice(&src[src_idx * cols..(src_idx + 1) * cols]);
        }
        out
    }

    #[test]
    fn approx_scpc_matches_release_and_is_thread_count_invariant() {
        let (y, x, coords) = load_bench_csv("bench_1000.csv");
        let n = y.len();
        let p = 3usize;
        let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);

        for metric in [Metric::Euclidean, Metric::GreatCircleFractionLonLat] {
            let (coords_norm, perm) = crate::nystrom::normalize_coords(&coords, n, 2, metric).unwrap();
            let wj_p = permute_rows(&wj, &perm, p);
            let xj_p = permute_rows(&xj, &perm, p);
            let xjs_p = permute_rows(&xjs, &perm, p);
            let q_p = permute_rows(&q_mat, &perm, proj_rank);

            let out_1 = run_approx_scpc(&ApproxScpcRequest {
                coords: &coords_norm,
                n,
                d: 2,
                metric,
                avc: 0.03,
                uncond: false,
                cvs: true,
                coef: &coef,
                wj: &wj_p,
                cluster_id: None,
                fe_id: None,
                xj: Some(&xj_p),
                xjs: Some(&xjs_p),
                proj_u: Some(&q_p),
                proj_v: Some(&q_p),
                proj_rank,
                cv_tolerance: 0.001,
                nthreads: Some(1),
                seed: 1,
                landmarks: 1000,
                samples: 20,
                pair_count: 1_000_000,
            })
            .unwrap();

            let out_4 = run_approx_scpc(&ApproxScpcRequest {
                coords: &coords_norm,
                n,
                d: 2,
                metric,
                avc: 0.03,
                uncond: false,
                cvs: true,
                coef: &coef,
                wj: &wj_p,
                cluster_id: None,
                fe_id: None,
                xj: Some(&xj_p),
                xjs: Some(&xjs_p),
                proj_u: Some(&q_p),
                proj_v: Some(&q_p),
                proj_rank,
                cv_tolerance: 0.001,
                nthreads: Some(4),
                seed: 1,
                landmarks: 1000,
                samples: 20,
                pair_count: 1_000_000,
            })
            .unwrap();

            assert!(!out_1.tied);
            assert!(out_1.q >= 1);
            assert!(out_1.c0 > 0.0 && out_1.cmax >= out_1.c0 && out_1.cvfin > 0.0);
            assert_eq!(out_1.q, out_4.q);
            assert!((out_1.c0 - out_4.c0).abs() < 1e-12);
            assert!((out_1.cmax - out_4.cmax).abs() < 1e-12);
            assert!((out_1.cvfin - out_4.cvfin).abs() < 1e-12);
            for (row1, row4) in out_1.scpcstats.chunks_exact(6).zip(out_4.scpcstats.chunks_exact(6)) {
                assert!(row1[1] > 0.0 && row1[4] < row1[0] && row1[0] < row1[5]);
                assert!(row1[3] >= 0.0 && row1[3] <= 1.0);
                for (&a, &b) in row1.iter().zip(row4.iter()) {
                    assert!((a - b).abs() < 1e-12, "approx stats mismatch a={a}, b={b}");
                }
            }
            for (&a, &b) in out_1
                .scpccvs
                .as_ref()
                .unwrap()
                .iter()
                .zip(out_4.scpccvs.as_ref().unwrap().iter())
            {
                assert!((a - b).abs() < 1e-12, "approx cvs mismatch a={a}, b={b}");
            }
            for (&a, &b) in out_1.wfin.iter().zip(out_4.wfin.iter()) {
                assert!((a - b).abs() < 1e-12, "approx wfin mismatch a={a}, b={b}");
            }
        }

        // Unconditional Nystrom run on Euclidean coords is also thread-count invariant.
        let (coords_norm, perm) =
            crate::nystrom::normalize_coords(&coords, n, 2, Metric::Euclidean).unwrap();
        let wj_p = permute_rows(&wj, &perm, p);
        let out_u1 = run_approx_scpc(&ApproxScpcRequest {
            coords: &coords_norm,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: true,
            cvs: true,
            coef: &coef,
            wj: &wj_p,
            cluster_id: None,
            fe_id: None,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
            proj_rank: 0,
            cv_tolerance: 0.001,
            nthreads: Some(1),
            seed: 1,
            landmarks: 1000,
            samples: 20,
            pair_count: 1_000_000,
        })
        .unwrap();
        let out_u4 = run_approx_scpc(&ApproxScpcRequest {
            coords: &coords_norm,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: true,
            cvs: true,
            coef: &coef,
            wj: &wj_p,
            cluster_id: None,
            fe_id: None,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
            proj_rank: 0,
            cv_tolerance: 0.001,
            nthreads: Some(4),
            seed: 1,
            landmarks: 1000,
            samples: 20,
            pair_count: 1_000_000,
        })
        .unwrap();
        assert_eq!(out_u1.q, out_u4.q);
        for (&a, &b) in out_u1.scpcstats.iter().zip(out_u4.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
        for (&a, &b) in out_u1
            .scpccvs
            .as_ref()
            .unwrap()
            .iter()
            .zip(out_u4.scpccvs.as_ref().unwrap().iter())
        {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn iv_conditional_weights_use_transposed_residual_maker_and_match_ols_when_z_equals_x() {
        // Reproducer 1 from rust_core_handoff.md: s=[0,1,2,4,6], x=[-1,-4,2,2,1], z=x, y=x+(s==0).
        let n = 5usize;
        let coords = [0.0, 1.0, 2.0, 4.0, 6.0];
        let x = [-1.0, -4.0, 2.0, 2.0, 1.0];
        let y: Vec<f64> = x
            .iter()
            .enumerate()
            .map(|(i, &xi)| xi + if i == 0 { 1.0 } else { 0.0 })
            .collect();

        // Design [1, x], fit OLS y ~ 1 + x
        let mut xtx = [0.0; 4];
        let mut xty = [0.0; 2];
        for i in 0..n {
            xtx[0] += 1.0;
            xtx[1] += x[i];
            xtx[2] += x[i];
            xtx[3] += x[i] * x[i];
            xty[0] += y[i];
            xty[1] += x[i] * y[i];
        }
        let det = xtx[0] * xtx[3] - xtx[1] * xtx[2];
        let inv = [xtx[3] / det, -xtx[1] / det, -xtx[2] / det, xtx[0] / det];
        let beta = [
            inv[0] * xty[0] + inv[1] * xty[1],
            inv[2] * xty[0] + inv[3] * xty[1],
        ];
        let bread = inv;
        let mut wj = vec![0.0; n * 2];
        let mut xj = vec![0.0; n * 2];
        let mut xjs = vec![0.0; n * 2];
        let mut x_full = vec![0.0; n * 2];
        for i in 0..n {
            let resid = y[i] - (beta[0] + beta[1] * x[i]);
            let row = [1.0, x[i]];
            x_full[i * 2] = 1.0;
            x_full[i * 2 + 1] = x[i];
            for j in 0..2 {
                let d_ij = (n as f64) * (bread[j * 2] * row[0] + bread[j * 2 + 1] * row[1]);
                wj[i * 2 + j] = d_ij * resid + beta[j];
                xj[i * 2 + j] = d_ij;
                xjs[i * 2 + j] = if d_ij > 0.0 {
                    1.0
                } else if d_ij < 0.0 {
                    -1.0
                } else {
                    0.0
                };
            }
        }
        let q_mat = modified_gram_schmidt(&x_full, n, 2);
        let out_ols = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 1,
            metric: Metric::Euclidean,
            avc: 0.5,
            uncond: false,
            cvs: true,
            coef: &beta,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank: 2,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();

        // For IV with z=x: V = P_Z X (X' P_Z X)^{-1} = X (X' X)^{-1}, U = X.
        let mut v_iv = vec![0.0; n * 2];
        for i in 0..n {
            v_iv[i * 2] = inv[0] + inv[1] * x[i];
            v_iv[i * 2 + 1] = inv[2] + inv[3] * x[i];
        }
        let out_iv = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 1,
            metric: Metric::Euclidean,
            avc: 0.5,
            uncond: false,
            cvs: true,
            coef: &beta,
            wj: &wj,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&v_iv),
            proj_v: Some(&x_full),
            proj_rank: 2,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();

        assert!((out_ols.scpcstats[6 + 3] - 0.0241534748).abs() < 1e-8);
        assert!((out_ols.scpccvs.as_ref().unwrap()[4 + 2] - 17.62678224).abs() < 1e-4);
        for (&a, &b) in out_ols.scpcstats.iter().zip(&out_iv.scpcstats) {
            assert!((a - b).abs() < 1e-11, "OLS vs IV stat mismatch: {a} vs {b}");
        }
        for (&a, &b) in out_ols
            .scpccvs
            .as_ref()
            .unwrap()
            .iter()
            .zip(out_iv.scpccvs.as_ref().unwrap())
        {
            assert!((a - b).abs() < 1e-11, "OLS vs IV cv mismatch: {a} vs {b}");
        }

        // Independent quadrature for this covariance gives 0.07143375445 at
        // the previous IV cutoff. The corrected cutoff must control size.
        let conditional = build_conditional_basis(
            &out_iv.wfin, n, out_iv.q, None, None, &xj, &xjs, 2, 1,
            Some(&v_iv), Some(&x_full), 2,
        ).unwrap();
        let dist = distance_matrix(&coords, n, 1, Metric::Euclidean).unwrap();
        let omega = covariance_grids(&dist, n, &[&conditional], &[out_iv.q + 1], &[out_iv.c0])
            .unwrap().remove(0);
        let mut probability = ProbabilityWorkspace::new(out_iv.q);
        let root_q = (out_iv.q as f64).sqrt();
        let old_p = probability.evaluate(&omega, out_iv.q + 1, 14.50861818 / root_q).unwrap();
        assert!((old_p - 0.07143375445).abs() < 1e-8, "old rejection probability {old_p}");
        let corrected_p = probability.evaluate(
            &omega, out_iv.q + 1, out_iv.scpccvs.as_ref().unwrap()[6] / root_q,
        ).unwrap();
        assert!(corrected_p <= 0.05, "corrected rejection probability {corrected_p}");

        // Genuinely oblique IV projection: M = I - U V', check build_conditional_basis == D_xjs M' D_xj B_q
        // U = [1, x], V = P_Z X (X' P_Z X)^{-1} where X is oblique to Z
        let z_ob = [-2.0, -1.0, 0.5, 3.0, -0.5];
        // Center z_ob and x to keep exact orthogonality with [1]
        let mx = x.iter().sum::<f64>() / (n as f64);
        let mz = z_ob.iter().sum::<f64>() / (n as f64);
        let xc: Vec<f64> = x.iter().map(|v| v - mx).collect();
        let zc: Vec<f64> = z_ob.iter().map(|v| v - mz).collect();
        let ztz: f64 = zc.iter().map(|v| v * v).sum();
        let ztx: f64 = zc.iter().zip(&xc).map(|(a, b)| a * b).sum();
        let pzx1: Vec<f64> = zc.iter().map(|&zi| zi * (ztx / ztz)).collect();
        let xt_pzx: f64 = xc.iter().zip(&pzx1).map(|(a, b)| a * b).sum();
        let mut u_ob = vec![0.0; n * 2];
        let mut v_ob = vec![0.0; n * 2];
        let mut w_raw = vec![0.0; n * 2];
        for i in 0..n {
            u_ob[i * 2] = 1.0;
            u_ob[i * 2 + 1] = xc[i];
            v_ob[i * 2] = 1.0 / (n as f64);
            v_ob[i * 2 + 1] = pzx1[i] / xt_pzx;
            w_raw[i * 2] = 1.0 / (n as f64).sqrt();
            w_raw[i * 2 + 1] = ((i as f64) * 1.1).sin();
        }
        // Build expected D_xjs (I - V U') D_xj b_1 directly
        let basis_got = build_conditional_basis(
            &w_raw,
            n,
            1,
            None,
            None,
            &xj,
            &xjs,
            2,
            1,
            Some(&v_ob),
            Some(&u_ob),
            2,
        )
        .unwrap();
        let mut u_t_rx = [0.0; 2];
        for i in 0..n {
            let rx_i = w_raw[i * 2 + 1] * xj[i * 2 + 1];
            u_t_rx[0] += u_ob[i * 2] * rx_i;
            u_t_rx[1] += u_ob[i * 2 + 1] * rx_i;
        }
        for i in 0..n {
            let rx_i = w_raw[i * 2 + 1] * xj[i * 2 + 1];
            let m_t_rx_i = rx_i - (v_ob[i * 2] * u_t_rx[0] + v_ob[i * 2 + 1] * u_t_rx[1]);
            let exp_i = m_t_rx_i * xjs[i * 2 + 1];
            assert!((basis_got[i * 2 + 1] - exp_i).abs() < 1e-12);
        }
    }
}
