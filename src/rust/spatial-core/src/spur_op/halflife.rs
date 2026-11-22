//! SPUR half-life likelihoods, simulation, and confidence-set inversion.

use super::{
    center_columns, covariance_shapes_differ, factor_omega, prepare_y_projections,
    select_leading_basis, validate_emat, OmegaPairFactor, SpurGeometry, TILE_COLS,
};
use crate::distance::Metric;
use crate::threads::{is_serial, run_with_threads};
use crate::{linalg, CoreError};
use rayon::prelude::*;

pub(super) const MAX_HALFLIFE_TILE_BYTES: usize = 4 * 1024 * 1024;
pub(super) const HL_HA_LEN: usize = 50;
const HALF_LIFE_SEARCH_STEPS: usize = 128;
pub const HALF_LIFE_MAX_COMPONENTS: usize = HALF_LIFE_SEARCH_STEPS / 2 + 1;
const HALF_LIFE_BOUND_TOL: f64 = 1e-7;

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
    pub components: Vec<[f64; 2]>,
    pub max_dist: f64,
    pub pv_mat: Vec<f64>,
    pub hl_grid_ho: Vec<f64>,
    pub tied: bool,
}

#[inline(always)]
pub(super) fn half_power(x: f64, q: usize) -> f64 {
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
    acc
}

pub(super) fn halflife_tile_size(q: usize, p_pairs: usize, nrep: usize) -> (usize, bool) {
    // Bound both the shared packed coefficient bank and per-job draw scratch.
    let base_per_draw = 2 * q;
    let budget = MAX_HALFLIFE_TILE_BYTES / size_of::<f64>();
    let packed_width = base_per_draw + p_pairs + HL_HA_LEN + 2;
    let packed = p_pairs <= budget / HL_HA_LEN && packed_width <= budget;
    let per_draw = if packed { packed_width } else { base_per_draw };
    ((budget / per_draw).clamp(1, TILE_COLS).min(nrep), packed)
}

struct HalfLifeScratch {
    tile_e: Vec<f64>,
    e_tile: Vec<f64>,
    u_tile: Vec<f64>,
    q_tile: Vec<f64>,
    sum_ha_over_ho: Vec<f64>,
    col_sq: Vec<f64>,
    xc_obs: Vec<f64>,
}

impl HalfLifeScratch {
    fn new(q: usize, max_tile: usize, packed_pairs: Option<usize>) -> Self {
        Self {
            tile_e: vec![0.0; q * max_tile],
            e_tile: vec![0.0; q * max_tile],
            u_tile: vec![0.0; packed_pairs.unwrap_or(0) * max_tile],
            q_tile: vec![0.0; packed_pairs.map_or(0, |_| HL_HA_LEN * max_tile)],
            sum_ha_over_ho: vec![0.0; packed_pairs.map_or(0, |_| max_tile)],
            col_sq: vec![0.0; packed_pairs.map_or(0, |_| max_tile)],
            xc_obs: vec![0.0; q],
        }
    }
}

pub(super) fn halflife_cancellation_bound(matrix: &[f64], p_pairs: usize) -> f64 {
    // Sum of absolute quadratic-form terms is bounded by ||A||_F * ||z||².
    let scale = matrix.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let scaled_frobenius = matrix
        .iter()
        .map(|v| (v / scale).powi(2))
        .sum::<f64>()
        .sqrt();
    scale * (8.0 * f64::EPSILON * (p_pairs as f64) * scaled_frobenius)
}

#[inline(always)]
pub(super) fn add_halflife_likelihood_terms(
    sum: &mut [f64],
    forms: &[f64],
    inv_norms: &[f64],
    determinant_ratio: f64,
    power: impl Fn(f64) -> f64,
) {
    for ((total, &form), &inv_norm) in sum.iter_mut().zip(forms).zip(inv_norms) {
        let denominator = power(form * inv_norm);
        // Check before reciprocal. Non-normal powers send this draw to the
        // log-space calculation, even if the other alternatives are normal.
        let term = if denominator <= f64::MIN_POSITIVE.recip() {
            determinant_ratio * denominator.recip()
        } else {
            f64::NAN
        };
        *total += term;
    }
}

#[inline(never)]
fn halflife_log_draw_sum(
    alternatives: &[OmegaPairFactor],
    log_weights: &[f64],
    packed_forms: &[f64],
    draws: &[f64],
    cols: usize,
    c: usize,
    q: f64,
    log_null_density: f64,
) -> Result<f64, CoreError> {
    let mut log_ha = [0.0; HL_HA_LEN];
    for (j, ha) in alternatives.iter().enumerate() {
        let form = if packed_forms.is_empty() {
            let value = quadratic_form_from_upper(&ha.ch_omi.upper, draws, cols, c);
            if !value.is_finite() || value <= 0.0 {
                return Err(CoreError::degenerate(
                    "half-life simulation draw has an invalid covariance quadratic form",
                ));
            }
            value
        } else {
            packed_forms[j * cols + c]
        };
        log_ha[j] = log_weights[j] + ha.ch_omi.log_sqrt_det - 0.5 * q * form.ln();
    }
    Ok(log_sum_exp(&log_ha) - log_null_density)
}

pub(super) fn quadratic_form_from_upper(
    upper: &[f64],
    draws: &[f64],
    cols: usize,
    c: usize,
) -> f64 {
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

pub(super) fn log_sum_exp(terms: &[f64]) -> f64 {
    let max = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let scaled_sum = terms.iter().map(|&term| (term - max).exp()).sum::<f64>();
    max + scaled_sum.ln()
}

pub(super) fn accepts_halflife_null(exceedances: usize, nrep: usize, level: f64) -> bool {
    // Compare the complement from the integer count so an exact grid boundary
    // is not accepted due to rounding 1 - level.
    (nrep - exceedances) as f64 / ((nrep + 1) as f64) < level
}

/// Composite five-point Gauss--Legendre quadrature on ten segments of [0, 1].
/// The segments get shorter near zero, where small inter-site distances matter.
pub(super) fn half_life_mixture() -> ([f64; HL_HA_LEN], [f64; HL_HA_LEN]) {
    let nodes = [
        -0.906179845938664,
        -0.5384693101056831,
        0.0,
        0.5384693101056831,
        0.906179845938664,
    ];
    let weights = [
        0.2369268850561891,
        0.4786286704993665,
        128.0 / 225.0,
        0.4786286704993665,
        0.2369268850561891,
    ];
    let mut rates = [0.0; HL_HA_LEN];
    let mut masses = [0.0; HL_HA_LEN];
    let mut left = 0.0;
    for segment in 0..10 {
        let right = (2_u32.pow(segment) as f64) / 512.0;
        let half_width = (right - left) / 2.0;
        for j in 0..5 {
            let hl = left + half_width * (1.0 + nodes[j]);
            let index = segment as usize * 5 + j;
            rates[index] = std::f64::consts::LN_2 / hl;
            masses[index] = half_width * weights[j];
        }
        left = right;
    }
    (rates, masses)
}

fn from_compact_half_life(t: f64) -> f64 {
    t / (1.0 - t)
}

pub(super) fn half_life_covariances(
    geom: &SpurGeometry,
    mr: &[f64],
    q: usize,
    points: &[f64],
) -> Vec<Vec<f64>> {
    let rates: Vec<_> = points
        .iter()
        .map(|&t| std::f64::consts::LN_2 * (1.0 - t) / t)
        .collect();
    let regular: Vec<_> = rates
        .iter()
        .copied()
        .filter(|rate| rate.is_finite() && *rate >= 0.01)
        .collect();
    let mut regular = geom.projected_omega_grid(&regular, mr, q).into_iter();
    rates
        .into_iter()
        .map(|rate| {
            if rate.is_finite() && rate >= 0.01 {
                return regular.next().expect("regular half-life covariance");
            }
            let kernel: Vec<_> = geom
                .norm_distmat
                .iter()
                .map(|&d| {
                    if rate.is_infinite() {
                        if d == 0.0 {
                            1.0
                        } else {
                            0.0
                        }
                    } else if rate == 0.0 {
                        -0.5 * d
                    } else {
                        0.5 * (-rate * d).exp_m1() / rate
                    }
                })
                .collect();
            let mut kr = vec![0.0; geom.n * q];
            linalg::multiply(&kernel, mr, &mut kr, geom.n, geom.n, q, false);
            let mut omega = vec![0.0; q * q];
            linalg::multiply(mr, &kr, &mut omega, q, geom.n, q, true);
            super::symmetrize(&mut omega, q);
            omega
        })
        .collect()
}

pub(super) fn invert_half_life(
    grid: &[f64],
    accepted: &[bool],
    evaluate: impl Fn(&[f64]) -> Result<Vec<bool>, CoreError>,
) -> Result<Vec<[f64; 2]>, CoreError> {
    let mut brackets: Vec<_> = grid
        .windows(2)
        .zip(accepted.windows(2))
        .filter(|(_, flags)| flags[0] != flags[1])
        .map(|(points, flags)| (points[0], points[1], flags[0]))
        .collect();
    while brackets
        .first()
        .is_some_and(|&(lo, hi, _)| hi - lo > HALF_LIFE_BOUND_TOL)
    {
        let midpoints: Vec<_> = brackets
            .iter()
            .map(|&(lo, hi, _)| (lo + hi) / 2.0)
            .collect();
        let status = evaluate(&midpoints)?;
        for ((lo, hi, left_accepted), (mid, mid_accepted)) in
            brackets.iter_mut().zip(midpoints.into_iter().zip(status))
        {
            if mid_accepted == *left_accepted {
                *lo = mid;
            } else {
                *hi = mid;
            }
        }
    }
    let mut components = Vec::new();
    let mut start = 0.0;
    // Use the outside bracket endpoints so boundary rounding includes the accepted side.
    for (lo, hi, left_accepted) in brackets {
        if left_accepted {
            components.push([from_compact_half_life(start), from_compact_half_life(hi)]);
        } else {
            start = lo;
        }
    }
    if accepted.last().copied().unwrap_or(false) {
        components.push([from_compact_half_life(start), f64::INFINITY]);
    }
    Ok(components)
}

pub fn run_spur_halflife(req: &SpurHalfLifeRequest<'_>) -> Result<SpurHalfLifeOutput, CoreError> {
    validate_emat(req.emat, req.n, req.q, req.nrep, 1)?;
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
        let geom = SpurGeometry::from_coords(req.coords, req.n, req.d, req.metric)?;
        let q = req.q;
        let n = req.n;
        let nrep = req.nrep;

        let sigdm_bm = geom.sigma_lbm_dm();
        let (basis_r, tied) = select_leading_basis(&sigdm_bm, n, q)?;
        let mut mr = basis_r.clone();
        center_columns(&mut mr, n, q);
        let p_vec = prepare_y_projections(req.y, n, 1, &basis_r, q)?;

        let grid: Vec<_> = (0..=HALF_LIFE_SEARCH_STEPS)
            .map(|i| i as f64 / HALF_LIFE_SEARCH_STEPS as f64)
            .collect();
        let (rates, weights) = half_life_mixture();
        let log_weights = weights.map(f64::ln);
        let om_ha = geom.projected_omega_grid(&rates, &mr, q);
        let om_ho = half_life_covariances(&geom, &mr, q, &grid);
        let mut informative = false;
        for candidate in om_ho[1..].iter().chain(&om_ha) {
            if covariance_shapes_differ(candidate, &om_ho[0], q)? {
                informative = true;
                break;
            }
        }
        if !informative {
            return Err(CoreError::degenerate(
                "half-life experiment is uninformative: retained covariance shapes are proportional",
            ));
        }
        let ha_factors: Vec<_> = om_ha
            .iter()
            .map(|omega| factor_omega(omega, q))
            .collect::<Result<_, _>>()?;

        let q_f = q as f64;
        let ha_x_log: Vec<f64> = ha_factors
            .iter()
            .zip(log_weights)
            .map(|(ha, log_weight)| {
                let mut xc = vec![0.0_f64; q];
                linalg::multiply(&ha.ch_omi.upper, &p_vec, &mut xc, q, q, 1, false);
                let qnorm: f64 = xc.iter().map(|&v| v * v).sum();
                if !qnorm.is_finite() || qnorm <= 0.0 {
                    return Err(CoreError::degenerate(
                        "projected SPUR vector has zero quadratic form",
                    ));
                }
                Ok(log_weight + ha.ch_omi.log_sqrt_det - 0.5 * q_f * qnorm.ln())
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
        let emat_inv: Vec<f64> = emat_sq.iter().map(|value| value.recip()).collect();

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
                cancellation_scale.push(halflife_cancellation_bound(&ha.om_inv, p_pairs));
            }
        }

        let make_scratch = || HalfLifeScratch::new(q, max_tile, packed.then_some(p_pairs));
        let eval_null_hl = |scratch: &mut HalfLifeScratch,
                            ho: &OmegaPairFactor|
         -> Result<usize, CoreError> {
            let HalfLifeScratch {
                tile_e,
                e_tile,
                u_tile,
                q_tile,
                sum_ha_over_ho,
                col_sq,
                xc_obs,
            } = scratch;
            linalg::multiply(&ho.ch_omi.upper, &p_vec, xc_obs, q, q, 1, false);
            let qnorm_obs: f64 = xc_obs.iter().map(|&v| v * v).sum();
            if !qnorm_obs.is_finite() || qnorm_obs <= 0.0 {
                return Err(CoreError::degenerate(
                    "projected SPUR vector has zero quadratic form",
                ));
            }
            let log_ho_obs = ho.ch_omi.log_sqrt_det - 0.5 * q_f * qnorm_obs.ln();
            let log_obs_sum = log_sum_exp(&ha_x_log) - log_ho_obs;
            let obs_sum = log_obs_sum.exp();

            let det_ratios: [f64; HL_HA_LEN] = std::array::from_fn(|j| {
                weights[j] * (ha_factors[j].ch_omi.log_sqrt_det - ho.ch_omi.log_sqrt_det).exp()
            });
            let use_fast_pow = packed
                && q <= 60
                && obs_sum.is_finite()
                && obs_sum >= f64::MIN_POSITIVE
                && det_ratios
                    .iter()
                    .all(|&v| v.is_finite() && v > 1e-150 && v < 1e150);

            let mut exceed = 0usize;
            for start in (0..nrep).step_by(max_tile) {
                let cols = (nrep - start).min(max_tile);
                let e_tile = &mut e_tile[..q * cols];
                for (dst, row) in tile_e[..q * cols]
                    .chunks_exact_mut(cols)
                    .zip(req.emat.chunks_exact(nrep))
                {
                    dst.copy_from_slice(&row[start..start + cols]);
                }

                // e_tile = L_ho * tile_e has covariance om_ho and exact null quadratic form ||tile_e||^2.
                linalg::multiply(
                    &ho.ch_om.lower,
                    &tile_e[..q * cols],
                    e_tile,
                    q,
                    q,
                    cols,
                    false,
                );

                let emat_inv = &emat_inv[start..start + cols];
                if packed {
                    col_sq[..cols].fill(0.0);
                    let mut product_rows = u_tile[..p_pairs * cols].chunks_exact_mut(cols);
                    for (b, row_b) in e_tile.chunks_exact(cols).enumerate() {
                        for (a, row_a) in e_tile.chunks_exact(cols).take(b + 1).enumerate() {
                            let dst = product_rows.next().unwrap();
                            for ((product, &left), &right) in dst.iter_mut().zip(row_a).zip(row_b) {
                                *product = left * right;
                            }
                            if a == b {
                                for (sq, &term) in col_sq[..cols].iter_mut().zip(dst.iter()) {
                                    *sq += term;
                                }
                            }
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
                    if use_fast_pow {
                        sum_ha_over_ho[..cols].fill(0.0);
                    }
                    let max_norm = col_sq[..cols].iter().fold(0.0_f64, |max, &norm| {
                        if norm > 0.0 {
                            max.max(norm)
                        } else {
                            f64::INFINITY
                        }
                    });
                    for (j, row) in q_tile[..HL_HA_LEN * cols]
                        .chunks_exact_mut(cols)
                        .enumerate()
                    {
                        let thresh = cancellation_scale[j];
                        // Screen against the largest draw norm. If this bound
                        // fails, use each draw's own norm before repairing it.
                        let bound = thresh * max_norm;
                        let needs_repair = row.iter().fold(false, |bad, &form| {
                            bad | !form.is_finite() | !(form > bound)
                        });
                        if needs_repair {
                            for (c, (form, &norm)) in row.iter_mut().zip(col_sq.iter()).enumerate()
                            {
                                if !form.is_finite() || !(norm > 0.0 && *form > thresh * norm) {
                                    *form = quadratic_form_from_upper(
                                        &ha_factors[j].ch_omi.upper,
                                        e_tile,
                                        cols,
                                        c,
                                    );
                                    if !form.is_finite() || *form <= 0.0 {
                                        return Err(CoreError::degenerate(
                                            "half-life simulation draw has an invalid covariance quadratic form",
                                        ));
                                    }
                                }
                            }
                        }
                        if use_fast_pow {
                            let sum = &mut sum_ha_over_ho[..cols];
                            let coeff = det_ratios[j];
                            // Select the common exponents outside the per-draw loop.
                            match q {
                                8 => {
                                    add_halflife_likelihood_terms(sum, row, emat_inv, coeff, |x| {
                                        let x2 = x * x;
                                        x2 * x2
                                    })
                                }
                                15 => {
                                    add_halflife_likelihood_terms(sum, row, emat_inv, coeff, |x| {
                                        let x2 = x * x;
                                        let x4 = x2 * x2;
                                        let x7 = x4 * x2 * x;
                                        x7 * x.sqrt()
                                    })
                                }
                                _ => {
                                    add_halflife_likelihood_terms(sum, row, emat_inv, coeff, |x| {
                                        half_power(x, q)
                                    })
                                }
                            }
                        }
                    }
                }

                let q_tile = &q_tile[..if packed { HL_HA_LEN * cols } else { 0 }];
                for (c, &norm) in emat_sq[start..start + cols].iter().enumerate() {
                    let exceeds = if use_fast_pow
                        && sum_ha_over_ho[c].is_finite()
                        && sum_ha_over_ho[c] >= f64::MIN_POSITIVE
                        && emat_inv[c].is_finite()
                    {
                        sum_ha_over_ho[c] >= obs_sum
                    } else {
                        let log_ho = ho.ch_omi.log_sqrt_det - 0.5 * q_f * norm.ln();
                        halflife_log_draw_sum(
                            &ha_factors,
                            &log_weights,
                            q_tile,
                            e_tile,
                            cols,
                            c,
                            q_f,
                            log_ho,
                        )? >= log_obs_sum
                    };
                    if exceeds {
                        exceed += 1;
                    }
                }
            }
            Ok(exceed)
        };

        let evaluate = |covariances: &[Vec<f64>]| -> Result<Vec<usize>, CoreError> {
            let at = |scratch: &mut HalfLifeScratch, omega: &Vec<f64>| {
                eval_null_hl(scratch, &factor_omega(omega, q)?)
            };
            if is_serial() {
                let mut scratch = make_scratch();
                covariances
                    .iter()
                    .map(|omega| at(&mut scratch, omega))
                    .collect()
            } else {
                covariances
                    .par_iter()
                    .with_min_len(4)
                    .map_init(make_scratch, at)
                    .collect()
            }
        };
        let counts = evaluate(&om_ho)?;
        let accepted: Vec<_> = counts
            .iter()
            .map(|&count| accepts_halflife_null(count, nrep, req.level))
            .collect();
        let components = invert_half_life(&grid, &accepted, |points| {
            let omegas = half_life_covariances(&geom, &mr, q, points);
            Ok(evaluate(&omegas)?
                .into_iter()
                .map(|count| accepts_halflife_null(count, nrep, req.level))
                .collect())
        })?;
        let ci_lower = components.first().map_or(f64::NAN, |interval| interval[0]);
        let ci_upper = components.last().map_or(f64::NAN, |interval| interval[1]);
        Ok(SpurHalfLifeOutput {
            ci_lower,
            ci_upper,
            components,
            max_dist: geom.max_dist,
            pv_mat: counts
                .into_iter()
                .map(|count| (count + 1) as f64 / (nrep + 1) as f64)
                .collect(),
            hl_grid_ho: grid.into_iter().map(from_compact_half_life).collect(),
            tied,
        })
    })?
}
