//! SPUR half-life likelihoods, simulation, and confidence-set inversion.

use super::{
    center_columns, covariance_shapes_differ, factor_omega, prepare_y_projections,
    select_leading_basis, validate_emat, OmegaPairFactor, SpurGeometry, RHO_BM, TILE_COLS,
};
use crate::distance::Metric;
use crate::threads::{is_serial, run_with_threads};
use crate::{linalg, CoreError};
use rayon::prelude::*;

pub(super) const MAX_HALFLIFE_TILE_BYTES: usize = 4 * 1024 * 1024;
pub(super) const HL_HA_LEN: usize = 50;
const HL_HO_LEN: usize = HL_HO_SMALL_LEN + HL_HO_LARGE_LEN + 1;
const HL_HO_LARGE_LEN: usize = 30;
const HL_HO_SMALL_LEN: usize = 100;

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
        log_ha[j] = ha.ch_omi.log_sqrt_det - 0.5 * q * form.ln();
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
        let sigdm_bm = geom.kernel_dm(c_bm);
        let (basis_r, tied) = select_leading_basis(&sigdm_bm, n, q)?;
        let mut mr = basis_r.clone();
        center_columns(&mut mr, n, q);
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
                (ha_factors[j].ch_omi.log_sqrt_det - ho.ch_omi.log_sqrt_det).exp()
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
                        sum_ha_over_ho[c] > obs_sum
                    } else {
                        let log_ho = ho.ch_omi.log_sqrt_det - 0.5 * q_f * norm.ln();
                        halflife_log_draw_sum(&ha_factors, q_tile, e_tile, cols, c, q_f, log_ho)?
                            > log_obs_sum
                    };
                    if exceeds {
                        exceed += 1;
                    }
                }
            }
            Ok(exceed)
        };

        let exceedances = if is_serial() {
            let mut scratch = make_scratch();
            ho_factors
                .iter()
                .map(|ho| eval_null_hl(&mut scratch, ho))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            ho_factors
                .par_iter()
                .with_min_len(4)
                .map_init(make_scratch, eval_null_hl)
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
