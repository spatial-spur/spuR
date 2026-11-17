//! Scale-invariant SCPC decay calibration with co-location floor handling.

use crate::{linalg, CoreError};
use rayon::prelude::*;

const CHUNK_SIZE: usize = 32_768;

/// Evaluate average pairwise correlation `mean(exp(-c * dist))`.
pub fn avg_correlation(c: f64, distances: &[f64]) -> f64 {
    if distances.is_empty() {
        return f64::NAN;
    }
    // Chunk identically across serial and parallel paths so sum order is thread-count invariant.
    let chunk_sums: Vec<f64> = if !crate::threads::is_serial() && distances.len() >= 4 * CHUNK_SIZE {
        distances
            .par_chunks(CHUNK_SIZE)
            .map_init(Vec::new, |scratch, chunk| {
                scratch.resize(chunk.len(), 0.0);
                linalg::exp_scaled(chunk, -c, scratch);
                scratch.iter().sum::<f64>()
            })
            .collect()
    } else {
        let mut scratch = Vec::with_capacity(distances.len().min(CHUNK_SIZE));
        distances
            .chunks(CHUNK_SIZE)
            .map(|chunk| {
                scratch.resize(chunk.len(), 0.0);
                linalg::exp_scaled(chunk, -c, &mut scratch);
                scratch.iter().sum::<f64>()
            })
            .collect()
    };
    chunk_sums.into_iter().sum::<f64>() / distances.len() as f64
}

/// Fraction of zero-distance pairs in `distances` (co-location correlation floor).
pub fn correlation_floor(distances: &[f64]) -> Result<f64, CoreError> {
    if distances.is_empty() {
        return Err(CoreError::invalid_input(
            "at least two locations are required",
        ));
    }
    let mut zeros = 0usize;
    for &d in distances {
        if !d.is_finite() || d < 0.0 {
            return Err(CoreError::invalid_input(
                "distances must be finite and nonnegative",
            ));
        }
        if d == 0.0 {
            zeros += 1;
        }
    }
    Ok(zeros as f64 / distances.len() as f64)
}

fn eval_residual_and_slope_chunk(
    distances: &[f64],
    c: f64,
    near_one: bool,
    scratch: &mut Vec<f64>,
) -> (f64, f64) {
    if near_one {
        return distances
            .iter()
            .fold((0.0, 0.0), |(delta_sum, slope_sum), &d| {
                if d == 0.0 {
                    (delta_sum, slope_sum)
                } else {
                    let m1 = (-c * d).exp_m1();
                    (delta_sum + m1, slope_sum + d * (1.0 + m1))
                }
            });
    }
    scratch.resize(distances.len(), 0.0);
    linalg::exp_scaled(distances, -c, scratch);
    let mut mass = 0.0_f64;
    let mut slope = 0.0_f64;
    for (&e, &d) in scratch.iter().zip(distances) {
        mass += e;
        slope += e * d;
    }
    (mass, slope)
}

fn eval_residual_and_slope(
    distances: &[f64],
    c: f64,
    target_avc: f64,
    near_one: bool,
) -> (f64, f64) {
    let n = distances.len() as f64;
    let chunk_pairs: Vec<(f64, f64)> =
        if !crate::threads::is_serial() && distances.len() >= 4 * CHUNK_SIZE {
            distances
                .par_chunks(CHUNK_SIZE)
                .map_init(Vec::new, |scratch, chunk| {
                    eval_residual_and_slope_chunk(chunk, c, near_one, scratch)
                })
                .collect()
        } else {
            let mut scratch = Vec::with_capacity(distances.len().min(CHUNK_SIZE));
            distances
                .chunks(CHUNK_SIZE)
                .map(|chunk| eval_residual_and_slope_chunk(chunk, c, near_one, &mut scratch))
                .collect()
        };
    let (mass_or_delta, pos_slope) = chunk_pairs
        .into_iter()
        .fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
    let residual = if near_one {
        (1.0 - target_avc) + mass_or_delta / n
    } else {
        mass_or_delta / n - target_avc
    };
    (residual, -pos_slope / n)
}

/// Solve for exponential decay parameter `c` matching `target_avc` using a
/// scale-invariant bracket and safeguarded Newton + bisection.
pub fn solve_c0_from_avc(distances: &[f64], target_avc: f64) -> Result<f64, CoreError> {
    if distances.is_empty() {
        return Err(CoreError::invalid_input(
            "at least two locations are required",
        ));
    }
    let mut zeros = 0usize;
    let mut sum_pos = 0.0_f64;
    let mut max_d = 0.0_f64;
    for &d in distances {
        if !d.is_finite() || d < 0.0 {
            return Err(CoreError::invalid_input(
                "distances must be finite and nonnegative",
            ));
        }
        if d == 0.0 {
            zeros += 1;
        } else {
            sum_pos += d;
            max_d = max_d.max(d);
        }
    }
    let n = distances.len();
    let floor = zeros as f64 / n as f64;
    if max_d == 0.0
        || !target_avc.is_finite()
        || target_avc <= floor
        || target_avc >= 1.0
    {
        return Err(CoreError::calibration(
            "target correlation must exceed the co-located-pair fraction and be below one",
        ));
    }

    let pos_count = (n - zeros) as f64;
    let mean_pos = sum_pos / pos_count;
    let target_pos = (target_avc - floor) / (1.0 - floor);
    let near_one = target_avc > 0.9999;

    // Jensen's inequality gives mean_{d>0} exp(-c*d) >= exp(-c*mean_pos), so
    // c_init is invariant to positive coordinate rescaling and sits near the root.
    let raw_init = if mean_pos.is_finite() && mean_pos > 0.0 {
        let neg_ln = -target_pos.ln();
        let cand = neg_ln / mean_pos;
        if cand.is_finite() && cand > 0.0 {
            cand
        } else {
            max_d.recip().min(f64::MAX * 0.5).max(f64::MIN_POSITIVE)
        }
    } else {
        max_d.recip().min(f64::MAX * 0.5).max(f64::MIN_POSITIVE)
    };

    let mut lo = raw_init;
    let mut hi = raw_init;
    let (mut f_lo, mut d_lo) = eval_residual_and_slope(distances, lo, target_avc, near_one);
    if f_lo == 0.0 {
        return Ok(lo);
    }

    if f_lo < 0.0 {
        // Initial guess is above the root; halve downwards until residual >= 0.
        hi = lo;
        for _ in 0..1024 {
            lo *= 0.5;
            if lo == 0.0 || !lo.is_finite() {
                return Err(CoreError::calibration(
                    "decay calibration lower bracket underflowed; distances may be degenerate",
                ));
            }
            let (f, d) = eval_residual_and_slope(distances, lo, target_avc, near_one);
            f_lo = f;
            d_lo = d;
            if f_lo == 0.0 {
                return Ok(lo);
            }
            if f_lo > 0.0 {
                break;
            }
            hi = lo;
        }
        if f_lo < 0.0 {
            return Err(CoreError::calibration(
                "decay calibration lower bracket did not enclose root",
            ));
        }
    } else {
        // Initial guess is below the root; double upwards until residual <= 0.
        let mut f_hi = f_lo;
        for _ in 0..1024 {
            lo = hi;
            f_lo = f_hi;
            if hi >= f64::MAX * 0.5 {
                hi = f64::MAX;
                let (f, _) = eval_residual_and_slope(distances, hi, target_avc, near_one);
                f_hi = f;
                break;
            }
            hi *= 2.0;
            let (f, d) = eval_residual_and_slope(distances, hi, target_avc, near_one);
            f_hi = f;
            if f_hi == 0.0 {
                return Ok(hi);
            }
            if f_hi < 0.0 {
                break;
            }
            d_lo = d;
        }
        if f_hi > 0.0 {
            return Err(CoreError::calibration(
                "spatial decay calibration failed to bracket upper decay scale",
            ));
        }
    }

    let residual_tol = 1e-12 * (target_avc - floor).min(1.0 - target_avc).max(1e-15);
    let mut c = if d_lo < 0.0 && f_lo.is_finite() {
        let newton = lo - f_lo / d_lo;
        if newton > lo && newton < hi && newton.is_finite() {
            newton
        } else if lo > 0.0 && hi / lo > 4.0 {
            (lo * hi).sqrt()
        } else {
            0.5 * (lo + hi)
        }
    } else {
        0.5 * (lo + hi)
    };

    for _ in 0..128 {
        let (f, slope) = eval_residual_and_slope(distances, c, target_avc, near_one);
        if f.abs() <= residual_tol {
            return Ok(c);
        }
        if f > 0.0 {
            lo = c;
        } else {
            hi = c;
        }
        if (hi - lo) <= 1e-13 * c.max(f64::MIN_POSITIVE) {
            c = 0.5 * (lo + hi);
            break;
        }

        let midpoint = if lo > 0.0 && hi / lo > 4.0 {
            let geom = (lo * hi).sqrt();
            if geom > lo && geom < hi {
                geom
            } else {
                lo + 0.5 * (hi - lo)
            }
        } else {
            lo + 0.5 * (hi - lo)
        };

        let next = if slope < 0.0 && slope.is_finite() && f.is_finite() {
            let cand = c - f / slope;
            if cand > lo && cand < hi && cand.is_finite() {
                cand
            } else {
                midpoint
            }
        } else {
            midpoint
        };

        if next == c {
            c = midpoint;
            if c == lo || c == hi {
                break;
            }
        } else {
            c = next;
        }
    }

    let (final_f, _) = eval_residual_and_slope(distances, c, target_avc, near_one);
    let accept_tol = 1e-8 * (target_avc - floor).min(1.0 - target_avc).max(1e-15);
    if final_f.abs() <= accept_tol || (hi - lo) <= 1e-12 * c.max(f64::MIN_POSITIVE) {
        Ok(c)
    } else {
        Err(CoreError::calibration(
            "spatial decay calibration did not converge",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ErrorKind;

    #[test]
    fn decay_solver_is_scale_invariant_across_twenty_four_orders_of_magnitude() {
        let base = [1.0, 2.5, 4.0, 7.2, 10.5, 13.0];
        let target = 0.03;
        let c_base = solve_c0_from_avc(&base, target).unwrap();
        assert!((avg_correlation(c_base, &base) - target).abs() < 1e-12);

        for scale in [1e-12, 1e-6, 1e6, 1e12] {
            let scaled: Vec<f64> = base.iter().map(|&d| d * scale).collect();
            let c_scaled = solve_c0_from_avc(&scaled, target).unwrap();
            let rel = (c_scaled * scale - c_base).abs() / c_base;
            assert!(rel < 1e-11, "scale={scale}, rel={rel}");
            assert!((avg_correlation(c_scaled, &scaled) - target).abs() < 1e-11);
        }
    }

    #[test]
    fn decay_solver_handles_near_one_targets_and_colocation_floor() {
        // 2 out of 6 pairs are co-located -> floor = 1/3
        let distances = [0.0, 0.0, 1.0, 2.0, 3.0, 4.0];
        let floor = correlation_floor(&distances).unwrap();
        assert!((floor - 1.0 / 3.0).abs() < 1e-15);

        // Above floor succeeds
        let c_hi = solve_c0_from_avc(&distances, 0.4).unwrap();
        assert!((avg_correlation(c_hi, &distances) - 0.4).abs() < 1e-12);

        // Near-one target succeeds accurately
        let c_near1 = solve_c0_from_avc(&distances, 0.99999).unwrap();
        assert!((avg_correlation(c_near1, &distances) - 0.99999).abs() < 1e-12);

        // At or below floor fails cleanly with CalibrationFailed
        for bad in [floor, floor - 0.05, 1.0, 0.0] {
            let err = solve_c0_from_avc(&distances, bad).unwrap_err();
            assert_eq!(err.kind(), ErrorKind::CalibrationFailed);
        }

        // All identical distances fail cleanly
        let identical = [0.0, 0.0, 0.0];
        let err = solve_c0_from_avc(&identical, 0.5).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::CalibrationFailed);
    }
}
