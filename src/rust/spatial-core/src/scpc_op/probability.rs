//! Gaussian rejection probabilities and critical-value calibration.

use crate::linalg::EigenWorkspace;
use crate::quadrature::RULE;
use crate::student_t::student_t_ppf;
use crate::CoreError;
use std::f64::consts::{FRAC_2_PI, FRAC_PI_2};
use std::sync::OnceLock;

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
    if !ratios.iter().any(|&v| !v.is_finite()) {
        let product: f64 = ratios.iter().map(|value| cos2 / (cos2 + value)).product();
        if product >= f64::MIN_POSITIVE {
            return product.sqrt();
        }
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

pub(super) struct ProbabilityWorkspace {
    eigensolver: EigenWorkspace,
    ratios: Vec<f64>,
    logs: Vec<f64>,
}

impl ProbabilityWorkspace {
    pub(super) fn new(q: usize) -> Self {
        Self {
            eigensolver: EigenWorkspace::new(q + 1),
            ratios: Vec::with_capacity(q),
            logs: Vec::with_capacity(q),
        }
    }

    /// Evaluate tail-stable adaptive angular/complement rejection probability.
    pub(super) fn evaluate(
        &mut self,
        omega: &[f64],
        width: usize,
        cv: f64,
    ) -> Result<f64, CoreError> {
        if let Some(p) = self.eigensolver.rejection_eigen_ratios(
            omega,
            width,
            cv,
            &mut self.ratios,
            &mut self.logs,
        )? {
            return Ok(p);
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

    pub(super) fn maximum(
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
    pub(super) fn critical(
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
        let mut converged = false;

        for _ in 0..=grid_count + 4 {
            if rp <= level {
                converged = true;
                break;
            }
            let omega = &omegas[active * matrix_size..(active + 1) * matrix_size];
            cv1 = cv0;
            let mut p_upper = self.evaluate(omega, width, cv1)?;
            for iter in 0..1024 {
                if p_upper > level {
                    cv0 = cv1;
                    cv1 += if iter < 32 { step } else { cv1 };
                    if !cv1.is_finite() {
                        return Err(CoreError::calibration(
                            "SCPC critical value could not be bracketed",
                        ));
                    }
                    p_upper = self.evaluate(omega, width, cv1)?;
                } else {
                    break;
                }
            }
            if p_upper > level {
                return Err(CoreError::calibration(
                    "SCPC critical value could not be bracketed",
                ));
            }

            while cv1 - cv0 > tol {
                let mid = cv0 + 0.5 * (cv1 - cv0);
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
                converged = true;
                break;
            }
            active = max_i;
            rp = max_p;
            cv0 = cv1;
        }
        if !converged {
            return Err(CoreError::calibration(
                "SCPC critical value calibration did not converge across scale grid",
            ));
        }
        Ok(cv1 * root_q)
    }
}
