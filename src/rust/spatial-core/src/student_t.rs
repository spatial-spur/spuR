//! Pure-Rust log-gamma and Student-t inverse CDF matching `scipy.stats.t.ppf` / `R::qt`.

use crate::CoreError;
use std::f64::consts::PI;

/// Lanczos approximation to `ln(Gamma(z))` for `z > 0`.
pub(crate) fn ln_gamma(z: f64) -> f64 {
    const G: f64 = 7.0;
    const COEFFS: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    if z < 0.5 {
        let sin_pi_z = (PI * z).sin();
        return PI.ln() - sin_pi_z.ln() - ln_gamma(1.0 - z);
    }
    let z = z - 1.0;
    let mut x = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        x += c / (z + i as f64);
    }
    let t = z + G + 0.5;
    0.5 * (2.0 * PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const MAX_IT: usize = 200;
    const EPS: f64 = f64::EPSILON;
    const TINY: f64 = 1e-30;

    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAX_IT {
        let m_f = m as f64;
        let m2 = 2.0 * m_f;
        let aa = m_f * (b - m_f) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;

        let aa = -(a + m_f) * (qab + m_f) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() <= EPS {
            break;
        }
    }
    h
}

pub(crate) fn betainc(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let ln_beta = ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b);
    let bt = (a * x.ln() + b * (1.0 - x).ln() - ln_beta).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

/// Student-t quantile `F^{-1}(p; df)` for `0 < p < 1` and `df >= 1`.
pub(crate) fn student_t_ppf(p: f64, df: usize) -> Result<f64, CoreError> {
    if !p.is_finite() || p <= 0.0 || p >= 1.0 || df == 0 {
        return Err(CoreError::invalid_input(
            "invalid Student-t quantile parameters",
        ));
    }
    if (p - 0.5).abs() == 0.0 {
        return Ok(0.0);
    }
    let upper = p > 0.5;
    let tail = if upper { 1.0 - p } else { p };
    let target = 2.0 * tail;
    let nu = df as f64;
    if df == 1 {
        let val = (PI * (p - 0.5)).tan();
        return Ok(val);
    }
    if df == 2 {
        let denom = (target * (2.0 - target)).sqrt();
        let mag = (1.0 - target) * 2.0_f64.sqrt() / denom;
        return Ok(if upper { mag } else { -mag });
    }

    let a = 0.5 * nu;
    let b = 0.5;
    let ln_beta = ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b);
    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;
    let mut x = 0.5_f64;
    for _ in 0..80 {
        let f = betainc(a, b, x) - target;
        if f.abs() <= 1e-15 * target.min(1.0) {
            break;
        }
        if f < 0.0 {
            lo = x;
        } else {
            hi = x;
        }
        let deriv = ((a - 1.0) * x.ln() + (b - 1.0) * (1.0 - x).ln() - ln_beta).exp();
        let next = if deriv.is_finite() && deriv > 0.0 {
            x - f / deriv
        } else {
            0.5 * (lo + hi)
        };
        x = if next > lo && next < hi {
            next
        } else {
            0.5 * (lo + hi)
        };
        if hi - lo <= 1e-15 * hi.max(1.0) {
            break;
        }
    }
    let mag = (nu * (1.0 - x) / x).sqrt();
    Ok(if upper { mag } else { -mag })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn student_t_quantiles_match_reference() {
        // Values checked against R qt(1 - level / 2, df)
        let cases = [
            (1, 0.05, 12.7062047361746924),
            (2, 0.05, 4.3026527297494619),
            (8, 0.05, 2.3060041352041658),
            (19, 0.05, 2.0930240544083092),
            (19, 0.01, 2.8609346064649790),
            (19, 0.10, 1.7291328115213689),
            (19, 0.32, 1.0211724139253480),
            (60, 0.05, 2.0002978220142595),
        ];
        for (df, level, expected) in cases {
            let got = student_t_ppf(1.0 - level * 0.5, df).unwrap();
            assert!(
                (got - expected).abs() < 1e-12,
                "df={df} level={level}: got={got}, expected={expected}"
            );
        }
    }
}
