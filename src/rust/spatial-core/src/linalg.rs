use crate::CoreError;

#[cfg(all(target_vendor = "apple", not(feature = "portable")))]
#[link(name = "Accelerate", kind = "framework")]
extern "C" {
    fn cblas_dgemm(
        order: i32,
        trans_a: i32,
        trans_b: i32,
        m: i32,
        n: i32,
        k: i32,
        alpha: f64,
        a: *const f64,
        lda: i32,
        b: *const f64,
        ldb: i32,
        beta: f64,
        c: *mut f64,
        ldc: i32,
    );
    fn vvexp(output: *mut f64, input: *const f64, count: *const i32);

}

#[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
pub(crate) fn faer_par() -> faer::Par {
    faer::Par::Seq
}

#[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
#[inline(always)]
fn fast_exp_block8(x: [f64; 8]) -> [f64; 8] {
    const INV_LN2: f64 = 1.4426950408889634;
    const LN2_HI: f64 = 0.6931471804855391;
    const LN2_LO: f64 = 7.440617110012397e-11;
    const SHIFT: f64 = 6755399441055744.0;
    const C2: f64 = 0.5;
    const C3: f64 = 1.0 / 6.0;
    const C4: f64 = 1.0 / 24.0;
    const C5: f64 = 1.0 / 120.0;
    const C6: f64 = 1.0 / 720.0;
    const C7: f64 = 1.0 / 5040.0;
    const C8: f64 = 1.0 / 40320.0;
    const C9: f64 = 1.0 / 362880.0;
    const C10: f64 = 1.0 / 3628800.0;
    const C11: f64 = 1.0 / 39916800.0;
    const C12: f64 = 1.0 / 479001600.0;

    let mut out = [0.0_f64; 8];
    for i in 0..8 {
        let z = x[i];
        let t = z * INV_LN2 + SHIFT;
        let kf = t - SHIFT;
        let r = (z - kf * LN2_HI) - kf * LN2_LO;
        let r2 = r * r;
        let r4 = r2 * r2;
        let e0 = 1.0 + r + r2 * (C2 + r * C3);
        let e1 = C4 + r * C5 + r2 * (C6 + r * C7);
        let e2 = C8 + r * C9 + r2 * (C10 + r * C11) + r4 * C12;
        let p = e0 + r4 * (e1 + r4 * e2);
        let scale_pow2 = f64::from_bits(t.to_bits().wrapping_add(1023) << 52);
        out[i] = scale_pow2 * p;
    }
    out
}

pub(crate) fn exp_scaled(input: &[f64], scale: f64, output: &mut [f64]) {
    debug_assert_eq!(input.len(), output.len());
    if scale == 0.0 {
        output.fill(1.0);
        return;
    }
    #[cfg(all(target_vendor = "apple", not(feature = "portable")))]
    {
        for (out_chunk, in_chunk) in output
            .chunks_mut(i32::MAX as usize)
            .zip(input.chunks(i32::MAX as usize))
        {
            for (value, distance) in out_chunk.iter_mut().zip(in_chunk) {
                *value = scale * distance;
            }
            let count = out_chunk.len() as i32;
            unsafe { vvexp(out_chunk.as_mut_ptr(), out_chunk.as_ptr(), &count) };
        }
    }
    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    {
        let mut in_chunks = input.chunks_exact(8);
        let mut out_chunks = output.chunks_exact_mut(8);
        for (ich, och) in (&mut in_chunks).zip(&mut out_chunks) {
            let z = [
                scale * ich[0],
                scale * ich[1],
                scale * ich[2],
                scale * ich[3],
                scale * ich[4],
                scale * ich[5],
                scale * ich[6],
                scale * ich[7],
            ];
            let fast_ok = (z[0].abs() <= 700.0)
                & (z[1].abs() <= 700.0)
                & (z[2].abs() <= 700.0)
                & (z[3].abs() <= 700.0)
                & (z[4].abs() <= 700.0)
                & (z[5].abs() <= 700.0)
                & (z[6].abs() <= 700.0)
                & (z[7].abs() <= 700.0);
            if fast_ok {
                och.copy_from_slice(&fast_exp_block8(z));
            } else {
                for j in 0..8 {
                    och[j] = z[j].exp();
                }
            }
        }
        for (value, distance) in out_chunks
            .into_remainder()
            .iter_mut()
            .zip(in_chunks.remainder())
        {
            *value = (scale * distance).exp();
        }
    }
}

/// Covariance and basis products use the same row-major convention.
pub(crate) fn multiply(
    a: &[f64],
    b: &[f64],
    c: &mut [f64],
    m: usize,
    k: usize,
    n: usize,
    transpose_a: bool,
) {
    debug_assert!(a.len() >= m * k && b.len() >= k * n && c.len() >= m * n);
    if m == 0 || n == 0 {
        return;
    }
    if k == 0 {
        c[..m * n].fill(0.0);
        return;
    }
    #[cfg(all(target_vendor = "apple", not(feature = "portable")))]
    {
        let lda = if transpose_a { m } else { k };
        if let (Ok(m_i32), Ok(n_i32), Ok(k_i32), Ok(lda_i32)) = (
            i32::try_from(m),
            i32::try_from(n),
            i32::try_from(k),
            i32::try_from(lda),
        ) {
            unsafe {
                cblas_dgemm(
                    101,
                    if transpose_a { 112 } else { 111 },
                    111,
                    m_i32,
                    n_i32,
                    k_i32,
                    1.0,
                    a.as_ptr(),
                    lda_i32,
                    b.as_ptr(),
                    n_i32,
                    0.0,
                    c.as_mut_ptr(),
                    n_i32,
                );
            }
        } else {
            c[..m * n].fill(0.0);
            if transpose_a {
                for l in 0..k {
                    let a_row = &a[l * m..(l + 1) * m];
                    let b_row = &b[l * n..(l + 1) * n];
                    for i in 0..m {
                        let a_val = a_row[i];
                        let c_row = &mut c[i * n..(i + 1) * n];
                        for j in 0..n {
                            c_row[j] += a_val * b_row[j];
                        }
                    }
                }
            } else {
                for i in 0..m {
                    let a_row = &a[i * k..(i + 1) * k];
                    let c_row = &mut c[i * n..(i + 1) * n];
                    for l in 0..k {
                        let a_val = a_row[l];
                        let b_row = &b[l * n..(l + 1) * n];
                        for j in 0..n {
                            c_row[j] += a_val * b_row[j];
                        }
                    }
                }
            }
        }
    }
    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    {
        let a_view = if transpose_a {
            faer::MatRef::from_row_major_slice(&a[..k * m], k, m).transpose()
        } else {
            faer::MatRef::from_row_major_slice(&a[..m * k], m, k)
        };
        let b_view = faer::MatRef::from_row_major_slice(&b[..k * n], k, n);
        let mut c_view = faer::MatMut::from_row_major_slice_mut(&mut c[..m * n], m, n);
        faer::linalg::matmul::matmul(
            &mut c_view,
            faer::Accum::Replace,
            a_view,
            b_view,
            1.0,
            faer_par(),
        );
    }
}

/// The q=1 strict rejection event is (X-cY)(X+cY)>0. The Gaussian
/// quadrant probability below also covers singular covariance, including ties.
/// Scaling standard deviations before squaring keeps unlike units representable.
fn bivariate_rejection(a: f64, b: f64, d: f64, cv: f64) -> Result<f64, CoreError> {
    let tol = 64.0 * f64::EPSILON;
    let error = || CoreError::invalid_input("SCPC covariance is not positive semidefinite");
    let magnitude = a.abs().max(d.abs()).max(b.abs());
    if magnitude == 0.0 {
        return Ok(0.0);
    }
    let (v0, v1) = (a / magnitude + tol, d / magnitude + tol);
    if v0 < 0.0 || v1 < 0.0 || b.abs() / magnitude > v0.sqrt() * v1.sqrt() {
        return Err(error());
    }
    let (a, d) = (a.max(0.0), d.max(0.0));
    if a == 0.0 || d == 0.0 {
        return Ok(if a > 0.0 { 1.0 } else { 0.0 });
    }
    // P = atan2(2 c sqrt(ad-b²), c²d-a) / pi. Independent power-of-two
    // scales keep both products representable; retain their rounding errors
    // before subtracting, since a one-ulp determinant is not an atom.
    let exponent = |v: f64| v.log2().floor().clamp(-1022.0, 1023.0) as i32;
    let (ea, ed) = (exponent(a), exponent(d));
    let eb = (ea + ed).div_euclid(2);
    let an = a / 2.0_f64.powi(ea);
    let dn = d / 2.0_f64.powi(ed) * 2.0_f64.powi((ea + ed).rem_euclid(2));
    let bn = b / 2.0_f64.powi(eb);
    let ad = an * dn;
    let bb = bn * bn;
    let det = if bb > 2.0 * ad {
        0.0
    } else {
        ((ad - bb) + (an.mul_add(dn, -ad) - bn.mul_add(bn, -bb))).max(0.0)
    };
    let t = cv * 2.0_f64.powi(eb - ea);
    let (spread, difference) = if t > 2.0_f64.powi(500) {
        (2.0 * det.sqrt() / t, dn)
    } else {
        let td = t * dn;
        let ttd = t * td;
        let delta = (ttd - an) + (t * t.mul_add(dn, -td) + td.mul_add(t, -ttd));
        (2.0 * t * det.sqrt(), delta)
    };
    if spread == 0.0 {
        return Ok(if difference < 0.0 { 1.0 } else { 0.0 });
    }
    Ok(spread.atan2(difference) / std::f64::consts::PI)
}

/// Certify zero 2x2 minors exactly. Comparing only the rounded products can
/// misclassify a genuinely positive but tiny variance as an atom.
fn equal_products(a: f64, b: f64, c: f64, d: f64) -> bool {
    let p = a * b;
    if p != c * d {
        return false;
    }
    if p == 0.0 {
        return (a == 0.0 || b == 0.0) && (c == 0.0 || d == 0.0);
    }
    p.abs() >= f64::MIN_POSITIVE / f64::EPSILON.powi(2) && a.mul_add(b, -p) == c.mul_add(d, -p)
}

fn one_factor_probability(omega: &[f64], n: usize, width: usize, cv: f64) -> Option<f64> {
    if omega[0] <= 0.0 || omega[0] * omega[width + 1] != omega[width] * omega[width] {
        return None;
    }
    let exp_half = |v: f64| -> i32 {
        if v <= 0.0 || !v.is_finite() {
            0
        } else {
            (v.log2().floor() as i32).div_euclid(2)
        }
    };
    let e0 = exp_half(omega[0]);
    let s0 = 2.0_f64.powi(-e0);
    let first = (omega[0] * s0) * s0;
    if !first.is_finite() || first <= 0.0 {
        return None;
    }
    let mut denom_max_e = i32::MIN;
    for i in 1..n {
        let diag_i = omega[i * width + i];
        if diag_i < 0.0 || !diag_i.is_finite() {
            return None;
        }
        let ei = exp_half(diag_i);
        denom_max_e = denom_max_e.max(ei);
        let si = 2.0_f64.powi(-ei);
        let left = (omega[i * width] * s0) * si;
        for j in 1..=i {
            let diag_j = omega[j * width + j];
            if diag_j < 0.0 || !diag_j.is_finite() {
                return None;
            }
            let ej = exp_half(diag_j);
            let sj = 2.0_f64.powi(-ej);
            let sample = (omega[i * width + j] * si) * sj;
            let right = (omega[j * width] * s0) * sj;
            if !sample.is_finite()
                || !left.is_finite()
                || !right.is_finite()
                || !equal_products(first, sample, left, right)
            {
                return None;
            }
        }
    }
    let sd = 2.0_f64.powi(-denom_max_e);
    let denominator: f64 = (1..n).map(|i| (omega[i * width + i] * sd) * sd).sum();
    if !denominator.is_finite() || denominator <= 0.0 {
        return None;
    }
    let log_thresh = 0.5 * (first.ln() - denominator.ln())
        + ((e0 - denom_max_e) as f64) * std::f64::consts::LN_2;
    Some(if cv == 0.0 || cv.ln() < log_thresh {
        1.0
    } else {
        0.0
    })
}

fn log_sum_exp(iter: impl Iterator<Item = f64> + Clone) -> f64 {
    let max_val = iter.clone().fold(f64::NEG_INFINITY, f64::max);
    if !max_val.is_finite() {
        return max_val;
    }
    max_val + iter.map(|x| (x - max_val).exp()).sum::<f64>().ln()
}

#[inline]
fn log1p_exp(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}

/// Solve `H(u) = s * e^{-u} + sum w_j / (e^u + a_j) - 1 = 0` for `u = ln(lambda_pos) <= 0`.
/// Returns `Ok(None)` if `s == 0` and `alpha^2 >= sum w_j / mu_j` (no positive eigenvalue).
fn solve_lambda_pos_ln(
    active: &[(f64, f64)],
    s: f64,
    two_ln_alpha: f64,
) -> Result<Option<f64>, CoreError> {
    if active.is_empty() {
        return Ok(if s > 0.0 { Some(s.ln()) } else { None });
    }
    if s == 0.0 {
        let ln_g0_over_alpha_sq = log_sum_exp(
            active
                .iter()
                .map(|&(mu, w)| w.ln() - two_ln_alpha - mu.ln()),
        );
        if !ln_g0_over_alpha_sq.is_finite() || ln_g0_over_alpha_sq <= 0.0 {
            return Ok(None);
        }
    }

    let eval_h = |u: f64| -> (f64, f64) {
        let s_term = if s > 0.0 { (s.ln() - u).exp() } else { 0.0 };
        let mut sum_t = 0.0_f64;
        let mut sum_tp = 0.0_f64;
        for &(mu, w) in active {
            let log_a = two_ln_alpha + mu.ln();
            let diff = log_a - u;
            let t = (w.ln() - u - log1p_exp(diff)).exp();
            let p = if diff > 700.0 {
                0.0
            } else if diff < -700.0 {
                1.0
            } else {
                (1.0 + diff.exp()).recip()
            };
            sum_t += t;
            sum_tp += t * p;
        }
        (s_term + sum_t - 1.0, -(s_term + sum_tp))
    };

    let mut u_hi = 0.0_f64;
    let (h_hi, _) = eval_h(u_hi);
    if h_hi >= 0.0 {
        return Ok(Some(0.0));
    }

    let mut u_lo = if s > 0.0 {
        let lo = s.ln();
        let (h_lo, _) = eval_h(lo);
        if h_lo <= 0.0 {
            return Ok(Some(lo));
        }
        lo
    } else {
        let min_log_a = active
            .iter()
            .map(|&(mu, _)| two_ln_alpha + mu.ln())
            .fold(f64::INFINITY, f64::min);
        let floor_u = min_log_a - 64.0;
        let mut u = min_log_a.min(0.0) - 1.0;
        let mut step = 1.0;
        loop {
            let (h_lo, _) = eval_h(u);
            if h_lo > 0.0 {
                break u;
            }
            if u <= floor_u || step > 1e12 {
                return Ok(None);
            }
            u -= step;
            step *= 2.0;
        }
    };

    let mut u = 0.5 * (u_lo + u_hi);
    for _ in 0..128 {
        if (u_hi - u_lo) <= 1e-15 {
            break;
        }
        let (h, dh) = eval_h(u);
        if h == 0.0 {
            return Ok(Some(u));
        }
        if h > 0.0 {
            u_lo = u;
        } else {
            u_hi = u;
        }
        let next = if dh != 0.0 && dh.is_finite() {
            u - h / dh
        } else {
            f64::NAN
        };
        if next.is_finite() && next > u_lo && next < u_hi {
            u = next;
        } else {
            u = 0.5 * (u_lo + u_hi);
        }
    }
    Ok(Some(0.5 * (u_lo + u_hi)))
}

/// Solve the negative root `eta in (a, b)` of
/// `-s / eta + sum_{j} w_j / (mu_j - eta) = alpha^2` with shifted coordinates
/// around `mid = 0.5 * (a + b)` so poles at `a` and `b` never cancel.
fn solve_eta_interval(
    active: &[(f64, f64)],
    s: f64,
    two_ln_alpha: f64,
    a: f64,
    k: usize,
) -> Result<f64, CoreError> {
    let b = active[k].0;
    let log_rhs = two_ln_alpha + b.ln();
    if log_rhs > 100.0 {
        return Ok(b);
    }
    let rhs = if log_rhs < -745.0 { 0.0 } else { log_rhs.exp() };
    let a_hat = a / b;
    let half = 0.5 * (1.0 - a_hat);
    if half <= 0.0 {
        return Ok(b);
    }
    let mid_hat = 0.5 * (a_hat + 1.0);

    let eval_direct = |u: f64| -> f64 {
        let mut val = if s > 0.0 { -s / u } else { 0.0 };
        for &(mu, w) in active {
            val += w / (mu / b - u);
        }
        val - rhs
    };

    if eval_direct(mid_hat) >= 0.0 {
        // Root lies in (a_hat, mid_hat]; parameterize by t = eta_hat - a_hat in (0, half].
        let pole_w = if a == 0.0 { s } else { active[k - 1].1 };
        let eval_near_a = |t: f64| -> (f64, f64, f64) {
            let mut non_pole = if a > 0.0 && s > 0.0 {
                -s / (a_hat + t)
            } else {
                0.0
            };
            let mut deriv = if s > 0.0 {
                let d = a_hat + t;
                s / (d * d)
            } else {
                0.0
            };
            for (j, &(mu, w)) in active.iter().enumerate() {
                if a > 0.0 && j + 1 == k {
                    continue;
                }
                let denom = if j < k {
                    -((a - mu) / b + t)
                } else {
                    (mu - a) / b - t
                };
                non_pole += w / denom;
                deriv += w / (denom * denom);
            }
            let pole = -pole_w / t;
            deriv += pole_w / (t * t);
            (pole + non_pole - rhs, non_pole - rhs, deriv)
        };

        let (_, r_half, _) = eval_near_a(half);
        let mut t_lo = if r_half > 0.0 {
            (0.5 * pole_w / r_half)
                .min(0.5 * half)
                .max(f64::MIN_POSITIVE)
        } else {
            half * 1e-6
        };
        for _ in 0..128 {
            let (f_lo, _, _) = eval_near_a(t_lo);
            if f_lo < 0.0 || t_lo <= f64::MIN_POSITIVE {
                break;
            }
            t_lo *= 0.25;
        }
        let mut t_hi = half;
        let mut t = (0.5 * (t_lo + t_hi)).min(t_lo.sqrt() * t_hi.sqrt().max(t_lo));
        if !(t > t_lo && t < t_hi) {
            t = 0.5 * (t_lo + t_hi);
        }
        for _ in 0..128 {
            if (t_hi - t_lo) <= 1e-15 * (a_hat + t_hi) {
                break;
            }
            let (f, _, df) = eval_near_a(t);
            if f == 0.0 {
                return Ok(b * (a_hat + t));
            }
            if f < 0.0 {
                t_lo = t;
            } else {
                t_hi = t;
            }
            let next = if df > 0.0 && df.is_finite() {
                t - f / df
            } else {
                f64::NAN
            };
            if next.is_finite() && next > t_lo && next < t_hi {
                t = next;
            } else if t_lo > 0.0 && t_hi / t_lo > 16.0 {
                t = t_lo.sqrt() * t_hi.sqrt();
            } else {
                t = 0.5 * (t_lo + t_hi);
            }
        }
        let eta_hat = a_hat + 0.5 * (t_lo + t_hi);
        Ok(b * eta_hat)
    } else {
        // Root lies in (mid_hat, 1.0); parameterize by t = 1.0 - eta_hat in (0, half).
        let pole_w = active[k].1;
        let eval_near_b = |t: f64| -> (f64, f64) {
            let mut val = if s > 0.0 { -s / (1.0 - t) } else { 0.0 };
            let mut deriv = if s > 0.0 {
                let d = 1.0 - t;
                -s / (d * d)
            } else {
                0.0
            };
            for (j, &(mu, w)) in active.iter().enumerate() {
                if j == k {
                    val += pole_w / t;
                    deriv -= pole_w / (t * t);
                } else if j < k {
                    let denom = (1.0 - mu / b) - t;
                    val -= w / denom;
                    deriv -= w / (denom * denom);
                } else {
                    let denom = (mu / b - 1.0) + t;
                    val += w / denom;
                    deriv -= w / (denom * denom);
                }
            }
            (val - rhs, deriv)
        };

        let mut t_lo = 0.0_f64;
        let mut t_hi = half;
        let mut t = if rhs > 0.0 {
            let est = (pole_w / rhs).clamp(1e-16, 0.5 * half);
            est
        } else {
            0.5 * half
        };
        for _ in 0..128 {
            if (t_hi - t_lo) <= 1e-15 {
                break;
            }
            let (f, df) = eval_near_b(t);
            if f == 0.0 {
                return Ok(b * (1.0 - t));
            }
            if f > 0.0 {
                t_lo = t;
            } else {
                t_hi = t;
            }
            let next = if df < 0.0 && df.is_finite() {
                t - f / df
            } else {
                f64::NAN
            };
            if next.is_finite() && next > t_lo && next < t_hi {
                t = next;
            } else {
                t = 0.5 * (t_lo + t_hi);
            }
        }
        Ok(b * (1.0 - 0.5 * (t_lo + t_hi)))
    }
}

pub(crate) struct EigenWorkspace {
    n: usize,
    matrix: Vec<f64>,
    auxiliary: Vec<f64>,
    scratch: Vec<f64>,
    coupling: Vec<f64>,
    terms: Vec<(f64, f64)>,
    active: Vec<(f64, f64)>,
    deflated: Vec<f64>,
}

impl EigenWorkspace {
    pub(crate) fn new(n: usize) -> Self {
        let q = n.saturating_sub(1);
        Self {
            n,
            matrix: vec![0.0; q * q],
            auxiliary: vec![0.0; q],
            scratch: vec![0.0; q * q],
            coupling: vec![0.0; q],
            terms: Vec::with_capacity(q),
            active: Vec::with_capacity(q),
            deflated: Vec::with_capacity(q),
        }
    }

    pub(crate) fn singular_or_bivariate_probability(
        &self,
        omega: &[f64],
        width: usize,
        cv: f64,
    ) -> Result<Option<f64>, CoreError> {
        let n = self.n;
        let mut active = (1..n).filter(|&i| omega[i * width + i] != 0.0);
        let j = active.next().unwrap_or(1);
        if n > 2
            && (active.next().is_some()
                || (1..n).any(|i| {
                    i != j
                        && (0..n)
                            .any(|k| omega[i * width + k] != 0.0 || omega[k * width + i] != 0.0)
                }))
        {
            // A strict-boundary atom requires a single common Gaussian factor.
            let probability = one_factor_probability(omega, n, width, cv);
            if probability.is_some() {
                for i in 0..n {
                    for j in 0..=i {
                        if !omega[i * width + j].is_finite()
                            || omega[i * width + j] != omega[j * width + i]
                        {
                            return Ok(None);
                        }
                    }
                }
            }
            return Ok(probability);
        }
        // An identically zero denominator coordinate does not change the statistic.
        let (a, b, transpose, d) = (omega[0], omega[j], omega[j * width], omega[j * width + j]);
        if ![a, b, transpose, d].iter().all(|v| v.is_finite()) {
            return Err(CoreError::invalid_input(
                "SCPC covariance contains non-finite values",
            ));
        }
        if (b - transpose).abs() > 64.0 * f64::EPSILON * a.abs().max(d.abs()).max(b.abs()) {
            return Err(CoreError::invalid_input("SCPC covariance is not symmetric"));
        }
        bivariate_rejection(a, b.midpoint(transpose), d, cv).map(Some)
    }

    /// Equilibrate Omega into a 2-block Schur complement / DPR1 problem and return
    /// either an exact closed-form probability (`Ok(Some(p))`) or the dimensionless
    /// eigenvalue ratios and log-ratios needed for tail-stable integration (`Ok(None)`).
    pub(crate) fn rejection_eigen_ratios(
        &mut self,
        omega: &[f64],
        width: usize,
        cv: f64,
        ratios: &mut Vec<f64>,
        logs: &mut Vec<f64>,
    ) -> Result<Option<f64>, CoreError> {
        ratios.clear();
        logs.clear();
        if !cv.is_finite() {
            return Err(CoreError::invalid_input(
                "SCPC critical value must be finite",
            ));
        }
        let n = self.n;
        if n == 0 {
            return Ok(Some(0.0));
        }
        for i in 0..n {
            for j in 0..n {
                if !omega[i * width + j].is_finite() {
                    return Err(CoreError::invalid_input(
                        "SCPC covariance contains non-finite values",
                    ));
                }
            }
        }

        let raw_pos_diag = omega[0];
        let raw_neg_max_diag = (1..n).map(|i| omega[i * width + i]).fold(0.0_f64, f64::max);
        let overall_max_diag = raw_pos_diag.max(raw_neg_max_diag);
        let global_tol = 32.0 * (n as f64) * f64::EPSILON * overall_max_diag;

        for i in 0..n {
            for j in 0..i {
                if (omega[i * width + j] - omega[j * width + i]).abs() > global_tol {
                    return Err(CoreError::invalid_input("SCPC covariance is not symmetric"));
                }
            }
        }
        if overall_max_diag <= 0.0 {
            for i in 0..n {
                for j in 0..=i {
                    if omega[i * width + j] != 0.0 {
                        return Err(CoreError::invalid_input(
                            "SCPC covariance is not positive semidefinite",
                        ));
                    }
                }
            }
            return Ok(Some(0.0));
        }

        if raw_pos_diag < -global_tol || (1..n).any(|i| omega[i * width + i] < -global_tol) {
            return Err(CoreError::invalid_input(
                "SCPC covariance is not positive semidefinite",
            ));
        }

        let pos_diag = raw_pos_diag.max(0.0);
        let neg_max_diag = raw_neg_max_diag.max(0.0);
        let base_tol = 32.0 * (n as f64) * f64::EPSILON;

        if neg_max_diag == 0.0 {
            for i in 1..n {
                for j in 0..=i {
                    let sym = omega[i * width + j].midpoint(omega[j * width + i]);
                    if sym.abs() > global_tol {
                        return Err(CoreError::invalid_input(
                            "SCPC covariance is not positive semidefinite",
                        ));
                    }
                }
            }
            if cv < 0.0 {
                return Ok(Some(1.0));
            }
            return Ok(Some(if pos_diag > 0.0 { 1.0 } else { 0.0 }));
        }

        let q = n - 1;
        if n == 2 {
            if cv < 0.0 {
                return Ok(Some(1.0));
            }
            let sym01 = omega[1].midpoint(omega[width]);
            return bivariate_rejection(pos_diag, sym01, omega[width + 1].max(0.0), cv).map(Some);
        }

        // Fast exact-atom checks when cv >= 0 and all entries are strictly valid.
        if cv < 0.0 {
            // Still validate PSD below when n >= 3 before short-circuiting.
        } else if let Some(p) = self.singular_or_bivariate_probability(omega, width, cv)? {
            return Ok(Some(p));
        }

        let s_neg = neg_max_diag.sqrt();
        for i in 0..q {
            for j in 0..=i {
                let sym =
                    omega[(i + 1) * width + (j + 1)].midpoint(omega[(j + 1) * width + (i + 1)]);
                let val = (sym / s_neg) / s_neg;
                if !val.is_finite() {
                    return Err(CoreError::invalid_input(
                        "SCPC covariance normalization overflowed",
                    ));
                }
                self.matrix[i * q + j] = val;
                self.matrix[j * q + i] = val;
            }
        }
        symmetric_eigen(
            &mut self.matrix[..q * q],
            &mut self.auxiliary[..q],
            Some(&mut self.scratch[..q * q]),
        )?;

        for &mu in &self.auxiliary[..q] {
            if mu < -base_tol && (mu * s_neg) * s_neg < -global_tol {
                return Err(CoreError::invalid_input(
                    "SCPC covariance is not positive semidefinite",
                ));
            }
        }

        if pos_diag == 0.0 {
            for j in 0..q {
                let sym0 = omega[j + 1].midpoint(omega[(j + 1) * width]);
                if sym0.abs() > global_tol {
                    return Err(CoreError::invalid_input(
                        "SCPC covariance is not positive semidefinite",
                    ));
                }
            }
            return Ok(Some(if cv < 0.0 { 1.0 } else { 0.0 }));
        }

        let s0 = pos_diag.sqrt();
        for i in 0..q {
            let sym = omega[(i + 1) * width].midpoint(omega[i + 1]);
            let r_i = (sym / s0) / s_neg;
            if !r_i.is_finite() {
                return Err(CoreError::invalid_input(
                    "SCPC covariance normalization overflowed",
                ));
            }
            self.coupling[i] = r_i;
        }

        self.terms.clear();
        let f1 = 1.0 + global_tol / pos_diag;
        for j in 0..q {
            let mu = self.auxiliary[j];
            let proj: f64 = (0..q)
                .map(|k| self.scratch[j * q + k] * self.coupling[k])
                .sum();
            let proj_sq = proj * proj;
            if mu <= base_tol {
                let f2 = mu.max(0.0) + global_tol / neg_max_diag;
                if f1.is_finite() && f2.is_finite() && proj_sq > f1 * f2 {
                    return Err(CoreError::invalid_input(
                        "SCPC covariance is not positive semidefinite",
                    ));
                }
                continue;
            }
            let w = proj_sq / mu;
            if w > 1.0 + base_tol {
                return Err(CoreError::invalid_input(
                    "SCPC covariance is not positive semidefinite",
                ));
            }
            self.terms.push((mu, w));
        }

        let w_sum: f64 = self.terms.iter().map(|&(_, w)| w).sum();
        if w_sum > 1.0 + base_tol {
            return Err(CoreError::invalid_input(
                "SCPC covariance is not positive semidefinite",
            ));
        }
        let mut s = (1.0 - w_sum).max(0.0);
        if w_sum > 1.0 - 1e-12 && w_sum <= 1.0 + base_tol {
            // Recompute Schur complement without sqrt(pos_diag) rounding when w_sum ~ 1.0
            let mut b_quad = 0.0_f64;
            for j in 0..q {
                let mu = self.auxiliary[j];
                if mu <= base_tol {
                    continue;
                }
                let raw_proj: f64 = (0..q)
                    .map(|k| {
                        let sym = omega[(k + 1) * width].midpoint(omega[k + 1]);
                        self.scratch[j * q + k] * (sym / s_neg)
                    })
                    .sum();
                b_quad += (raw_proj * raw_proj) / mu;
            }
            if b_quad.is_finite() {
                let diff = pos_diag - b_quad;
                if diff > 0.0 {
                    s = (diff / pos_diag).min(1.0);
                } else if diff / pos_diag >= -base_tol {
                    s = 0.0;
                }
            }
        }
        if w_sum > 0.0 && (w_sum > 1.0 || w_sum > 1.0 - s) {
            let target = if w_sum > 1.0 { 1.0 } else { 1.0 - s };
            let scale = target / w_sum;
            for (_, w) in &mut self.terms {
                *w *= scale;
            }
        }

        if cv <= 0.0 {
            return Ok(Some(1.0));
        }

        self.terms.sort_by(|a, b| a.0.total_cmp(&b.0));
        self.deflated.clear();
        self.active.clear();
        for &(mu, w) in &self.terms {
            if w <= 1e-30 {
                self.deflated.push(mu);
            } else if let Some(last) = self.active.last_mut() {
                if mu - last.0 <= 32.0 * f64::EPSILON * mu {
                    let w_tot = last.1 + w;
                    let mu_merged = (last.1 * last.0 + w * mu) / w_tot;
                    let mu_defl = (w * last.0 + last.1 * mu) / w_tot;
                    last.0 = mu_merged;
                    last.1 = w_tot;
                    self.deflated.push(mu_defl);
                } else {
                    self.active.push((mu, w));
                }
            } else {
                self.active.push((mu, w));
            }
        }

        let two_ln_alpha = 2.0 * (cv.ln() + 0.5 * (neg_max_diag.ln() - pos_diag.ln()));
        let Some(ln_lambda_pos) = solve_lambda_pos_ln(&self.active, s, two_ln_alpha)? else {
            return Ok(Some(0.0));
        };

        if !self.active.is_empty() {
            let m = self.active.len();
            if s > 0.0 {
                let eta0 = solve_eta_interval(&self.active, s, two_ln_alpha, 0.0, 0)?;
                if eta0 > 0.0 {
                    let log_r = two_ln_alpha + eta0.ln() - ln_lambda_pos;
                    logs.push(log_r);
                    ratios.push(log_r.exp());
                }
            }
            for k in 1..m {
                let a = self.active[k - 1].0;
                let eta_k = solve_eta_interval(&self.active, s, two_ln_alpha, a, k)?;
                if eta_k > 0.0 {
                    let log_r = two_ln_alpha + eta_k.ln() - ln_lambda_pos;
                    logs.push(log_r);
                    ratios.push(log_r.exp());
                }
            }
        }
        for &eta in &self.deflated {
            if eta > 0.0 {
                let log_r = two_ln_alpha + eta.ln() - ln_lambda_pos;
                logs.push(log_r);
                ratios.push(log_r.exp());
            }
        }

        if logs.is_empty() {
            return Ok(Some(1.0));
        }

        Ok(None)
    }
}

pub(crate) fn symmetric_eigen(
    a: &mut [f64],
    values: &mut [f64],
    mut vectors: Option<&mut [f64]>,
) -> Result<(), CoreError> {
    let n = values.len();
    if let Some(vectors) = vectors.as_deref_mut() {
        vectors.fill(0.0);
        for i in 0..n {
            vectors[i * n + i] = 1.0;
        }
    }
    for _ in 0..64 {
        let mut changed = false;
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p * n + q];
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                if apq == 0.0
                    || apq.abs() <= 4.0 * f64::EPSILON * (app.abs().sqrt() * aqq.abs().sqrt())
                {
                    continue;
                }
                changed = true;
                let difference = (aqq - app) * 0.5;
                let tangent = if difference == 0.0 {
                    1.0
                } else {
                    apq / (difference + difference.signum() * difference.hypot(apq))
                };
                let cosine = 1.0 / (1.0 + tangent * tangent).sqrt();
                let sine = tangent * cosine;
                a[p * n + p] = app - tangent * apq;
                a[q * n + q] = aqq + tangent * apq;
                a[p * n + q] = 0.0;
                a[q * n + p] = 0.0;
                for k in (0..n).filter(|&k| k != p && k != q) {
                    let left = a[k * n + p];
                    let right = a[k * n + q];
                    a[k * n + p] = cosine * left - sine * right;
                    a[k * n + q] = sine * left + cosine * right;
                    a[p * n + k] = a[k * n + p];
                    a[q * n + k] = a[k * n + q];
                }
                if let Some(vectors) = vectors.as_deref_mut() {
                    for k in 0..n {
                        let left = vectors[p * n + k];
                        let right = vectors[q * n + k];
                        vectors[p * n + k] = cosine * left - sine * right;
                        vectors[q * n + k] = sine * left + cosine * right;
                    }
                }
            }
        }
        if !changed {
            for (i, value) in values.iter_mut().enumerate() {
                *value = a[i * n + i];
            }
            return Ok(());
        }
    }
    Err(CoreError::internal(
        "SCPC covariance eigensolver failed to converge",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiply_respects_row_major_layout_and_transpose() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let b = [7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
        let mut product = [0.0; 4];
        multiply(&a, &b, &mut product, 2, 3, 2, false);
        assert_eq!(product, [58.0, 64.0, 139.0, 154.0]);

        let mut transposed = [0.0; 9];
        multiply(&a, &b, &mut transposed, 3, 2, 3, true);
        assert_eq!(
            transposed,
            [47.0, 52.0, 57.0, 64.0, 71.0, 78.0, 81.0, 90.0, 99.0]
        );

        let mut product = [1.0; 6];
        multiply(&[], &[], &mut product, 2, 0, 3, false);
        assert_eq!(product, [0.0; 6]);
    }

    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    #[test]
    fn exp_scaled_matches_libm_across_range_and_boundaries() {
        let mut xs = Vec::new();
        for i in 0..=2000 {
            xs.push(-700.0 + 1400.0 * (i as f64 / 2000.0));
        }
        xs.extend_from_slice(&[
            -700.0,
            -700.1,
            -708.0,
            -750.0,
            0.0,
            -0.0,
            700.0,
            700.1,
            710.0,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ]);
        let mut out = vec![0.0; xs.len()];
        exp_scaled(&xs, 1.0, &mut out);
        for (i, (&x, &got)) in xs.iter().zip(&out).enumerate() {
            let expected = x.exp();
            if expected.is_nan() {
                assert!(got.is_nan(), "idx {i}");
            } else if expected.is_infinite() {
                assert_eq!(got, expected, "idx {i}");
            } else if expected == 0.0 {
                assert_eq!(got, 0.0, "idx {i}");
            } else {
                let rel = (got - expected).abs() / expected;
                assert!(
                    rel < 2e-15,
                    "idx {i} x={x} got={got} expected={expected} rel={rel}"
                );
            }
        }
        let mut scaled = vec![0.0; xs.len()];
        exp_scaled(&xs, -0.5, &mut scaled);
        for (i, (&x, &got)) in xs.iter().zip(&scaled).enumerate() {
            let expected = (-0.5 * x).exp();
            if expected.is_finite() && expected > 0.0 {
                let rel = (got - expected).abs() / expected;
                assert!(rel < 2e-15, "scaled idx {i} x={x} rel={rel}");
            }
        }
    }

    #[test]
    fn rejection_eigen_ratios_handles_zero_cv_scale_disparity_and_indefiniteness() {
        let mut ratios = Vec::new();
        let mut logs = Vec::new();

        // 1. Scale disparity IV and FE blocks produce non-trivial pos and neg modes.
        let cov_iv = [1.0, 5e6, 5e6, 1e14];
        let mut ws2 = EigenWorkspace::new(2);
        assert!(matches!(
            ws2.rejection_eigen_ratios(&cov_iv, 2, 1.0, &mut ratios, &mut logs),
            Ok(Some(p)) if p > 0.0 && p < 1.0
        ));
        assert!(matches!(
            ws2.rejection_eigen_ratios(&[1.0, 0.0, 0.0, 1e14], 2, 1.0, &mut ratios, &mut logs),
            Ok(Some(p)) if p > 0.0 && p < 1e-6
        ));
        assert!(matches!(
            ws2.rejection_eigen_ratios(&[1e14, 0.0, 0.0, 1.0], 2, 1.0, &mut ratios, &mut logs),
            Ok(Some(p)) if p > 0.999999 && p < 1.0
        ));

        // 2. Correlated 10x10 SPD matrix at cv == 0.0 must return exact 1.0 (zero Jacobi roundoff).
        let n10 = 10usize;
        let mut omega10 = vec![0.0_f64; n10 * n10];
        for i in 0..n10 {
            for j in 0..n10 {
                omega10[i * n10 + j] = (-0.35 * ((i as f64) - (j as f64)).abs()).exp();
            }
        }
        let mut ws10 = EigenWorkspace::new(n10);
        let p0 = ws10
            .rejection_eigen_ratios(&omega10, n10, 0.0, &mut ratios, &mut logs)
            .unwrap();
        assert_eq!(p0, Some(1.0));

        // 3. Extreme scale disparity diagonal Diag(1e200, 1e-120, 1e-120) at cv = 1e160 -> ratios ~ 1.0.
        let mut ws3 = EigenWorkspace::new(3);
        let res_extreme = ws3
            .rejection_eigen_ratios(
                &[1e200, 0.0, 0.0, 0.0, 1e-120, 0.0, 0.0, 0.0, 1e-120],
                3,
                1e160,
                &mut ratios,
                &mut logs,
            )
            .unwrap();
        assert!(res_extreme.is_none());
        assert_eq!(ratios.len(), 2);
        for &r in &ratios {
            assert!((r - 1.0).abs() < 1e-12, "unexpected ratio {r}");
        }

        // 4. Rank-1 extreme scale disparity v = (2^400, 2^-400, 2^-400).
        let hi = 2.0_f64.powi(400);
        let lo = 2.0_f64.powi(-400);
        let v = [hi, lo, lo];
        let mut r1 = [0.0_f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                r1[i * 3 + j] = v[i] * v[j];
            }
        }
        let cv_thresh = (hi / lo) / 2.0_f64.sqrt();
        assert_eq!(
            ws3.rejection_eigen_ratios(&r1, 3, 0.5 * cv_thresh, &mut ratios, &mut logs)
                .unwrap(),
            Some(1.0)
        );
        assert_eq!(
            ws3.rejection_eigen_ratios(&r1, 3, 2.0 * cv_thresh, &mut ratios, &mut logs)
                .unwrap(),
            Some(0.0)
        );

        // 5. Validation guards: non-finite, asymmetric, indefinite, roundoff.
        assert!(ws2
            .rejection_eigen_ratios(
                &[f64::INFINITY, 0.0, 0.0, 1.0],
                2,
                1.0,
                &mut ratios,
                &mut logs
            )
            .is_err());
        assert!(ws2
            .rejection_eigen_ratios(&[1.0, 0.5, 0.0, 1.0], 2, 1.0, &mut ratios, &mut logs)
            .is_err());
        assert!(ws2
            .rejection_eigen_ratios(&[0.0, 1.0, 1.0, 0.0], 2, 1.0, &mut ratios, &mut logs)
            .is_err());

        let eps = f64::EPSILON;
        let roundoff = [1.0, 0.0, 0.0, 0.0, eps, 16.0 * eps, 0.0, 16.0 * eps, 0.0];
        let ro_res = ws3
            .rejection_eigen_ratios(&roundoff, 3, 1.0, &mut ratios, &mut logs)
            .unwrap();
        assert!(ro_res.is_none() || ro_res == Some(1.0));
        assert!(ratios.iter().all(|&r| r < 1e-12));

        let indefinite = [1.0, 0.0, 0.0, 0.0, eps, 1e-3, 0.0, 1e-3, 0.0];
        assert!(ws3
            .rejection_eigen_ratios(&indefinite, 3, 1.0, &mut ratios, &mut logs)
            .is_err());
        assert!(ws3
            .rejection_eigen_ratios(
                &[1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 2.0, 1.0],
                3,
                0.0,
                &mut ratios,
                &mut logs,
            )
            .is_err());
    }

    #[test]
    fn dense_omega_eigen_ratios_match_direct_decomposition() {
        // Test non-trivial 4x4 dense covariance against direct Omega^{1/2} D Omega^{1/2} eigen.
        let n = 4usize;
        let coords = [0.0_f64, 0.3, 1.1, 2.4];
        let mut omega = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                omega[i * n + j] = (-0.7 * (coords[i] - coords[j]).abs()).exp();
            }
        }
        let mut a = omega.clone();
        let mut vals = vec![0.0; n];
        let mut vecs = vec![0.0; n * n];
        symmetric_eigen(&mut a, &mut vals, Some(&mut vecs)).unwrap();
        let cv = 1.5_f64;
        let mut mat = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let inner: f64 = (1..n).map(|k| vecs[i * n + k] * vecs[j * n + k]).sum();
                let v =
                    (vecs[i * n] * vecs[j * n] - cv * cv * inner) * vals[i].sqrt() * vals[j].sqrt();
                mat[i * n + j] = v;
                mat[j * n + i] = v;
            }
        }
        let mut evs = vec![0.0; n];
        symmetric_eigen(&mut mat, &mut evs, None).unwrap();
        let max_p = evs.iter().copied().fold(0.0_f64, f64::max);
        let mut ref_r: Vec<f64> = evs
            .iter()
            .filter(|&&v| v < 0.0)
            .map(|&v| -v / max_p)
            .collect();
        ref_r.sort_by(|a, b| a.total_cmp(b));

        let mut ws = EigenWorkspace::new(n);
        let mut ratios = Vec::new();
        let mut logs = Vec::new();
        ws.rejection_eigen_ratios(&omega, n, cv, &mut ratios, &mut logs)
            .unwrap();
        ratios.sort_by(|a, b| a.total_cmp(b));
        assert_eq!(ratios.len(), ref_r.len());
        for (a, b) in ratios.iter().zip(&ref_r) {
            assert!((a - b).abs() < 1e-11, "got={ratios:?} vs ref={ref_r:?}");
        }
    }
}
