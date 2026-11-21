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
        for (value, distance) in output.iter_mut().zip(input) {
            *value = scale * distance;
        }
        unsafe { vvexp(output.as_mut_ptr(), output.as_ptr(), &(output.len() as i32)) };
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
        unsafe {
            cblas_dgemm(
                101,
                if transpose_a { 112 } else { 111 },
                111,
                m as i32,
                n as i32,
                k as i32,
                1.0,
                a.as_ptr(),
                lda as i32,
                b.as_ptr(),
                n as i32,
                0.0,
                c.as_mut_ptr(),
                n as i32,
            );
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
    let magnitude = (0..n).map(|i| omega[i * width + i]).fold(0.0, f64::max);
    let exponent = magnitude.log2().floor().clamp(-1022.0, 1023.0) as i32;
    let scale = 2.0_f64.powi(exponent);
    let first = omega[0] / scale;
    for i in 1..n {
        for j in 1..=i {
            let sample = omega[i * width + j] / scale;
            let left = omega[i * width] / scale;
            let right = omega[j * width] / scale;
            if sample * scale != omega[i * width + j]
                || left * scale != omega[i * width]
                || right * scale != omega[j * width]
                || first * scale != omega[0]
                || !equal_products(first, sample, left, right)
            {
                return None;
            }
        }
    }
    let denominator: f64 = (1..n).map(|i| omega[i * width + i] / scale).sum();
    let positive = if cv > 1.0 {
        (first / cv) / cv > denominator
    } else {
        first > (denominator * cv) * cv
    };
    Some(if positive { 1.0 } else { 0.0 })
}

pub(crate) struct EigenWorkspace {
    matrix: Vec<f64>,
    values: Vec<f64>,
    auxiliary: Vec<f64>,
    scratch: Vec<f64>,
    permutation: Vec<usize>,
    diag_scale: Vec<f64>,
    corr: Vec<f64>,
}

impl EigenWorkspace {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            matrix: vec![0.0; n * n],
            values: vec![0.0; n],
            auxiliary: vec![0.0; n],
            scratch: vec![0.0; n * n],
            permutation: (0..n).collect(),
            diag_scale: vec![0.0; n],
            corr: vec![0.0; n * n],
        }
    }

    pub(crate) fn singular_or_bivariate_probability(
        &self,
        omega: &[f64],
        width: usize,
        cv: f64,
    ) -> Result<Option<f64>, CoreError> {
        let n = self.values.len();
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
        bivariate_rejection(a, b + 0.5 * (transpose - b), d, cv).map(Some)
    }

    /// Factor Omega = L L' before diagonalizing the symmetric L' D L.
    pub(crate) fn weighted_eigenvalues(
        &mut self,
        omega: &[f64],
        width: usize,
        positive: f64,
        negative: f64,
    ) -> Result<&[f64], CoreError> {
        let n = self.values.len();
        if n == 0 {
            return Ok(&self.values);
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

        let pos_diag = omega[0];
        let neg_max_diag = (1..n).map(|i| omega[i * width + i]).fold(0.0_f64, f64::max);
        let overall_max_diag = pos_diag.max(neg_max_diag);
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
            self.values.fill(0.0);
            return Ok(&self.values);
        }

        if pos_diag < -global_tol || (1..n).any(|i| omega[i * width + i] < -global_tol) {
            return Err(CoreError::invalid_input(
                "SCPC covariance is not positive semidefinite",
            ));
        }

        let s0 = if pos_diag > 0.0 { pos_diag.sqrt() } else { 0.0 };
        let s_neg = if neg_max_diag > 0.0 {
            neg_max_diag.sqrt()
        } else {
            0.0
        };

        self.diag_scale[0] = s0;
        for i in 1..n {
            self.diag_scale[i] = s_neg;
        }

        let neg_tol = 32.0 * (n as f64) * f64::EPSILON;
        self.corr.fill(0.0);
        let mut build_ok = true;
        for i in 0..n {
            let si = self.diag_scale[i];
            if si == 0.0 {
                for j in 0..=i {
                    if omega[i * width + j].abs() > global_tol {
                        build_ok = false;
                        break;
                    }
                }
                continue;
            }
            for j in 0..=i {
                let sj = self.diag_scale[j];
                let entry = omega[i * width + j];
                if sj == 0.0 {
                    if entry.abs() > global_tol {
                        build_ok = false;
                        break;
                    }
                } else {
                    let c_ij = (entry / si) / sj;
                    if !c_ij.is_finite() {
                        build_ok = false;
                        break;
                    }
                    self.corr[i * n + j] = c_ij;
                    self.corr[j * n + i] = c_ij;
                }
            }
        }
        if !build_ok {
            return self.repair_roundoff(omega, width, positive, negative, overall_max_diag);
        }

        let run_pivoted_chol = |corr: &[f64],
                                scratch: &mut [f64],
                                permutation: &mut [usize],
                                auxiliary: &mut [f64]|
         -> Result<usize, ()> {
            scratch.fill(0.0);
            for i in 0..n {
                permutation[i] = i;
                auxiliary[i] = corr[i * n + i];
            }
            let mut rank = 0;
            for k in 0..n {
                let p = (k..n)
                    .max_by(|&a, &b| auxiliary[a].total_cmp(&auxiliary[b]))
                    .expect("pivot range is nonempty");
                if auxiliary[p] <= 0.0 {
                    let negligible = auxiliary[k..].iter().all(|value| value.abs() <= neg_tol)
                        && (k..n).all(|i| {
                            (k..i).all(|j| {
                                let product: f64 = (0..k)
                                    .map(|l| scratch[i * n + l] * scratch[j * n + l])
                                    .sum();
                                let residual = corr[permutation[i] * n + permutation[j]] - product;
                                residual.abs() <= neg_tol
                            })
                        });
                    if !negligible {
                        return Err(());
                    }
                    break;
                }
                permutation.swap(k, p);
                auxiliary.swap(k, p);
                for j in 0..k {
                    scratch.swap(k * n + j, p * n + j);
                }
                let root = auxiliary[k].sqrt();
                scratch[k * n + k] = root;
                for i in k + 1..n {
                    let product: f64 = (0..k)
                        .map(|j| scratch[i * n + j] * scratch[k * n + j])
                        .sum();
                    let value = (corr[permutation[i] * n + permutation[k]] - product) / root;
                    scratch[i * n + k] = value;
                    auxiliary[i] -= value * value;
                }
                rank += 1;
            }
            Ok(rank)
        };

        let rank = match run_pivoted_chol(
            &self.corr,
            &mut self.scratch,
            &mut self.permutation,
            &mut self.auxiliary,
        ) {
            Ok(r) => r,
            Err(()) => {
                return self.repair_roundoff(omega, width, positive, negative, overall_max_diag)
            }
        };

        if rank == 1 && n > 2 {
            // The exact one-factor boundary was checked above. A second
            // representable direction may be lost in a rounded Cholesky pivot.
            return self.repair_roundoff(omega, width, positive, negative, overall_max_diag);
        }

        let s0_w = if self.diag_scale[0] > 0.0 {
            self.diag_scale[0] * positive.abs().sqrt()
        } else {
            0.0
        };
        let sneg_w = if s_neg > 0.0 {
            s_neg * negative.abs().sqrt()
        } else {
            0.0
        };
        let w_scale = s0_w.max(sneg_w);
        if w_scale == 0.0 {
            self.values.fill(0.0);
            return Ok(&self.values);
        }
        let pos_norm = if self.diag_scale[0] > 0.0 {
            let r = s0_w / w_scale;
            r * r * positive.signum()
        } else {
            0.0
        };
        let neg_norm = if s_neg > 0.0 {
            let r = sneg_w / w_scale;
            r * r * negative.signum()
        } else {
            0.0
        };

        self.matrix.fill(0.0);
        for i in 0..rank {
            for j in i..rank {
                let value = (j..n)
                    .map(|k| {
                        let weight = if self.permutation[k] == 0 {
                            pos_norm
                        } else {
                            neg_norm
                        };
                        self.scratch[k * n + i] * weight * self.scratch[k * n + j]
                    })
                    .sum();
                self.matrix[i * n + j] = value;
                self.matrix[j * n + i] = value;
            }
        }
        symmetric_eigen(&mut self.matrix, &mut self.values, None)?;
        Ok(&self.values)
    }

    // A PSD block can acquire a negative eigenvalue at roundoff scale in the
    // covariance multiplication. Judge that residual in the original units.
    fn repair_roundoff(
        &mut self,
        omega: &[f64],
        width: usize,
        positive: f64,
        negative: f64,
        scale: f64,
    ) -> Result<&[f64], CoreError> {
        let n = self.values.len();
        for i in 0..n {
            for j in 0..n {
                self.corr[i * n + j] = omega[i * width + j] / scale;
            }
        }
        symmetric_eigen(&mut self.corr, &mut self.values, Some(&mut self.scratch))?;
        if self
            .values
            .iter()
            .any(|&v| v < -32.0 * n as f64 * f64::EPSILON)
        {
            return Err(CoreError::invalid_input(
                "SCPC covariance is not positive semidefinite",
            ));
        }
        for value in &mut self.values {
            *value = value.max(0.0).sqrt();
        }
        // The Jacobi solver stores one eigenvector per row.
        for i in 0..n {
            for j in i..n {
                let inner: f64 = (1..n)
                    .map(|k| self.scratch[i * n + k] * self.scratch[j * n + k])
                    .sum();
                let v = (positive * self.scratch[i * n] * self.scratch[j * n] + negative * inner)
                    * self.values[i]
                    * self.values[j];
                self.matrix[i * n + j] = v;
                self.matrix[j * n + i] = v;
            }
        }
        symmetric_eigen(&mut self.matrix, &mut self.values, None)?;
        Ok(&self.values)
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
    fn covariance_eigenfactor_preserves_weighted_quadratic_form() {
        let omega = [1.0, 0.3, 0.2, 0.3, 0.6, 0.4, 0.2, 0.4, 0.9];
        let weights = [0.7, -0.4, -0.4];
        let mut ws = EigenWorkspace::new(3);
        let eig = ws.repair_roundoff(&omega, 3, weights[0], weights[1], 1.0).unwrap();
        let mut trace = 0.0;
        let mut trace_squared = 0.0;
        for i in 0..3 {
            trace += weights[i] * omega[i * 3 + i];
            for j in 0..3 {
                trace_squared += weights[i] * weights[j] * omega[i * 3 + j].powi(2);
            }
        }
        assert!((eig.iter().sum::<f64>() - trace).abs() < 1e-13);
        assert!((eig.iter().map(|v| v * v).sum::<f64>() - trace_squared).abs() < 1e-13);
    }

    #[test]
    fn weighted_eigenvalues_handles_scale_disparity_roundoff_and_indefiniteness() {
        let cov_iv = [1.0, 5e6, 5e6, 1e14];
        let mut ws2 = EigenWorkspace::new(2);
        let eigs = ws2
            .weighted_eigenvalues(&cov_iv, 2, 1.0, -1.0)
            .unwrap()
            .to_vec();
        let max_pos = eigs.iter().copied().fold(0.0, f64::max);
        assert!(max_pos > 0.0, "scaled IV dropped positive block: {eigs:?}");
        assert!(
            eigs.iter().any(|&v| v < 0.0),
            "scaled IV dropped negative block: {eigs:?}"
        );

        let cov_iv_diag = [1.0, 0.0, 0.0, 1e14];
        let eigs_d = ws2
            .weighted_eigenvalues(&cov_iv_diag, 2, 1.0, -1.0)
            .unwrap()
            .to_vec();
        assert!(eigs_d.iter().any(|&v| v > 0.0) && eigs_d.iter().any(|&v| v < 0.0));

        let cov_fe = [1e14, 0.0, 0.0, 1.0];
        let eigs_fe = ws2
            .weighted_eigenvalues(&cov_fe, 2, 1.0, -1.0)
            .unwrap()
            .to_vec();
        assert!(eigs_fe.iter().any(|&v| v > 0.0) && eigs_fe.iter().any(|&v| v < 0.0));

        assert!(ws2
            .weighted_eigenvalues(&[f64::INFINITY, 0.0, 0.0, -1e100], 2, 1.0, -1.0)
            .is_err());
        assert!(ws2
            .weighted_eigenvalues(&[0.0, 1.0, 1.0, 0.0], 2, 1.0, -1.0)
            .is_err());

        let eps = f64::EPSILON;
        let mut ws3 = EigenWorkspace::new(3);
        let roundoff = [1.0, 0.0, 0.0, 0.0, eps, 16.0 * eps, 0.0, 16.0 * eps, 0.0];
        let r_eigs = ws3
            .weighted_eigenvalues(&roundoff, 3, 1.0, -1.0)
            .unwrap()
            .to_vec();
        assert!(
            (r_eigs[0] - 1.0).abs() < 1e-12 || r_eigs.iter().any(|&v| (v - 1.0).abs() < 1e-12)
        );
        assert_eq!(r_eigs.iter().filter(|&&v| v.abs() < 1e-12).count(), 2);

        let indefinite = [1.0, 0.0, 0.0, 0.0, eps, 1e-3, 0.0, 1e-3, 0.0];
        assert!(ws3
            .weighted_eigenvalues(&indefinite, 3, 1.0, -1.0)
            .is_err());
    }
}
