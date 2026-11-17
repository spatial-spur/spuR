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
    fn dgeev_(
        job_left: *const u8,
        job_right: *const u8,
        n: *const i32,
        matrix: *mut f64,
        lda: *const i32,
        real: *mut f64,
        imaginary: *mut f64,
        left: *mut f64,
        ld_left: *const i32,
        right: *mut f64,
        ld_right: *const i32,
        workspace: *mut f64,
        work_len: *const i32,
        info: *mut i32,
    );
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
        // Accelerate explicitly supports in-place vector exponential.
        unsafe { vvexp(output.as_mut_ptr(), output.as_ptr(), &(output.len() as i32)) };
    }
    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    for (value, distance) in output.iter_mut().zip(input) {
        *value = (scale * distance).exp();
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
    let lda = if transpose_a { m } else { k };
    #[cfg(all(target_vendor = "apple", not(feature = "portable")))]
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
    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    unsafe {
        let (row_a, col_a) = if transpose_a {
            (1, lda as isize)
        } else {
            (lda as isize, 1)
        };
        matrixmultiply::dgemm(
            m,
            k,
            n,
            1.0,
            a.as_ptr(),
            row_a,
            col_a,
            b.as_ptr(),
            n as isize,
            1,
            0.0,
            c.as_mut_ptr(),
            n as isize,
            1,
        );
    }
}

pub(crate) struct EigenWorkspace {
    matrix: Vec<f64>,
    values: Vec<f64>,
    auxiliary: Vec<f64>,
    scratch: Vec<f64>,
    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    permutation: Vec<usize>,
}

impl EigenWorkspace {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            matrix: vec![0.0; n * n],
            values: vec![0.0; n],
            auxiliary: vec![0.0; n],
            scratch: vec![
                0.0;
                if cfg!(all(target_vendor = "apple", not(feature = "portable"))) {
                    8 * n
                } else {
                    n * n
                }
            ],
            #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
            permutation: (0..n).collect(),
        }
    }

    /// Eigenvalues of diag(positive, negative, ...) * Omega.
    #[cfg(all(target_vendor = "apple", not(feature = "portable")))]
    pub(crate) fn weighted_eigenvalues(
        &mut self,
        omega: &[f64],
        width: usize,
        positive: f64,
        negative: f64,
    ) -> Result<&[f64], CoreError> {
        let n = self.values.len();
        for j in 0..n {
            for i in 0..n {
                let scale = if i == 0 { positive } else { negative };
                self.matrix[j * n + i] = scale * omega[i * width + j];
            }
        }
        let n = n as i32;
        let work_len = self.scratch.len() as i32;
        let mut dummy_left = 0.0;
        let mut dummy_right = 0.0;
        let mut status = 0;
        unsafe {
            dgeev_(
                b"N".as_ptr(),
                b"N".as_ptr(),
                &n,
                self.matrix.as_mut_ptr(),
                &n,
                self.values.as_mut_ptr(),
                self.auxiliary.as_mut_ptr(),
                &mut dummy_left,
                &1,
                &mut dummy_right,
                &1,
                self.scratch.as_mut_ptr(),
                &work_len,
                &mut status,
            );
        }
        if status != 0 {
            return Err(CoreError::internal(
                "SCPC covariance eigensolver failed to converge",
            ));
        }
        Ok(&self.values)
    }

    /// A pivoted square root turns the quadratic form into a symmetric
    /// eigenproblem on platforms without Accelerate. This includes singular PSD covariances.
    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    pub(crate) fn weighted_eigenvalues(
        &mut self,
        omega: &[f64],
        width: usize,
        positive: f64,
        negative: f64,
    ) -> Result<&[f64], CoreError> {
        let n = self.values.len();
        self.scratch.fill(0.0);
        for i in 0..n {
            self.permutation[i] = i;
            self.auxiliary[i] = omega[i * width + i];
        }
        let tolerance =
            32.0 * n as f64 * f64::EPSILON * self.auxiliary.iter().copied().fold(0.0, f64::max);
        let mut rank = 0;
        for k in 0..n {
            let p = (k..n)
                .max_by(|&a, &b| self.auxiliary[a].total_cmp(&self.auxiliary[b]))
                .expect("pivot range is nonempty");
            if self.auxiliary[p] <= 0.0 {
                if self.auxiliary[k..].iter().any(|&value| value < -tolerance) {
                    return Err(CoreError::invalid_input(
                        "SCPC covariance is not positive semidefinite",
                    ));
                }
                break;
            }
            self.permutation.swap(k, p);
            self.auxiliary.swap(k, p);
            for j in 0..k {
                self.scratch.swap(k * n + j, p * n + j);
            }
            let root = self.auxiliary[k].sqrt();
            self.scratch[k * n + k] = root;
            for i in k + 1..n {
                let product: f64 = (0..k)
                    .map(|j| self.scratch[i * n + j] * self.scratch[k * n + j])
                    .sum();
                let value =
                    (omega[self.permutation[i] * width + self.permutation[k]] - product) / root;
                self.scratch[i * n + k] = value;
                self.auxiliary[i] -= value * value;
            }
            rank += 1;
        }
        self.matrix.fill(0.0);
        for i in 0..rank {
            for j in i..rank {
                let value = (j..n)
                    .map(|k| {
                        let sign = if self.permutation[k] == 0 {
                            positive
                        } else {
                            negative
                        };
                        self.scratch[k * n + i] * sign * self.scratch[k * n + j]
                    })
                    .sum();
                self.matrix[i * n + j] = value;
                self.matrix[j * n + i] = value;
            }
        }
        symmetric_eigen(&mut self.matrix, &mut self.values, None)?;
        Ok(&self.values)
    }
}

#[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
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
