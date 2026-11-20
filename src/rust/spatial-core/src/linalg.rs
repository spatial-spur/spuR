use crate::CoreError;
use oxiblas_blas::level3::{gemm::gemm_transposed, Trans};
use oxiblas_matrix::{MatMut, MatRef};

pub(crate) fn exp_scaled(input: &[f64], scale: f64, output: &mut [f64]) {
    debug_assert_eq!(input.len(), output.len());
    if scale == 0.0 {
        output.fill(1.0);
        return;
    }
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
    // A row-major buffer is a column-major view of its transpose. Compute
    // Cᵀ = Bᵀ Aᵀ directly into the output buffer without copying either input.
    let b_t = MatRef::from_column_major(b, n, k).expect("right matrix fits its buffer");
    let c_t = MatMut::from_column_major(c, n, m).expect("output matrix fits its buffer");
    let (rows, cols, transpose) = if transpose_a {
        (m, k, Trans::Trans)
    } else {
        (k, m, Trans::NoTrans)
    };
    let a_view = MatRef::from_column_major(a, rows, cols).expect("left matrix fits its buffer");
    gemm_transposed(Trans::NoTrans, transpose, 1.0, b_t, a_view, 0.0, c_t);
}

pub(crate) struct EigenWorkspace {
    matrix: Vec<f64>,
    values: Vec<f64>,
    auxiliary: Vec<f64>,
    scratch: Vec<f64>,
    permutation: Vec<usize>,
}

impl EigenWorkspace {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            matrix: vec![0.0; n * n],
            values: vec![0.0; n],
            auxiliary: vec![0.0; n],
            scratch: vec![0.0; n * n],
            permutation: (0..n).collect(),
        }
    }

    /// A pivoted square root turns the quadratic form into a symmetric
    /// eigenproblem, including when the covariance is singular and PSD.
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
            if self.auxiliary[p] <= tolerance {
                // A roundoff-sized positive pivot can amplify noise when used
                // as a divisor. Stop only if the whole remaining Schur block
                // is within tolerance, including its off-diagonal entries.
                let negligible = self.auxiliary[k..]
                    .iter()
                    .all(|value| value.abs() <= tolerance)
                    && (k..n).all(|i| {
                        (k..i).all(|j| {
                            let product: f64 = (0..k)
                                .map(|l| self.scratch[i * n + l] * self.scratch[j * n + l])
                                .sum();
                            let residual =
                                omega[self.permutation[i] * width + self.permutation[j]] - product;
                            residual.abs() <= tolerance
                        })
                    });
                if !negligible {
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
        jacobi_eigenvalues(&mut self.matrix, &mut self.values)?;
        Ok(&self.values)
    }
}

// The SCPC probability evaluator solves many small eigenvalues-only problems
// and reuses its workspace across critical-value trials.
fn jacobi_eigenvalues(a: &mut [f64], values: &mut [f64]) -> Result<(), CoreError> {
    let n = values.len();
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
    use super::{multiply, EigenWorkspace};

    #[test]
    fn multiply_respects_row_major_layout_and_transpose() {
        let right = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut product = [0.0; 6];
        let expected = [9.0, 12.0, 15.0, 19.0, 26.0, 33.0];

        multiply(&[1.0, 2.0, 3.0, 4.0], &right, &mut product, 2, 2, 3, false);
        assert_eq!(product, expected);
        multiply(&[1.0, 3.0, 2.0, 4.0], &right, &mut product, 2, 2, 3, true);
        assert_eq!(product, expected);
        multiply(&[], &[], &mut product, 2, 0, 3, false);
        assert_eq!(product, [0.0; 6]);
    }

    #[test]
    fn weighted_eigenvalues_stop_at_roundoff_pivots() {
        let eps = f64::EPSILON;
        let mut covariance = [1.0, 0.0, 0.0, 0.0, eps, 16.0 * eps, 0.0, 16.0 * eps, 0.0];
        let mut workspace = EigenWorkspace::new(3);
        let values = workspace
            .weighted_eigenvalues(&covariance, 3, 1.0, -1.0)
            .unwrap();
        assert_eq!(values, [1.0, 0.0, 0.0]);

        // Roundoff-sized diagonal entries can conceal an indefinite block.
        covariance[5] = 1e-3;
        covariance[7] = 1e-3;
        assert!(workspace
            .weighted_eigenvalues(&covariance, 3, 1.0, -1.0)
            .is_err());
    }
}
