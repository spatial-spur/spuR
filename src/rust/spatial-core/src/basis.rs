//! Symmetric eigensolvers and the centered SCPC working basis.

use crate::{checked_len, linalg, CoreError};
use oxiblas_lapack::evd::SymmetricEvdDc;
use oxiblas_matrix::MatRef;
use rayon::prelude::*;

pub struct ScpcBasis {
    /// Row-major weights `n x (q + 1)`, with a normalized constant in column 0.
    pub weights: Vec<f64>,
    /// True when leading eigenvalues are numerically tied or degenerate across the cutoff.
    pub tied: bool,
}

/// The leading algebraic eigenvectors of the exponential covariance restricted
/// to the centered subspace, with the normalized constant prepended.
pub fn scpc_basis(
    distances: &[f64],
    n: usize,
    scale: f64,
    q: usize,
) -> Result<ScpcBasis, CoreError> {
    if q == 0 || q >= n || distances.len() != checked_len(&[n, n])? {
        return Err(CoreError::invalid_input(
            "SCPC basis dimensions do not match",
        ));
    }
    if !scale.is_finite() || scale < 0.0 {
        return Err(CoreError::invalid_input(
            "SCPC basis scale must be finite and nonnegative",
        ));
    }
    if distances.iter().any(|&d| !d.is_finite() || d < 0.0) {
        return Err(CoreError::invalid_input(
            "distances must be finite and nonnegative",
        ));
    }

    let dimension = n - 1;
    let covariance = contrast_covariance(distances, n, scale);
    let eigen = leading_eigenpairs(&covariance, dimension, q)?;

    let root_n = (n as f64).sqrt();
    let mut result = vec![root_n.recip(); checked_len(&[n, q + 1])?];
    for (j, vector) in eigen.vectors.chunks_exact(dimension).enumerate() {
        let sum = vector.iter().sum::<f64>();
        let correction = sum / (root_n * (root_n + 1.0));
        for (row, &value) in result.chunks_exact_mut(q + 1).zip(vector) {
            row[j + 1] = value - correction;
        }
        result[dimension * (q + 1) + j + 1] = -sum / root_n;
    }
    Ok(ScpcBasis {
        weights: result,
        tied: eigen.tied,
    })
}

pub struct SymmetricEigen {
    pub values: Vec<f64>,
    pub vectors: Vec<f64>,
    pub tied: bool,
}

pub(crate) fn symmetric_eigenpairs(
    matrix: &[f64],
    n: usize,
    q: usize,
) -> Result<SymmetricEigen, CoreError> {
    if n == 0 || q == 0 || q > n || matrix.len() != checked_len(&[n, n])? {
        return Err(CoreError::invalid_input(
            "symmetric eigenproblem dimensions do not match",
        ));
    }
    leading_eigenpairs(matrix, n, q)
}

fn leading_eigenpairs(covariance: &[f64], n: usize, q: usize) -> Result<SymmetricEigen, CoreError> {
    let steps = (2 * q + 16).max(64);
    if n >= 256 && steps < n / 2 {
        if let Some(pairs) = partial_eigenvectors(covariance, n, q, steps)? {
            return Ok(SymmetricEigen {
                values: pairs.values,
                vectors: pairs.vectors,
                tied: false,
            });
        }
    }
    let mut pairs = dense_eigenpairs(covariance, n, (q + 1).min(n))?;
    let tied = has_tied_values(&pairs.values);
    pairs.values.truncate(q);
    pairs.vectors.truncate(q * n);
    Ok(SymmetricEigen {
        values: pairs.values,
        vectors: pairs.vectors,
        tied,
    })
}

fn contrast_covariance(distances: &[f64], n: usize, scale: f64) -> Vec<f64> {
    // The first n-1 columns of the Householder reflection taking the last
    // coordinate vector to -1/sqrt(n) span the centered subspace.
    let m = n - 1;
    let root_n = (n as f64).sqrt();
    let mut last_row = vec![0.0; n];
    linalg::exp_scaled(&distances[m * n..], -scale, &mut last_row);
    let projection = |row_sum: f64, last: f64| (row_sum / root_n + last) / (root_n + 1.0);
    let mut projected = vec![0.0; n];
    let mut covariance = vec![0.0; m * m];
    let fill_row = |i, (row, result): (&mut [f64], &mut f64)| {
        let start = i * n;
        linalg::exp_scaled(&distances[start..start + m], -scale, row);
        *result = projection(row.iter().sum::<f64>() + last_row[i], last_row[i]);
    };
    let parallel = !crate::threads::is_serial() && m >= 256;
    if parallel {
        covariance
            .par_chunks_exact_mut(m)
            .zip(projected[..m].par_iter_mut())
            .enumerate()
            .for_each(|(i, row)| fill_row(i, row));
    } else {
        for (i, row) in covariance
            .chunks_exact_mut(m)
            .zip(&mut projected)
            .enumerate()
        {
            fill_row(i, row);
        }
    }
    projected[m] = projection(last_row.iter().sum(), last_row[m]);
    let common = projection(projected.iter().sum(), projected[m]);
    let center_row = |i, row: &mut [f64]| {
        for (j, value) in row.iter_mut().enumerate() {
            *value += common - (projected[i] + projected[j]);
        }
    };
    if parallel {
        covariance
            .par_chunks_exact_mut(m)
            .enumerate()
            .for_each(|(i, row)| center_row(i, row));
    } else {
        for (i, row) in covariance.chunks_exact_mut(m).enumerate() {
            center_row(i, row);
        }
    }
    covariance
}

pub(crate) struct Eigenpairs {
    pub(crate) values: Vec<f64>,
    pub(crate) vectors: Vec<f64>,
}

fn eigen_tolerance(leading: f64) -> f64 {
    2e-12 * leading.abs().max(1.0)
}

fn has_tied_values(values: &[f64]) -> bool {
    let tolerance = 8.0 * eigen_tolerance(values[0]);
    values.windows(2).any(|pair| pair[0] - pair[1] <= tolerance)
}

pub(crate) fn dense_eigenpairs(a: &[f64], n: usize, q: usize) -> Result<Eigenpairs, CoreError> {
    // For this symmetric matrix, its row-major buffer is also a valid
    // column-major view. Eigenvectors are stored by column in descending order.
    let matrix = MatRef::from_column_major(a, n, n).ok_or_else(|| {
        CoreError::invalid_input("symmetric eigenproblem dimensions do not match")
    })?;
    let eigen = SymmetricEvdDc::compute(matrix)
        .map_err(|_| CoreError::internal("SCPC basis eigensolver failed to converge"))?;
    let vectors = eigen.eigenvectors();
    Ok(Eigenpairs {
        values: eigen.eigenvalues().iter().rev().take(q).copied().collect(),
        vectors: (n - q..n)
            .rev()
            .flat_map(|j| (0..n).map(move |i| vectors[(i, j)]))
            .collect(),
    })
}

fn orthogonalize(vectors: &[f64], work: &mut [f64], coefficients: &mut [f64]) -> f64 {
    coefficients.fill(0.0);
    for _ in 0..2 {
        for (vector, coefficient) in vectors.chunks_exact(work.len()).zip(&mut *coefficients) {
            let dot = vector.iter().zip(&*work).map(|(x, y)| x * y).sum::<f64>();
            for (value, &direction) in work.iter_mut().zip(vector) {
                *value -= dot * direction;
            }
            *coefficient += dot;
        }
    }
    work.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn seed_vector(state: &mut u64, values: &mut [f64]) {
    for value in values {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *value = ((*state >> 11) as f64) * (1.0 / (1_u64 << 53) as f64) - 0.5;
    }
}

fn partial_eigenvectors(
    a: &[f64],
    n: usize,
    q: usize,
    steps: usize,
) -> Result<Option<Eigenpairs>, CoreError> {
    const STARTS: usize = 2;
    let capacity = steps + STARTS;
    let retained = q + 4;
    let mut vectors = vec![0.0; n * capacity];
    let mut projected = vec![0.0; steps * steps];
    let mut coupling = vec![0.0; STARTS * steps];
    let mut work = vec![0.0; n];
    let mut coefficients = vec![0.0; capacity];
    let mut basis = vec![0.0; n * retained];
    let mut product = vec![0.0; n * retained];
    let mut state = 0xa0761d6478bd642f_u64;
    let mut available = 0;
    let mut processed = 0;

    for _ in 0..8 {
        for j in processed..steps {
            while available < STARTS || j == available {
                seed_vector(&mut state, &mut work);
                let length = orthogonalize(&vectors[..available * n], &mut work, &mut coefficients);
                if length <= 1e-12 {
                    return Ok(None);
                }
                for (out, value) in vectors[available * n..(available + 1) * n]
                    .iter_mut()
                    .zip(&work)
                {
                    *out = value / length;
                }
                available += 1;
            }
            linalg::multiply(a, &vectors[j * n..(j + 1) * n], &mut work, n, n, 1, false);
            let length = orthogonalize(&vectors[..available * n], &mut work, &mut coefficients);
            for (i, &value) in coefficients[..available].iter().enumerate() {
                if i < steps {
                    projected[j * steps + i] = value;
                    projected[i * steps + j] = value;
                } else {
                    coupling[(i - steps) * steps + j] = value;
                }
            }
            if available < capacity && length > 1e-13 * coefficients[j].abs().max(1.0) {
                for (out, value) in vectors[available * n..(available + 1) * n]
                    .iter_mut()
                    .zip(&work)
                {
                    *out = value / length;
                }
                if available >= steps {
                    coupling[(available - steps) * steps + j] = length;
                }
                available += 1;
            }
        }

        let ritz = dense_eigenpairs(&projected, steps, retained)?;
        let row_major_coefficients: Vec<_> = (0..steps)
            .flat_map(|i| {
                ritz.vectors
                    .chunks_exact(steps)
                    .map(move |vector| vector[i])
            })
            .collect();
        linalg::multiply(
            &vectors[..n * steps],
            &row_major_coefficients,
            &mut basis,
            n,
            steps,
            retained,
            true,
        );
        let tolerance = eigen_tolerance(ritz.values[0]);
        let residuals_small = ritz.vectors.chunks_exact(steps).take(q + 1).all(|vector| {
            coupling
                .chunks_exact(steps)
                .map(|row| {
                    row.iter()
                        .zip(vector)
                        .map(|(x, y)| x * y)
                        .sum::<f64>()
                        .powi(2)
                })
                .sum::<f64>()
                .sqrt()
                <= tolerance
        });
        let transpose = |j| basis.chunks_exact(retained).map(move |row| row[j]);
        if residuals_small {
            linalg::multiply(a, &basis, &mut product, n, n, retained, false);
            let certified = (0..=q).all(|j| {
                basis
                    .chunks_exact(retained)
                    .zip(product.chunks_exact(retained))
                    .map(|(x, ax)| (ax[j] - ritz.values[j] * x[j]).powi(2))
                    .sum::<f64>()
                    .sqrt()
                    <= tolerance
            });
            return Ok(
                (certified && !has_tied_values(&ritz.values[..=q])).then(|| Eigenpairs {
                    values: ritz.values[..q].to_vec(),
                    vectors: (0..q).flat_map(transpose).collect(),
                }),
            );
        }
        for (out, value) in vectors[..n * retained]
            .iter_mut()
            .zip((0..retained).flat_map(transpose))
        {
            *out = value;
        }
        vectors.copy_within(steps * n..available * n, retained * n);
        available = retained + (available - steps);
        processed = retained;
        projected.fill(0.0);
        for (i, &value) in ritz.values.iter().enumerate() {
            projected[i * steps + i] = value;
        }
        coupling.fill(0.0);
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partial_solver_checks_interior_and_cutoff_ties() {
        // n and q enter the partial-solver branch. Well-separated low
        // eigenvalues make its q+1 residual certification deterministic.
        let n = 260;
        let q = 2;
        let steps = (2 * q + 16).max(64);
        let mut matrix = vec![0.0; n * n];
        for i in 0..n {
            matrix[i * n + i] = 1.0 / (i as f64 + 1.0);
        }
        for (i, &v) in [10.0, 9.0, 8.0].iter().enumerate() {
            matrix[i * n + i] = v;
        }
        assert!(partial_eigenvectors(&matrix, n, q, steps)
            .unwrap()
            .is_some());
        assert!(!symmetric_eigenpairs(&matrix, n, q).unwrap().tied);

        // Ties within the retained prefix and across its cutoff both force
        // dense fallback, which reports ambiguity to the host.
        for leading in [[10.0, 10.0, 8.0], [10.0, 9.0, 9.0]] {
            for (i, &v) in leading.iter().enumerate() {
                matrix[i * n + i] = v;
            }
            assert!(partial_eigenvectors(&matrix, n, q, steps)
                .unwrap()
                .is_none());
            assert!(symmetric_eigenpairs(&matrix, n, q).unwrap().tied);
        }
    }

    #[test]
    fn dense_full_eigenpairs_reconstructs_and_is_orthonormal() {
        let n = 24usize;
        let mut a = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..=i {
                let d = ((i as f64) - (j as f64)).abs() / (n as f64);
                let v = (-2.0 * d).exp();
                a[i * n + j] = v;
                a[j * n + i] = v;
            }
        }
        let pairs = dense_eigenpairs(&a, n, n).unwrap();
        assert_eq!(pairs.values.len(), n);
        assert_eq!(pairs.vectors.len(), n * n);

        // Descending eigenvalues
        for w in pairs.values.windows(2) {
            assert!(w[0] >= w[1] - 1e-14);
        }

        // Orthonormality V^T V = I_n
        for i in 0..n {
            for j in 0..n {
                let dot: f64 = (0..n)
                    .map(|k| pairs.vectors[i * n + k] * pairs.vectors[j * n + k])
                    .sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dot - expected).abs() < 1e-12, "V'V[{i},{j}]={dot}");
            }
        }

        // Reconstruction A = V diag(lambda) V^T
        for i in 0..n {
            for j in 0..n {
                let rec: f64 = (0..n)
                    .map(|k| pairs.values[k] * pairs.vectors[k * n + i] * pairs.vectors[k * n + j])
                    .sum();
                assert!(
                    (rec - a[i * n + j]).abs() < 1e-12,
                    "recon[{i},{j}]={rec} vs {}",
                    a[i * n + j]
                );
            }
        }
    }
}
