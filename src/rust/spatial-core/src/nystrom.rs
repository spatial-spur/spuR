//! Landmark Nyström basis construction and sampled covariance grids for large-n SCPC.

use crate::distance::{Distances, Metric};
use crate::{basis, calibration, checked_len, linalg, CoreError};
use rayon::prelude::*;

const PAIR_TILE: usize = 4096;
const ROW_TILE: usize = 128;

#[derive(Clone, Copy)]
struct Lcg(u32);

impl Lcg {
    #[inline]
    fn next_u(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(64389).wrapping_add(1);
        f64::from(self.0) / 4_294_967_296.0
    }

    #[inline]
    fn index(&mut self, count: usize) -> usize {
        if count <= 1 {
            self.0 = self.0.wrapping_mul(64389).wrapping_add(1);
            return 0;
        }
        let raw = (self.next_u() * count as f64) as usize;
        raw.min(count - 1)
    }
}

/// Normalize coordinates in the exact host convention so randomized subsampling is invariant.
///
/// For `Metric::Euclidean`, centers columns, rotates by PCA eigenvectors (descending eigenvalues)
/// with deterministic sign on column 0, sorts rows lexicographically, shifts by column minimum,
/// and scales by overall maximum.
/// For geodesic metrics, centers and wraps longitude to `[-180, 180)` and sorts rows by `(lat, lon)`.
pub fn normalize_coords(
    coords: &[f64],
    n: usize,
    d: usize,
    metric: Metric,
) -> Result<(Vec<f64>, Vec<usize>), CoreError> {
    if n == 0 || d == 0 || coords.len() != checked_len(&[n, d])? {
        return Err(CoreError::invalid_input(
            "invalid coordinate dimensions for normalization",
        ));
    }
    if coords.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input("coordinates must be finite"));
    }

    match metric {
        Metric::Euclidean => {
            let mut bounds = vec![(f64::INFINITY, f64::NEG_INFINITY); d];
            for row in coords.chunks_exact(d) {
                for (&value, (low, high)) in row.iter().zip(&mut bounds) {
                    *low = low.min(value);
                    *high = high.max(value);
                }
            }
            let anchors: Vec<_> = bounds.iter().map(|(lo, hi)| lo.midpoint(*hi)).collect();
            let mut centered = vec![0.0_f64; n * d];
            let mut max_abs = 0.0_f64;
            for (dst, src) in centered.chunks_exact_mut(d).zip(coords.chunks_exact(d)) {
                for c in 0..d {
                    let v = src[c] - anchors[c];
                    dst[c] = v;
                    max_abs = max_abs.max(v.abs());
                }
            }
            let mut rotated = vec![0.0_f64; n * d];
            if max_abs > 0.0 && max_abs.is_finite() {
                let inv_scale = max_abs.recip();
                if inv_scale.is_finite() {
                    for v in &mut centered {
                        *v *= inv_scale;
                    }
                } else {
                    for v in &mut centered {
                        *v /= max_abs;
                    }
                }
                let mut means = vec![0.0; d];
                for row in centered.chunks_exact(d) {
                    for (mean, value) in means.iter_mut().zip(row) {
                        *mean += value;
                    }
                }
                for mean in &mut means {
                    *mean /= n as f64;
                }
                for row in centered.chunks_exact_mut(d) {
                    for (value, mean) in row.iter_mut().zip(&means) {
                        *value -= mean;
                    }
                }
                let mut cov = vec![0.0_f64; d * d];
                linalg::multiply(&centered, &centered, &mut cov, d, n, d, true);
                let eigen = basis::symmetric_eigenpairs(&cov, d, d)?;

                // eigen.vectors stores each eigenvector r in [r * d .. (r + 1) * d].
                // Build rot_t (d x d) where column r is eigenvector r: rot[c * d + r] = vectors[r * d + c].
                let mut rot = vec![0.0_f64; d * d];
                for r in 0..d {
                    for c in 0..d {
                        rot[c * d + r] = eigen.vectors[r * d + c];
                    }
                }
                linalg::multiply(&centered, &rot, &mut rotated, n, d, d, false);
            }

            let mut max_col0 = f64::NEG_INFINITY;
            let mut max_abs_col0 = 0.0_f64;
            for row in rotated.chunks_exact(d) {
                max_col0 = max_col0.max(row[0]);
                max_abs_col0 = max_abs_col0.max(row[0].abs());
            }
            if max_col0 != max_abs_col0 {
                for v in &mut rotated {
                    *v = -*v;
                }
            }

            let mut perm: Vec<usize> = (0..n).collect();
            perm.sort_by(|&a, &b| {
                let ra = &rotated[a * d..(a + 1) * d];
                let rb = &rotated[b * d..(b + 1) * d];
                for c in 0..d {
                    let ord = ra[c].total_cmp(&rb[c]);
                    if ord != std::cmp::Ordering::Equal {
                        return ord;
                    }
                }
                a.cmp(&b)
            });

            let mut out = vec![0.0_f64; n * d];
            let mut col_min = vec![f64::INFINITY; d];
            for (dst_idx, &src_idx) in perm.iter().enumerate() {
                let src_row = &rotated[src_idx * d..(src_idx + 1) * d];
                let dst_row = &mut out[dst_idx * d..(dst_idx + 1) * d];
                for c in 0..d {
                    dst_row[c] = src_row[c];
                    col_min[c] = col_min[c].min(src_row[c]);
                }
            }
            let mut global_max = 0.0_f64;
            for row in out.chunks_exact_mut(d) {
                for c in 0..d {
                    row[c] -= col_min[c];
                    global_max = global_max.max(row[c]);
                }
            }
            if global_max.is_finite() && global_max > 0.0 {
                let inv = global_max.recip();
                for v in &mut out {
                    *v *= inv;
                }
            }
            Ok((out, perm))
        }
        Metric::GreatCircleFractionLatLon
        | Metric::GreatCircle { .. }
        | Metric::GreatCircleFractionLonLat => {
            if d != 2 {
                return Err(CoreError::invalid_input(
                    "geographic normalization requires 2 coordinate columns",
                ));
            }
            let (lat_col, lon_col) = match metric {
                Metric::GreatCircleFractionLonLat => (1usize, 0usize),
                _ => (0usize, 1usize),
            };
            let mean_lon = coords.chunks_exact(2).map(|r| r[lon_col]).sum::<f64>() / n as f64;
            let mut wrapped = coords.to_vec();
            for row in wrapped.chunks_exact_mut(2) {
                let shifted = row[lon_col] - mean_lon;
                row[lon_col] = (shifted + 180.0).rem_euclid(360.0) - 180.0;
            }
            let mut perm: Vec<usize> = (0..n).collect();
            perm.sort_by(|&a, &b| {
                let ra = &wrapped[a * 2..a * 2 + 2];
                let rb = &wrapped[b * 2..b * 2 + 2];
                ra[lat_col]
                    .total_cmp(&rb[lat_col])
                    .then_with(|| ra[lon_col].total_cmp(&rb[lon_col]))
                    .then_with(|| a.cmp(&b))
            });
            let mut out = vec![0.0_f64; n * 2];
            for (dst_idx, &src_idx) in perm.iter().enumerate() {
                out[dst_idx * 2] = wrapped[src_idx * 2];
                out[dst_idx * 2 + 1] = wrapped[src_idx * 2 + 1];
            }
            Ok((out, perm))
        }
    }
}

// Consecutive links are sampled with replacement, excluding self-pairs.
struct PairChain {
    n: usize,
    current: usize,
    rng: Lcg,
}

impl PairChain {
    fn new(n: usize, count: usize, seed: u32) -> Result<Self, CoreError> {
        if n < 2 || n > i32::MAX as usize || count == 0 || count >= isize::MAX as usize / 8 {
            return Err(CoreError::invalid_input(
                "invalid SCPC pair-chain dimensions",
            ));
        }
        let mut rng = Lcg(seed);
        Ok(Self {
            n,
            current: rng.index(n),
            rng,
        })
    }

    fn next(&mut self) -> usize {
        self.current = (self.current + 1 + self.rng.index(self.n - 1)) % self.n;
        self.current
    }

    fn final_state(mut self) -> u32 {
        // Released raninds() takes one extra step after writing its last index.
        self.rng.index(self.n - 1);
        self.rng.0
    }
}

/// Match the released SCPC pair chain, returning indices and final LCG state.
pub fn pair_chain(n: usize, count: usize, seed: u32) -> Result<(Vec<usize>, u32), CoreError> {
    let mut chain = PairChain::new(n, count, seed)?;
    let mut indices = Vec::with_capacity(count + 1);
    indices.push(chain.current);
    indices.extend((0..count).map(|_| chain.next()));
    Ok((indices, chain.final_state()))
}

fn landmarks(n: usize, size: usize, samples: usize, seed: u32) -> (Vec<Vec<usize>>, u32) {
    let mut rng = Lcg(seed);
    let mut order: Vec<usize> = (0..n).collect();
    let subsets = (0..samples)
        .map(|_| {
            for i in 0..size {
                let j = i + rng.index(n - i);
                order.swap(i, j);
            }
            order[..size].to_vec()
        })
        .collect();
    (subsets, rng.0)
}

/// Lower and upper SCPC decay scales with co-location floor adjustment.
pub fn decay_limits(distances: &[f64], avc: f64, min_avc: f64) -> Result<[f64; 2], CoreError> {
    let floor = calibration::correlation_floor(distances)?;
    if avc <= floor {
        return Err(CoreError::calibration("unattainable sampled SCPC decay").with_message(format!(
            "Target correlation avc={avc} must exceed the pooled landmark co-located-pair fraction {floor:.6}."
        )));
    }
    let c0 = calibration::solve_c0_from_avc(distances, avc)?;
    let cmax = calibration::solve_c0_from_avc(distances, floor + (1.0 - floor) * min_avc)?;
    Ok([c0, cmax])
}

/// Controls the landmark Nyström approximation.
pub struct LandmarkRequest {
    pub avc: f64,
    pub min_avc: f64,
    pub components: usize,
    pub landmarks: usize,
    pub samples: usize,
    pub seed: u32,
}

/// Row-major Nyström basis (`n x (q + 1)`), decay limits `[c0, cmax]`, and terminal LCG state.
pub struct LandmarkBasis {
    pub weights: Vec<f64>,
    pub decay: [f64; 2],
    pub state: u32,
}

/// Build the ensemble Nyström basis. Returns `Err(TiedEigenspace)` when a subset or combined
/// eigenspace has a split tie so `auto` can preserve host orientation via fallback.
pub fn landmark_basis(
    coords: &[f64],
    n: usize,
    width: usize,
    metric: Metric,
    request: &LandmarkRequest,
) -> Result<LandmarkBasis, CoreError> {
    let distance = Distances::new(coords, n, width, metric)?;
    let q = request.components;
    let m = request.landmarks.min(n);
    if q == 0 || q >= m || request.samples == 0 || n > i32::MAX as usize {
        return Err(CoreError::invalid_input("invalid SCPC landmark dimensions"));
    }
    let total_columns = checked_len(&[q, request.samples])?;
    let wall_len = checked_len(&[n, total_columns])?;
    let pairs = checked_len(&[m, m - 1])? / 2;
    let (subsets, state) = landmarks(n, m, request.samples, request.seed);

    let mut distances = vec![0.0; checked_len(&[request.samples, pairs])?];
    if crate::threads::is_serial() || request.samples == 1 {
        for (block, subset) in distances.chunks_exact_mut(pairs).zip(&subsets) {
            let mut k = 0;
            for j in 0..m - 1 {
                for i in j + 1..m {
                    block[k] = distance.pair(subset[i], subset[j]);
                    k += 1;
                }
            }
        }
    } else {
        distances
            .par_chunks_exact_mut(pairs)
            .zip(&subsets)
            .for_each(|(block, subset)| {
                let mut k = 0;
                for j in 0..m - 1 {
                    for i in j + 1..m {
                        block[k] = distance.pair(subset[i], subset[j]);
                        k += 1;
                    }
                }
            });
    }

    let decay = decay_limits(&distances, request.avc, request.min_avc)?;

    let extend_one = |(triangle, subset): (&[f64], &Vec<usize>)| -> Result<Vec<f64>, CoreError> {
        let mut dense = vec![0.0; m * m];
        let mut k = 0;
        for j in 0..m - 1 {
            for i in j + 1..m {
                dense[i * m + j] = triangle[k];
                dense[j * m + i] = triangle[k];
                k += 1;
            }
        }
        let local = basis::scpc_basis(&dense, m, decay[0], q)?;
        if local.tied {
            return Err(CoreError::tied(
                "tied eigenspace in Nystrom landmark subset",
            ));
        }
        let weights: Vec<f64> = local
            .weights
            .chunks_exact(q + 1)
            .flat_map(|row| row[1..].iter().copied())
            .collect();
        Ok(project_landmarks(
            &distance, n, subset, &weights, q, decay[0],
        )?)
    };

    let projected: Vec<Vec<f64>> = if crate::threads::is_serial() || request.samples == 1 {
        distances
            .chunks_exact(pairs)
            .zip(&subsets)
            .map(extend_one)
            .collect::<Result<Vec<_>, _>>()?
    } else {
        distances
            .par_chunks_exact(pairs)
            .zip(&subsets)
            .map(extend_one)
            .collect::<Result<Vec<_>, _>>()?
    };
    drop(distances);

    let mut wall = vec![0.0; wall_len];
    for (sample, projection) in projected.into_iter().enumerate() {
        for (destination, source) in wall
            .chunks_exact_mut(total_columns)
            .zip(projection.chunks_exact(q))
        {
            for (j, &value) in source.iter().enumerate() {
                destination[j * request.samples + sample] = value;
            }
        }
    }
    let weights = combine_landmarks(wall, n, q, request.samples)?;
    Ok(LandmarkBasis {
        weights,
        decay,
        state,
    })
}

fn project_landmarks(
    distance: &Distances<'_>,
    n: usize,
    subset: &[usize],
    weights: &[f64],
    q: usize,
    decay: f64,
) -> Result<Vec<f64>, CoreError> {
    let m = subset.len();
    let mut projected = vec![0.0; n * q];
    let mut d = vec![0.0; ROW_TILE.min(n) * m];
    let mut kernel = d.clone();
    for (block_index, output) in projected.chunks_mut(ROW_TILE * q).enumerate() {
        let rows = output.len() / q;
        for (row, distances) in d[..rows * m].chunks_exact_mut(m).enumerate() {
            let i = block_index * ROW_TILE + row;
            for (value, &j) in distances.iter_mut().zip(subset) {
                let dij = distance.pair(i, j);
                if !dij.is_finite() {
                    return Err(CoreError::invalid_input(
                        "non-finite pairwise distance in Nystrom landmark projection",
                    ));
                }
                *value = dij;
            }
        }
        linalg::exp_scaled(&d[..rows * m], -decay, &mut kernel[..rows * m]);
        linalg::multiply(&kernel, weights, output, rows, m, q, false);
    }
    for j in 0..q {
        let mean = projected.chunks_exact(q).map(|row| row[j]).sum::<f64>() / n as f64;
        let mut norm = 0.0;
        for row in projected.chunks_exact_mut(q) {
            row[j] -= mean;
            norm += row[j] * row[j];
        }
        if norm.is_finite() && norm > 0.0 {
            let inverse_norm = norm.sqrt().recip();
            for row in projected.chunks_exact_mut(q) {
                row[j] *= inverse_norm;
            }
        }
    }
    Ok(projected)
}

fn combine_landmarks(
    mut wall: Vec<f64>,
    n: usize,
    q: usize,
    samples: usize,
) -> Result<Vec<f64>, CoreError> {
    let total_columns = q * samples;
    let mut weights = vec![(n as f64).sqrt().recip(); checked_len(&[n, q + 1])?];
    let mut projection = vec![0.0; total_columns];
    let mut component = vec![0.0; n];

    for j in 1..=q {
        let columns = j * samples;
        let mut prefix = Vec::with_capacity(n * columns);
        for row in wall.chunks_exact(total_columns) {
            prefix.extend_from_slice(&row[..columns]);
        }
        let mut gram = vec![0.0; checked_len(&[columns, columns])?];
        linalg::multiply(&prefix, &prefix, &mut gram, columns, n, columns, true);
        let eigen = basis::symmetric_eigenpairs(&gram, columns, 1)?;
        if eigen.tied {
            return Err(CoreError::tied(
                "tied eigenspace while consolidating Nystrom landmark basis",
            ));
        }
        linalg::multiply(
            &prefix,
            &eigen.vectors,
            &mut component,
            n,
            columns,
            1,
            false,
        );
        let norm = component.iter().map(|v| v * v).sum::<f64>().sqrt();
        if !norm.is_finite() || norm <= f64::EPSILON.sqrt() {
            return Err(CoreError::degenerate(
                "degenerate principal direction in Nystrom landmark consolidation",
            ));
        }
        let inv = norm.recip();
        for (row, comp) in weights.chunks_exact_mut(q + 1).zip(&mut component) {
            *comp *= inv;
            row[j] = *comp;
        }
        if j < q {
            linalg::multiply(
                &component,
                &wall,
                &mut projection,
                1,
                n,
                total_columns,
                true,
            );
            for (row, &weight) in wall.chunks_exact_mut(total_columns).zip(&component) {
                for (value, &proj) in row.iter_mut().zip(&projection) {
                    *value -= weight * proj;
                }
            }
        }
    }
    Ok(weights)
}

/// Check whether every scale matrix's leading `(q + 1) x (q + 1)` submatrix is positive semidefinite
/// via pivoted Cholesky factorization.
pub(crate) fn is_psd_prefix_grid(omegas: &[f64], width: usize, q: usize) -> bool {
    let Some(k) = q.checked_add(1) else {
        return false;
    };
    let Ok(matrix_len) = checked_len(&[width, width]) else {
        return false;
    };
    if k > width || omegas.is_empty() || omegas.len() % matrix_len != 0 {
        return false;
    }
    let mut diag = vec![0.0_f64; k];
    let mut perm = vec![0usize; k];
    let mut factors = vec![0.0_f64; k * k];

    for omega in omegas.chunks_exact(matrix_len) {
        let mut scale = 0.0_f64;
        for i in 0..k {
            for &value in &omega[i * width..i * width + k] {
                if !value.is_finite() {
                    return false;
                }
                scale = scale.max(value.abs());
            }
        }
        if scale == 0.0 {
            continue;
        }
        let tol = 32.0 * (k as f64) * f64::EPSILON;
        for i in 0..k {
            perm[i] = i;
            diag[i] = omega[i * width + i] / scale;
            for j in 0..i {
                if (omega[i * width + j] / scale - omega[j * width + i] / scale).abs() > tol {
                    return false;
                }
            }
        }
        if diag.iter().any(|&d| d < -tol) {
            return false;
        }
        factors.fill(0.0);
        for step in 0..k {
            let mut p = step;
            let mut best = diag[step];
            for cand in step + 1..k {
                if diag[cand] > best {
                    best = diag[cand];
                    p = cand;
                }
            }
            if best <= tol {
                if diag[step..].iter().any(|&v| v < -tol) {
                    return false;
                }
                // A zero Schur diagonal can only have zero covariances in its
                // remaining row/column. Pivot diagonals alone miss e.g. [0, 1; 1, 0].
                for i in step..k {
                    for j in step..i {
                        let dot: f64 = (0..step)
                            .map(|l| factors[i * k + l] * factors[j * k + l])
                            .sum();
                        let residual = omega[perm[i] * width + perm[j]] / scale - dot;
                        if !residual.is_finite() || residual.abs() > tol {
                            return false;
                        }
                    }
                }
                break;
            }
            perm.swap(step, p);
            diag.swap(step, p);
            for j in 0..step {
                factors.swap(step * k + j, p * k + j);
            }
            let root = best.sqrt();
            factors[step * k + step] = root;
            for i in step + 1..k {
                let dot: f64 = (0..step)
                    .map(|j| factors[i * k + j] * factors[step * k + j])
                    .sum();
                let val = (omega[perm[i] * width + perm[step]] / scale - dot) / root;
                factors[i * k + step] = val;
                diag[i] -= val * val;
                if diag[i] < -tol {
                    return false;
                }
            }
        }
    }
    true
}

// The dense representation coalesces repeated pairs before evaluating any
// coefficient. The tiled path retains one shared chain and distance per link.
enum PairGeometry {
    Dense {
        frequencies: Vec<f64>,
        distances: Vec<f64>,
    },
    Tiled {
        chain: Vec<usize>,
        distances: Vec<f64>,
    },
}

pub(crate) struct SampledPairs {
    n: usize,
    count: usize,
    state: u32,
    geometry: PairGeometry,
}

impl SampledPairs {
    pub(crate) fn new(
        coordinates: &[f64],
        n: usize,
        width: usize,
        metric: Metric,
        count: usize,
        seed: u32,
    ) -> Result<Self, CoreError> {
        let distance = Distances::new(coordinates, n, width, metric)?;
        if n <= 512 && count >= n * n.saturating_sub(1) {
            let mut chain = PairChain::new(n, count, seed)?;
            let mut frequencies = vec![0.0; n * n];
            for _ in 0..count {
                let i = chain.current;
                let j = chain.next();
                frequencies[i * n + j] += 1.0;
                frequencies[j * n + i] += 1.0;
            }
            let mut distances = vec![0.0_f64; n * n];
            for i in 0..n {
                for j in 0..i {
                    let d = distance.pair(i, j);
                    if !d.is_finite() {
                        return Err(CoreError::invalid_input(
                            "non-finite pairwise distance in sampled SCPC pairs",
                        ));
                    }
                    distances[i * n + j] = d;
                    distances[j * n + i] = d;
                }
            }

            return Ok(Self {
                n,
                count,
                state: chain.final_state(),
                geometry: PairGeometry::Dense {
                    frequencies,
                    distances,
                },
            });
        }
        let (chain, state) = pair_chain(n, count, seed)?;
        let distances: Vec<f64> = if crate::threads::is_serial() || count < 8192 {
            chain
                .windows(2)
                .map(|pair| {
                    let d = distance.pair(pair[0], pair[1]);
                    if d.is_finite() {
                        Ok(d)
                    } else {
                        Err(CoreError::invalid_input(
                            "non-finite pairwise distance in sampled SCPC pairs",
                        ))
                    }
                })
                .collect::<Result<Vec<_>, _>>()?
        } else {
            chain
                .par_windows(2)
                .map(|pair| {
                    let d = distance.pair(pair[0], pair[1]);
                    if d.is_finite() {
                        Ok(d)
                    } else {
                        Err(CoreError::invalid_input(
                            "non-finite pairwise distance in sampled SCPC pairs",
                        ))
                    }
                })
                .collect::<Result<Vec<_>, _>>()?
        };

        Ok(Self {
            n,
            count,
            state,
            geometry: PairGeometry::Tiled { chain, distances },
        })
    }

    pub(crate) fn covariance_grid(
        &self,
        weights: &[f64],
        q: usize,
        scales: &[f64],
    ) -> Result<Vec<f64>, CoreError> {
        let (n, count) = (self.n, self.count);
        if q == 0
            || weights.len() != checked_len(&[n, q])?
            || weights.iter().any(|v| !v.is_finite())
            || scales.is_empty()
            || scales.iter().any(|s| !s.is_finite() || *s < 0.0)
        {
            return Err(CoreError::invalid_input(
                "invalid sampled SCPC covariance inputs",
            ));
        }

        let matrix_size = checked_len(&[q, q])?;
        let mut gram = vec![0.0; matrix_size];
        linalg::multiply(weights, weights, &mut gram, q, n, q, true);
        let mut grid = vec![0.0; checked_len(&[scales.len(), matrix_size])?];
        let half_multiplier = 0.5 * (n as f64) * ((n - 1) as f64) / (count as f64);

        let (chain, distances) = match &self.geometry {
            PairGeometry::Dense {
                frequencies,
                distances,
            } => {
                let eval_scale = |output: &mut [f64], scale: f64| {
                    let mut kernel = vec![0.0_f64; n * n];
                    let mut work = vec![0.0_f64; n * q];
                    linalg::exp_scaled(distances, -scale, &mut kernel);
                    for (val, &freq) in kernel.iter_mut().zip(frequencies) {
                        *val *= freq;
                    }
                    linalg::multiply(&kernel, weights, &mut work, n, n, q, false);
                    linalg::multiply(weights, &work, output, q, n, q, true);
                    for a in 0..q {
                        for b in 0..=a {
                            let ab = a * q + b;
                            let ba = b * q + a;
                            let value =
                                gram[ab] + 0.5 * half_multiplier * (output[ab] + output[ba]);
                            output[ab] = value;
                            output[ba] = value;
                        }
                    }
                };
                if crate::threads::is_serial() || scales.len() <= 1 {
                    for (output, &scale) in grid.chunks_exact_mut(matrix_size).zip(scales) {
                        eval_scale(output, scale);
                    }
                } else {
                    grid.par_chunks_exact_mut(matrix_size)
                        .zip(scales.par_iter())
                        .for_each(|(output, &scale)| eval_scale(output, scale));
                }
                return Ok(grid);
            }
            PairGeometry::Tiled { chain, distances } => (chain, distances),
        };
        let tile_len = PAIR_TILE.min(count);
        let eval_scale_tiled = |output: &mut [f64], scale: f64| {
            let mut w1_buf = vec![0.0_f64; tile_len * q];
            let mut w2_buf = vec![0.0_f64; tile_len * q];
            let mut cd_buf = vec![0.0_f64; tile_len];
            let mut temp = vec![0.0_f64; matrix_size];
            output.copy_from_slice(&gram);

            for start in (0..count).step_by(tile_len) {
                let len = (count - start).min(tile_len);
                linalg::exp_scaled(&distances[start..start + len], -scale, &mut cd_buf[..len]);
                for idx in 0..len {
                    let i1 = chain[start + idx];
                    let i2 = chain[start + idx + 1];
                    let cd = cd_buf[idx];
                    w1_buf[idx * q..(idx + 1) * q].copy_from_slice(&weights[i1 * q..(i1 + 1) * q]);
                    let src2 = &weights[i2 * q..(i2 + 1) * q];
                    let dst2 = &mut w2_buf[idx * q..(idx + 1) * q];
                    for c in 0..q {
                        dst2[c] = src2[c] * cd;
                    }
                }
                linalg::multiply(
                    &w1_buf[..len * q],
                    &w2_buf[..len * q],
                    &mut temp,
                    q,
                    len,
                    q,
                    true,
                );
                for a in 0..q {
                    for b in 0..=a {
                        let ab = a * q + b;
                        let ba = b * q + a;
                        let val = half_multiplier * (temp[ab] + temp[ba]);
                        output[ab] += val;
                        if ab != ba {
                            output[ba] += val;
                        }
                    }
                }
            }
        };

        if crate::threads::is_serial() || scales.len() <= 1 {
            for (output, &scale) in grid.chunks_exact_mut(matrix_size).zip(scales) {
                eval_scale_tiled(output, scale);
            }
        } else {
            grid.par_chunks_exact_mut(matrix_size)
                .zip(scales.par_iter())
                .for_each(|(output, &scale)| eval_scale_tiled(output, scale));
        }

        Ok(grid)
    }
}

/// Sampled covariance with exact W'W diagonal and linked off-diagonal pairs.
#[allow(clippy::too_many_arguments)]
pub fn sampled_covariance_grid(
    coordinates: &[f64],
    n: usize,
    width: usize,
    metric: Metric,
    weights: &[f64],
    q: usize,
    scales: &[f64],
    count: usize,
    seed: u32,
) -> Result<(Vec<f64>, u32), CoreError> {
    let pairs = SampledPairs::new(coordinates, n, width, metric, count, seed)?;
    Ok((pairs.covariance_grid(weights, q, scales)?, pairs.state))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::distance_matrix;

    #[test]
    fn pair_chain_and_landmarks_match_released_lcg() {
        let (indices, state) = pair_chain(5, 6, 1).unwrap();
        assert_eq!(indices, vec![0, 4, 3, 1, 2, 4, 0]);
        assert_eq!(state, 3_259_150_681);

        let (subsets, lm_state) = landmarks(4, 3, 1, 1);
        assert_eq!(subsets[0], vec![0, 3, 1]);
        assert_eq!(lm_state, 3_598_220_700);
    }

    #[test]
    fn sampled_covariance_preserves_nonorthonormal_diagonal() {
        let n = 40;
        let coords: Vec<f64> = (0..n).flat_map(|i| [i as f64, 2.0 * i as f64]).collect();
        let w: Vec<f64> = (0..n)
            .flat_map(|i| [2.0, (i + 1) as f64 / n as f64, ((i + 1) as f64).sin()])
            .collect();
        let mut gram = vec![0.0; 9];
        linalg::multiply(&w, &w, &mut gram, 3, n, 3, true);
        let (grid, _) = sampled_covariance_grid(
            &coords,
            n,
            2,
            Metric::Euclidean,
            &w,
            3,
            &[1e6, 2e6],
            2000,
            7,
        )
        .unwrap();
        for mat in grid.chunks_exact(9) {
            for (a, b) in mat.iter().zip(&gram) {
                assert!((a - b).abs() < 1e-11, "a={a}, b={b}");
            }
        }
        assert!(is_psd_prefix_grid(&grid, 3, 2));
    }

    #[test]
    fn sampled_psd_check_is_scale_relative_and_checks_the_full_schur_remainder() {
        let singular = [1.0, 1.0, 1.0, 1.0];
        let indefinite = [0.0, 1.0, 1.0, 0.0];
        // Eliminating one pivot leaves [0, 1; 1, 0] rather than a zero block.
        let indefinite_schur = [1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 2.0, 1.0];
        for scale in [1e-200, 1.0, 1e200] {
            assert!(is_psd_prefix_grid(&singular.map(|v| v * scale), 2, 1));
            assert!(!is_psd_prefix_grid(&indefinite.map(|v| v * scale), 2, 1));
            assert!(!is_psd_prefix_grid(
                &indefinite_schur.map(|v| v * scale),
                3,
                2
            ));
        }
        assert!(is_psd_prefix_grid(&[0.0; 4], 2, 1));
        assert!(!is_psd_prefix_grid(&[1.0, 0.1, 0.2, 1.0], 2, 1));
        assert!(!is_psd_prefix_grid(&[1.0, 0.0, 0.0, f64::NAN], 2, 1));
    }

    #[test]
    fn full_geographic_landmark_extension_spans_exact_subspace() {
        let n = 36;
        let mut coords = Vec::with_capacity(n * 2);
        for i in 0..n {
            let fi = i as f64;
            coords.push(-170.0 + (fi * 37.3).rem_euclid(340.0));
            coords.push(-65.0 + (fi * 19.7).rem_euclid(130.0));
        }
        let lb = landmark_basis(
            &coords,
            n,
            2,
            Metric::GreatCircleFractionLonLat,
            &LandmarkRequest {
                avc: 0.1,
                min_avc: 1e-4,
                components: 3,
                landmarks: n,
                samples: 1,
                seed: 7,
            },
        )
        .unwrap();
        let dm = distance_matrix(&coords, n, 2, Metric::GreatCircleFractionLonLat).unwrap();
        let exact = basis::scpc_basis(&dm, n, lb.decay[0], 3).unwrap();
        // Check column-space projector P_approx == P_exact on columns 1..=3
        let mut p_approx = vec![0.0; n * n];
        let mut p_exact = vec![0.0; n * n];
        let a_sub: Vec<f64> = lb
            .weights
            .chunks_exact(4)
            .flat_map(|r| r[1..].iter().copied())
            .collect();
        let e_sub: Vec<f64> = exact
            .weights
            .chunks_exact(4)
            .flat_map(|r| r[1..].iter().copied())
            .collect();
        let mut a_sub_t = vec![0.0; 3 * n];
        let mut e_sub_t = vec![0.0; 3 * n];
        for i in 0..n {
            for c in 0..3 {
                a_sub_t[c * n + i] = a_sub[i * 3 + c];
                e_sub_t[c * n + i] = e_sub[i * 3 + c];
            }
        }
        linalg::multiply(&a_sub, &a_sub_t, &mut p_approx, n, 3, n, false);
        linalg::multiply(&e_sub, &e_sub_t, &mut p_exact, n, 3, n, false);
        for (a, b) in p_approx.iter().zip(&p_exact) {
            assert!((a - b).abs() < 1e-7, "projector diff a={a}, b={b}");
        }
    }

    #[test]
    fn normalize_coords_is_scale_invariant_at_extreme_magnitudes() {
        let base = vec![-1.5, 0.25, 0.0, -0.75, 2.0, 1.5, -0.5, 2.25, 1.25, -1.0];
        let (norm_ref, perm_ref) = normalize_coords(&base, 5, 2, Metric::Euclidean).unwrap();
        for factor in [1e-200_f64, 1e200_f64] {
            let scaled: Vec<f64> = base.iter().map(|&v| v * factor).collect();
            let (norm, perm) = normalize_coords(&scaled, 5, 2, Metric::Euclidean).unwrap();
            assert_eq!(perm, perm_ref);
            for (a, b) in norm.iter().zip(&norm_ref) {
                assert!(a.is_finite() && (a - b).abs() < 1e-12, "a={a}, b={b}");
            }
        }
        let shifted: Vec<_> = base.iter().map(|&v| 1e308 + v * 1e306).collect();
        let (norm, perm) = normalize_coords(&shifted, 5, 2, Metric::Euclidean).unwrap();
        assert_eq!(perm, perm_ref);
        for (a, b) in norm.iter().zip(&norm_ref) {
            assert!((a - b).abs() < 1e-12, "a={a}, b={b}");
        }
    }
}
