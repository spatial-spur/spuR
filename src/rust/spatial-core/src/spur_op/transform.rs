//! Spatial transformations by nearest neighbors, radius, clusters, or LBM GLS.

use super::{col_demean_matrix, symmetrize, SpurGeometry};
use crate::distance::{distance_matrix, Distances, Metric};
use crate::threads::{is_serial, run_with_threads};
use crate::{basis, checked_len, linalg, CoreError};
use rayon::prelude::*;

const LBM_GLS_EIG_TOL: f64 = 1e-10;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpurTransformKind {
    Nn { exclude_zero_distances: bool },
    Iso { radius: f64, strict_lt: bool },
    Cluster,
    LbmGls,
}

#[derive(Clone, Copy, Debug)]
pub struct SpurTransformRequest<'a> {
    pub coords: Option<&'a [f64]>,
    pub n: usize,
    pub d: usize,
    pub metric: Metric,
    pub kind: SpurTransformKind,
    pub cluster_id: Option<&'a [usize]>,
    pub y: Option<&'a [f64]>,
    pub n_y: usize,
    pub return_matrix: bool,
    pub demean_lbmgls: bool,
    pub nthreads: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct SpurTransformOutput {
    pub matrix: Option<Vec<f64>>,
    pub transformed: Option<Vec<f64>>,
}

fn transform_by_neighbors(
    req: &SpurTransformRequest<'_>,
    select: impl Fn(usize, &mut Vec<usize>) -> Result<(), CoreError> + Sync,
) -> Result<SpurTransformOutput, CoreError> {
    let n = req.n;
    let (width, y) = if req.return_matrix {
        (n, None)
    } else {
        (req.y.map_or(0, |_| req.n_y), req.y)
    };
    let mut values = vec![0.0; checked_len(&[n, width])?];
    let write_row = |neighbors: &mut Vec<usize>, (i, row): (usize, &mut [f64])| {
        neighbors.clear();
        select(i, neighbors)?;
        if neighbors.is_empty() || row.is_empty() {
            return Ok(());
        }
        let weight = (neighbors.len() as f64).recip();
        if let Some(y) = y {
            row.copy_from_slice(&y[i * width..(i + 1) * width]);
            for &j in neighbors.iter() {
                for (out, &value) in row.iter_mut().zip(&y[j * width..(j + 1) * width]) {
                    *out -= weight * value;
                }
            }
        } else {
            row[i] = 1.0;
            for &j in neighbors.iter() {
                row[j] = -weight;
            }
        }
        Ok(())
    };
    if width == 0 {
        let mut neighbors = Vec::new();
        for i in 0..n {
            write_row(&mut neighbors, (i, &mut []))?;
        }
    } else if is_serial() || n < 64 {
        let mut neighbors = Vec::new();
        for row in values.chunks_exact_mut(width).enumerate() {
            write_row(&mut neighbors, row)?;
        }
    } else {
        values
            .par_chunks_mut(width)
            .enumerate()
            .try_for_each_init(Vec::new, write_row)?;
    }

    if req.return_matrix {
        let transformed = req.y.map(|y| {
            let mut out = vec![0.0; y.len()];
            linalg::multiply(&values, y, &mut out, n, n, req.n_y, false);
            out
        });
        Ok(SpurTransformOutput {
            matrix: Some(values),
            transformed,
        })
    } else {
        Ok(SpurTransformOutput {
            matrix: None,
            transformed: y.map(|_| values),
        })
    }
}

fn neighbor_transform(
    req: &SpurTransformRequest<'_>,
    coords: &[f64],
) -> Result<SpurTransformOutput, CoreError> {
    let n = req.n;
    let distances = Distances::new(coords, n, req.d, req.metric)?;
    // Cache pairwise distances when the caller already needs a dense matrix.
    // A transformed-only request keeps only the current row's neighbor indices.
    let dense = req
        .return_matrix
        .then(|| distance_matrix(coords, n, req.d, req.metric))
        .transpose()?;
    let distance = |i: usize, j: usize| {
        let value = match &dense {
            Some(matrix) => matrix[i * n + j],
            None => distances.pair(i.max(j), i.min(j)),
        };
        if !value.is_finite() {
            return Err(CoreError::invalid_input(
                "pairwise distance exceeds the floating-point range",
            ));
        }
        Ok(value)
    };
    match req.kind {
        SpurTransformKind::Nn {
            exclude_zero_distances,
        } => {
            transform_by_neighbors(req, |i, neighbors| {
                let mut nearest = f64::INFINITY;
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    // Comparing raw distances preserves distinctions that can
                    // disappear when dividing by a much larger diameter.
                    let d = distance(i, j)?;
                    if exclude_zero_distances && d == 0.0 {
                        continue;
                    }
                    if d < nearest {
                        nearest = d;
                        neighbors.clear();
                    }
                    if d == nearest {
                        neighbors.push(j);
                    }
                }
                if neighbors.is_empty() {
                    return Err(CoreError::invalid_input(
                        "nearest-neighbor transform requires distinct locations when excluding zero distances",
                    ));
                }
                Ok(())
            })
        }
        SpurTransformKind::Iso { radius, strict_lt } => {
            transform_by_neighbors(req, |i, neighbors| {
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let d = distance(i, j)?;
                    if d > 0.0 && if strict_lt { d < radius } else { d <= radius } {
                        neighbors.push(j);
                    }
                }
                Ok(())
            })
        }
        _ => unreachable!("neighbor transform requires NN or ISO"),
    }
}

pub fn run_spur_transform(
    req: &SpurTransformRequest<'_>,
) -> Result<SpurTransformOutput, CoreError> {
    let n = req.n;
    let min_n = if matches!(req.kind, SpurTransformKind::Cluster) {
        1
    } else {
        2
    };
    if n < min_n {
        return Err(CoreError::invalid_input(if min_n == 1 {
            "cluster transform requires at least one observation"
        } else {
            "spatial transform requires at least 2 observations"
        }));
    }
    if let Some(y) = req.y {
        if req.n_y == 0 || y.len() != checked_len(&[n, req.n_y])? {
            return Err(CoreError::invalid_input(
                "invalid Y dimensions for transform",
            ));
        }
        if y.iter().any(|v| !v.is_finite()) {
            return Err(CoreError::invalid_input("Y values must be finite"));
        }
    }

    run_with_threads(req.nthreads, || match req.kind {
        SpurTransformKind::Cluster => {
            let cid = req.cluster_id.ok_or_else(|| {
                CoreError::invalid_input("cluster transformation requires cluster_id")
            })?;
            if cid.len() != n {
                return Err(CoreError::invalid_input("cluster_id length mismatch"));
            }
            let mut group_map = std::collections::HashMap::new();
            let mut group_idx = Vec::with_capacity(n);
            let mut counts = Vec::new();
            for &raw_id in cid {
                let next_idx = group_map.len();
                let idx = *group_map.entry(raw_id).or_insert(next_idx);
                if idx == counts.len() {
                    counts.push(1usize);
                } else {
                    counts[idx] += 1;
                }
                group_idx.push(idx);
            }
            let n_groups = counts.len();
            let transformed = if let Some(y) = req.y {
                let ny = req.n_y;
                let sums_len = checked_len(&[n_groups, ny])?;
                let out_len = checked_len(&[n, ny])?;
                let mut sums = vec![0.0_f64; sums_len];
                for (i, &g) in group_idx.iter().enumerate() {
                    for col in 0..ny {
                        sums[g * ny + col] += y[i * ny + col];
                    }
                }
                let mut out = vec![0.0_f64; out_len];
                for (i, &g) in group_idx.iter().enumerate() {
                    let inv_cnt = (counts[g] as f64).recip();
                    for col in 0..ny {
                        out[i * ny + col] = y[i * ny + col] - sums[g * ny + col] * inv_cnt;
                    }
                }
                Some(out)
            } else {
                None
            };
            let matrix = if req.return_matrix {
                let mut mat = vec![0.0_f64; checked_len(&[n, n])?];
                for i in 0..n {
                    let gi = group_idx[i];
                    let w = (counts[gi] as f64).recip();
                    for j in 0..n {
                        if group_idx[j] == gi {
                            mat[i * n + j] = -w;
                        }
                    }
                    mat[i * n + i] += 1.0;
                }
                Some(mat)
            } else {
                None
            };
            Ok(SpurTransformOutput {
                matrix,
                transformed,
            })
        }
        SpurTransformKind::Nn { .. } => {
            let coords = req
                .coords
                .ok_or_else(|| CoreError::invalid_input("nn transform requires coordinates"))?;
            neighbor_transform(req, coords)
        }
        SpurTransformKind::Iso { radius, .. } => {
            if !radius.is_finite() || radius <= 0.0 {
                return Err(CoreError::invalid_input(
                    "iso transform radius must be positive and finite",
                ));
            }
            let coords = req
                .coords
                .ok_or_else(|| CoreError::invalid_input("iso transform requires coordinates"))?;
            neighbor_transform(req, coords)
        }
        SpurTransformKind::LbmGls => {
            let coords = req
                .coords
                .ok_or_else(|| CoreError::invalid_input("lbmgls transform requires coordinates"))?;
            let geom = SpurGeometry::from_coords(coords, n, req.d, req.metric)?;
            apply_lbmgls_with_geom(&geom, req.y, req.n_y, req.return_matrix, req.demean_lbmgls)
        }
    })?
}

pub(super) fn apply_lbmgls_with_geom(
    geom: &SpurGeometry,
    y: Option<&[f64]>,
    n_y: usize,
    return_matrix: bool,
    demean_lbmgls: bool,
) -> Result<SpurTransformOutput, CoreError> {
    let n = geom.n;
    let mut sigdm = geom.sigma_lbm_dm();
    let pairs = basis::dense_eigenpairs(&mut sigdm, n, n)?;
    let keep: Vec<usize> = pairs
        .values
        .iter()
        .enumerate()
        .filter(|(_, &v)| v > LBM_GLS_EIG_TOL)
        .map(|(idx, _)| idx)
        .collect();
    let rank = keep.len();
    if rank == 0 {
        return Err(CoreError::degenerate(
            "LBM-GLS covariance has no eigenvalues above 1e-10",
        ));
    }

    // Build V_keep (n x rank) and V_scaled = V_keep * diag(eval^{-1/2}) (n x rank).
    let mut v_keep = vec![0.0_f64; n * rank];
    let mut v_scaled = vec![0.0_f64; n * rank];
    for (col_idx, &eig_idx) in keep.iter().enumerate() {
        let inv_sqrt = pairs.values[eig_idx].sqrt().recip();
        let vec_slice = &pairs.vectors[eig_idx * n..(eig_idx + 1) * n];
        for i in 0..n {
            let v = vec_slice[i];
            v_keep[i * rank + col_idx] = v;
            v_scaled[i * rank + col_idx] = v * inv_sqrt;
        }
    }

    let matrix = if return_matrix {
        let mut v_keep_t = vec![0.0_f64; rank * n];
        for i in 0..n {
            for r in 0..rank {
                v_keep_t[r * n + i] = v_keep[i * rank + r];
            }
        }
        let mut mat = vec![0.0_f64; n * n];
        linalg::multiply(&v_scaled, &v_keep_t, &mut mat, n, rank, n, false);
        symmetrize(&mut mat, n);
        Some(mat)
    } else {
        None
    };

    let transformed = y.map(|y_mat| {
        let mut vt_y = vec![0.0_f64; rank * n_y];
        linalg::multiply(&v_keep, y_mat, &mut vt_y, rank, n, n_y, true);
        let mut out = vec![0.0_f64; n * n_y];
        linalg::multiply(&v_scaled, &vt_y, &mut out, n, rank, n_y, false);
        if demean_lbmgls {
            col_demean_matrix(&out, n, n_y)
        } else {
            out
        }
    });

    Ok(SpurTransformOutput {
        matrix,
        transformed,
    })
}
