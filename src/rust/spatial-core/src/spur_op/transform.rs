//! Spatial transformations by nearest neighbors, radius, clusters, or LBM GLS.

use super::{center_columns, SpurGeometry};
use crate::distance::{Distances, Metric};
use crate::threads::{is_serial, run_with_threads};
use crate::{basis, checked_len, linalg, CoreError};
use rayon::prelude::*;
use std::collections::HashMap;

const LBM_GLS_EIG_TOL: f64 = 1e-10;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpurTransformKind {
    Nn,
    Iso { radius: f64 },
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
    pub y: &'a [f64],
    pub n_y: usize,
    pub nthreads: Option<usize>,
}

fn transform_by_neighbors(
    req: &SpurTransformRequest<'_>,
    select: impl Fn(usize, &mut Vec<usize>) -> Result<(), CoreError> + Sync,
) -> Result<Vec<f64>, CoreError> {
    let width = req.n_y;
    let mut values = vec![0.0; req.y.len()];
    let write_row = |neighbors: &mut Vec<usize>, (i, row): (usize, &mut [f64])| {
        neighbors.clear();
        select(i, neighbors)?;
        if neighbors.is_empty() {
            return Ok(());
        }
        let weight = (neighbors.len() as f64).recip();
        row.copy_from_slice(&req.y[i * width..(i + 1) * width]);
        for &j in neighbors.iter() {
            for (out, &value) in row.iter_mut().zip(&req.y[j * width..(j + 1) * width]) {
                *out -= weight * value;
            }
        }
        Ok(())
    };

    // Keep only the neighbor indices for the current row.
    if is_serial() || req.n < 64 {
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
    Ok(values)
}

fn neighbor_transform(
    req: &SpurTransformRequest<'_>,
    coords: &[f64],
) -> Result<Vec<f64>, CoreError> {
    let n = req.n;
    let distances = Distances::new(coords, n, req.d, req.metric)?;
    let distance = |i: usize, j: usize| {
        let value = distances.pair(i.max(j), i.min(j));
        if !value.is_finite() {
            return Err(CoreError::invalid_input(
                "pairwise distance exceeds the floating-point range",
            ));
        }
        Ok(value)
    };
    match req.kind {
        SpurTransformKind::Nn => transform_by_neighbors(req, |i, neighbors| {
            let mut nearest = f64::INFINITY;
            for j in 0..n {
                if i == j {
                    continue;
                }
                // Comparing raw distances preserves distinctions that can
                // disappear when dividing by a much larger diameter.
                let d = distance(i, j)?;
                if d == 0.0 {
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
                    "nearest-neighbor transform requires at least two distinct locations",
                ));
            }
            Ok(())
        }),
        SpurTransformKind::Iso { radius } => transform_by_neighbors(req, |i, neighbors| {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let d = distance(i, j)?;
                if d < radius {
                    neighbors.push(j);
                }
            }
            Ok(())
        }),
        _ => unreachable!("neighbor transform requires NN or ISO"),
    }
}

fn cluster_transform(req: &SpurTransformRequest<'_>) -> Result<Vec<f64>, CoreError> {
    let cid = req
        .cluster_id
        .ok_or_else(|| CoreError::invalid_input("cluster transformation requires cluster_id"))?;
    if cid.len() != req.n {
        return Err(CoreError::invalid_input("cluster_id length mismatch"));
    }
    let mut group_map = HashMap::new();
    let mut group_idx = Vec::with_capacity(req.n);
    for &raw_id in cid {
        let next_idx = group_map.len();
        let idx = *group_map.entry(raw_id).or_insert(next_idx);
        group_idx.push(idx);
    }
    let mut out = req.y.to_vec();
    crate::preparation::within(&mut out, req.n_y, &group_idx, false)?;
    Ok(out)
}

/// Transform the response columns and return them in row-major order.
pub fn run_spur_transform(req: &SpurTransformRequest<'_>) -> Result<Vec<f64>, CoreError> {
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
    if req.n_y == 0 || req.y.len() != checked_len(&[n, req.n_y])? {
        return Err(CoreError::invalid_input(
            "invalid Y dimensions for transform",
        ));
    }
    if req.y.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input("Y values must be finite"));
    }

    run_with_threads(req.nthreads, || match req.kind {
        SpurTransformKind::Cluster => cluster_transform(req),
        SpurTransformKind::Nn => {
            let coords = req
                .coords
                .ok_or_else(|| CoreError::invalid_input("nn transform requires coordinates"))?;
            neighbor_transform(req, coords)
        }
        SpurTransformKind::Iso { radius } => {
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
            apply_lbmgls(geom.sigma_lbm_dm(), n, req.y, req.n_y)
        }
    })?
}

pub(super) fn apply_lbmgls(
    mut sigdm: Vec<f64>,
    n: usize,
    y: &[f64],
    n_y: usize,
) -> Result<Vec<f64>, CoreError> {
    let pairs = basis::dense_eigenpairs(&mut sigdm, n, n)?;
    drop(sigdm);
    // Eigenvalues are in descending order and each eigenvector is contiguous.
    let rank = pairs
        .values
        .iter()
        .take_while(|&&value| value > LBM_GLS_EIG_TOL)
        .count();
    if rank == 0 {
        return Err(CoreError::degenerate(
            "LBM-GLS covariance has no eigenvalues above 1e-10",
        ));
    }

    // Apply V diag(lambda^-1/2) V' to Y without materializing the transform.
    let mut vt_y = vec![0.0_f64; rank * n_y];
    let mut centered = y.to_vec();
    crate::preparation::within(&mut centered, n_y, &vec![0; n], false)?;
    linalg::multiply(&pairs.vectors, &centered, &mut vt_y, rank, n, n_y, false);
    let mut v_scaled = vec![0.0_f64; n * rank];
    for (j, &value) in pairs.values[..rank].iter().enumerate() {
        let inv_sqrt = value.sqrt().recip();
        for i in 0..n {
            v_scaled[i * rank + j] = pairs.vectors[j * n + i] * inv_sqrt;
        }
    }
    let mut out = vec![0.0_f64; y.len()];
    linalg::multiply(&v_scaled, &vt_y, &mut out, n, rank, n_y, false);
    center_columns(&mut out, n, n_y);
    Ok(out)
}
