//! Shared `extendr` R FFI adapter for `scpcR` and `spuR`.

use extendr_api::prelude::*;
use spatial_core::{
    run_approx_scpc, run_exact_scpc, run_spur_halflife, run_spur_i0, run_spur_i0_resid,
    run_spur_i1, run_spur_i1_resid, run_spur_pipeline, run_spur_transform, ApproxScpcRequest,
    CoreError, ErrorKind, ExactScpcOutput, ExactScpcRequest, Metric, SpurHalfLifeOutput,
    SpurHalfLifeRequest, SpurI0Output, SpurI1Output, SpurPipelineRequest, SpurTestRequest,
    SpurTransformKind, SpurTransformRequest,
};

fn core_err(err: CoreError) -> Error {
    match err.kind() {
        ErrorKind::InvalidInput | ErrorKind::TiedEigenspace | ErrorKind::InvalidSampledCovariance => {
            Error::Other(err.message().to_string())
        }
        ErrorKind::UnsupportedCase => Error::Other(format!("unsupported case: {}", err.message())),
        ErrorKind::CalibrationFailed => {
            Error::Other(format!("calibration failed: {}", err.message()))
        }
        ErrorKind::DegenerateVariance => {
            Error::Other(format!("degenerate variance: {}", err.message()))
        }
        ErrorKind::InternalFailure => Error::Other(format!("internal failure: {}", err.message())),
    }
}

/// The extendr result envelope carries a classed R condition for approximate
/// SCPC outcomes that the host needs to distinguish at dispatch.
pub struct ApproxScpcError {
    message: String,
    class: Option<&'static str>,
}

impl From<Error> for ApproxScpcError {
    fn from(error: Error) -> Self {
        Self {
            message: error.to_string(),
            class: None,
        }
    }
}

impl From<CoreError> for ApproxScpcError {
    fn from(error: CoreError) -> Self {
        let class = match error.kind() {
            ErrorKind::TiedEigenspace => Some("spatial_tied_eigenspace"),
            ErrorKind::InvalidSampledCovariance => Some("spatial_invalid_sampled_covariance"),
            _ => None,
        };
        Self {
            message: core_err(error).to_string(),
            class,
        }
    }
}

impl From<ApproxScpcError> for Robj {
    fn from(error: ApproxScpcError) -> Self {
        if let Some(class) = error.class {
            let mut condition = list!(message = error.message, call = NULL);
            condition
                .set_class([class, "error", "condition"])
                .expect("failed to set the SCPC error condition class");
            condition.into()
        } else {
            error.message.into()
        }
    }
}

/// Convert an R column-major matrix buffer (`nrows x ncols`) into a contiguous row-major `Vec<f64>`.
pub fn col_major_to_row_major(
    data: &[f64],
    nrows: usize,
    ncols: usize,
) -> extendr_api::Result<Vec<f64>> {
    let expected = nrows
        .checked_mul(ncols)
        .ok_or_else(|| Error::Other("matrix dimensions overflow usize".into()))?;
    if data.len() != expected {
        return Err(Error::Other(format!(
            "column-major buffer length {} does not match {} x {}",
            data.len(),
            nrows,
            ncols
        )));
    }
    let mut row_major = vec![0.0; expected];
    for col in 0..ncols {
        let src_col = &data[col * nrows..(col + 1) * nrows];
        for row in 0..nrows {
            row_major[row * ncols + col] = src_col[row];
        }
    }
    Ok(row_major)
}

/// Convert a contiguous row-major slice (`nrows x ncols`) into an R column-major `RMatrix<f64>`.
pub fn row_major_to_rmatrix(
    row_major: &[f64],
    nrows: usize,
    ncols: usize,
) -> extendr_api::Result<RMatrix<f64>> {
    let expected = nrows
        .checked_mul(ncols)
        .ok_or_else(|| Error::Other("matrix dimensions overflow usize".into()))?;
    if row_major.len() != expected {
        return Err(Error::Other(format!(
            "row-major buffer length {} does not match {} x {}",
            row_major.len(),
            nrows,
            ncols
        )));
    }
    if i32::try_from(nrows).is_err() || i32::try_from(ncols).is_err() {
        return Err(Error::Other(
            "matrix dimensions exceed R 32-bit dim limit".into(),
        ));
    }
    let mut mat = RMatrix::<f64>::new(nrows, ncols);
    let dst = mat.data_mut();
    for col in 0..ncols {
        let dst_col = &mut dst[col * nrows..(col + 1) * nrows];
        for row in 0..nrows {
            dst_col[row] = row_major[row * ncols + col];
        }
    }
    Ok(mat)
}

/// Bridge one whole exact SCPC inference call from R `SEXP` wrappers to `spatial_core::run_exact_scpc`.
///
/// All R column-major matrices are copied into row-major Rust buffers on the R main thread
/// prior to entering Rayon-parallel numerical routines, and R output matrices are allocated
/// back on the R main thread after `run_exact_scpc` completes.
#[allow(clippy::too_many_arguments)]
pub fn run_exact_scpc_r(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_exact_scpc_r_impl(
        coords,
        latlong,
        avc,
        uncond,
        cvs,
        coef,
        wj,
        xj,
        xjs,
        proj_u,
        proj_v,
        cv_tolerance,
        nthreads,
        Nullable::Null,
        Nullable::Null,
    )
}

#[cfg_attr(feature = "exports", extendr)]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_model_scpc(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    cvs: bool,
    coef: Doubles,
    scores: RMatrix<f64>,
    bread: RMatrix<f64>,
    design: Nullable<RMatrix<f64>>,
    instruments: Nullable<RMatrix<f64>>,
    direction_design: Nullable<RMatrix<f64>>,
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
    nthreads: Nullable<i32>,
    seed: Nullable<f64>,
) -> std::result::Result<List, ApproxScpcError> {
    use spatial_core::preparation::{run_model_scpc, ScpcPreparation};

    let n = coords.nrows();
    let d = coords.ncols();
    let p = scores.ncols();
    let k = coef.len();
    if bread.nrows() != k || bread.ncols() != p {
        return Err(Error::Other("score or bread dimensions do not match the model".into()).into());
    }
    let cluster = parse_r_index_vec(cluster_id, n, "cluster_id", "clusters")?;
    let n_obs = cluster.as_ref().map_or(n, Vec::len);
    if scores.nrows() != n_obs {
        return Err(Error::Other("scores must have one row per observation".into()).into());
    }
    let fe = parse_r_index_vec(fe_id, n_obs, "fe_id", "observations")?;
    let read_design = |matrix: Nullable<RMatrix<f64>>| -> extendr_api::Result<(Option<Vec<f64>>, usize)> {
        match matrix {
            Nullable::Null => Ok((None, 0)),
            Nullable::NotNull(m) => {
                if m.nrows() != n_obs {
                    return Err(Error::Other("a model matrix has the wrong number of observations".into()));
                }
                Ok((Some(col_major_to_row_major(m.data(), n_obs, m.ncols())?), m.ncols()))
            }
        }
    };
    let (x, design_cols) = read_design(design)?;
    let (z, instrument_cols) = read_design(instruments)?;
    let (directions, direction_cols) = read_design(direction_design)?;
    if directions.is_some() && direction_cols != p {
        return Err(Error::Other("the fitted direction matrix must have one column per score".into()).into());
    }
    let seed = match seed {
        Nullable::Null => None,
        Nullable::NotNull(s) if (0.0..4_294_967_296.0).contains(&s) && s.fract() == 0.0 => Some(s as u32),
        _ => return Err(Error::Other("seed must be an integer in [0, 2^32)".into()).into()),
    };
    let coord_data = col_major_to_row_major(coords.data(), n, d)?;
    let obs_scores = col_major_to_row_major(scores.data(), n_obs, p)?;
    let bread_rows = col_major_to_row_major(bread.data(), k, p)?;
    let coefficients = coef.as_robj().as_real_slice()
        .ok_or_else(|| Error::ExpectedReal(coef.as_robj().clone()))?;
    let model = ScpcPreparation {
        n_sites: n, n_obs, n_scores: p, coefficients, bread: &bread_rows,
        scores: &obs_scores, design: x.as_deref(), design_cols,
        instruments: z.as_deref(), instrument_cols, direction_design: directions.as_deref(),
        cluster_ids: cluster.as_deref(), fixed_effects: fe.as_deref(),
    };
    let out = run_model_scpc(model, &coord_data, d,
        if latlong { Metric::GreatCircleFractionLonLat } else { Metric::Euclidean },
        avc, cvs, parse_r_nthreads(nthreads)?, seed)?;
    Ok(build_r_scpc_list(out, n, k)?)
}

struct ParsedRScpcInputs {
    coords: Vec<f64>,
    n: usize,
    d: usize,
    metric: Metric,
    coef: Vec<f64>,
    k_use: usize,
    wj: Vec<f64>,
    cluster_id: Option<Vec<usize>>,
    fe_id: Option<Vec<usize>>,
    xj: Option<Vec<f64>>,
    xjs: Option<Vec<f64>>,
    proj_u: Option<Vec<f64>>,
    proj_v: Option<Vec<f64>>,
    proj_rank: usize,
    nthreads: Option<usize>,
}

#[allow(clippy::too_many_arguments)]
fn parse_r_scpc_inputs(
    coords: RMatrix<f64>,
    latlong: bool,
    uncond: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    nthreads: Nullable<i32>,
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
) -> extendr_api::Result<ParsedRScpcInputs> {
    let n = coords.nrows();
    let d = coords.ncols();
    let k_use = coef.len();

    if wj.nrows() != n || wj.ncols() != k_use {
        return Err(Error::Other(format!(
            "wj dimensions ({} x {}) do not match n = {} and k_use = {}",
            wj.nrows(),
            wj.ncols(),
            n,
            k_use
        )));
    }

    let nthreads = parse_r_nthreads(nthreads)?;
    let cluster_vec = parse_r_index_vec(cluster_id, n, "cluster_id", "clusters")?;
    let m_obs = cluster_vec.as_ref().map_or(n, Vec::len);
    let fe_vec = parse_r_index_vec(fe_id, m_obs, "fe_id", "observations")?;

    let coef_slice = coef
        .as_robj()
        .as_real_slice()
        .ok_or_else(|| Error::ExpectedReal(coef.as_robj().clone()))?
        .to_vec();

    let coords_rm = col_major_to_row_major(coords.data(), n, d)?;
    let wj_rm = col_major_to_row_major(wj.data(), n, k_use)?;

    let (xj_rm, xjs_rm, proj_u_rm, proj_v_rm, proj_rank) = if uncond {
        (None, None, None, None, 0usize)
    } else {
        let xj_mat = match xj {
            Nullable::NotNull(m) => m,
            Nullable::Null => return Err(Error::Other("conditional SCPC requires `xj` matrix".into())),
        };
        let xjs_mat = match xjs {
            Nullable::NotNull(m) => m,
            Nullable::Null => {
                return Err(Error::Other("conditional SCPC requires `xjs` matrix".into()))
            }
        };
        if xj_mat.nrows() != m_obs || xj_mat.ncols() != k_use {
            return Err(Error::Other(format!(
                "xj dimensions ({} x {}) do not match m_obs = {} and k_use = {}",
                xj_mat.nrows(),
                xj_mat.ncols(),
                m_obs,
                k_use
            )));
        }
        if xjs_mat.nrows() != m_obs || xjs_mat.ncols() != k_use {
            return Err(Error::Other(format!(
                "xjs dimensions ({} x {}) do not match m_obs = {} and k_use = {}",
                xjs_mat.nrows(),
                xjs_mat.ncols(),
                m_obs,
                k_use
            )));
        }
        let xj_rm = col_major_to_row_major(xj_mat.data(), m_obs, k_use)?;
        let xjs_rm = col_major_to_row_major(xjs_mat.data(), m_obs, k_use)?;

        match (proj_u, proj_v) {
            (Nullable::Null, Nullable::Null) => (Some(xj_rm), Some(xjs_rm), None, None, 0usize),
            (Nullable::NotNull(u_mat), Nullable::NotNull(v_mat)) => {
                let rank = u_mat.ncols();
                if u_mat.nrows() != m_obs || v_mat.nrows() != m_obs || v_mat.ncols() != rank {
                    return Err(Error::Other(
                        "conditional projection U/V dimensions do not match".into(),
                    ));
                }
                if rank == 0 {
                    (Some(xj_rm), Some(xjs_rm), None, None, 0usize)
                } else {
                    let u_slice = u_mat.data();
                    let v_slice = v_mat.data();
                    let u_rm = col_major_to_row_major(u_slice, m_obs, rank)?;
                    let v_rm = if std::ptr::eq(u_slice.as_ptr(), v_slice.as_ptr()) {
                        None
                    } else {
                        Some(col_major_to_row_major(v_slice, m_obs, rank)?)
                    };
                    (Some(xj_rm), Some(xjs_rm), Some(u_rm), v_rm, rank)
                }
            }
            _ => {
                return Err(Error::Other(
                    "proj_u and proj_v must both be provided or both be NULL".into(),
                ))
            }
        }
    };

    let metric = if latlong {
        Metric::GreatCircleFractionLonLat
    } else {
        Metric::Euclidean
    };

    Ok(ParsedRScpcInputs {
        coords: coords_rm,
        n,
        d,
        metric,
        coef: coef_slice,
        k_use,
        wj: wj_rm,
        cluster_id: cluster_vec,
        fe_id: fe_vec,
        xj: xj_rm,
        xjs: xjs_rm,
        proj_u: proj_u_rm,
        proj_v: proj_v_rm,
        proj_rank,
        nthreads,
    })
}

fn build_r_scpc_list(out: ExactScpcOutput, n: usize, k_use: usize) -> extendr_api::Result<List> {
    let scpcstats = row_major_to_rmatrix(&out.scpcstats, k_use, 6)?;
    let scpccvs: Robj = match out.scpccvs {
        Some(ref cvs_rm) => row_major_to_rmatrix(cvs_rm, k_use, 4)?.into(),
        None => r!(NULL),
    };
    let wfin = row_major_to_rmatrix(&out.wfin, n, out.q + 1)?;
    let q_i32 = i32::try_from(out.q)
        .map_err(|_| Error::Other("selected q exceeds 32-bit R integer".into()))?;

    List::from_names_and_values(
        [
            "scpcstats",
            "scpccvs",
            "wfin",
            "c0",
            "cmax",
            "cvfin",
            "q",
            "tied",
        ],
        [
            Robj::from(scpcstats),
            scpccvs,
            Robj::from(wfin),
            Robj::from(out.c0),
            Robj::from(out.cmax),
            Robj::from(out.cvfin),
            Robj::from(q_i32),
            Robj::from(out.tied),
        ],
    )
}

#[allow(clippy::too_many_arguments)]
pub fn run_exact_scpc_r_impl(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
) -> extendr_api::Result<List> {
    let inputs = parse_r_scpc_inputs(
        coords, latlong, uncond, coef, wj, xj, xjs, proj_u, proj_v, nthreads, cluster_id, fe_id,
    )?;
    let proj_u_ref = inputs.proj_u.as_deref();
    let proj_v_ref = inputs.proj_v.as_deref().or(proj_u_ref);

    let req = ExactScpcRequest {
        coords: &inputs.coords,
        n: inputs.n,
        d: inputs.d,
        metric: inputs.metric,
        avc,
        uncond,
        cvs,
        coef: &inputs.coef,
        wj: &inputs.wj,
        cluster_id: inputs.cluster_id.as_deref(),
        fe_id: inputs.fe_id.as_deref(),
        xj: inputs.xj.as_deref(),
        xjs: inputs.xjs.as_deref(),
        proj_u: proj_u_ref,
        proj_v: proj_v_ref,
        proj_rank: inputs.proj_rank,
        cv_tolerance,
        nthreads: inputs.nthreads,
    };

    let out = run_exact_scpc(&req).map_err(core_err)?;
    build_r_scpc_list(out, inputs.n, inputs.k_use)
}

#[allow(clippy::too_many_arguments)]
pub fn run_approx_scpc_r_impl(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
    seed: f64,
    landmarks: i32,
    samples: i32,
    pair_count: i32,
) -> std::result::Result<List, ApproxScpcError> {
    if !seed.is_finite() || seed < 0.0 || seed >= 4_294_967_296.0 || seed.fract() != 0.0 {
        return Err(
            Error::Other("`seed` must be an integer-valued number in [0, 2^32).".into()).into(),
        );
    }
    if landmarks < 2 || samples < 1 || pair_count < 1 {
        return Err(Error::Other("invalid Nystrom landmark or sampling parameters.".into()).into());
    }
    let inputs = parse_r_scpc_inputs(
        coords, latlong, uncond, coef, wj, xj, xjs, proj_u, proj_v, nthreads, cluster_id, fe_id,
    )?;
    let proj_u_ref = inputs.proj_u.as_deref();
    let proj_v_ref = inputs.proj_v.as_deref().or(proj_u_ref);

    let req = ApproxScpcRequest {
        coords: &inputs.coords,
        n: inputs.n,
        d: inputs.d,
        metric: inputs.metric,
        avc,
        uncond,
        cvs,
        coef: &inputs.coef,
        wj: &inputs.wj,
        cluster_id: inputs.cluster_id.as_deref(),
        fe_id: inputs.fe_id.as_deref(),
        xj: inputs.xj.as_deref(),
        xjs: inputs.xjs.as_deref(),
        proj_u: proj_u_ref,
        proj_v: proj_v_ref,
        proj_rank: inputs.proj_rank,
        cv_tolerance,
        nthreads: inputs.nthreads,
        seed: seed as u32,
        landmarks: landmarks as usize,
        samples: samples as usize,
        pair_count: pair_count as usize,
    };

    let out = run_approx_scpc(&req)?;
    Ok(build_r_scpc_list(out, inputs.n, inputs.k_use)?)
}

#[cfg(feature = "exports")]
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_exact_scpc(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_exact_scpc_r(
        coords,
        latlong,
        avc,
        uncond,
        cvs,
        coef,
        wj,
        xj,
        xjs,
        proj_u,
        proj_v,
        cv_tolerance,
        nthreads,
    )
}

#[cfg(feature = "exports")]
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_exact_scpc_clustered(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
    cluster_id: Nullable<Integers>,
) -> extendr_api::Result<List> {
    run_exact_scpc_r_impl(
        coords,
        latlong,
        avc,
        uncond,
        cvs,
        coef,
        wj,
        xj,
        xjs,
        proj_u,
        proj_v,
        cv_tolerance,
        nthreads,
        cluster_id,
        Nullable::Null,
    )
}

#[cfg(feature = "exports")]
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_exact_scpc_full(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
) -> extendr_api::Result<List> {
    run_exact_scpc_r_impl(
        coords,
        latlong,
        avc,
        uncond,
        cvs,
        coef,
        wj,
        xj,
        xjs,
        proj_u,
        proj_v,
        cv_tolerance,
        nthreads,
        cluster_id,
        fe_id,
    )
}

#[cfg(feature = "exports")]
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_approx_scpc(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
    seed: f64,
    landmarks: i32,
    samples: i32,
    pair_count: i32,
) -> std::result::Result<List, ApproxScpcError> {
    run_approx_scpc_r_impl(
        coords,
        latlong,
        avc,
        uncond,
        cvs,
        coef,
        wj,
        xj,
        xjs,
        proj_u,
        proj_v,
        cv_tolerance,
        nthreads,
        cluster_id,
        fe_id,
        seed,
        landmarks,
        samples,
        pair_count,
    )
}

fn parse_r_nthreads(nthreads: Nullable<i32>) -> extendr_api::Result<Option<usize>> {
    match nthreads {
        Nullable::Null => Ok(None),
        Nullable::NotNull(t) => {
            if t == i32::MIN || t <= 0 {
                return Err(Error::Other(
                    "`nthreads` must be a positive integer.".into(),
                ));
            }
            Ok(Some(t as usize))
        }
    }
}

fn parse_r_index_vec(
    ids: Nullable<Integers>,
    upper_bound: usize,
    name: &str,
    label: &str,
) -> extendr_api::Result<Option<Vec<usize>>> {
    match ids {
        Nullable::Null => Ok(None),
        Nullable::NotNull(ints) => {
            let slice = ints
                .as_robj()
                .as_integer_slice()
                .ok_or_else(|| Error::ExpectedInteger(ints.as_robj().clone()))?;
            let mut out = Vec::with_capacity(slice.len());
            for &val in slice {
                if val == i32::MIN || val < 0 || (val as usize) >= upper_bound {
                    return Err(Error::Other(format!(
                        "`{name}` entries must be non-NA integers in 0..{} for {label} (got {val})",
                        upper_bound
                    )));
                }
                out.push(val as usize);
            }
            Ok(Some(out))
        }
    }
}

fn build_r_spur_i0_list(out: SpurI0Output, n_y: usize) -> extendr_api::Result<List> {
    let n_rho = out.rho_grid.len();
    let cvalue_mat = row_major_to_rmatrix(&out.cvalue_mat, n_rho, 3)?;
    let pvalue_mat = row_major_to_rmatrix(&out.pvalue_mat, n_rho, n_y)?;
    List::from_names_and_values(
        [
            "LR",
            "pvalue",
            "cvalue",
            "ha_parm",
            "cvalue_mat",
            "pvalue_mat",
            "rho_grid",
            "tied",
        ],
        [
            Robj::from(out.lr),
            Robj::from(out.pvalue),
            Robj::from(out.cvalue),
            Robj::from(out.ha_param),
            Robj::from(cvalue_mat),
            Robj::from(pvalue_mat),
            Robj::from(out.rho_grid),
            Robj::from(out.tied),
        ],
    )
}

fn build_r_spur_i1_list(out: SpurI1Output) -> extendr_api::Result<List> {
    List::from_names_and_values(
        ["LR", "pvalue", "cv_vec", "ha_parm", "tied"],
        [
            Robj::from(out.lr),
            Robj::from(out.pvalue),
            Robj::from(out.cvalue),
            Robj::from(out.ha_param),
            Robj::from(out.tied),
        ],
    )
}

fn build_r_spur_hl_list(out: SpurHalfLifeOutput) -> extendr_api::Result<List> {
    List::from_names_and_values(
        [
            "hl_ci_lower",
            "hl_ci_upper",
            "max_dist",
            "pv_mat",
            "hl_grid_ho",
            "tied",
        ],
        [
            Robj::from(out.ci_lower),
            Robj::from(out.ci_upper),
            Robj::from(out.max_dist),
            Robj::from(out.pv_mat),
            Robj::from(out.hl_grid_ho),
            Robj::from(out.tied),
        ],
    )
}

fn run_spur_test_r_impl(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: Nullable<RMatrix<f64>>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
    is_i0: bool,
) -> extendr_api::Result<List> {
    let nthreads_usize = parse_r_nthreads(nthreads)?;
    let n = coords.nrows();
    let d = coords.ncols();
    if y.nrows() != n || y.ncols() == 0 {
        return Err(Error::Other(format!(
            "Y dimensions ({} x {}) incompatible with coords rows n = {n}",
            y.nrows(),
            y.ncols()
        )));
    }
    let n_y = y.ncols();
    let q = emat.nrows();
    let nrep = emat.ncols();
    let metric = if latlong {
        Metric::GreatCircleFractionLonLat
    } else {
        Metric::Euclidean
    };

    let coords_rm = col_major_to_row_major(coords.data(), n, d)?;
    let y_rm = col_major_to_row_major(y.data(), n, n_y)?;
    let emat_rm = col_major_to_row_major(emat.data(), q, nrep)?;
    let (x_rm, k_x) = match x_in {
        Nullable::Null => (None, 0usize),
        Nullable::NotNull(xm) => {
            if xm.nrows() != n || xm.ncols() == 0 {
                return Err(Error::Other(format!(
                    "X_in dimensions ({} x {}) incompatible with n = {n}",
                    xm.nrows(),
                    xm.ncols()
                )));
            }
            let kx = xm.ncols();
            (Some(col_major_to_row_major(xm.data(), n, kx)?), kx)
        }
    };

    let req = SpurTestRequest {
        coords: &coords_rm,
        n,
        d,
        metric,
        y: &y_rm,
        n_y,
        x_in: x_rm.as_deref(),
        k_x,
        emat: &emat_rm,
        q,
        nrep,
        nthreads: nthreads_usize,
    };

    if is_i0 {
        let out = if req.x_in.is_some() {
            run_spur_i0_resid(&req)
        } else {
            run_spur_i0(&req)
        }
        .map_err(core_err)?;
        build_r_spur_i0_list(out, n_y)
    } else {
        let out = if req.x_in.is_some() {
            run_spur_i1_resid(&req)
        } else {
            run_spur_i1(&req)
        }
        .map_err(core_err)?;
        build_r_spur_i1_list(out)
    }
}

pub fn run_spur_i0_r(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_test_r_impl(coords, latlong, y, Nullable::Null, emat, nthreads, true)
}

pub fn run_spur_i1_r(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_test_r_impl(coords, latlong, y, Nullable::Null, emat, nthreads, false)
}

pub fn run_spur_i0_resid_r(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_test_r_impl(coords, latlong, y, Nullable::NotNull(x_in), emat, nthreads, true)
}

pub fn run_spur_i1_resid_r(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_test_r_impl(coords, latlong, y, Nullable::NotNull(x_in), emat, nthreads, false)
}

pub fn run_spur_halflife_r(
    coords: RMatrix<f64>,
    latlong: bool,
    y: Doubles,
    emat: RMatrix<f64>,
    level: f64,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    let nthreads_usize = parse_r_nthreads(nthreads)?;
    let n = coords.nrows();
    let d = coords.ncols();
    let y_slice = y
        .as_robj()
        .as_real_slice()
        .ok_or_else(|| Error::ExpectedReal(y.as_robj().clone()))?;
    if y_slice.len() != n {
        return Err(Error::Other(format!(
            "y length {} does not match coords rows n = {n}",
            y_slice.len()
        )));
    }
    let q = emat.nrows();
    let nrep = emat.ncols();
    let metric = if latlong {
        Metric::GreatCircleFractionLonLat
    } else {
        Metric::Euclidean
    };
    let coords_rm = col_major_to_row_major(coords.data(), n, d)?;
    let emat_rm = col_major_to_row_major(emat.data(), q, nrep)?;
    let out = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords_rm,
        n,
        d,
        metric,
        y: y_slice,
        emat: &emat_rm,
        q,
        nrep,
        level,
        nthreads: nthreads_usize,
    })
    .map_err(core_err)?;
    build_r_spur_hl_list(out)
}

#[allow(clippy::too_many_arguments)]
pub fn run_spur_transform_r(
    kind: &str,
    n: i32,
    y: Nullable<RMatrix<f64>>,
    coords: Nullable<RMatrix<f64>>,
    latlong: bool,
    radius: Nullable<f64>,
    strict_lt: bool,
    exclude_zero_distances: bool,
    cluster_id: Nullable<Integers>,
    return_matrix: bool,
    demean_lbmgls: bool,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    let n_usize = usize::try_from(n)
        .map_err(|_| Error::Other("transform requires a nonnegative sample size".into()))?;
    let nthreads_usize = parse_r_nthreads(nthreads)?;
    let tr_kind = match kind {
        "nn" => SpurTransformKind::Nn {
            exclude_zero_distances,
        },
        "iso" => {
            let r = match radius {
                Nullable::NotNull(v) => v,
                Nullable::Null => return Err(Error::Other("iso transform requires radius".into())),
            };
            SpurTransformKind::Iso {
                radius: r,
                strict_lt,
            }
        }
        "cluster" => SpurTransformKind::Cluster,
        "lbmgls" => SpurTransformKind::LbmGls,
        other => {
            return Err(Error::Other(format!(
                "unknown transform kind '{other}' (expected nn, iso, cluster, or lbmgls)"
            )));
        }
    };

    let (coords_rm, d, metric) = match coords {
        Nullable::Null => (None, 0usize, Metric::Euclidean),
        Nullable::NotNull(cm) => {
            if cm.nrows() != n_usize {
                return Err(Error::Other(format!(
                    "coords rows {} does not match expected n = {n_usize}",
                    cm.nrows()
                )));
            }
            let dc = cm.ncols();
            let met = if latlong {
                Metric::GreatCircleFractionLonLat
            } else {
                Metric::Euclidean
            };
            (Some(col_major_to_row_major(cm.data(), n_usize, dc)?), dc, met)
        }
    };

    let cluster_vec = parse_r_index_vec(cluster_id, n_usize, "cluster_id", "observations")?;
    let (y_rm, n_y) = match y {
        Nullable::Null => (None, 0usize),
        Nullable::NotNull(ym) => {
            if ym.nrows() != n_usize || ym.ncols() == 0 {
                return Err(Error::Other(format!(
                    "Y dimensions ({} x {}) incompatible with n = {n_usize}",
                    ym.nrows(),
                    ym.ncols()
                )));
            }
            let ny = ym.ncols();
            (Some(col_major_to_row_major(ym.data(), n_usize, ny)?), ny)
        }
    };

    let out = run_spur_transform(&SpurTransformRequest {
        coords: coords_rm.as_deref(),
        n: n_usize,
        d,
        metric,
        kind: tr_kind,
        cluster_id: cluster_vec.as_deref(),
        y: y_rm.as_deref(),
        n_y,
        return_matrix,
        demean_lbmgls,
        nthreads: nthreads_usize,
    })
    .map_err(core_err)?;

    let mat_obj: Robj = match out.matrix {
        Some(ref m) => row_major_to_rmatrix(m, n_usize, n_usize)?.into(),
        None => r!(NULL),
    };
    let tr_obj: Robj = match out.transformed {
        Some(ref t) => row_major_to_rmatrix(t, n_usize, n_y)?.into(),
        None => r!(NULL),
    };
    List::from_names_and_values(["matrix", "transformed"], [mat_obj, tr_obj])
}

#[allow(clippy::too_many_arguments)]
pub fn run_spur_pipeline_r(
    coords: RMatrix<f64>,
    latlong: bool,
    y: Doubles,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    transform_vars: RMatrix<f64>,
    demean_lbmgls: bool,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    let nthreads_usize = parse_r_nthreads(nthreads)?;
    let n = coords.nrows();
    let d = coords.ncols();
    let y_slice = y
        .as_robj()
        .as_real_slice()
        .ok_or_else(|| Error::ExpectedReal(y.as_robj().clone()))?;
    if y_slice.len() != n {
        return Err(Error::Other(format!(
            "y length {} does not match coords rows n = {n}",
            y_slice.len()
        )));
    }
    if x_in.nrows() != n || x_in.ncols() == 0 {
        return Err(Error::Other(format!(
            "X_in dimensions ({} x {}) incompatible with n = {n}",
            x_in.nrows(),
            x_in.ncols()
        )));
    }
    if transform_vars.nrows() != n || transform_vars.ncols() == 0 {
        return Err(Error::Other(format!(
            "transform_vars dimensions ({} x {}) incompatible with n = {n}",
            transform_vars.nrows(),
            transform_vars.ncols()
        )));
    }
    let k_x = x_in.ncols();
    let n_tr = transform_vars.ncols();
    let q = emat.nrows();
    let nrep = emat.ncols();
    let metric = if latlong {
        Metric::GreatCircleFractionLonLat
    } else {
        Metric::Euclidean
    };

    let coords_rm = col_major_to_row_major(coords.data(), n, d)?;
    let x_rm = col_major_to_row_major(x_in.data(), n, k_x)?;
    let emat_rm = col_major_to_row_major(emat.data(), q, nrep)?;
    let tr_rm = col_major_to_row_major(transform_vars.data(), n, n_tr)?;

    let out = run_spur_pipeline(&SpurPipelineRequest {
        coords: &coords_rm,
        n,
        d,
        metric,
        y: y_slice,
        x_in: &x_rm,
        k_x,
        emat: &emat_rm,
        q,
        nrep,
        transform_vars: &tr_rm,
        n_transform_vars: n_tr,
        demean_lbmgls,
        nthreads: nthreads_usize,
    })
    .map_err(core_err)?;

    let tr_mat = row_major_to_rmatrix(&out.transformed_vars, n, n_tr)?;
    List::from_names_and_values(
        ["i0", "i1", "i0resid", "i1resid", "transformed_vars"],
        [
            Robj::from(build_r_spur_i0_list(out.i0, 1)?),
            Robj::from(build_r_spur_i1_list(out.i1)?),
            Robj::from(build_r_spur_i0_list(out.i0resid, 1)?),
            Robj::from(build_r_spur_i1_list(out.i1resid)?),
            Robj::from(tr_mat),
        ],
    )
}

#[cfg(feature = "exports")]
#[extendr]
pub fn rust_run_spur_i0(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i0_r(coords, latlong, y, emat, nthreads)
}

#[cfg(feature = "exports")]
#[extendr]
pub fn rust_run_spur_i1(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i1_r(coords, latlong, y, emat, nthreads)
}

#[cfg(feature = "exports")]
#[extendr]
pub fn rust_run_spur_i0_resid(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i0_resid_r(coords, latlong, y, x_in, emat, nthreads)
}

#[cfg(feature = "exports")]
#[extendr]
pub fn rust_run_spur_i1_resid(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i1_resid_r(coords, latlong, y, x_in, emat, nthreads)
}

#[cfg(feature = "exports")]
#[extendr]
pub fn rust_run_spur_halflife(
    coords: RMatrix<f64>,
    latlong: bool,
    y: Doubles,
    emat: RMatrix<f64>,
    level: f64,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_halflife_r(coords, latlong, y, emat, level, nthreads)
}

#[cfg(feature = "exports")]
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_spur_transform(
    kind: &str,
    n: i32,
    y: Nullable<RMatrix<f64>>,
    coords: Nullable<RMatrix<f64>>,
    latlong: bool,
    radius: Nullable<f64>,
    strict_lt: bool,
    exclude_zero_distances: bool,
    cluster_id: Nullable<Integers>,
    return_matrix: bool,
    demean_lbmgls: bool,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_transform_r(
        kind,
        n,
        y,
        coords,
        latlong,
        radius,
        strict_lt,
        exclude_zero_distances,
        cluster_id,
        return_matrix,
        demean_lbmgls,
        nthreads,
    )
}

#[cfg(feature = "exports")]
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_spur_pipeline(
    coords: RMatrix<f64>,
    latlong: bool,
    y: Doubles,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    transform_vars: RMatrix<f64>,
    demean_lbmgls: bool,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_pipeline_r(
        coords,
        latlong,
        y,
        x_in,
        emat,
        transform_vars,
        demean_lbmgls,
        nthreads,
    )
}

/// Return compiled `spatial-core` version, active linear algebra provider, and ABI versions.
pub fn build_info_r() -> List {
    list!(
        version = spatial_core::VERSION,
        spec_version = spatial_core::SPEC_VERSION,
        os = std::env::consts::OS,
        arch = std::env::consts::ARCH,
        linalg_provider = spatial_core::linalg_provider(),
        scpc_api_version = spatial_core::SCPC_API_VERSION as i32,
        spur_api_version = spatial_core::SPUR_API_VERSION as i32
    )
}

#[cfg(feature = "exports")]
#[extendr]
pub fn rust_build_info() -> List {
    build_info_r()
}

#[cfg(feature = "exports")]
extendr_module! {
    mod spatial_r;
    fn rust_run_model_scpc;
    fn rust_build_info;
    fn rust_run_exact_scpc;
    fn rust_run_exact_scpc_clustered;
    fn rust_run_exact_scpc_full;
    fn rust_run_approx_scpc;
    fn rust_run_spur_i0;
    fn rust_run_spur_i1;
    fn rust_run_spur_i0_resid;
    fn rust_run_spur_i1_resid;
    fn rust_run_spur_halflife;
    fn rust_run_spur_transform;
    fn rust_run_spur_pipeline;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn col_major_to_row_major_transposes_correctly() {
        // R matrix 2 x 3:
        // [1 3 5]
        // [2 4 6]
        let col_major = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let row_major = col_major_to_row_major(&col_major, 2, 3).unwrap();
        assert_eq!(row_major, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);
    }

    #[test]
    fn col_major_to_row_major_rejects_mismatched_length() {
        assert!(col_major_to_row_major(&[1.0, 2.0], 2, 2).is_err());
    }

    #[test]
    fn core_err_maps_all_error_kinds() {
        let cases = [
            (CoreError::invalid_input("bad arg"), "bad arg"),
            (CoreError::unsupported("approx"), "unsupported case: approx"),
            (CoreError::tied("ambiguous"), "ambiguous"),
            (CoreError::invalid_sampled_covariance("not psd"), "not psd"),
            (CoreError::calibration("no root"), "calibration failed: no root"),
            (CoreError::degenerate("zero var"), "degenerate variance: zero var"),
            (CoreError::internal("lapack"), "internal failure: lapack"),
        ];
        for (err, expected) in cases {
            match core_err(err) {
                Error::Other(msg) => assert_eq!(msg, expected),
                other => panic!("unexpected Error variant: {other:?}"),
            }
        }
    }

    #[test]
    fn approx_scpc_error_sets_expected_classes_and_messages() {
        let tied = ApproxScpcError::from(CoreError::tied("tied basis"));
        assert_eq!(tied.class, Some("spatial_tied_eigenspace"));
        assert_eq!(tied.message, "tied basis");

        let invalid_cov =
            ApproxScpcError::from(CoreError::invalid_sampled_covariance("indefinite grid"));
        assert_eq!(invalid_cov.class, Some("spatial_invalid_sampled_covariance"));
        assert_eq!(invalid_cov.message, "indefinite grid");

        let other_core = ApproxScpcError::from(CoreError::degenerate("zero var"));
        assert_eq!(other_core.class, None);
        assert_eq!(other_core.message, "degenerate variance: zero var");

        let extendr_err = ApproxScpcError::from(Error::Other("direct msg".to_string()));
        assert_eq!(extendr_err.class, None);
        assert_eq!(extendr_err.message, "direct msg");
    }
}
