use extendr_api::prelude::*;
use spatial_r::{
    run_spur_halflife_r, run_spur_i0_r, run_spur_i0_resid_r, run_spur_i1_r, run_spur_i1_resid_r,
    run_spur_pipeline_r, run_spur_transform_r,
};

#[extendr]
fn rust_run_spur_i0(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i0_r(coords, latlong, y, emat, nthreads)
}

#[extendr]
fn rust_run_spur_i1(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i1_r(coords, latlong, y, emat, nthreads)
}

#[extendr]
fn rust_run_spur_i0_resid(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i0_resid_r(coords, latlong, y, x_in, emat, nthreads)
}

#[extendr]
fn rust_run_spur_i1_resid(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i1_resid_r(coords, latlong, y, x_in, emat, nthreads)
}

#[extendr]
fn rust_run_spur_halflife(
    coords: RMatrix<f64>,
    latlong: bool,
    y: Doubles,
    emat: RMatrix<f64>,
    level: f64,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_halflife_r(coords, latlong, y, emat, level, nthreads)
}

#[extendr]
#[allow(clippy::too_many_arguments)]
fn rust_run_spur_transform(
    kind: &str,
    n: i32,
    y: RMatrix<f64>,
    coords: Nullable<RMatrix<f64>>,
    latlong: bool,
    radius: Nullable<f64>,
    cluster_id: Nullable<Integers>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_transform_r(kind, n, y, coords, latlong, radius, cluster_id, nthreads)
}

#[extendr]
fn rust_run_spur_pipeline(
    coords: RMatrix<f64>,
    latlong: bool,
    y: Doubles,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    transform_vars: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_pipeline_r(coords, latlong, y, x_in, emat, transform_vars, nthreads)
}

#[extendr]
fn rust_build_info() -> List {
    spatial_r::build_info_r()
}
