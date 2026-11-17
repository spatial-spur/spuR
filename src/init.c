#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifdef SPATIAL_USE_RUST
extern SEXP wrap__rust_run_spur_i0(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP emat,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_i1(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP emat,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_i0_resid(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP x_in,
    SEXP emat,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_i1_resid(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP x_in,
    SEXP emat,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_halflife(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP emat,
    SEXP level,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_transform(
    SEXP kind,
    SEXP n,
    SEXP y,
    SEXP coords,
    SEXP latlong,
    SEXP radius,
    SEXP strict_lt,
    SEXP exclude_zero_distances,
    SEXP cluster_id,
    SEXP return_matrix,
    SEXP demean_lbmgls,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_pipeline(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP x_in,
    SEXP emat,
    SEXP transform_vars,
    SEXP demean_lbmgls,
    SEXP nthreads
);

static const R_CallMethodDef CallEntries[] = {
    {"wrap__rust_run_spur_i0", (DL_FUNC) &wrap__rust_run_spur_i0, 5},
    {"wrap__rust_run_spur_i1", (DL_FUNC) &wrap__rust_run_spur_i1, 5},
    {"wrap__rust_run_spur_i0_resid", (DL_FUNC) &wrap__rust_run_spur_i0_resid, 6},
    {"wrap__rust_run_spur_i1_resid", (DL_FUNC) &wrap__rust_run_spur_i1_resid, 6},
    {"wrap__rust_run_spur_halflife", (DL_FUNC) &wrap__rust_run_spur_halflife, 6},
    {"wrap__rust_run_spur_transform", (DL_FUNC) &wrap__rust_run_spur_transform, 12},
    {"wrap__rust_run_spur_pipeline", (DL_FUNC) &wrap__rust_run_spur_pipeline, 8},
    {NULL, NULL, 0}
};
#else
static const R_CallMethodDef CallEntries[] = {
    {NULL, NULL, 0}
};
#endif

void R_init_spuR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
