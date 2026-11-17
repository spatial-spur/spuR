#' @keywords internal
.resolve_coords_input_scpc <- function(data,
                                       obs_index,
                                       lon = NULL,
                                       lat = NULL,
                                       coords_euclidean = NULL) {
  use_geodesic <- !is.null(lon) || !is.null(lat)
  use_euclidean <- !is.null(coords_euclidean)

  if (use_geodesic && use_euclidean) {
    stop("Specify either `lon`/`lat` or `coords_euclidean`, not both.")
  }
  if (!use_geodesic && !use_euclidean) {
    stop("Specify coordinates via `lon`/`lat` or `coords_euclidean`.")
  }

  if (use_geodesic) {
    if (is.null(lon) || is.null(lat)) {
      stop("For geodesic coordinates, provide both `lon` and `lat`.")
    }
    if (!is.character(lon) || length(lon) != 1L || !nzchar(lon) ||
      !is.character(lat) || length(lat) != 1L || !nzchar(lat)) {
      stop("`lon` and `lat` must each be a single column name.")
    }

    miss <- setdiff(c(lon, lat), names(data))
    if (length(miss) > 0L) {
      stop("Coordinate variables not found in data: ", paste(miss, collapse = ", "))
    }

    coords_lon_lat <- data[obs_index, c(lon, lat), drop = FALSE]
    if (!all(vapply(coords_lon_lat, is.numeric, logical(1)))) {
      stop("`lon` and `lat` must reference numeric columns.")
    }
    if (any(!is.finite(as.matrix(coords_lon_lat)))) {
      stop("Geodesic coordinates must be finite.")
    }
    if (any(coords_lon_lat[[lon]] < -180 | coords_lon_lat[[lon]] > 180)) {
      stop("Longitude values must be in [-180, 180].")
    }
    if (any(coords_lon_lat[[lat]] < -90 | coords_lon_lat[[lat]] > 90)) {
      stop("Latitude values must be in [-90, 90].")
    }

    # Internal SPUR distance code expects columns ordered as (lat, lon).
    coords <- as.matrix(coords_lon_lat[, c(lat, lon), drop = FALSE])
    return(list(coords = coords, latlong = TRUE))
  }

  if (!is.character(coords_euclidean) || length(coords_euclidean) < 1L) {
    stop("`coords_euclidean` must be a character vector with at least one column name.")
  }

  miss <- setdiff(coords_euclidean, names(data))
  if (length(miss) > 0L) {
    stop("Coordinate variables not found in data: ", paste(miss, collapse = ", "))
  }

  coords <- data[obs_index, coords_euclidean, drop = FALSE]
  if (!all(vapply(coords, is.numeric, logical(1)))) {
    stop("`coords_euclidean` columns must be numeric.")
  }
  if (any(!is.finite(as.matrix(coords)))) {
    stop("Euclidean coordinates must be finite.")
  }

  list(coords = as.matrix(coords), latlong = FALSE)
}

#' @keywords internal
.resolve_spur_coords <- function(data,
                                 use_rows,
                                 lon = NULL,
                                 lat = NULL,
                                 coords_euclidean = NULL,
                                 verbose = FALSE) {
  coord_info <- .resolve_coords_input_scpc(
    data = data,
    obs_index = which(use_rows),
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean
  )

  coords <- get_s_matrix(
    coords = coord_info$coords,
    latlong = coord_info$latlong,
    verbose = verbose
  )

  list(coords = coords, latlong = coord_info$latlong)
}

#' @keywords internal
.validate_common_args <- function(data, q, nrep) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  if (length(q) != 1L || !is.finite(q) || q < 1 || q != as.integer(q)) {
    stop("`q` must be a positive integer.")
  }
  if (q < 2L) {
    stop("`q` must be at least 2 for informative SPUR inference.")
  }
  if (length(nrep) != 1L || !is.finite(nrep) || nrep < 1 || nrep != as.integer(nrep)) {
    stop("`nrep` must be a positive integer.")
  }
  invisible(NULL)
}

#' @keywords internal
.parse_single_var_formula <- function(formula, data, fn_name) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.")
  }

  if (length(formula) == 2L) {
    rhs_vars <- all.vars(formula[[2L]])
    if (length(rhs_vars) != 1L) {
      stop("`formula` must reference exactly one variable.")
    }
    var <- rhs_vars[[1L]]
  } else if (length(formula) == 3L) {
    lhs_vars <- all.vars(formula[[2L]])
    rhs_vars <- all.vars(formula[[3L]])
    if (length(lhs_vars) != 1L || length(rhs_vars) != 0L) {
      stop(sprintf(
        "For %s use a single-variable formula like `y ~ 1` or `~ y`.",
        fn_name
      ))
    }
    var <- lhs_vars[[1L]]
  } else {
    stop("`formula` must be one-sided (`~ y`) or two-sided (`y ~ 1`).")
  }

  if (!var %in% names(data)) {
    stop("The variable specified in `formula` does not exist in data.")
  }
  if (!is.numeric(data[[var]])) {
    stop("The variable specified in `formula` must be numeric.")
  }

  list(var = var, use_rows = !is.na(data[[var]]))
}

#' @keywords internal
.parse_residual_formula <- function(formula, data, fn_name) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop(sprintf("`formula` for %s must be two-sided, e.g. `y ~ x1 + x2`.", fn_name))
  }

  mf <- stats::model.frame(formula = formula, data = data, na.action = stats::na.pass)
  use_rows <- stats::complete.cases(mf)
  if (!any(use_rows)) {
    stop("No complete observations after applying `formula`.")
  }

  mf_use <- mf[use_rows, , drop = FALSE]
  y <- stats::model.response(mf_use)
  if (is.null(y)) {
    stop("`formula` must have a dependent variable on the left-hand side.")
  }
  if (is.matrix(y) && ncol(y) != 1L) {
    stop("Dependent variable in `formula` must be univariate.")
  }
  if (!is.numeric(y)) {
    stop("Dependent variable in `formula` must be numeric.")
  }

  terms_obj <- stats::terms(formula, data = data)
  x_in <- stats::model.matrix(stats::delete.response(terms_obj), data = mf_use)
  if (!"(Intercept)" %in% colnames(x_in)) {
    x_in <- cbind("(Intercept)" = 1, x_in)
  }

  list(Y = as.matrix(y), X_in = as.matrix(x_in), use_rows = use_rows)
}

#' @keywords internal
.parse_transform_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.")
  }

  vars <- unique(all.vars(formula))
  if (!length(vars)) {
    stop("`formula` must reference at least one variable to transform.")
  }
  if (!all(vars %in% names(data))) {
    stop("At least one variable in `formula` was not found in data.")
  }
  if (!all(vapply(data[vars], is.numeric, logical(1)))) {
    stop("All variables in `formula` must be numeric.")
  }

  vars
}

#' SPUR I(0) Test
#'
#' Formula interface for the SPUR spatial I(0) test.
#'
#' @param formula Formula identifying a single variable (`y ~ 1` or `~ y`).
#' @param data Data frame containing test and coordinate variables.
#' @param q Number of low-frequency weighted averages (at least 2).
#' @param nrep Number of Monte Carlo draws.
#' @param lon Optional longitude column name (must be supplied with `lat`).
#' @param lat Optional latitude column name (must be supplied with `lon`).
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param seed Optional integer seed for reproducibility of Monte Carlo draws.
#' @param verbose Logical; if `TRUE`, print test summary and diagnostics.
#' @param backend Backend selector: `"auto"` (default), `"rust"`, or `"fallback"`.
#' @param nthreads Optional positive integer thread count for parallel Rust execution.
#' @return A list of class `"spur_test"` with `pvalue`, `teststat`, `ha_param`,
#'   `cv`, and `full`.
#' @examples
#' \donttest{
#' data(spur_example)
#' res <- spurtest_i0(am ~ 1, data = spur_example,
#'                    lon = "lon", lat = "lat", nrep = 999, seed = 42)
#' res
#' }
#' @export
spurtest_i0 <- function(formula,
                        data,
                        q = 15,
                        nrep = 100000,
                        lon = NULL,
                        lat = NULL,
                        coords_euclidean = NULL,
                        seed = 42L,
                        verbose = FALSE,
                        backend = c("auto", "rust", "fallback"),
                        nthreads = NULL) {
  backend_resolved <- .resolve_spur_backend(backend, missing(backend))
  nthreads_resolved <- .resolve_spur_nthreads(nthreads)
  .validate_common_args(data = data, q = q, nrep = nrep)
  parsed <- .parse_single_var_formula(formula = formula, data = data, fn_name = "spurtest_i0()")
  coord_info <- .resolve_spur_coords(
    data = data,
    use_rows = parsed$use_rows,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    verbose = verbose
  )

  y <- as.matrix(data[[parsed$var]][parsed$use_rows])
  if (q >= nrow(y)) {
    stop("`q` must be less than the number of observations.")
  }
  .validate_spur_y(y, name = parsed$var)
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = q, nrep = nrep)

  rust_loaded <- .require_spur_rust_if_requested(backend_resolved, "wrap__rust_run_spur_i0")
  tied_from_rust <- FALSE
  if (!identical(backend_resolved, "fallback") && rust_loaded) {
    rust_attempt <- .try_spur_test_rust(
      "wrap__rust_run_spur_i0",
      coords = coord_info$coords,
      latlong = coord_info$latlong,
      y_mat = y,
      emat = emat,
      is_i0 = TRUE,
      nthreads_resolved = nthreads_resolved,
      backend_resolved = backend_resolved
    )
    tied_from_rust <- isTRUE(rust_attempt$tied)
    if (!is.null(rust_attempt$res)) {
      return(.format_spur_test_from_full(
        rust_attempt$res,
        test_type = "I(0)",
        title = "Spatial I(0) Test Results",
        stat_label = "LFST",
        cv = rust_attempt$cv,
        verbose = verbose,
        backend = "rust",
        tied = tied_from_rust
      ))
    }
  }

  distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
  res <- spur_i0_test(Y = y, distmat = distmat, emat = emat)
  .format_spur_test_from_full(
    res,
    test_type = "I(0)",
    title = "Spatial I(0) Test Results",
    stat_label = "LFST",
    cv = res$cvalue,
    verbose = verbose,
    backend = "fallback",
    tied = tied_from_rust
  )
}

#' SPUR I(1) Test
#'
#' Formula interface for the SPUR spatial I(1) test.
#'
#' @param formula Formula identifying a single variable (`y ~ 1` or `~ y`).
#' @param data Data frame containing test and coordinate variables.
#' @param q Number of low-frequency weighted averages (at least 2).
#' @param nrep Number of Monte Carlo draws.
#' @param lon Optional longitude column name (must be supplied with `lat`).
#' @param lat Optional latitude column name (must be supplied with `lon`).
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param seed Optional integer seed for reproducibility of Monte Carlo draws.
#' @param verbose Logical; if `TRUE`, print test summary and diagnostics.
#' @param backend Backend selector: `"auto"` (default), `"rust"`, or `"fallback"`.
#' @param nthreads Optional positive integer thread count for parallel Rust execution.
#' @return A list of class `"spur_test"` with `pvalue`, `teststat`, `ha_param`,
#'   `cv`, and `full`.
#' @examples
#' \donttest{
#' data(spur_example)
#' res <- spurtest_i1(am ~ 1, data = spur_example,
#'                    lon = "lon", lat = "lat", nrep = 999, seed = 42)
#' res
#' }
#' @export
spurtest_i1 <- function(formula,
                        data,
                        q = 15,
                        nrep = 100000,
                        lon = NULL,
                        lat = NULL,
                        coords_euclidean = NULL,
                        seed = 42L,
                        verbose = FALSE,
                        backend = c("auto", "rust", "fallback"),
                        nthreads = NULL) {
  backend_resolved <- .resolve_spur_backend(backend, missing(backend))
  nthreads_resolved <- .resolve_spur_nthreads(nthreads)
  .validate_common_args(data = data, q = q, nrep = nrep)
  parsed <- .parse_single_var_formula(formula = formula, data = data, fn_name = "spurtest_i1()")
  coord_info <- .resolve_spur_coords(
    data = data,
    use_rows = parsed$use_rows,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    verbose = verbose
  )

  y <- as.matrix(data[[parsed$var]][parsed$use_rows])
  if (q >= nrow(y)) {
    stop("`q` must be less than the number of observations.")
  }
  .validate_spur_y(y, name = parsed$var)
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = q, nrep = nrep)

  rust_loaded <- .require_spur_rust_if_requested(backend_resolved, "wrap__rust_run_spur_i1")
  tied_from_rust <- FALSE
  if (!identical(backend_resolved, "fallback") && rust_loaded) {
    rust_attempt <- .try_spur_test_rust(
      "wrap__rust_run_spur_i1",
      coords = coord_info$coords,
      latlong = coord_info$latlong,
      y_mat = y,
      emat = emat,
      is_i0 = FALSE,
      nthreads_resolved = nthreads_resolved,
      backend_resolved = backend_resolved
    )
    tied_from_rust <- isTRUE(rust_attempt$tied)
    if (!is.null(rust_attempt$res)) {
      return(.format_spur_test_from_full(
        rust_attempt$res,
        test_type = "I(1)",
        title = "Spatial I(1) Test Results",
        stat_label = "LFUR",
        cv = rust_attempt$cv,
        verbose = verbose,
        backend = "rust",
        tied = tied_from_rust
      ))
    }
  }

  distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
  res <- spur_i1_test(Y = y, distmat = distmat, emat = emat)
  .format_spur_test_from_full(
    res,
    test_type = "I(1)",
    title = "Spatial I(1) Test Results",
    stat_label = "LFUR",
    cv = res$cv_vec,
    verbose = verbose,
    backend = "fallback",
    tied = tied_from_rust
  )
}

#' SPUR I(0) Residual Test
#'
#' Formula interface for the SPUR spatial I(0) residual test.
#'
#' @param formula Two-sided model formula (`y ~ x1 + x2`).
#' @param data Data frame containing model and coordinate variables.
#' @param q Number of low-frequency weighted averages (at least 2).
#' @param nrep Number of Monte Carlo draws.
#' @param lon Optional longitude column name (must be supplied with `lat`).
#' @param lat Optional latitude column name (must be supplied with `lon`).
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param seed Optional integer seed for reproducibility of Monte Carlo draws.
#' @param verbose Logical; if `TRUE`, print test summary and diagnostics.
#' @param backend Backend selector: `"auto"` (default), `"rust"`, or `"fallback"`.
#' @param nthreads Optional positive integer thread count for parallel Rust execution.
#' @return A list of class `"spur_test"` with `pvalue`, `teststat`, `ha_param`,
#'   `cv`, and `full`.
#' @examples
#' \donttest{
#' data(spur_example)
#' res <- spurtest_i0resid(am ~ gini, data = spur_example,
#'                         lon = "lon", lat = "lat", nrep = 999, seed = 42)
#' res
#' }
#' @export
spurtest_i0resid <- function(formula,
                             data,
                             q = 15,
                             nrep = 100000,
                             lon = NULL,
                             lat = NULL,
                             coords_euclidean = NULL,
                             seed = 42L,
                             verbose = FALSE,
                             backend = c("auto", "rust", "fallback"),
                             nthreads = NULL) {
  backend_resolved <- .resolve_spur_backend(backend, missing(backend))
  nthreads_resolved <- .resolve_spur_nthreads(nthreads)
  .validate_common_args(data = data, q = q, nrep = nrep)
  parsed <- .parse_residual_formula(formula = formula, data = data, fn_name = "spurtest_i0resid()")
  coord_info <- .resolve_spur_coords(
    data = data,
    use_rows = parsed$use_rows,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    verbose = verbose
  )

  .validate_spur_y(parsed$Y, name = "Y")
  q_x <- .orthonormal_design_basis(parsed$X_in)
  resid_dof <- nrow(parsed$Y) - ncol(q_x)
  if (q > resid_dof) {
    stop(sprintf(
      "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
      q, resid_dof
    ))
  }
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = q, nrep = nrep)

  rust_loaded <- .require_spur_rust_if_requested(backend_resolved, "wrap__rust_run_spur_i0_resid")
  tied_from_rust <- FALSE
  if (!identical(backend_resolved, "fallback") && rust_loaded) {
    rust_attempt <- .try_spur_test_rust(
      "wrap__rust_run_spur_i0_resid",
      coords = coord_info$coords,
      latlong = coord_info$latlong,
      y_mat = parsed$Y,
      emat = emat,
      x_in = parsed$X_in,
      is_i0 = TRUE,
      nthreads_resolved = nthreads_resolved,
      backend_resolved = backend_resolved
    )
    tied_from_rust <- isTRUE(rust_attempt$tied)
    if (!is.null(rust_attempt$res)) {
      return(.format_spur_test_from_full(
        rust_attempt$res,
        test_type = "I(0) Residual",
        title = "Spatial I(0) Test Results for Residuals",
        stat_label = "LFST",
        cv = rust_attempt$cv,
        verbose = verbose,
        backend = "rust",
        tied = tied_from_rust
      ))
    }
  }

  distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
  res <- spur_i0_resid_test(Y = parsed$Y, X_in = parsed$X_in, distmat = distmat, emat = emat)
  .format_spur_test_from_full(
    res,
    test_type = "I(0) Residual",
    title = "Spatial I(0) Test Results for Residuals",
    stat_label = "LFST",
    cv = res$cvalue,
    verbose = verbose,
    backend = "fallback",
    tied = tied_from_rust
  )
}

#' SPUR I(1) Residual Test
#'
#' Formula interface for the SPUR spatial I(1) residual test.
#'
#' @param formula Two-sided model formula (`y ~ x1 + x2`).
#' @param data Data frame containing model and coordinate variables.
#' @param q Number of low-frequency weighted averages (at least 2).
#' @param nrep Number of Monte Carlo draws.
#' @param lon Optional longitude column name (must be supplied with `lat`).
#' @param lat Optional latitude column name (must be supplied with `lon`).
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param seed Optional integer seed for reproducibility of Monte Carlo draws.
#' @param verbose Logical; if `TRUE`, print test summary and diagnostics.
#' @param backend Backend selector: `"auto"` (default), `"rust"`, or `"fallback"`.
#' @param nthreads Optional positive integer thread count for parallel Rust execution.
#' @return A list of class `"spur_test"` with `pvalue`, `teststat`, `ha_param`,
#'   `cv`, and `full`.
#' @examples
#' \donttest{
#' data(spur_example)
#' res <- spurtest_i1resid(am ~ gini, data = spur_example,
#'                         lon = "lon", lat = "lat", nrep = 999, seed = 42)
#' res
#' }
#' @export
spurtest_i1resid <- function(formula,
                             data,
                             q = 15,
                             nrep = 100000,
                             lon = NULL,
                             lat = NULL,
                             coords_euclidean = NULL,
                             seed = 42L,
                             verbose = FALSE,
                             backend = c("auto", "rust", "fallback"),
                             nthreads = NULL) {
  backend_resolved <- .resolve_spur_backend(backend, missing(backend))
  nthreads_resolved <- .resolve_spur_nthreads(nthreads)
  .validate_common_args(data = data, q = q, nrep = nrep)
  parsed <- .parse_residual_formula(formula = formula, data = data, fn_name = "spurtest_i1resid()")
  coord_info <- .resolve_spur_coords(
    data = data,
    use_rows = parsed$use_rows,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    verbose = verbose
  )

  .validate_spur_y(parsed$Y, name = "Y")
  q_x <- .orthonormal_design_basis(parsed$X_in)
  resid_dof <- nrow(parsed$Y) - ncol(q_x)
  if (q > resid_dof) {
    stop(sprintf(
      "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
      q, resid_dof
    ))
  }
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = q, nrep = nrep)

  rust_loaded <- .require_spur_rust_if_requested(backend_resolved, "wrap__rust_run_spur_i1_resid")
  tied_from_rust <- FALSE
  if (!identical(backend_resolved, "fallback") && rust_loaded) {
    rust_attempt <- .try_spur_test_rust(
      "wrap__rust_run_spur_i1_resid",
      coords = coord_info$coords,
      latlong = coord_info$latlong,
      y_mat = parsed$Y,
      emat = emat,
      x_in = parsed$X_in,
      is_i0 = FALSE,
      nthreads_resolved = nthreads_resolved,
      backend_resolved = backend_resolved
    )
    tied_from_rust <- isTRUE(rust_attempt$tied)
    if (!is.null(rust_attempt$res)) {
      return(.format_spur_test_from_full(
        rust_attempt$res,
        test_type = "I(1) Residual",
        title = "Spatial I(1) Test Results for Residuals",
        stat_label = "LFUR",
        cv = rust_attempt$cv,
        verbose = verbose,
        backend = "rust",
        tied = tied_from_rust
      ))
    }
  }

  distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
  res <- spur_i1_resid_test(Y = parsed$Y, X_in = parsed$X_in, distmat = distmat, emat = emat)
  .format_spur_test_from_full(
    res,
    test_type = "I(1) Residual",
    title = "Spatial I(1) Test Results for Residuals",
    stat_label = "LFUR",
    cv = res$cv_vec,
    verbose = verbose,
    backend = "fallback",
    tied = tied_from_rust
  )
}

#' SPUR Half-Life Confidence Interval
#'
#' Formula interface for SPUR half-life confidence intervals.
#'
#' @param formula Formula identifying a single variable (`y ~ 1` or `~ y`).
#' @param data Data frame containing test and coordinate variables.
#' @param q Number of low-frequency weighted averages (at least 2).
#' @param nrep Number of Monte Carlo draws.
#' @param level Confidence level in percent.
#' @param normdist Logical; return CI as fractions of maximum distance.
#' @param lon Optional longitude column name (must be supplied with `lat`).
#' @param lat Optional latitude column name (must be supplied with `lon`).
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param seed Optional integer seed for reproducibility of Monte Carlo draws.
#' @param verbose Logical; if `TRUE`, print CI summary and diagnostics.
#' @param backend Backend selector: `"auto"` (default), `"rust"`, or `"fallback"`.
#' @param nthreads Optional positive integer thread count for parallel Rust execution.
#' @return A list of class `"spur_halflife"` with `ci_l`, `ci_u`, `max_dist`,
#'   and `full`.
#' @examples
#' \donttest{
#' data(spur_example)
#' hl <- spurhalflife(am ~ 1, data = spur_example,
#'                    lon = "lon", lat = "lat", nrep = 999, seed = 42)
#' hl
#' }
#' @export
spurhalflife <- function(formula,
                         data,
                         q = 15,
                         nrep = 100000,
                         level = 95,
                         normdist = FALSE,
                         lon = NULL,
                         lat = NULL,
                         coords_euclidean = NULL,
                         seed = 42L,
                         verbose = FALSE,
                         backend = c("auto", "rust", "fallback"),
                         nthreads = NULL) {
  backend_resolved <- .resolve_spur_backend(backend, missing(backend))
  nthreads_resolved <- .resolve_spur_nthreads(nthreads)
  .validate_common_args(data = data, q = q, nrep = nrep)
  if (q < 2L) {
    stop("`q` must be at least 2 for half-life inference.")
  }
  if (length(level) != 1L || !is.finite(level)) {
    stop("`level` must be a finite scalar.")
  }
  if (level >= 100 || level <= 0) {
    stop("Invalid level.")
  }

  parsed <- .parse_single_var_formula(formula = formula, data = data, fn_name = "spurhalflife()")
  coord_info <- .resolve_spur_coords(
    data = data,
    use_rows = parsed$use_rows,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    verbose = verbose
  )

  z <- as.matrix(data[[parsed$var]][parsed$use_rows])
  if (q >= nrow(z)) {
    stop("`q` must be less than the number of observations.")
  }
  .validate_spur_y(z, name = parsed$var)
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = q, nrep = nrep)

  rust_loaded <- .require_spur_rust_if_requested(backend_resolved, "wrap__rust_run_spur_halflife")
  tied_from_rust <- FALSE
  used_rust <- FALSE
  ci <- NULL
  max_dist <- NA_real_

  if (!identical(backend_resolved, "fallback") && rust_loaded) {
    e_m <- as.matrix(emat)
    storage.mode(e_m) <- "double"
    rust_res <- .call_spur_rust(
      "wrap__rust_run_spur_halflife",
      .swap_coords_for_rust(coord_info$coords, coord_info$latlong),
      isTRUE(coord_info$latlong),
      as.double(z),
      e_m,
      as.double(level / 100),
      nthreads_resolved
    )
    tied_from_rust <- isTRUE(rust_res$tied)
    if (!tied_from_rust || identical(backend_resolved, "rust")) {
      used_rust <- TRUE
      ci <- list(
        hl_ci_lower = as.numeric(rust_res$hl_ci_lower),
        hl_ci_upper = as.numeric(rust_res$hl_ci_upper),
        pv_mat = as.numeric(rust_res$pv_mat),
        hl_grid_ho = as.numeric(rust_res$hl_grid_ho),
        backend = "rust",
        tied = tied_from_rust
      )
      max_dist <- as.numeric(rust_res$max_dist)
      if (isTRUE(coord_info$latlong)) {
        max_dist <- max_dist * .SPUR_HALFLIFE_SCALE_M
      }
    }
  }

  if (!used_rust) {
    distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
    ci <- spur_persistence(Z = z, distmat = distmat, emat = emat, level = level / 100)
    ci$backend <- "fallback"
    ci$tied <- tied_from_rust
    max_dist <- if (isTRUE(coord_info$latlong)) {
      max(get_distmat_lat_lon(coord_info$coords)) * .SPUR_HALFLIFE_SCALE_M
    } else {
      max(get_distmat_euclidean(coord_info$coords))
    }
  }

  ci_l <- ci$hl_ci_lower
  ci_u <- ci$hl_ci_upper
  if (!is.na(ci_u) && ci_u >= 100) {
    ci_u <- NA_real_
  }

  if (!isTRUE(normdist)) {
    ci_l <- ci_l * max_dist
    ci_u <- ci_u * max_dist
  }

  units <- if (isTRUE(normdist)) {
    paste("fractions of maximum distance", format(max_dist, digits = 4))
  } else if (isTRUE(coord_info$latlong)) {
    "metres"
  } else {
    "unit of original coordinates"
  }

  if (isTRUE(verbose)) {
    cat(sprintf("Spatial half life %s%% confidence interval (in %s)\n", level, units))
    cat("---------------------------------------\n")
    cat(sprintf("Lower bound: %9.4f\n", ci_l))
    if (is.na(ci_u)) {
      cat("Upper bound:    inf\n")
    } else {
      cat(sprintf("Upper bound: %9.4f\n", ci_u))
    }
    cat("---------------------------------------\n")
  }

  structure(
    list(
      ci_l = ci_l,
      ci_u = ci_u,
      max_dist = max_dist,
      full = ci,
      backend = if (used_rust) "rust" else "fallback",
      tied = tied_from_rust
    ),
    class = "spur_halflife",
    level = level
  )
}

#' Transform Variables Using SPUR Spatial Transformations
#'
#' Formula interface for `spurtransform`.
#'
#' @param formula Formula listing variables to transform. All variables
#'   in the formula are transformed. Both one-sided (`~ y + x1 + x2`)
#'   and two-sided (`y ~ x1 + x2`) formulas work, so you can pass your
#'   analysis formula directly.
#' @param data A data frame containing variables and coordinates.
#' @param prefix Prefix for transformed variable names.
#' @param transformation One of `"lbmgls"`, `"nn"`, `"iso"`, or `"cluster"`.
#' @param radius Radius used by `transformation = "iso"`.
#' @param clustvar Cluster variable used by `transformation = "cluster"`.
#' @param lon Optional longitude column name (must be supplied with `lat`).
#' @param lat Optional latitude column name (must be supplied with `lon`).
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param separately Logical; handle missingness separately by variable.
#' @param replace Logical; allow overwriting existing output columns.
#' @param verbose Logical; if `TRUE`, print coordinate diagnostics.
#' @param backend Backend selector: `"auto"` (default), `"rust"`, or `"fallback"`.
#' @param nthreads Optional positive integer thread count for parallel Rust execution.
#' @return A data frame with transformed variables appended.
#' @examples
#' data(spur_example)
#'
#' # One-sided formula: transform listed variables
#' out <- spurtransform(~ am + rm, data = spur_example,
#'                      lon = "lon", lat = "lat", transformation = "lbmgls")
#' head(out[, c("am", "h_am", "rm", "h_rm")])
#'
#' # Two-sided formula: pass your analysis formula to transform all variables
#' out2 <- spurtransform(am ~ gini + fracblack, data = spur_example,
#'                       lon = "lon", lat = "lat", transformation = "lbmgls")
#' head(out2[, c("h_am", "h_gini", "h_fracblack")])
#' @export
spurtransform <- function(formula,
                          data,
                          prefix = "h_",
                          transformation = "lbmgls",
                          radius = NULL,
                          clustvar = NULL,
                          lon = NULL,
                          lat = NULL,
                          coords_euclidean = NULL,
                          separately = FALSE,
                          replace = FALSE,
                          verbose = FALSE,
                          backend = c("auto", "rust", "fallback"),
                          nthreads = NULL) {
  backend_resolved <- .resolve_spur_backend(backend, missing(backend))
  nthreads_resolved <- .resolve_spur_nthreads(nthreads)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  vars <- .parse_transform_formula(formula = formula, data = data)

  allowed <- c("lbmgls", "nn", "iso", "cluster")
  if (!transformation %in% allowed) {
    stop("Invalid transformation.")
  }

  if (transformation != "iso" && !is.null(radius)) {
    stop("Option radius only allowed with transformation(iso).")
  }
  if (transformation == "iso" && is.null(radius)) {
    stop("Radius missing.")
  }
  if (transformation == "iso" && radius <= 0) {
    stop("Radius must be positive.")
  }

  if (transformation != "cluster" && !is.null(clustvar)) {
    stop("Option clustvar only allowed with transformation(cluster).")
  }
  if (transformation == "cluster" && is.null(clustvar)) {
    stop("Clustvar missing.")
  }
  if (!is.null(clustvar) && !clustvar %in% names(data)) {
    stop("clustvar was not found in data.")
  }

  if (!is.character(prefix) || length(prefix) != 1L || !nzchar(prefix) || make.names(prefix) != prefix) {
    stop("prefix must be a valid name prefix.")
  }

  rust_loaded <- .require_spur_rust_if_requested(backend_resolved, "wrap__rust_run_spur_transform")
  use_rust_path <- !identical(backend_resolved, "fallback") && rust_loaded

  base_sample <- if (isTRUE(separately)) {
    rep(TRUE, nrow(data))
  } else {
    stats::complete.cases(data[vars])
  }

  cluster_group <- if (transformation == "cluster") {
    as.numeric(factor(data[[clustvar]]))
  } else {
    NULL
  }

  result <- data

  if (!isTRUE(separately)) {
    use_rows <- base_sample
    if (transformation == "cluster") {
      use_rows <- use_rows & !is.na(data[[clustvar]])
    }

    if (use_rust_path && any(use_rows) && (sum(use_rows) >= 2L || identical(backend_resolved, "rust"))) {
      for (var in vars) {
        out_var <- paste0(prefix, var)
        if (out_var %in% names(result) && !isTRUE(replace)) {
          stop(sprintf("%s already defined or invalid name", out_var))
        }
      }
      tr_mat <- .transform_vars_rust(
        data = data,
        vars = vars,
        use_rows = use_rows,
        transformation = transformation,
        radius = radius,
        cluster_group = cluster_group,
        lon = lon,
        lat = lat,
        coords_euclidean = coords_euclidean,
        verbose = verbose,
        nthreads_resolved = nthreads_resolved
      )
      for (j in seq_along(vars)) {
        out_var <- paste0(prefix, vars[[j]])
        result_vec <- rep(NA_real_, nrow(data))
        result_vec[use_rows] <- as.numeric(tr_mat[, j])
        result[[out_var]] <- result_vec
      }
      return(result)
    }

    shared <- NULL
    if (any(use_rows)) {
      coord_info <- .resolve_spur_coords(
        data = data,
        use_rows = use_rows,
        lon = lon,
        lat = lat,
        coords_euclidean = coords_euclidean,
        verbose = verbose
      )

      cluster_vals <- if (transformation == "cluster") {
        get_cluster_matrix(cluster_values = cluster_group, use_rows = use_rows)
      } else {
        NULL
      }

      shared <- list(
        use_rows = use_rows,
        H = make_transform(
          s = coord_info$coords,
          transformation = transformation,
          radius = radius,
          cluster = cluster_vals,
          latlong = coord_info$latlong
        )
      )
    } else {
      shared <- list(use_rows = use_rows, H = NULL)
    }

    for (var in vars) {
      out_var <- paste0(prefix, var)
      if (out_var %in% names(result) && !isTRUE(replace)) {
        stop(sprintf("%s already defined or invalid name", out_var))
      }
      if (!any(shared$use_rows)) {
        warning(sprintf("No valid observations for variable: %s", var))
        next
      }
      transformed <- spur_transform_vector(
        data = data,
        varname = var,
        H = shared$H,
        transformation = transformation,
        use_rows = shared$use_rows
      )
      result_vec <- rep(NA_real_, nrow(data))
      result_vec[shared$use_rows] <- transformed
      result[[out_var]] <- result_vec
    }
    return(result)
  }

  for (var in vars) {
    out_var <- paste0(prefix, var)
    if (out_var %in% names(result) && !isTRUE(replace)) {
      stop(sprintf("%s already defined or invalid name", out_var))
    }

    use_rows <- base_sample & !is.na(data[[var]])
    if (transformation == "cluster") {
      use_rows <- use_rows & !is.na(data[[clustvar]])
    }

    if (!any(use_rows)) {
      warning(sprintf("No valid observations for variable: %s", var))
      next
    }

    if (use_rust_path && (sum(use_rows) >= 2L || identical(backend_resolved, "rust"))) {
      tr_mat <- .transform_vars_rust(
        data = data,
        vars = var,
        use_rows = use_rows,
        transformation = transformation,
        radius = radius,
        cluster_group = cluster_group,
        lon = lon,
        lat = lat,
        coords_euclidean = coords_euclidean,
        verbose = verbose,
        nthreads_resolved = nthreads_resolved
      )
      result_vec <- rep(NA_real_, nrow(data))
      result_vec[use_rows] <- as.numeric(tr_mat[, 1L])
      result[[out_var]] <- result_vec
      next
    }

    coord_info <- .resolve_spur_coords(
      data = data,
      use_rows = use_rows,
      lon = lon,
      lat = lat,
      coords_euclidean = coords_euclidean,
      verbose = verbose
    )

    cluster_vals <- if (transformation == "cluster") {
      get_cluster_matrix(cluster_values = cluster_group, use_rows = use_rows)
    } else {
      NULL
    }

    H <- make_transform(
      s = coord_info$coords,
      transformation = transformation,
      radius = radius,
      cluster = cluster_vals,
      latlong = coord_info$latlong
    )

    transformed <- spur_transform_vector(
      data = data,
      varname = var,
      H = H,
      transformation = transformation,
      use_rows = use_rows
    )

    result_vec <- rep(NA_real_, nrow(data))
    result_vec[use_rows] <- transformed
    result[[out_var]] <- result_vec
  }

  result
}
