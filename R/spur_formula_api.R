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
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q < 1 || q != trunc(q)) {
    stop("`q` must be a positive integer.")
  }
  if (q < 2L) {
    stop(sprintf("`q` must be at least 2 for informative SPUR inference (requested q=%g).", q))
  }
  if (!is.numeric(nrep) || length(nrep) != 1L || !is.finite(nrep) || nrep < 1 || nrep != trunc(nrep)) {
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
    eval_formula <- stats::as.formula(call("~", formula[[2L]]), env = environment(formula))
  } else if (length(formula) == 3L) {
    if (length(all.vars(formula[[3L]])) != 0L ||
        !identical(deparse(formula[[3L]]), "1")) {
      stop(sprintf(
        "For %s use a single-variable formula like `y ~ 1` or `~ y`.",
        fn_name
      ))
    }
    eval_formula <- stats::as.formula(call("~", formula[[2L]]), env = environment(formula))
  } else {
    stop("`formula` must be one-sided (`~ y`) or two-sided (`y ~ 1`).")
  }

  mf <- stats::model.frame(formula = eval_formula, data = data, na.action = stats::na.pass)
  if (ncol(mf) != 1L) {
    stop("`formula` must reference a single evaluated response.")
  }
  values <- mf[[1L]]
  if (!is.numeric(values) || (!is.null(dim(values)) && ncol(values) != 1L)) {
    stop("The variable specified in `formula` must be numeric.")
  }
  values <- as.numeric(values)

  list(values = values, var = names(mf)[[1L]], use_rows = !is.na(values))
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
  if (ncol(x_in) == 0L) {
    stop("Residual test formula must include at least one term (e.g., an intercept or regressor).")
  }

  list(Y = as.matrix(y), X_in = x_in, use_rows = use_rows)
}

#' @keywords internal
.parse_transform_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.")
  }
  frame <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  if (!is.null(stats::model.offset(frame))) {
    stop("Offsets are not supported by spurtransform().")
  }

  values <- stats::model.matrix(attr(frame, "terms"), data = frame)
  values <- values[, colnames(values) != "(Intercept)", drop = FALSE]
  if (length(formula) == 3L) {
    response <- stats::model.response(frame)
    if (!is.numeric(response) || NCOL(response) != 1L) {
      stop("The formula response must be univariate and numeric.")
    }
    labels <- c(paste(deparse(formula[[2L]]), collapse = ""), colnames(values))
    values <- cbind(response, values)
    colnames(values) <- labels
  }
  if (ncol(values) == 0L) {
    stop("`formula` must reference at least one variable to transform.")
  }
  colnames(values) <- make.names(colnames(values), unique = TRUE)
  values
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
#'   `cv`, and `full`. The `tied` flag indicates a retained tied block and
#'   `all_directions` indicates exhausted positive rank. If tuning misses 50%,
#'   `best_tuning_power` gives the largest estimated power on the independent draws.
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

  y <- as.matrix(parsed$values[parsed$use_rows])
  if (q >= nrow(y)) {
    stop(sprintf("`q` must be less than the number of observations (q=%g, n=%d).", q, nrow(y)))
  }
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = 2L * q, nrep = nrep)

  reason <- .spur_rust_fallback_reason(backend_resolved, "wrap__rust_run_spur_i0")
  if (is.null(reason)) {
    rust_attempt <- .try_spur_test_rust(
      coords = coord_info$coords,
      latlong = coord_info$latlong,
      y_mat = y,
      emat = emat,
      is_i0 = TRUE,
      nthreads_resolved = nthreads_resolved,
      backend_resolved = backend_resolved
    )
    if (!is.null(rust_attempt$res)) {
      return(.format_spur_test_from_full(
        rust_attempt$res,
        test_type = "I(0)",
        title = "Spatial I(0) Test Results",
        stat_label = "LFST",
        cv = rust_attempt$cv,
        verbose = verbose,
        backend = "rust",
        q = q, nrep = nrep
      ))
    }
  }

  distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
  res <- .run_spur_test(y, distmat, emat, is_i0 = TRUE)
  .format_spur_test_from_full(
    res,
    test_type = "I(0)",
    title = "Spatial I(0) Test Results",
    stat_label = "LFST",
    cv = res$cvalue,
    verbose = verbose,
    backend = "fallback",
    q = q, nrep = nrep
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
#'   `cv`, and `full`. The `tied` flag indicates a retained tied block and
#'   `all_directions` indicates exhausted positive rank. If tuning misses 50%,
#'   `best_tuning_power` gives the largest estimated power on the independent draws.
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

  y <- as.matrix(parsed$values[parsed$use_rows])
  if (q >= nrow(y)) {
    stop(sprintf("`q` must be less than the number of observations (q=%g, n=%d).", q, nrow(y)))
  }
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = 2L * q, nrep = nrep)

  reason <- .spur_rust_fallback_reason(backend_resolved, "wrap__rust_run_spur_i1")
  if (is.null(reason)) {
    rust_attempt <- .try_spur_test_rust(
      coords = coord_info$coords,
      latlong = coord_info$latlong,
      y_mat = y,
      emat = emat,
      is_i0 = FALSE,
      nthreads_resolved = nthreads_resolved,
      backend_resolved = backend_resolved
    )
    if (!is.null(rust_attempt$res)) {
      return(.format_spur_test_from_full(
        rust_attempt$res,
        test_type = "I(1)",
        title = "Spatial I(1) Test Results",
        stat_label = "LFUR",
        cv = rust_attempt$cv,
        verbose = verbose,
        backend = "rust",
        q = q, nrep = nrep
      ))
    }
  }

  distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
  res <- .run_spur_test(y, distmat, emat, is_i0 = FALSE)
  .format_spur_test_from_full(
    res,
    test_type = "I(1)",
    title = "Spatial I(1) Test Results",
    stat_label = "LFUR",
    cv = res$cv_vec,
    verbose = verbose,
    backend = "fallback",
    q = q, nrep = nrep
  )
}

#' SPUR I(0) Residual Test
#'
#' Formula interface for the SPUR spatial I(0) residual test.
#' The fitted design must span a constant.
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
#'   `cv`, and `full`. The `tied` flag indicates a retained tied block and
#'   `all_directions` indicates exhausted positive rank. If tuning misses 50%,
#'   `best_tuning_power` gives the largest estimated power on the independent draws.
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

  # Numerical rank and response checks belong to the selected operation.
  resid_dof <- nrow(parsed$Y) - ncol(parsed$X_in)
  if (q > resid_dof) {
    stop(sprintf(
      "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
      q, resid_dof
    ))
  }
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = 2L * q, nrep = nrep)

  reason <- .spur_rust_fallback_reason(backend_resolved, "wrap__rust_run_spur_i0_resid")
  if (is.null(reason)) {
    rust_attempt <- .try_spur_test_rust(
      coords = coord_info$coords,
      latlong = coord_info$latlong,
      y_mat = parsed$Y,
      emat = emat,
      x_in = parsed$X_in,
      is_i0 = TRUE,
      nthreads_resolved = nthreads_resolved,
      backend_resolved = backend_resolved
    )
    if (!is.null(rust_attempt$res)) {
      return(.format_spur_test_from_full(
        rust_attempt$res,
        test_type = "I(0) Residual",
        title = "Spatial I(0) Test Results for Residuals",
        stat_label = "LFST",
        cv = rust_attempt$cv,
        verbose = verbose,
        backend = "rust",
        q = q, nrep = nrep
      ))
    }
  }

  distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
  res <- .run_spur_test(parsed$Y, distmat, emat, is_i0 = TRUE,
                        X_in = parsed$X_in)
  .format_spur_test_from_full(
    res,
    test_type = "I(0) Residual",
    title = "Spatial I(0) Test Results for Residuals",
    stat_label = "LFST",
    cv = res$cvalue,
    verbose = verbose,
    backend = "fallback",
    q = q, nrep = nrep
  )
}

#' SPUR I(1) Residual Test
#'
#' Formula interface for the SPUR spatial I(1) residual test.
#' The fitted design must span a constant.
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
#'   `cv`, and `full`. The `tied` flag indicates a retained tied block and
#'   `all_directions` indicates exhausted positive rank. If tuning misses 50%,
#'   `best_tuning_power` gives the largest estimated power on the independent draws.
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

  # Numerical rank and response checks belong to the selected operation.
  resid_dof <- nrow(parsed$Y) - ncol(parsed$X_in)
  if (q > resid_dof) {
    stop(sprintf(
      "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
      q, resid_dof
    ))
  }
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = 2L * q, nrep = nrep)

  reason <- .spur_rust_fallback_reason(backend_resolved, "wrap__rust_run_spur_i1_resid")
  if (is.null(reason)) {
    rust_attempt <- .try_spur_test_rust(
      coords = coord_info$coords,
      latlong = coord_info$latlong,
      y_mat = parsed$Y,
      emat = emat,
      x_in = parsed$X_in,
      is_i0 = FALSE,
      nthreads_resolved = nthreads_resolved,
      backend_resolved = backend_resolved
    )
    if (!is.null(rust_attempt$res)) {
      return(.format_spur_test_from_full(
        rust_attempt$res,
        test_type = "I(1) Residual",
        title = "Spatial I(1) Test Results for Residuals",
        stat_label = "LFUR",
        cv = rust_attempt$cv,
        verbose = verbose,
        backend = "rust",
        q = q, nrep = nrep
      ))
    }
  }

  distmat <- .normalized_distmat(coords = coord_info$coords, latlong = coord_info$latlong)
  res <- .run_spur_test(parsed$Y, distmat, emat, is_i0 = FALSE,
                        X_in = parsed$X_in)
  .format_spur_test_from_full(
    res,
    test_type = "I(1) Residual",
    title = "Spatial I(1) Test Results for Residuals",
    stat_label = "LFUR",
    cv = res$cv_vec,
    verbose = verbose,
    backend = "fallback",
    q = q, nrep = nrep
  )
}

#' SPUR Half-Life Confidence Set
#'
#' Inverts the continuous half-life tests of Mueller and Watson (2024).
#'
#' @param formula Formula identifying a single variable (`y ~ 1` or `~ y`).
#' @param data Data frame containing test and coordinate variables.
#' @param q Number of low-frequency weighted averages, from 2 through 100 for half-life quadrature.
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
#'   `components`, `full`, and the basis flags `tied` and `all_directions`.
#' @details The calculation evaluates both limits and adds short-half-life
#'   probes based on the smallest positive inter-site distance. It refines
#'   detected changes in acceptance. Inspect `components` for detected components
#'   and `ci_l` and `ci_u` for their hull. Narrow components or gaps can escape
#'   detection. Initial probes and p-values are in
#'   `full$hl_grid_ho` and `full$pv_mat`.
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
  if (q > 100L) stop(sprintf("q=%d exceeds the validated half-life quadrature limit of 100.", q))
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

  z <- as.matrix(parsed$values[parsed$use_rows])
  if (q >= nrow(z)) {
    stop(sprintf("`q` must be less than the number of observations (q=%g, n=%d).", q, nrow(z)))
  }
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = q, nrep = nrep)

  reason <- .spur_rust_fallback_reason(backend_resolved, "wrap__rust_run_spur_halflife")
  used_rust <- FALSE
  ci <- NULL
  max_dist <- NA_real_

  if (is.null(reason)) {
    e_m <- as.matrix(emat)
    storage.mode(e_m) <- "double"
    .mark_spur_accelerate_used()
    rust_res <- .Call(
      "wrap__rust_run_spur_halflife",
      .swap_coords_for_rust(coord_info$coords, coord_info$latlong),
      isTRUE(coord_info$latlong),
      as.double(z),
      e_m,
      as.double(level / 100),
      nthreads_resolved,
      PACKAGE = "spuR"
    )
    rust_res <- .check_spur_rust_result(rust_res)
    used_rust <- TRUE
    ci <- list(
      hl_ci_lower = as.numeric(rust_res$hl_ci_lower),
      hl_ci_upper = as.numeric(rust_res$hl_ci_upper),
      pv_mat = as.numeric(rust_res$pv_mat),
      hl_grid_ho = as.numeric(rust_res$hl_grid_ho),
      components = rust_res$components,
      backend = "rust",
      tied = isTRUE(rust_res$tied), all_directions = isTRUE(rust_res$all_directions)
    )
    max_dist <- as.numeric(rust_res$max_dist)
    if (isTRUE(coord_info$latlong)) {
      max_dist <- max_dist * .SPUR_EARTH_CIRCUMFERENCE_M
    }
  }

  if (!used_rust) {
    distances <- .spur_distance_setup(coord_info$coords, coord_info$latlong)
    ci <- spur_persistence(Z = z, distmat = distances$normalized, emat = emat, level = level / 100)
    ci$backend <- "fallback"
    max_dist <- distances$max_dist
    if (isTRUE(coord_info$latlong)) max_dist <- max_dist * .SPUR_EARTH_CIRCUMFERENCE_M
  }

  .warn_spur_inference(q, nrep, level = level, all_directions = ci$all_directions)
  ci_l <- ci$hl_ci_lower
  ci_u <- ci$hl_ci_upper
  components <- ci$components
  colnames(components) <- c("lower", "upper")
  if (!isTRUE(normdist)) {
    components <- components * max_dist
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
    cat(sprintf("Approximate spatial half-life %s%% confidence set (%s)\n", level, units))
    cat("---------------------------------------\n")
    .print_halflife_components(components)
    cat("---------------------------------------\n")
  }

  structure(
    list(
      ci_l = ci_l,
      ci_u = ci_u,
      max_dist = max_dist,
      components = components,
      full = ci,
      backend = if (used_rust) "rust" else "fallback",
      tied = ci$tied, all_directions = ci$all_directions
    ),
    class = "spur_halflife",
    level = level
  )
}

#' Transform Variables Using SPUR Spatial Transformations
#'
#' Formula interface for `spurtransform`.
#'
#' @param formula One-sided (`~ y + x1 + x2`) or two-sided (`y ~ x1 + x2`)
#'   formula. The evaluated response and design columns are transformed,
#'   including interactions and expressions such as `log(y)`. The
#'   intercept is excluded. Evaluated column names are made syntactic.
#' @param data A data frame containing variables and coordinates.
#' @param prefix Prefix for transformed variable names.
#' @param transformation One of `"lbmgls"`, `"nn"`, `"iso"`, or `"cluster"`.
#' @param radius Radius used by `transformation = "iso"`. Use meters for
#'   longitude and latitude or the units of the Euclidean coordinates.
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
#' @details Nearest-neighbor differencing averages the observations tied at the
#'   smallest positive distance. It requires at least two distinct locations.
#'   Isotropic differencing averages other observations strictly within the radius,
#'   including observations at the same site. An isolated observation transforms to zero.
#'   LBM-GLS centers each transformed column.
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
  values <- .parse_transform_formula(formula = formula, data = data)

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
  if (transformation == "iso" &&
      (!is.numeric(radius) || length(radius) != 1L || !is.finite(radius) || radius <= 0)) {
    stop("Radius must be positive and finite.")
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

  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !nzchar(prefix) || make.names(prefix) != prefix) {
    stop("prefix must be a valid name prefix.")
  }

  use_rust_path <- is.null(.spur_rust_fallback_reason(
    backend_resolved, "wrap__rust_run_spur_transform"
  ))

  cluster_group <- if (transformation == "cluster") {
    factor(data[[clustvar]])
  } else {
    NULL
  }

  output_names <- paste0(prefix, colnames(values))
  clashes <- intersect(output_names, names(data))
  if (!isTRUE(replace) && length(clashes)) {
    stop(sprintf("%s already defined or invalid name", clashes[[1L]]))
  }

  result <- data
  all_columns <- seq_len(ncol(values))
  groups <- if (isTRUE(separately)) as.list(all_columns) else list(all_columns)
  for (columns in groups) {
    y_mat <- values[, columns, drop = FALSE]
    use_rows <- stats::complete.cases(y_mat)
    if (!is.null(cluster_group)) use_rows <- use_rows & !is.na(cluster_group)
    if (!any(use_rows)) {
      for (term in colnames(y_mat)) warning(sprintf("No valid observations for variable: %s", term))
      result[output_names[columns]] <- rep(list(rep(NA_real_, nrow(data))), length(columns))
      next
    }

    y_mat <- y_mat[use_rows, , drop = FALSE]
    if (any(!is.finite(y_mat))) {
      stop("Evaluated transformation variables must be finite on the selected sample.")
    }
    coord_info <- if (is.null(cluster_group)) {
      .resolve_spur_coords(data, use_rows, lon, lat, coords_euclidean, verbose)
    } else {
      list(coords = NULL, latlong = FALSE)
    }
    cluster <- cluster_group[use_rows]
    transformed <- if (use_rust_path &&
        (transformation == "cluster" || nrow(y_mat) >= 2L || identical(backend_resolved, "rust"))) {
      .transform_vars_rust(y_mat, coord_info, transformation, radius, cluster, nthreads_resolved)
    } else if (transformation == "cluster") {
      .demean_spur_clusters(y_mat, cluster)
    } else if (transformation == "lbmgls") {
      distances <- .normalized_distmat(coord_info$coords, coord_info$latlong)
      .apply_lbm_transform(y_mat, eigen(get_sigma_lbm_dm(distances), symmetric = TRUE))
    } else {
      make_transform(coord_info$coords, transformation, radius, cluster, coord_info$latlong) %*% y_mat
    }

    out <- matrix(NA_real_, nrow(data), length(columns))
    out[use_rows, ] <- transformed
    result[output_names[columns]] <- as.data.frame(out)
  }

  result
}

.demean_spur_clusters <- function(values, groups) {
  for (rows in split(seq_len(nrow(values)), groups, drop = TRUE)) {
    block <- values[rows, , drop = FALSE]
    lo <- apply(block, 2L, min)
    hi <- apply(block, 2L, max)
    anchor <- 0.5 * lo + 0.5 * hi
    scale <- pmax(hi - anchor, anchor - lo)
    scale[scale == 0] <- 1
    scaled <- sweep(sweep(block, 2L, anchor, "-"), 2L, scale, "/")
    values[rows, ] <- sweep(sweep(scaled, 2L, colMeans(scaled), "-"), 2L, scale, "*")
  }
  if (any(!is.finite(values))) stop("Cluster demeaning produced a non-finite entry.")
  values
}
