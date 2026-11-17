# ---------------------------------------------------------------------------
# Backend, thread configuration, and native Rust dispatch helpers for spuR
# ---------------------------------------------------------------------------

.SPUR_EARTH_CIRCUMFERENCE_M <- 3.14159265359 * 6371000.009 * 2
.SPUR_HALFLIFE_SCALE_M <- pi * 6371000.009 * 2

.validate_spur_backend_value <- function(val, source = "`backend`") {
  valid <- c("auto", "rust", "fallback")
  if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(val) || !val %in% valid) {
    stop(source, " must be one of \"auto\", \"rust\", or \"fallback\".", call. = FALSE)
  }
  val
}

.resolve_spur_backend <- function(backend, backend_missing) {
  if (!backend_missing && !is.null(backend) && !identical(backend, c("auto", "rust", "fallback"))) {
    return(.validate_spur_backend_value(backend, "`backend`"))
  }

  opt <- getOption("spatial.backend")
  if (!is.null(opt)) {
    return(.validate_spur_backend_value(opt, "Option `spatial.backend`"))
  }

  env <- Sys.getenv("SPATIAL_BACKEND", unset = "")
  if (nzchar(env)) {
    return(.validate_spur_backend_value(env, "Environment variable `SPATIAL_BACKEND`"))
  }

  "auto"
}

.parse_spur_positive_int_threads <- function(val, source = "`nthreads`") {
  if (is.character(val)) {
    if (length(val) != 1L || is.na(val) || !nzchar(trimws(val)) || !grepl("^[0-9]+$", trimws(val))) {
      stop(source, " must be NULL or a single positive integer.", call. = FALSE)
    }
    val <- suppressWarnings(as.numeric(trimws(val)))
  }
  if (
    !is.numeric(val) ||
      length(val) != 1L ||
      is.na(val) ||
      !is.finite(val) ||
      val < 1 ||
      abs(val - round(val)) > 0 ||
      val > .Machine$integer.max
  ) {
    stop(source, " must be NULL or a single positive integer.", call. = FALSE)
  }
  as.integer(round(val))
}

.is_spur_r_check_limit_cores <- function() {
  val <- trimws(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))
  nzchar(val) && !tolower(val) %in% c("false", "0", "no")
}

.resolve_spur_nthreads <- function(nthreads = NULL) {
  resolved <- NULL
  if (!is.null(nthreads)) {
    resolved <- .parse_spur_positive_int_threads(nthreads, "`nthreads`")
  } else {
    opt_threads <- getOption("spuR.threads")
    if (!is.null(opt_threads)) {
      resolved <- .parse_spur_positive_int_threads(opt_threads, "Option `spuR.threads`")
    } else {
      env_threads <- Sys.getenv("SPUR_NUM_THREADS", unset = "")
      if (nzchar(trimws(env_threads))) {
        resolved <- .parse_spur_positive_int_threads(env_threads, "Environment variable `SPUR_NUM_THREADS`")
      } else {
        opt_mc <- getOption("mc.cores")
        if (!is.null(opt_mc)) {
          resolved <- .parse_spur_positive_int_threads(opt_mc, "Option `mc.cores`")
        }
      }
    }
  }

  if (.is_spur_r_check_limit_cores()) {
    if (is.null(resolved) || resolved > 2L) {
      resolved <- 2L
    }
  }

  resolved
}

.is_spur_rust_loaded <- function(sym = "wrap__rust_run_spur_i0") {
  is.loaded(sym, PACKAGE = "spuR")
}

.call_spur_rust <- function(sym, ...) {
  res <- .Call(sym, ..., PACKAGE = "spuR")
  if (inherits(res, "extendr_error")) {
    stop(res$value, call. = FALSE)
  }
  res
}

.swap_coords_for_rust <- function(coords, latlong) {
  coords_mat <- as.matrix(coords)
  if (isTRUE(latlong)) {
    coords_mat <- coords_mat[, c(2L, 1L), drop = FALSE]
  }
  storage.mode(coords_mat) <- "double"
  coords_mat
}

.require_spur_rust_if_requested <- function(backend_resolved, sym = "wrap__rust_run_spur_i0") {
  loaded <- .is_spur_rust_loaded(sym)
  if (identical(backend_resolved, "rust") && !loaded) {
    stop("`backend = \"rust\"` was requested, but the native Rust backend is not loaded in `spuR`.", call. = FALSE)
  }
  loaded
}

.try_spur_test_rust <- function(sym,
                                coords,
                                latlong,
                                y_mat,
                                emat,
                                x_in = NULL,
                                is_i0 = FALSE,
                                nthreads_resolved = NULL,
                                backend_resolved = "auto") {
  coords_rust <- .swap_coords_for_rust(coords, latlong)
  y_m <- as.matrix(y_mat)
  storage.mode(y_m) <- "double"
  e_m <- as.matrix(emat)
  storage.mode(e_m) <- "double"

  res <- if (is.null(x_in)) {
    .call_spur_rust(sym, coords_rust, isTRUE(latlong), y_m, e_m, nthreads_resolved)
  } else {
    x_m <- as.matrix(x_in)
    storage.mode(x_m) <- "double"
    .call_spur_rust(sym, coords_rust, isTRUE(latlong), y_m, x_m, e_m, nthreads_resolved)
  }

  tied <- isTRUE(res$tied)
  if (tied && identical(backend_resolved, "auto")) {
    return(list(res = NULL, tied = TRUE))
  }

  if (isTRUE(is_i0)) {
    full <- list(
      LR = as.numeric(res$LR),
      pvalue = as.numeric(res$pvalue),
      cvalue = as.numeric(res$cvalue),
      ha_parm = as.numeric(res$ha_parm),
      cvalue_mat = res$cvalue_mat,
      pvalue_mat = res$pvalue_mat,
      rho_grid = as.numeric(res$rho_grid),
      backend = "rust",
      tied = tied
    )
    list(res = full, tied = tied, cv = as.numeric(res$cvalue))
  } else {
    full <- list(
      LR = as.numeric(res$LR),
      pvalue = as.numeric(res$pvalue),
      cv_vec = as.numeric(res$cv_vec),
      ha_parm = as.numeric(res$ha_parm),
      backend = "rust",
      tied = tied
    )
    list(res = full, tied = tied, cv = as.numeric(res$cv_vec))
  }
}

.format_spur_test_from_full <- function(res,
                                        test_type,
                                        title,
                                        stat_label,
                                        cv,
                                        verbose,
                                        backend,
                                        tied) {
  res$backend <- backend
  res$tied <- tied
  teststat <- res$LR[1]
  pvalue <- res$pvalue[1]
  ha_param <- res$ha_parm

  .format_spur_test_output(
    title = title,
    stat_label = stat_label,
    teststat = teststat,
    pvalue = pvalue,
    verbose = verbose
  )

  structure(
    list(
      pvalue = pvalue,
      teststat = teststat,
      ha_param = ha_param,
      cv = cv,
      full = res,
      backend = backend,
      tied = tied
    ),
    class = "spur_test",
    test_type = test_type
  )
}

.transform_vars_rust <- function(data,
                                 vars,
                                 use_rows,
                                 transformation,
                                 radius,
                                 cluster_group,
                                 lon,
                                 lat,
                                 coords_euclidean,
                                 verbose,
                                 nthreads_resolved) {
  sub_n <- sum(use_rows)
  y_mat <- as.matrix(data[use_rows, vars, drop = FALSE])
  storage.mode(y_mat) <- "double"
  y_names <- colnames(y_mat)

  if (identical(transformation, "cluster")) {
    cl_raw <- factor(cluster_group[use_rows])
    cl_id <- as.integer(cl_raw) - 1L
    storage.mode(cl_id) <- "integer"
    res <- .call_spur_rust(
      "wrap__rust_run_spur_transform",
      "cluster",
      as.integer(sub_n),
      y_mat,
      NULL,
      FALSE,
      NULL,
      FALSE,
      TRUE,
      cl_id,
      FALSE,
      TRUE,
      nthreads_resolved
    )
    out_mat <- res$transformed
    colnames(out_mat) <- y_names
    return(out_mat)
  }

  coord_info <- .resolve_spur_coords(
    data = data,
    use_rows = use_rows,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    verbose = verbose
  )
  coords_rust <- .swap_coords_for_rust(coord_info$coords, coord_info$latlong)
  radius_arg <- if (identical(transformation, "iso")) {
    if (isTRUE(coord_info$latlong)) radius / .SPUR_EARTH_CIRCUMFERENCE_M else as.double(radius)
  } else {
    NULL
  }
  res <- .call_spur_rust(
    "wrap__rust_run_spur_transform",
    transformation,
    as.integer(sub_n),
    y_mat,
    coords_rust,
    isTRUE(coord_info$latlong),
    radius_arg,
    FALSE,
    TRUE,
    NULL,
    FALSE,
    TRUE,
    nthreads_resolved
  )
  out_mat <- res$transformed
  colnames(out_mat) <- y_names
  out_mat
}
