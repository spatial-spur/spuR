#' SPUR I(1) Residual Test
#'
#' R equivalent of Stata `spurtest i1resid`.
#'
#' @param dep_var Dependent variable name.
#' @param indep_vars Optional character vector of independent variable names.
#' @param q Number of low-frequency weighted averages.
#' @param nrep Number of Monte Carlo draws.
#' @param latlong Logical; if `TRUE`, use spherical distances.
#' @param data A data frame containing variables and coordinates.
#' @param coord_cols Optional explicit coordinate column names.
#' @param lat Name of latitude column when `coord_cols` is `NULL`.
#' @param lon Name of longitude column when `coord_cols` is `NULL`.
#' @return A list with `pvalue`, `teststat`, `ha_param`, `cv`, and `full`.
#' @export
spur_i1_resid <- function(dep_var,
                          indep_vars = character(0),
                          q = 15,
                          nrep = 100000,
                          latlong = TRUE,
                          data,
                          coord_cols = NULL,
                          lat = "lat",
                          lon = "lon") {
  vars <- c(dep_var, indep_vars)
  if (!all(vars %in% names(data))) {
    stop("At least one variable was not found in data.")
  }
  if (!all(vapply(data[vars], is.numeric, logical(1)))) {
    stop("dep_var and indep_vars must be numeric.")
  }

  use_rows <- stats::complete.cases(data[vars])
  Y <- as.matrix(data[[dep_var]][use_rows])

  if (length(indep_vars)) {
    X <- as.matrix(data[use_rows, indep_vars, drop = FALSE])
  } else {
    X <- matrix(numeric(0), nrow = sum(use_rows), ncol = 0)
  }
  X_in <- cbind(1, X)

  coord_cols <- .resolve_coord_cols(
    data = data,
    latlong = latlong,
    coord_cols = coord_cols,
    lat = lat,
    lon = lon
  )
  s <- .extract_coords(
    data = data,
    use_rows = use_rows,
    coord_cols = coord_cols,
    latlong = latlong
  )
  distmat <- .normalized_distmat(coords = s, latlong = latlong)

  emat <- .draw_emat(q = q, nrep = nrep)
  res <- spur_i1_resid_test(Y = Y, X_in = X_in, distmat = distmat, emat = emat)

  teststat <- res$LR[1]
  pvalue <- res$pvalue[1]
  ha_param <- res$ha_parm
  cv <- res$cv_vec

  .format_spur_test_output(
    title = "Spatial I(1) Test Results for Residuals",
    stat_label = "LFUR",
    teststat = teststat,
    pvalue = pvalue
  )

  list(
    pvalue = pvalue,
    teststat = teststat,
    ha_param = ha_param,
    cv = cv,
    full = res
  )
}
