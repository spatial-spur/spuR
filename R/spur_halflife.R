#' SPUR Half-Life Confidence Interval
#'
#' R equivalent of Stata `spurhalflife`.
#'
#' @param var Name of numeric variable.
#' @param q Number of low-frequency weighted averages.
#' @param nrep Number of Monte Carlo draws.
#' @param level Confidence level in percent.
#' @param latlong Logical; if `TRUE`, use spherical distances.
#' @param normdist Logical; return CI as fractions of maximum distance.
#' @param data A data frame containing variables and coordinates.
#' @param coord_cols Optional explicit coordinate column names.
#' @param lat Name of latitude column when `coord_cols` is `NULL`.
#' @param lon Name of longitude column when `coord_cols` is `NULL`.
#' @return A list with `ci_l`, `ci_u`, `max_dist`, and `full`.
#' @export
spur_halflife <- function(var,
                          q = 15,
                          nrep = 100000,
                          level = 95,
                          latlong = TRUE,
                          normdist = FALSE,
                          data,
                          coord_cols = NULL,
                          lat = "lat",
                          lon = "lon") {
  if (!var %in% names(data)) {
    stop("The variable specified does not exist in data.")
  }
  if (!is.numeric(data[[var]])) {
    stop("The variable specified must be numeric.")
  }
  if (level >= 100 || level <= 0) {
    stop("Invalid level.")
  }

  use_rows <- !is.na(data[[var]])
  Z <- as.matrix(data[[var]][use_rows])

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
  ci <- spur_persistence(Z = Z, distmat = distmat, emat = emat, level = level / 100)

  ci_l <- ci$hl_ci_lower
  ci_u <- ci$hl_ci_upper
  if (ci_u >= 100) {
    ci_u <- NA_real_
  }

  # Keep exact Stata scaling behavior.
  max_dist <- if (isTRUE(latlong)) {
    max(get_distmat_lat_lon(s))
  } else {
    max(get_distmat_euclidean(s))
  }
  max_dist <- max_dist * 3.14159265359 * 6371000.009 * 2

  if (!isTRUE(normdist)) {
    ci_l <- ci_l * max_dist
    ci_u <- ci_u * max_dist
  }

  units <- if (isTRUE(normdist)) {
    paste("fractions of maximum distance", format(max_dist, digits = 4))
  } else if (isTRUE(latlong)) {
    "metres"
  } else {
    "unit of original coordinates"
  }

  cat(sprintf("Spatial half life %s%% confidence interval (in %s)\n", level, units))
  cat("---------------------------------------\n")
  cat(sprintf("Lower bound: %9.4f\n", ci_l))
  if (is.na(ci_u)) {
    cat("Upper bound:    inf\n")
  } else {
    cat(sprintf("Upper bound: %9.4f\n", ci_u))
  }
  cat("---------------------------------------\n")

  list(
    ci_l = ci_l,
    ci_u = ci_u,
    max_dist = max_dist,
    full = ci
  )
}
