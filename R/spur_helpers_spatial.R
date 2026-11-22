#' Get Spatial Coordinates Matrix
#'
#' Validates coordinate data and reports the distance calculation mode.
#' R translation of the Stata Mata `get_s_matrix` function.
#'
#' @param coords A matrix or data frame containing coordinates.
#' @param latlong Logical; if `TRUE`, coordinates are latitude/longitude.
#' @param verbose Logical; if `TRUE`, print diagnostic messages.
#' @return The validated coordinates matrix.
#' @keywords internal
get_s_matrix <- function(coords, latlong = FALSE, verbose = FALSE) {
  if (any(is.na(coords))) {
    stop("Missing value(s) in coordinate variables; aborting")
  }
  if (any(!is.finite(as.matrix(coords)))) {
    stop("Non-finite value(s) in coordinate variables; aborting")
  }

  ns <- nrow(coords)
  if (ns < 2) {
    stop("Too few locations found; aborting")
  }

  if (isTRUE(verbose)) {
    cat("Found", ns, "observations and", ncol(coords), "-dimensional locations\n")
  }

  if (!latlong) {
    if (isTRUE(verbose)) {
      cat("Using Euclidean norm to compute distance between locations\n")
    }
  } else {
    if (ncol(coords) != 2) {
      stop("With latlong option, there must only be 2 coordinate columns")
    } else if (isTRUE(verbose)) {
      cat("Computing distances on surface of sphere treating coordinates as latitude and longitude\n")
    }
  }

  coords
}

#' Calculate Euclidean Distance Matrix
#'
#' Computes the matrix of Euclidean distances between all pairs of locations
#' in the input matrix.
#'
#' @param s A numeric matrix where rows represent locations and columns represent dimensions.
#' @return A numeric matrix of pairwise Euclidean distances between locations.
#' @keywords internal
get_distmat_euclidean <- function(s) {
  s <- as.matrix(s)
  distances <- stats::dist(s)
  if (all(is.finite(distances) & distances >= 1e-150)) return(as.matrix(distances))

  # Recompute only pairs that can lose precision while squaring.
  out <- as.matrix(distances)
  n <- nrow(s)
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    j <- j[!is.finite(out[i, j]) | out[i, j] < 1e-150]
    if (!length(j)) next
    delta <- sweep(s[j, , drop = FALSE], 2L, s[i, ], "-")
    scale <- apply(abs(delta), 1L, max)
    if (any(!is.finite(scale))) stop("Pairwise distance exceeds the floating-point range.")
    d <- scale * sqrt(rowSums((delta / ifelse(scale > 0, scale, 1))^2))
    if (any(!is.finite(d))) stop("Pairwise distance exceeds the floating-point range.")
    out[i, j] <- out[j, i] <- d
  }
  out
}

#' Calculate Distance Matrix for Latitude/Longitude Points
#'
#' Computes the matrix of great-circle distances between all pairs of locations
#' on a sphere using latitude and longitude coordinates.
#'
#' @param s A numeric matrix where rows represent locations and columns are
#'   `lat` and `lon`.
#' @return A numeric matrix of pairwise great-circle distances as fractions of
#'   the circumference.
#' @keywords internal
get_distmat_lat_lon <- function(s) {
  n <- nrow(s)
  distmat <- matrix(0, nrow = n, ncol = n)

  radians <- as.matrix(s) * (pi / 180)
  cos_lat <- cos(radians[, 1])

  for (i in seq_len(n)) {
    lat_sin <- sin((radians[i, 1] - radians[, 1]) * 0.5)
    lon_sin <- sin((radians[i, 2] - radians[, 2]) * 0.5)
    h <- lat_sin * lat_sin + cos_lat[i] * cos_lat * lon_sin * lon_sin
    distmat[, i] <- (1 / pi) * asin(sqrt(pmin(pmax(h, 0), 1)))
  }

  distmat
}

#' @keywords internal
.normalized_distmat <- function(coords, latlong) {
  .spur_distance_setup(coords, latlong)$normalized
}

.spur_distance_setup <- function(coords, latlong) {
  distmat <- if (isTRUE(latlong)) {
    get_distmat_lat_lon(coords)
  } else {
    get_distmat_euclidean(coords)
  }
  max_dist <- max(distmat)
  if (!is.finite(max_dist) || max_dist <= 0) {
    stop("Distance matrix must contain positive finite distances.")
  }
  list(normalized = distmat / max_dist, max_dist = max_dist)
}

#' @keywords internal
.format_spur_test_output <- function(title, stat_label, teststat, pvalue, verbose = FALSE) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat(title, "\n", sep = "")
  cat("---------------------------------------\n")
  cat(sprintf("Test Statistic (%s) : %9.4f\n", stat_label, teststat))
  cat(sprintf("P-value               : %9.4f\n", pvalue))
  cat("---------------------------------------\n")
  invisible(NULL)
}

#' @keywords internal
.draw_emat <- function(q, nrep) {
  matrix(stats::rnorm(q * nrep), nrow = q, ncol = nrep)
}
