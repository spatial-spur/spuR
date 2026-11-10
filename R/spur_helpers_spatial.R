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
  if (!is.matrix(coords)) {
    coords <- as.matrix(coords)
  }
  if (anyNA(coords)) {
    stop("Missing value(s) in coordinate variables; aborting")
  }
  if (!all(is.finite(coords))) {
    stop("Non-finite value(s) in coordinate variables; aborting")
  }

  ns <- nrow(coords)
  if (ns < 5) {
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
  storage.mode(s) <- "double"
  .Call(cpp_distmat_euclidean, s)
}

#' Calculate Distance Matrix for Latitude/Longitude Points
#'
#' Computes the matrix of great-circle distances between all pairs of locations
#' on a sphere using latitude and longitude coordinates.
#'
#' @param s A numeric matrix where rows represent locations and columns are
#'   `lat` and `lon`.
#' @return A numeric matrix of pairwise great-circle distances in normalized
#'   units, matching the SPUR Mata `getdistmat()` implementation.
#' @keywords internal
get_distmat_lat_lon <- function(s) {
  s <- as.matrix(s)
  storage.mode(s) <- "double"
  .Call(cpp_distmat_lat_lon, s)
}

#' Cluster Vector as Matrix
#'
#' Internal helper to convert cluster IDs into a column matrix and optionally
#' subset rows using a logical index.
#'
#' @param cluster_values Vector of cluster IDs.
#' @param use_rows Optional logical index vector.
#' @return A matrix (column vector) of cluster values for selected observations.
#' @keywords internal
get_cluster_matrix <- function(cluster_values, use_rows = NULL) {
  if (!is.null(use_rows)) {
    cluster_values <- cluster_values[use_rows]
  }
  as.matrix(cluster_values)
}

#' @keywords internal
.normalized_distmat <- function(coords, latlong) {
  distmat <- if (isTRUE(latlong)) {
    get_distmat_lat_lon(coords)
  } else {
    get_distmat_euclidean(coords)
  }
  max_dist <- max(distmat)
  if (!is.finite(max_dist) || max_dist <= 0) {
    stop("Distance matrix must contain positive finite distances.")
  }
  out <- distmat / max_dist
  attr(out, "max_dist") <- max_dist
  out
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
