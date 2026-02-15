#' Compute Demeaned Sigma (c)
#'
#' Computes a demeaned version of \eqn{\sigma(c)}, where
#' \eqn{\sigma(c) = \exp(-c \times \text{distmat})}, and then demeans the resulting matrix.
#'
#' @param distmat A numeric matrix of distances.
#' @param c A numeric constant used in the exponential transformation.
#' @return A numeric matrix representing the demeaned sigma.
#' @keywords internal
get_sigma_dm <- function(distmat, c) {
  n <- nrow(distmat)
  sigma <- exp(-c * distmat)
  sigma_dm <- sigma - matrix(colMeans(sigma), n, n, byrow = TRUE)
  sigma_dm <- sigma_dm - matrix(rowMeans(sigma_dm), n, n)
  sigma_dm
}

#' Compute LBM Covariance Matrix
#'
#' Computes the LBM covariance matrix from the given distance matrix, using the
#' first location as the origin.
#'
#' @param distmat A numeric matrix of distances.
#' @return A numeric matrix representing the LBM covariance matrix.
#' @keywords internal
get_sigma_lbm <- function(distmat) {
  n <- nrow(distmat)
  first_col <- matrix(distmat[, 1], n, n)
  first_row <- matrix(distmat[1, ], n, n, byrow = TRUE)
  0.5 * (first_col + first_row - distmat)
}

#' Compute Demeaned Sigma LBM
#'
#' Computes the sigma_lbm matrix using the provided distance matrix and then
#' demeans it.
#'
#' @param distmat A numeric matrix of distances.
#' @return A numeric matrix representing the demeaned sigma_lbm.
#' @keywords internal
get_sigma_lbm_dm <- function(distmat) {
  demean_sigma(get_sigma_lbm(distmat))
}

#' Residual Covariance Matrix
#'
#' Translation of SPUR Mata `get_sigma_residual.mata`.
#'
#' @param distmat A distance matrix.
#' @param c Spatial persistence parameter.
#' @param M Residual-maker matrix.
#' @return Residual covariance matrix.
#' @keywords internal
get_sigma_residual <- function(distmat, c, M) {
  sigma <- exp(-c * distmat)
  M %*% sigma %*% t(M)
}
