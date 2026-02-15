#' Demean Covariance Matrix
#'
#' Given a covariance matrix \code{sigma}, this function computes the variance of
#' \eqn{X - \bar{X}}, effectively demeaning the matrix by removing the influence
#' of the mean.
#'
#' @param sigma A numeric covariance matrix.
#' @return A numeric matrix that has been demeaned.
#' @keywords internal
demean_sigma <- function(sigma) {
  n <- nrow(sigma)
  sigma_dm <- sigma - matrix(colMeans(sigma), n, n, byrow = TRUE)
  sigma_dm <- sigma_dm - matrix(rowMeans(sigma_dm), n, n)
  sigma_dm
}

#' Lower Triangular Vectorization
#'
#' Extracts the elements of a matrix that are strictly below the main diagonal
#' and returns them as a vector.
#'
#' @param S A numeric matrix.
#' @return A numeric vector containing the lower-triangular elements.
#' @keywords internal
lvech <- function(S) {
  S[lower.tri(S)]
}
