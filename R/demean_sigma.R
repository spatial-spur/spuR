#' Demean a Covariance Matrix
#'
#' Given a covariance matrix \code{sigma}, this function computes the variance of
#' \eqn{X - \bar{X}}, effectively demeaning the matrix by removing the influence of the mean.
#'
#' @param sigma A numeric covariance matrix.
#' @return A numeric matrix that has been demeaned.
#' @export
demean_sigma <- function(sigma) {
  n <- nrow(sigma)
  # Subtract the column means replicated for each row
  sigma_dm <- sigma - matrix(colMeans(sigma), n, n, byrow = TRUE)
  # Subtract the row means replicated for each column
  sigma_dm <- sigma_dm - matrix(rowMeans(sigma_dm), n, n)
  return(sigma_dm)
}
