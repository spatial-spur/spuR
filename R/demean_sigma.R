#' Demean a Covariance Matrix
#'
#' Given a covariance matrix \code{sigma}, this function computes the variance of
#' \eqn{X - \bar{X}}, effectively demeaning the matrix by removing the influence of the mean.
#'
#' @param sigma A numeric covariance matrix.
#' @return A numeric matrix that has been demeaned.
demean_sigma <- function(sigma) {  
  n <- nrow(sigma)
  # Subtract column means - equivalent to repmat(mean(sigma,1),n,1)
  sigma_dm <- sigma - matrix(colMeans(sigma), n, n, byrow = TRUE)
  # Subtract row means - equivalent to repmat(mean(sigma_dm,2),1,n)
  sigma_dm <- sigma_dm - matrix(rowMeans(sigma_dm), n, n)
  return(sigma_dm)
}