#' Compute Demeaned Sigma (c)
#'
#' Computes a demeaned version of \eqn{\sigma(c)}, where
#' \eqn{\sigma(c) = \exp(-c \times \text{distmat})}, and then demeans the resulting matrix.
#'
#' @param distmat A numeric matrix of distances.
#' @param c A numeric constant used in the exponential transformation.
#' @return A numeric matrix representing the demeaned sigma.
#' @export
get_sigma_dm <- function(distmat, c) {
  n <- nrow(distmat)
  sigma <- exp(-c * distmat)
  sigma_dm <- sigma - matrix(colMeans(sigma), n, n, byrow = TRUE)
  sigma_dm <- sigma_dm - matrix(rowMeans(sigma_dm), n, n)
  return(sigma_dm)
}
