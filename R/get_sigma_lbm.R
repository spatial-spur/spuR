#' Compute LBM Covariance Matrix
#'
#' Computes the LBM covariance matrix from the given distance matrix, using the first location as the origin.
#'
#' @param distmat A numeric matrix of distances.
#' @return A numeric matrix representing the LBM covariance matrix.
get_sigma_lbm <- function(distmat) {
  n <- nrow(distmat)
  # Replicate the first column and first row of the distance matrix
  first_col <- matrix(distmat[, 1], n, n)
  first_row <- matrix(distmat[1, ], n, n, byrow = TRUE)
  
  sigma_lbm <- 0.5 * (first_col + first_row - distmat)
  return(sigma_lbm)
}
