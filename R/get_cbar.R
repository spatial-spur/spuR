#' Compute cbar Value
#'
#' Determines the value of \code{cbar} such that \eqn{mean(exp(-cbar \times vd))} is approximately equal to
#' a given \code{rhobar}, where \code{vd} is the vectorized lower-triangular portion (including the diagonal)
#' of the distance matrix \code{distmat}. This function uses iterative adjustments and bisection to find
#' an appropriate value.
#'
#' @param rhobar A numeric scalar representing the target value.
#' @param distmat A numeric matrix of distances.
#' @return A numeric scalar \code{cbar} satisfying the desired property.
#' @keywords internal
get_cbar <- function(rhobar, distmat) {
  storage.mode(distmat) <- "double"
  rcpp_get_cbar_vec(as.numeric(rhobar), distmat)
}
