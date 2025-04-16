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
#' @export
getcbar <- function(rhobar, distmat) {
  # Extract the lower-triangular part (excluding the diagonal) of distmat
  vd <- lvech(distmat)
  
  # Initialize c0 and c1
  c0 <- 10
  c1 <- 10
  
  # Adjust c0 so that mean(exp(-c0 * vd)) > rhobar
  i1 <- FALSE
  jj <- 0
  while (!i1) {
    v <- mean(exp(-c0 * vd))
    i1 <- (v > rhobar)
    if (!i1) {
      c1 <- c0
      c0 <- c0 / 2
      jj <- jj + 1
    }
    if (jj > 500) {
      stop("rhobar too large")
    }
  }
  
  # Adjust c1 so that mean(exp(-c1 * vd)) < rhobar
  i1 <- FALSE
  jj <- 0
  while (!i1) {
    v <- mean(exp(-c1 * vd))
    i1 <- (v < rhobar)
    if (!i1) {
      c0 <- c1
      c1 <- 2 * c1
      jj <- jj + 1
    }
    if (c1 > 10000) {  # safeguard to break if c1 becomes unreasonably large
      i1 <- TRUE
    }
    if (jj > 500) {
      stop("rhobar too small")
    }
  }
  
  # Bisection method to determine cbar
  while ((c1 - c0) > 0.001) {
    cm <- sqrt(c0 * c1)
    v <- mean(exp(-cm * vd))
    if (v < rhobar) {
      c1 <- cm
    } else {
      c0 <- cm
    }
  }
  
  cbar <- sqrt(c0 * c1)
  return(cbar)
}
