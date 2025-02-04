#' Compute Eigenvectors for Low-Frequency Weights
#'
#' Constructs eigenvectors corresponding to the largest \code{qmax} eigenvalues
#' of the \code{n x n} matrix \code{sig}. The eigenvectors are normalized such that
#' \code{t(R) \%*\% R = I} (if the optional normalization is applied).
#'
#' @param sig A numeric matrix.
#' @param qmax An integer specifying the number of largest eigenvalues/eigenvectors to retain.
#' @return A list with components:
#' \describe{
#'   \item{R}{A matrix whose columns are the eigenvectors corresponding to the largest \code{qmax} eigenvalues.}
#'   \item{DS}{A numeric vector containing the largest \code{qmax} eigenvalues.}
#' }
#' @export
Get_R <- function(sig, qmax) {
  n <- nrow(sig)
  eig_out <- eigen(sig)
  values <- eig_out$values
  vectors <- eig_out$vectors
  
  # Sort eigenvalues in descending order and get the sorting indices
  idx <- order(values, decreasing = TRUE)
  sorted_values <- values[idx]
  
  DS <- sorted_values[1:qmax]
  R <- vectors[, idx[1:qmax], drop = FALSE]
  
  # Optional normalization: Normalize so that t(R) %*% R/n = 1
  # Uncomment the following lines if such normalization is desired.
  # rr <- sum(R[, 1]^2)
  # scl <- sqrt(n / rr)
  # R <- scl * R
  
  return(list(R = R, DS = DS))
}
