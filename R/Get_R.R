#' Compute Eigenvectors for Low-Frequency Weights
#'
#' Constructs the eigenvectors corresponding to the largest \code{qmax}
#' eigenvalues of \code{sig}. This follows the ordering in the Stata/Mata
#' implementation (\code{order(d, -1)}).
#'
#' @param sig A numeric matrix.
#' @param qmax An integer specifying the number of largest eigenvalues/eigenvectors to retain.
#' @return A list with components:
#' \describe{
#'   \item{R}{A matrix whose columns are the eigenvectors corresponding to the largest \code{qmax} eigenvalues.}
#'   \item{DS}{A numeric vector containing the largest \code{qmax} eigenvalues.}
#' }
#' @keywords internal
get_r <- function(sig, qmax) {
  n <- nrow(sig)
  if (length(qmax) != 1L || !is.finite(qmax) || qmax < 1 || qmax != as.integer(qmax)) {
    stop("qmax must be a positive integer.")
  }
  qmax <- as.integer(qmax)
  if (qmax > n) {
    stop("qmax must be <= nrow(sig).")
  }

  # Use iterative partial eigendecomposition only when matrix is large.
  # For small n, full eigen() is more stable for strict parity checks.
  use_rspectra <- qmax < n && n >= 500L
  eig_out <- if (use_rspectra) {
    RSpectra::eigs_sym(sig, k = qmax, which = "LM")
  } else {
    eigen(sig, symmetric = TRUE)
  }
  idx <- order(eig_out$values, decreasing = TRUE)
  values <- eig_out$values[idx]
  vectors <- eig_out$vectors[, idx, drop = FALSE]

  list(
    R = vectors[, seq_len(qmax), drop = FALSE],
    DS = values[seq_len(qmax)]
  )
}
