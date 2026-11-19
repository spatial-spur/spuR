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
  storage.mode(sig) <- "double"
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
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    res <- RSpectra::eigs_sym(
      sig,
      k = qmax,
      which = "LA",
      opts = list(initvec = cos(seq_len(n) * 0.6180339887498949))
    )
    if (length(res$values) < qmax || ncol(res$vectors) < qmax) {
      eigen(sig, symmetric = TRUE)
    } else {
      res
    }
  } else {
    eigen(sig, symmetric = TRUE)
  }
  idx <- order(eig_out$values, decreasing = TRUE)
  values <- eig_out$values[idx][seq_len(qmax)]
  vectors <- eig_out$vectors[, idx, drop = FALSE][, seq_len(qmax), drop = FALSE]

  # Canonicalize eigenvector signs to match spatial-core: largest-magnitude
  # entry positive (earliest index within a roundoff tolerance).
  for (j in seq_len(qmax)) {
    col <- vectors[, j]
    magnitudes <- abs(col)
    peak <- max(magnitudes)
    max_idx <- which(magnitudes >= peak - 64 * .Machine$double.eps * peak)[1L]
    if (length(max_idx) == 1L && col[max_idx] < 0) {
      vectors[, j] <- -col
    }
  }

  list(
    R = vectors,
    DS = values
  )
}
