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
  use_rspectra <- .use_partial_spur_basis(n, qmax)
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
      k = qmax + 1L,
      which = "LA",
      opts = list(initvec = cos(seq_len(n) * 0.6180339887498949))
    )
    if (length(res$values) < qmax + 1L || ncol(res$vectors) < qmax + 1L) {
      eigen(sig, symmetric = TRUE)
    } else {
      res
    }
  } else {
    eigen(sig, symmetric = TRUE)
  }
  .select_spur_eigenvectors(eig_out, qmax)
}

.use_partial_spur_basis <- function(n, q) q + 1L < n && n >= 500L

.select_spur_eigenvectors <- function(eig, q) {
  order <- order(eig$values, decreasing = TRUE)
  sorted <- eig$values[order]
  gap_tol <- max(1.6e-11 * sorted[1L], 8 * .Machine$double.eps * nrow(eig$vectors))
  if (sorted[q] <= 1e-12 * sorted[1L] || sorted[1L] <= 0) {
    stop(sprintf("SPUR covariance basis has insufficient positive rank for requested q=%d. Positive rank is %d.", q, sum(sorted > 1e-12 * sorted[1L])))
  }
  if (length(sorted) > q && sorted[q] - sorted[q + 1L] <= gap_tol) {
    stop(sprintf("Requested q splits a tied eigenspace (q=%d). Choose q at an eigenvalue gap.", q))
  }
  idx <- order[seq_len(q)]
  values <- eig$values[idx]
  vectors <- eig$vectors[, idx, drop = FALSE]

  # Canonicalize eigenvector signs to match spatial-core: largest-magnitude
  # entry positive (earliest index within a roundoff tolerance).
  for (j in seq_len(q)) {
    col <- vectors[, j]
    magnitudes <- abs(col)
    peak <- max(magnitudes)
    max_idx <- which(peak - magnitudes <= 1e-12 * peak)[1L]
    if (length(max_idx) == 1L && col[max_idx] < 0) {
      vectors[, j] <- -col
    }
  }

  list(
    R = vectors,
    DS = values,
    tied = any(-diff(values) <= gap_tol),
    all_directions = q == length(sorted) || sorted[q + 1L] <= 1e-12 * sorted[1L]
  )
}
