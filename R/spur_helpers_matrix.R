#' Demean Covariance Matrix
#'
#' Given a covariance matrix \code{sigma}, this function computes the variance of
#' \eqn{X - \bar{X}}, effectively demeaning the matrix by removing the influence
#' of the mean.
#'
#' @param sigma A numeric covariance matrix.
#' @return A numeric matrix that has been demeaned.
#' @keywords internal
demean_sigma <- function(sigma) {
  n <- nrow(sigma)
  sigma_dm <- sigma - matrix(colMeans(sigma), n, n, byrow = TRUE)
  sigma_dm <- sigma_dm - matrix(rowMeans(sigma_dm), n, n)
  sigma_dm
}

#' Lower Triangular Vectorization
#'
#' Extracts the elements of a matrix that are strictly below the main diagonal
#' and returns them as a vector.
#'
#' @param S A numeric matrix.
#' @return A numeric vector containing the lower-triangular elements.
#' @keywords internal
lvech <- function(S) {
  S[lower.tri(S)]
}

#' @keywords internal
.safe_vec_norm <- function(x) {
  x <- as.double(x)
  if (!length(x)) {
    return(0)
  }
  max_abs <- max(abs(x))
  if (!is.finite(max_abs)) {
    return(NA_real_)
  }
  if (max_abs == 0) {
    return(0)
  }
  max_abs * sqrt(sum((x / max_abs)^2))
}

#' @keywords internal
.orthonormal_design_basis <- function(X_in, tol = 1e-10) {
  if (!is.matrix(X_in)) {
    X_in <- as.matrix(X_in)
  }
  if (!is.numeric(X_in)) {
    stop("Design matrix must be numeric.")
  }
  n <- nrow(X_in)
  k_x <- ncol(X_in)
  if (is.null(n) || is.null(k_x) || n < 2L || k_x < 1L || k_x >= n) {
    stop("Invalid residual test design matrix dimensions; residual tests require 1 <= ncol(X) < n.")
  }
  if (any(!is.finite(X_in))) {
    stop("Design matrix must contain only finite values.")
  }

  # Scale each column to unit L2 norm first so rank detection and M = I - Q_x Q_x'
  # are invariant to regressor units across arbitrarily small or large column scales.
  q_mat <- matrix(0, nrow = n, ncol = k_x)
  for (j in seq_len(k_x)) {
    col <- as.double(X_in[, j])
    max_abs <- max(abs(col))
    if (!is.finite(max_abs) || max_abs == 0) {
      stop("Design matrix is rank-deficient; residual tests require full column rank.")
    }
    col_scaled <- col / max_abs
    scaled_norm <- sqrt(sum(col_scaled^2))
    if (!is.finite(scaled_norm) || scaled_norm == 0) {
      stop("Design matrix is rank-deficient; residual tests require full column rank.")
    }
    q_mat[, j] <- col_scaled / scaled_norm
  }

  for (j in seq_len(k_x)) {
    if (j > 1L) {
      for (pass in seq_len(2L)) {
        for (prev in seq_len(j - 1L)) {
          dot <- sum(q_mat[, prev] * q_mat[, j])
          q_mat[, j] <- q_mat[, j] - dot * q_mat[, prev]
        }
      }
    }
    norm_j <- sqrt(sum(q_mat[, j]^2))
    if (!is.finite(norm_j) || norm_j <= tol) {
      stop("Design matrix is rank-deficient; residual tests require full column rank.")
    }
    q_mat[, j] <- q_mat[, j] / norm_j
  }

  q_mat
}

#' @keywords internal
.extract_spur_basis <- function(sigdm, q) {
  r_out <- get_r(sigdm, q)
  ds <- r_out$DS
  lead <- ds[1L]
  tail <- ds[q]
  if (!is.finite(lead) || lead <= 0 || !is.finite(tail) || tail <= 1e-12 * lead) {
    stop("SPUR covariance basis has insufficient positive rank for requested `q`.")
  }
  r_out$R
}

#' @keywords internal
.validate_spur_y <- function(Y, name = "Y") {
  if (is.null(dim(Y))) {
    Y <- matrix(Y, ncol = 1L)
  } else {
    Y <- as.matrix(Y)
  }
  if (!is.numeric(Y) || nrow(Y) < 2L || ncol(Y) < 1L) {
    stop(sprintf("Invalid response `%s` dimensions.", name))
  }
  if (any(!is.finite(Y))) {
    stop(sprintf("Response `%s` must contain only finite values.", name))
  }

  n <- nrow(Y)
  n_y <- ncol(Y)
  y_dm <- matrix(0, nrow = n, ncol = n_y)
  y_dm_norm <- numeric(n_y)

  for (j in seq_len(n_y)) {
    col <- as.double(Y[, j])
    if (max(col) - min(col) == 0) {
      stop(sprintf("Response `%s` has zero variation; SPUR inference requires non-constant `%s`.", name, name))
    }
    col_dm <- col - mean(col)
    dm_max <- max(abs(col_dm))
    if (!is.finite(dm_max) || dm_max == 0) {
      stop(sprintf("Response `%s` is effectively constant; SPUR inference requires variation in `%s`.", name, name))
    }
    y_dm[, j] <- col_dm / dm_max
    y_dm_norm[j] <- sqrt(sum(y_dm[, j]^2))
  }

  list(Y_dm = y_dm, y_dm_norm = y_dm_norm)
}

#' @keywords internal
.project_spur_y <- function(y_info, R) {
  p_mat <- crossprod(R, y_info$Y_dm)
  for (j in seq_len(ncol(p_mat))) {
    p_norm <- .safe_vec_norm(p_mat[, j])
    if (!is.finite(p_norm) || p_norm <= 1e-12 * y_info$y_dm_norm[j]) {
      stop("Projected SPUR vector has zero norm.")
    }
  }
  p_mat
}

#' @keywords internal
.spur_shapes_differ <- function(a, b) {
  ta <- sum(diag(a))
  tb <- sum(diag(b))
  if (!is.finite(ta) || !is.finite(tb) || ta <= 0 || tb <= 0) {
    stop("Projected SPUR covariance has non-positive trace.")
  }
  max(abs(a / ta - b / tb)) > 1e-12
}

#' @keywords internal
.validate_spur_i1_shapes <- function(distmat, R, M = NULL) {
  mr <- if (is.null(M)) sweep(R, 2L, colMeans(R), "-") else M %*% R
  shapes <- lapply(c(0.01, 0.5, 0.95), function(rho) {
    c_val <- get_cbar(rho, distmat)
    crossprod(mr, exp(-c_val * distmat) %*% mr)
  })
  if (!.spur_shapes_differ(shapes[[1L]], shapes[[2L]]) &&
      !.spur_shapes_differ(shapes[[1L]], shapes[[3L]])) {
    stop("SPUR experiment is uninformative: covariance shapes are proportional.")
  }
}

#' @keywords internal
.validate_core_spur_inputs <- function(distmat, emat, min_q = 2L) {
  if (!is.matrix(distmat) || !is.numeric(distmat) || nrow(distmat) != ncol(distmat) || nrow(distmat) < 2L) {
    stop("`distmat` must be a square numeric matrix with at least 2 rows.")
  }
  if (any(!is.finite(distmat)) || any(distmat < 0)) {
    stop("`distmat` must contain finite non-negative distances.")
  }
  if (!is.matrix(emat) || !is.numeric(emat)) {
    stop("`emat` must be a numeric matrix.")
  }
  q <- nrow(emat)
  nrep <- ncol(emat)
  n <- nrow(distmat)
  if (q < min_q) {
    if (min_q >= 2L) {
      stop("Informative SPUR inference requires `q >= 2`.")
    }
    stop("`q` must be at least 1.")
  }
  if (q >= n || nrep < 1L || any(!is.finite(emat))) {
    stop(sprintf(
      "`emat` dimensions (%d x %d) are invalid for n = %d (require %d <= q < n, nrep >= 1).",
      q, nrep, n, min_q
    ))
  }
  invisible(list(n = n, q = q, nrep = nrep))
}
