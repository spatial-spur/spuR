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
.extract_spur_basis <- function(sigdm, q, distmat, eig = NULL) {
  r_out <- tryCatch(if (is.null(eig)) get_r(sigdm, q) else .select_spur_eigenvectors(eig, q),
    error = function(e) stop(conditionMessage(e), " ", .spur_site_summary(distmat), call. = FALSE))
  r_out
}

.spur_site_summary <- function(d) {
  sites <- sum(vapply(seq_len(nrow(d)), function(i) !any(d[i, seq_len(i - 1L)] == 0), logical(1L)))
  sprintf("%d distinct locations, co-located-pair fraction %.6f.", sites, mean(lvech(d) == 0))
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
    if (max(col) == min(col)) {
      stop(sprintf("Response `%s` has zero variation; SPUR inference requires non-constant `%s`.", name, name))
    }
    centered <- col - (0.5 * min(col) + 0.5 * max(col))
    scaled <- centered / max(abs(centered))
    y_dm[, j] <- scaled - mean(scaled)
    y_dm_norm[j] <- sqrt(sum(y_dm[, j]^2))
  }

  list(Y_dm = y_dm, y_dm_norm = y_dm_norm)
}

#' @keywords internal
.residualize_spur <- function(values, q_x) {
  values - q_x %*% crossprod(q_x, values)
}

.spur_constant_residual <- function(q_x, design) {
  n <- nrow(q_x)
  constant <- rep(1 / sqrt(n), n)
  residual <- .residualize_spur(.residualize_spur(constant, q_x), q_x)
  # Squared norm bounds the fraction of the constant covariance left by projection.
  residual_norm <- sqrt(sum(residual^2))
  eps <- .Machine$double.eps
  if (residual_norm <= 8 * eps) return(NULL)
  scaled <- sweep(design, 2L, apply(abs(design), 2L, max), "/")
  fit <- svd(scaled, nv = 0L)
  coeff_norm <- sqrt(sum((colSums(fit$u) / (sqrt(n) * fit$d))^2))
  if (residual_norm > 64 * eps * (1 + fit$d[1L] * coeff_norm)) {
    stop("SPUR residual tests require a design that spans a constant.", call. = FALSE)
  }
  as.numeric(residual) * sqrt(n)
}

#' @keywords internal
.validate_spur_residual_y <- function(Y, q_x, constant_residual, name = "Y") {
  if (is.vector(Y)) {
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
  y_res <- matrix(0, nrow = n, ncol = n_y)
  y_res_norm <- numeric(n_y)

  for (j in seq_len(n_y)) {
    col <- as.double(Y[, j])
    if (max(col) == min(col)) {
      stop(sprintf("Response `%s` has zero variation; SPUR inference requires non-constant `%s`.", name, name))
    }
    anchor <- 0.5 * min(col) + 0.5 * max(col)
    centered <- col - anchor
    scale <- max(abs(centered))
    allowance <- 0
    if (is.null(constant_residual)) {
      col_scaled <- centered / scale
    } else {
      correction <- anchor * constant_residual
      scale <- max(scale, abs(correction))
      col_scaled <- centered / scale + correction / scale
      allowance <- (8 * .Machine$double.eps * abs(anchor) / scale) * sqrt(n)
    }
    col_norm <- .safe_vec_norm(col_scaled)
    u <- .residualize_spur(.residualize_spur(col_scaled, q_x), q_x)
    u_max <- max(abs(u))
    u_norm <- .safe_vec_norm(u)
    if (!is.finite(u_max) || u_max == 0 || !is.finite(u_norm) ||
        u_norm <= 1e-12 * col_norm + allowance) {
      stop(sprintf("Residualized response `%s` has zero variation; SPUR inference requires variation orthogonal to `X_in`.", name))
    }
    y_res[, j] <- u / u_max
    y_res_norm[j] <- u_norm / u_max
  }

  list(Y_dm = y_res, y_dm_norm = y_res_norm)
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
.validate_core_spur_inputs <- function(distmat, emat, min_q = 2L, blocks = 1L) {
  if (!is.matrix(distmat) || !is.numeric(distmat) || nrow(distmat) != ncol(distmat) || nrow(distmat) < 2L) {
    stop("`distmat` must be a square numeric matrix with at least 2 rows.")
  }
  if (any(!is.finite(distmat)) || any(distmat < 0)) {
    stop("`distmat` must contain finite non-negative distances.")
  }
  if (!is.matrix(emat) || !is.numeric(emat)) {
    stop("`emat` must be a numeric matrix.")
  }
  if (nrow(emat) %% blocks != 0L) {
    stop("Diagnostic `emat` must have two q-row blocks.")
  }
  q <- nrow(emat) %/% blocks
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
