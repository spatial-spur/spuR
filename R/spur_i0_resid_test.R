#' Core SPUR I(0) Residual Test
#'
#' Internal translation of SPUR Mata `spatial_i0_test_residual`.
#'
#' @param Y Numeric matrix (n x n_series).
#' @param X_in Numeric design matrix including constant.
#' @param distmat Numeric normalized distance matrix (n x n).
#' @param emat Numeric simulation matrix (q x nrep).
#' @param sigdm_bm Optional precomputed low-frequency covariance matrix.
#' @param R Optional precomputed low-frequency eigenvector basis.
#' @param R_M Optional basis after projecting out the regression design.
#' @param nthreads Optional positive worker-thread count; NULL uses package settings and applies core limits.
#' @return List with test outputs.
#' @keywords internal
spur_i0_resid_test <- function(Y, X_in, distmat, emat, sigdm_bm = NULL, R = NULL, R_M = NULL, nthreads = NULL) {
  nthreads <- .get_spur_threads(nthreads)
  q <- nrow(emat)
  if (q > nrow(X_in) - ncol(X_in)) {
    stop("q exceeds the residual degrees of freedom.")
  }
  storage.mode(distmat) <- "double"
  storage.mode(emat) <- "double"

  if (is.null(sigdm_bm) || is.null(R) || is.null(R_M)) {
    qr_x <- qr(X_in)
    if (qr_x$rank < ncol(X_in)) {
      stop("Design matrix is rank-deficient; residual tests require full column rank.")
    }
    storage.mode(X_in) <- "double"
    X_in <- qr.Q(qr_x)
    P_X <- X_in
    storage.mode(P_X) <- "double"

    rho_bm <- 0.999
    c_bm <- get_cbar(rho_bm, distmat)
    sigdm_bm <- rcpp_get_sigma_residual_proj(distmat, as.numeric(c_bm), X_in, P_X)
    R <- get_r(sigdm_bm, q)$R
    R_M <- R - P_X %*% crossprod(X_in, R)
  }
  storage.mode(R_M) <- "double"

  rho <- 0.001
  rho_grid <- seq(from = 0.0001, to = 0.03, length.out = 30)
  c_all <- get_cbar(c(rho, rho_grid), distmat)
  om_all <- rcpp_om_grid(distmat, R_M, c_all, as.integer(nthreads))

  om_rho <- om_all[, , 1L]
  om_bm <- crossprod(R, sigdm_bm %*% R)

  om_i0 <- om_rho
  om_ho <- om_rho
  ha_parm <- get_ha_param_i0(om_ho, om_i0, om_bm, emat)
  om_ha <- om_i0 + ha_parm * om_bm

  ch_omi_ho <- chol(solve(om_ho))
  ch_omi_ha <- chol(solve(om_ha))

  n_y <- ncol(Y)
  LR <- numeric(n_y)
  for (i in seq_len(n_y)) {
    X <- Y[, i] - mean(Y[, i])
    P <- crossprod(R, X)
    y_P_ho <- ch_omi_ho %*% P
    y_P_ha <- ch_omi_ha %*% P
    q_P_ho <- sum(y_P_ho^2)
    q_P_ha <- sum(y_P_ha^2)
    LR[i] <- q_P_ho / q_P_ha
  }

  ch_om_ho_mat <- array(NA_real_, dim = c(q, q, length(rho_grid)))
  for (i in seq_along(rho_grid)) {
    om_ho_i <- if (rho_grid[i] > 0) om_all[, , i + 1L] else diag(q)
    ch_om_ho_mat[, , i] <- chol(om_ho_i)
  }

  sz_vec <- c(0.01, 0.05, 0.10)
  mc_out <- rcpp_mc_lr_grid(
    ch_om_ho_mat,
    ch_omi_ho,
    ch_omi_ha,
    emat,
    LR,
    1 - sz_vec,
    as.integer(nthreads)
  )
  cvalue_mat <- mc_out[[1L]]
  pvalue_mat <- mc_out[[2L]]

  list(
    LR = LR,
    pvalue = apply(pvalue_mat, 2, max),
    cvalue = apply(cvalue_mat, 2, max),
    ha_parm = ha_parm,
    cvalue_mat = cvalue_mat,
    pvalue_mat = pvalue_mat,
    rho_grid = rho_grid
  )
}
