#' Core SPUR I(0) Residual Test
#'
#' Internal translation of SPUR Mata `spatial_i0_test_residual`.
#'
#' @param Y Numeric matrix (n x n_series).
#' @param X_in Numeric design matrix including constant.
#' @param distmat Numeric normalized distance matrix (n x n).
#' @param emat Numeric simulation matrix (q x nrep).
#' @return List with test outputs.
#' @keywords internal
spur_i0_resid_test <- function(Y, X_in, distmat, emat) {
  q <- nrow(emat)
  n <- nrow(distmat)

  M <- diag(n) - X_in %*% solve(t(X_in) %*% X_in) %*% t(X_in)

  rho_bm <- 0.999
  c_bm <- get_cbar(rho_bm, distmat)
  sigdm_bm <- get_sigma_residual(distmat, c_bm, M)

  R <- get_r(sigdm_bm, q)$R

  rho <- 0.001
  c_val <- get_cbar(rho, distmat)
  sigdm_rho <- get_sigma_residual(distmat, c_val, M)

  om_rho <- t(R) %*% sigdm_rho %*% R
  om_bm <- t(R) %*% sigdm_bm %*% R

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
    P <- t(R) %*% X
    y_P_ho <- ch_omi_ho %*% P
    y_P_ha <- ch_omi_ha %*% P
    q_P_ho <- sum(y_P_ho^2)
    q_P_ha <- sum(y_P_ha^2)
    LR[i] <- q_P_ho / q_P_ha
  }

  rho_grid <- seq(from = 0.0001, to = 0.03, length.out = 30)
  ch_om_ho_mat <- array(NA_real_, dim = c(q, q, length(rho_grid)))
  for (i in seq_along(rho_grid)) {
    om_ho_i <- diag(q)
    rho_i <- rho_grid[i]
    if (rho_i > 0) {
      c_i <- get_cbar(rho_i, distmat)
      sigdm_ho_i <- get_sigma_residual(distmat, c_i, M)
      om_ho_i <- t(R) %*% sigdm_ho_i %*% R
    }
    ch_om_ho_mat[, , i] <- chol(om_ho_i)
  }

  pvalue_mat <- matrix(NA_real_, nrow = length(rho_grid), ncol = n_y)
  sz_vec <- c(0.01, 0.05, 0.10)
  cvalue_mat <- matrix(NA_real_, nrow = length(rho_grid), ncol = length(sz_vec))

  for (ir in seq_along(rho_grid)) {
    ch_om_ho <- ch_om_ho_mat[, , ir]
    y_ho <- t(ch_om_ho) %*% emat
    y_ho_ho <- ch_omi_ho %*% y_ho
    y_ho_ha <- ch_omi_ha %*% y_ho
    q_ho_ho <- colSums(y_ho_ho^2)
    q_ho_ha <- colSums(y_ho_ha^2)
    lr_ho <- q_ho_ho / q_ho_ha

    cvalue_mat[ir, ] <- as.numeric(stats::quantile(lr_ho, probs = 1 - sz_vec))
    for (j in seq_len(n_y)) {
      pvalue_mat[ir, j] <- mean(lr_ho > LR[j])
    }
  }

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
