#' Core SPUR I(1) Residual Test
#'
#' Internal translation of SPUR Mata `spatial_i1_test_residual`.
#'
#' @param Y Numeric matrix (n x n_series).
#' @param X_in Numeric design matrix including constant.
#' @param distmat Numeric normalized distance matrix (n x n).
#' @param emat Numeric simulation matrix (q x nrep).
#' @return List with test outputs.
#' @keywords internal
spur_i1_resid_test <- function(Y, X_in, distmat, emat) {
  q <- nrow(emat)
  n <- nrow(distmat)
  n_y <- ncol(Y)

  xtx <- crossprod(X_in)
  if (qr(xtx)$rank < ncol(xtx)) {
    stop("Design matrix is rank-deficient; residual tests require full column rank.")
  }
  M <- diag(n) - X_in %*% solve(xtx, t(X_in))

  rho_bm <- 0.999
  c_bm <- get_cbar(rho_bm, distmat)
  sigdm_bm <- get_sigma_residual(distmat, c_bm, M)

  R <- get_r(sigdm_bm, q)$R
  om_ho <- t(R) %*% sigdm_bm %*% R

  ha_parm <- get_ha_param_i1_residual(om_ho, distmat, R, emat, M)
  sigdm_ha <- get_sigma_residual(distmat, ha_parm, M)
  om_ha <- t(R) %*% sigdm_ha %*% R

  ch_om_ho <- chol(om_ho)
  omi_ho <- solve(om_ho)
  omi_ha <- solve(om_ha)
  ch_omi_ho <- chol(omi_ho)
  ch_omi_ha <- chol(omi_ha)

  y_ho <- t(ch_om_ho) %*% emat
  y_ho_ho <- ch_omi_ho %*% y_ho
  y_ho_ha <- ch_omi_ha %*% y_ho
  q_ho_ho <- colSums(y_ho_ho^2)
  q_ho_ha <- colSums(y_ho_ha^2)
  lr_ho <- q_ho_ho / q_ho_ha

  sz_vec <- c(0.01, 0.05, 0.10)
  cv_vec <- as.numeric(stats::quantile(lr_ho, probs = 1 - sz_vec))

  LR <- numeric(n_y)
  pvalue <- numeric(n_y)
  for (i in seq_len(n_y)) {
    X <- Y[, i] - mean(Y[, i])
    P <- t(R) %*% X
    LR[i] <- (t(P) %*% omi_ho %*% P) / (t(P) %*% omi_ha %*% P)
    pvalue[i] <- mean(lr_ho > LR[i])
  }

  list(
    LR = as.vector(LR),
    pvalue = as.vector(pvalue),
    cv_vec = cv_vec,
    ha_parm = ha_parm
  )
}
