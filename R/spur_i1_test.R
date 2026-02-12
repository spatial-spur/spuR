#' Core SPUR I(1) Test Function
#'
#' This function implements the spatial I(1) test. It is called internally by \code{spur_i1()}.
#'
#' @param Y A numeric matrix (n x number_of_series).
#' @param distmat A numeric distance matrix (n x n).
#' @param emat A numeric simulation matrix (q x number_of_simulations).
#'
#' @return A list with elements: LR, pvalue, cvalue, ha_parm, cvalue_mat, pvalue_mat, and rho_grid.

spur_i1_test <- function(Y, distmat, emat) {
  q <- nrow(emat)

  sigdm_bm <- get_sigma_lbm_dm(distmat)
  R <- get_r(sigdm_bm, q)$R

  om_ho <- t(R) %*% sigdm_bm %*% R
  ha_parm <- get_ha_param_i1(om_ho, distmat, R, emat)
  sigdm_ha <- get_sigma_dm(distmat, ha_parm)
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
  cv_vec <- stats::quantile(lr_ho, probs = 1 - sz_vec, names = FALSE)

  n_y <- ncol(Y)
  LR <- numeric(n_y)
  pval <- numeric(n_y)

  for (i in seq_len(n_y)) {
    X <- Y[, i] - mean(Y[, i])
    P <- t(R) %*% X
    LR[i] <- (t(P) %*% omi_ho %*% P) / (t(P) %*% omi_ha %*% P)
    pval[i] <- mean(lr_ho > LR[i])
  }

  list(
    LR = as.vector(LR),
    pvalue = as.vector(pval),
    cv_vec = as.vector(cv_vec),
    ha_parm = ha_parm
  )
}
