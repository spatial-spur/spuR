#' Core SPUR I(1) Test Function
#'
#' This function implements the spatial I(1) test. It is called internally by \code{spurtest_i1()}.
#'
#' @param Y A numeric matrix (n x number_of_series).
#' @param distmat A numeric distance matrix (n x n).
#' @param emat A numeric simulation matrix (q x number_of_simulations).
#'
#' @return A list with elements: LR, pvalue, cv_vec, and ha_parm.
#' @keywords internal
spur_i1_test <- function(Y, distmat, emat) {
  dims <- .validate_core_spur_inputs(distmat, emat, min_q = 2L)
  q <- dims$q
  y_info <- .validate_spur_y(Y, name = "Y")
  if (nrow(y_info$Y_dm) != dims$n) {
    stop("`Y` row count must match `distmat`.")
  }

  sigdm_bm <- get_sigma_lbm_dm(distmat)
  R <- .extract_spur_basis(sigdm_bm, q)
  P_mat <- .project_spur_y(y_info, R)
  .validate_spur_i1_shapes(distmat, R)

  om_ho <- crossprod(R, sigdm_bm %*% R)
  om_ho <- 0.5 * (om_ho + t(om_ho))
  ha_parm <- get_ha_param_i1(om_ho, distmat, R, emat)
  sigdm_ha <- get_sigma_dm(distmat, ha_parm)
  om_ha <- crossprod(R, sigdm_ha %*% R)
  om_ha <- 0.5 * (om_ha + t(om_ha))

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

  n_y <- ncol(P_mat)
  LR <- numeric(n_y)
  pval <- numeric(n_y)

  for (i in seq_len(n_y)) {
    P <- P_mat[, i, drop = FALSE]
    q_P_ho <- sum((ch_omi_ho %*% P)^2)
    q_P_ha <- sum((ch_omi_ha %*% P)^2)
    if (!is.finite(q_P_ho) || !is.finite(q_P_ha) || q_P_ha <= 0) {
      stop("SPUR LR denominator is non-positive.")
    }
    LR[i] <- q_P_ho / q_P_ha
    pval[i] <- mean(lr_ho > LR[i])
  }

  list(
    LR = as.vector(LR),
    pvalue = as.vector(pval),
    cv_vec = as.vector(cv_vec),
    ha_parm = ha_parm
  )
}
