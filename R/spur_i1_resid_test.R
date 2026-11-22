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
  dims <- .validate_core_spur_inputs(distmat, emat, min_q = 2L)
  q <- dims$q
  n <- dims$n
  y_mat <- if (is.vector(Y)) matrix(Y, ncol = 1L) else as.matrix(Y)
  if (nrow(y_mat) != n) {
    stop("`Y` row count must match `distmat`.")
  }
  if (nrow(as.matrix(X_in)) != n) {
    stop("`X_in` row count must match `distmat`.")
  }

  q_x <- .orthonormal_design_basis(X_in)
  y_info <- .validate_spur_residual_y(y_mat, q_x, name = "Y")
  resid_dof <- n - ncol(q_x)
  if (q > resid_dof) {
    stop(sprintf(
      "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
      q, resid_dof
    ))
  }
  M <- diag(n) - tcrossprod(q_x)

  rho_bm <- 0.999
  c_bm <- get_cbar(rho_bm, distmat)
  sigdm_bm <- get_sigma_residual(distmat, c_bm, M)

  R <- .extract_spur_basis(sigdm_bm, q)
  P_mat <- .project_spur_y(y_info, R)
  .validate_spur_i1_shapes(distmat, R, M)

  om_ho <- crossprod(R, sigdm_bm %*% R)
  om_ho <- 0.5 * (om_ho + t(om_ho))

  ha_parm <- get_ha_param_i1_residual(om_ho, distmat, R, emat, M)
  sigdm_ha <- get_sigma_residual(distmat, ha_parm, M)
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
  cv_vec <- as.numeric(stats::quantile(lr_ho, probs = 1 - sz_vec))

  n_y <- ncol(P_mat)
  LR <- numeric(n_y)
  pvalue <- numeric(n_y)
  for (i in seq_len(n_y)) {
    P <- P_mat[, i, drop = FALSE]
    q_P_ho <- sum((ch_omi_ho %*% P)^2)
    q_P_ha <- sum((ch_omi_ha %*% P)^2)
    if (!is.finite(q_P_ho) || !is.finite(q_P_ha) || q_P_ha <= 0) {
      stop("SPUR LR denominator is non-positive.")
    }
    LR[i] <- q_P_ho / q_P_ha
    pvalue[i] <- mean(lr_ho > LR[i])
  }

  list(
    LR = as.vector(LR),
    pvalue = as.vector(pvalue),
    cv_vec = cv_vec,
    ha_parm = ha_parm
  )
}
