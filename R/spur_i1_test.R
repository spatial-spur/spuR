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
spur_i1_test <- function(Y, distmat, emat, sigdm_bm = NULL, R = NULL) {
  q <- nrow(emat)
  storage.mode(distmat) <- "double"
  storage.mode(emat) <- "double"

  if (is.null(sigdm_bm)) {
    sigdm_bm <- get_sigma_lbm_dm(distmat)
  }
  if (is.null(R)) {
    R <- get_r(sigdm_bm, q)$R
  }
  R_dm <- sweep(R, 2L, colMeans(R), "-")
  storage.mode(R_dm) <- "double"

  om_ho <- crossprod(R, sigdm_bm %*% R)
  ha_parm <- get_ha_param_i1(om_ho, distmat, R, emat)
  om_ha <- .Call(cpp_om_grid, distmat, R_dm, as.numeric(ha_parm))

  ch_om_ho <- chol(om_ho)
  omi_ho <- solve(om_ho)
  omi_ha <- solve(om_ha)
  ch_omi_ho <- chol(omi_ho)
  ch_omi_ha <- chol(omi_ha)

  n_y <- ncol(Y)
  LR <- numeric(n_y)

  for (i in seq_len(n_y)) {
    X <- Y[, i] - mean(Y[, i])
    P <- crossprod(R, X)
    LR[i] <- (t(P) %*% omi_ho %*% P) / (t(P) %*% omi_ha %*% P)
  }

  sz_vec <- c(0.01, 0.05, 0.10)
  mc_out <- .Call(
    cpp_mc_lr_grid,
    ch_om_ho,
    ch_omi_ho,
    ch_omi_ha,
    emat,
    LR,
    1 - sz_vec
  )
  cv_vec <- mc_out[[1L]][1L, ]
  pval <- mc_out[[2L]][1L, ]

  list(
    LR = as.vector(LR),
    pvalue = as.vector(pval),
    cv_vec = as.vector(cv_vec),
    ha_parm = ha_parm
  )
}
