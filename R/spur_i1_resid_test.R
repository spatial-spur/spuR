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
spur_i1_resid_test <- function(Y, X_in, distmat, emat, sigdm_bm = NULL, R = NULL, R_M = NULL) {
  q <- nrow(emat)
  n_y <- ncol(Y)
  storage.mode(distmat) <- "double"
  storage.mode(emat) <- "double"

  if (is.null(sigdm_bm) || is.null(R) || is.null(R_M)) {
    xtx <- crossprod(X_in)
    if (qr(xtx)$rank < ncol(xtx)) {
      stop("Design matrix is rank-deficient; residual tests require full column rank.")
    }
    storage.mode(X_in) <- "double"
    P_X <- X_in %*% solve(xtx)
    storage.mode(P_X) <- "double"

    rho_bm <- 0.999
    c_bm <- get_cbar(rho_bm, distmat)
    sigdm_bm <- .Call(cpp_get_sigma_residual_proj, distmat, as.numeric(c_bm), X_in, P_X)
    R <- get_r(sigdm_bm, q)$R
    R_M <- R - P_X %*% crossprod(X_in, R)
  }
  storage.mode(R_M) <- "double"

  om_ho <- crossprod(R, sigdm_bm %*% R)

  ha_parm <- get_ha_param_i1_residual(om_ho, distmat, R, emat, M = NULL, R_M = R_M)
  om_ha <- .Call(cpp_om_grid, distmat, R_M, as.numeric(ha_parm))

  ch_om_ho <- chol(om_ho)
  omi_ho <- solve(om_ho)
  omi_ha <- solve(om_ha)
  ch_omi_ho <- chol(omi_ho)
  ch_omi_ha <- chol(omi_ha)

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
  pvalue <- mc_out[[2L]][1L, ]

  list(
    LR = as.vector(LR),
    pvalue = as.vector(pvalue),
    cv_vec = as.vector(cv_vec),
    ha_parm = ha_parm
  )
}
