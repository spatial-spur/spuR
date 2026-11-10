#' Core SPUR I(0) Test Function
#'
#' This function implements the spatial I(0) test. It is called internally by \code{spurtest_i0()}.
#'
#' @param Y A numeric matrix (n x number_of_series).
#' @param distmat A numeric distance matrix (n x n).
#' @param emat A numeric simulation matrix (q x number_of_simulations).
#'
#' @return A list with elements: LR, pvalue, cvalue, ha_parm, cvalue_mat, pvalue_mat, and rho_grid.
#' @keywords internal
spur_i0_test <- function(Y, distmat, emat, sigdm_bm = NULL, R = NULL) {
  #--- Determine sizes ---
  q <- nrow(emat)
  storage.mode(distmat) <- "double"
  storage.mode(emat) <- "double"

  # 1) BM covariance matrix (approximation for demeaned value)
  if (is.null(sigdm_bm)) {
    sigdm_bm <- get_sigma_lbm_dm(distmat)
  }

  # 2) Construct R and eigenvalues for low-frequency weights
  if (is.null(R)) {
    R <- get_r(sigdm_bm, q)$R
  }
  R_dm <- sweep(R, 2L, colMeans(R), "-")
  storage.mode(R_dm) <- "double"

  # 3) Compute covariance matrices for rho = 0.001 and rho_grid
  rho <- 0.001
  rho_min <- 0.0001
  rho_max <- 0.03
  n_rho <- 30
  rho_grid <- seq(from = rho_min, to = rho_max, length.out = n_rho)

  c_all <- get_cbar(c(rho, rho_grid), distmat)
  om_all <- .Call(cpp_om_grid, distmat, R_dm, c_all)

  # 4) Matrices used in analysis
  om_rho <- om_all[, , 1L]
  om_bm <- crossprod(R, sigdm_bm %*% R)

  # 5) Find Ha parameter yielding ~50% power
  om_i0 <- om_rho
  om_ho <- om_rho
  ha_parm <- get_ha_param_i0(om_ho, om_i0, om_bm, emat)
  om_ha <- om_i0 + ha_parm * om_bm

  # 6) Compute Cholesky decompositions
  ch_omi_ho <- chol(solve(om_ho))
  ch_omi_ha <- chol(solve(om_ha))

  # 7) Get LR for data
  n_y <- ncol(Y)
  LR <- numeric(n_y)
  for(i in seq_len(n_y)) {
    X <- Y[, i] - mean(Y[, i])
    P <- crossprod(R, X)
    y_P_ho <- ch_omi_ho %*% P
    y_P_ha <- ch_omi_ha %*% P
    q_P_ho <- sum(y_P_ho^2)
    q_P_ha <- sum(y_P_ha^2)
    LR[i] <- q_P_ho / q_P_ha
  }

  # 8) Construct om_ho Cholesky factors for the grid of rho values
  ch_om_ho_mat <- array(NA_real_, dim = c(q, q, n_rho))
  for(i in seq_len(n_rho)) {
    om_ho_i <- if (rho_grid[i] > 0) om_all[, , i + 1L] else diag(q)
    ch_om_ho_mat[, , i] <- chol(om_ho_i)
  }

  # 9) Compute p-values and critical values
  sz_vec <- c(0.01, 0.05, 0.10)
  mc_out <- .Call(
    cpp_mc_lr_grid,
    ch_om_ho_mat,
    ch_omi_ho,
    ch_omi_ha,
    emat,
    LR,
    1 - sz_vec
  )
  cvalue_mat <- mc_out[[1L]]
  pvalue_mat <- mc_out[[2L]]

  # 10) Compute maximum critical values and p-values across the rho grid
  cvalue <- apply(cvalue_mat, 2, max)
  pvalue <- apply(pvalue_mat, 2, max)

  # 11) Return output as a list
  SP <- list(
    LR = LR,
    pvalue = pvalue,
    cvalue = cvalue,
    ha_parm = ha_parm,
    cvalue_mat = cvalue_mat,
    pvalue_mat = pvalue_mat,
    rho_grid = rho_grid
  )
  return(SP)
}
