#' Core SPUR I(0) Test Function
#'
#' This function implements the spatial I(0) test. It is called internally by \code{spur_i0()}.
#'
#' @param Y A numeric matrix (n x number_of_series).
#' @param distmat A numeric distance matrix (n x n).
#' @param emat A numeric simulation matrix (q x number_of_simulations).
#' @param q An integer specifying the number of eigenvectors.
#'
#' @return A list with elements: LR, pvalue, cvalue, ha_parm, cvalue_mat, pvalue_mat, and rho_grid.
#' @export
spur_i0_test <- function(Y, distmat, emat, q) {
  #--- Determine sizes ---
  n <- nrow(distmat)
  
  # 1) BM covariance matrix (approximation for demeaned value)
  sigdm_bm <- get_sigma_lbm_dm(distmat)
  
  # 2) Construct R and eigenvalues for low-frequency weights
  R_out <- Get_R(sigdm_bm, q)
  R <- R_out$R
  # (We assume the eigenvalues might be stored as lam if needed.)
  
  # 3) Compute covariance matrix for rho = 0.001
  rho <- 0.001
  c_val <- getcbar(rho, distmat)
  sigdm_rho <- get_sigma_dm(distmat, c_val)
  
  sigma <- diag(n)
  # Demean sigma: subtract column and row means
  sigma_dm <- sweep(sigma, 2, colMeans(sigma), "-")
  sigma_dm <- sweep(sigma_dm, 1, rowMeans(sigma_dm), "-")
  sigdm_wn <- sigma_dm
  
  # 4) Matrices used in analysis
  om_rho <- t(R) %*% sigdm_rho %*% R
  om_wn  <- t(R) %*% sigdm_wn  %*% R
  om_bm  <- t(R) %*% sigdm_bm %*% R
  
  # 5) Find Ha parameter yielding ~50% power
  om_i0   <- om_rho
  om_ho   <- om_rho
  ha_parm <- get_ha_parm_I0(om_ho, om_i0, om_bm, emat)
  om_ha   <- om_i0 + ha_parm * om_bm
  
  # 6) Compute Cholesky decompositions
  ch_omi_ho <- chol(solve(om_ho))
  ch_omi_ha <- chol(solve(om_ha))
  
  # 7) Get LR for data
  n_y <- ncol(Y)
  LR <- numeric(n_y)
  for(i in seq_len(n_y)) {
    X <- Y[, i] - mean(Y[, i])
    P <- t(R) %*% X
    y_P_ho <- ch_omi_ho %*% P
    y_P_ha <- ch_omi_ha %*% P
    q_P_ho <- sum(y_P_ho^2)
    q_P_ha <- sum(y_P_ha^2)
    LR[i] <- q_P_ho / q_P_ha
  }
  
  # 8) Construct om_ho for a grid of rho values
  rho_min <- 0.0001
  rho_max <- 0.03
  n_rho <- 30
  rho_grid <- seq(from = rho_min, to = rho_max, length.out = n_rho)
  ch_om_ho_mat <- array(NA_real_, dim = c(q, q, n_rho))
  
  for(i in seq_len(n_rho)) {
    om_ho_i <- diag(q)
    rho_i <- rho_grid[i]
    if (rho_i > 0) {
      c_i <- getcbar(rho_i, distmat)
      sigdm_ho_i <- get_sigma_dm(distmat, c_i)
      om_ho_i <- t(R) %*% sigdm_ho_i %*% R
    }
    ch_om_ho_mat[, , i] <- chol(om_ho_i)
  }
  
  # 9) Compute p-values and critical values
  pvalue_mat <- matrix(NA_real_, nrow = n_rho, ncol = n_y)
  sz_vec <- c(0.01, 0.05, 0.10)
  cvalue_mat <- matrix(NA_real_, nrow = n_rho, ncol = length(sz_vec))
  
  for(ir in seq_len(n_rho)) {
    ch_om_ho <- ch_om_ho_mat[, , ir]
    y_ho <- ch_om_ho %*% emat
    y_ho_ho <- ch_omi_ho %*% y_ho
    y_ho_ha <- ch_omi_ha %*% y_ho
    q_ho_ho <- colSums(y_ho_ho^2)
    q_ho_ha <- colSums(y_ho_ha^2)
    lr_ho <- q_ho_ho / q_ho_ha
    cv_vec <- quantile(lr_ho, probs = 1 - sz_vec)
    cvalue_mat[ir, ] <- cv_vec
    for(j in seq_len(n_y)) {
      pvalue_mat[ir, j] <- mean(lr_ho > LR[j])
    }
  }
  
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
