#' Core SPUR I(1) Test Function
#'
#' This function implements the spatial I(1) test. It is called internally by \code{spur_i1()}.
#'
#' @param Y A numeric matrix (n x number_of_series).
#' @param distmat A numeric distance matrix (n x n).
#' @param emat A numeric simulation matrix (q x number_of_simulations).
#' @param q An integer specifying the number of eigenvectors.
#'
#' @return A list with elements: LR, pvalue, cvalue, ha_parm, cvalue_mat, pvalue_mat, and rho_grid.
#' @export

spur_i1_test <- function(Y, distmat, emat) {
  ## dimensions -------------------------------------------------
  q <- nrow(emat)
  n <- nrow(distmat)

  ## Identity matrix and double‑demeaned version ----------------
  sigma     <- diag(n)
  sigma_dm  <- sigma - matrix(colMeans(sigma), n, n, byrow = TRUE)
  sigma_dm  <- sigma_dm - matrix(rowMeans(sigma_dm), n, n)

  ## Brownian‑motion covariance (demeaned) ----------------------
  sigdm_bm <- get_sigma_lbm_dm(distmat)        # << user‑supplied >>

  ## Low‑frequency eigenvectors --------------------------------
  R_lam <- Get_R(sigdm_bm, q)                  # << user‑supplied >>
  R     <- R_lam$R            # eigenvectors (cols)
  # lam <- R_lam$lam          # eigenvalues (if you need them)

  om_ho <- t(R) %*% sigdm_bm %*% R            # Ω₀
  ## Parameter under H₁ that yields ≈50 % power -----------------
  ha_parm  <- get_ha_parm_I1(om_ho, distmat, R, emat)  # << user‑supplied >>
  sigdm_ha <- get_sigma_dm(distmat, ha_parm)           # << user‑supplied >>
  om_ha    <- t(R) %*% sigdm_ha %*% R                  # Ω₁

  ## Draws of LR under H₀ --------------------------------------
  ch_om_ho   <- chol(om_ho)                     # upper‑triangular
  omi_ho     <- solve(om_ho)
  omi_ha     <- solve(om_ha)
  ch_omi_ho  <- chol(omi_ho)
  ch_omi_ha  <- chol(omi_ha)

  y_ho       <- t(ch_om_ho) %*% emat
  y_ho_ho    <- ch_omi_ho %*% y_ho
  y_ho_ha    <- ch_omi_ha %*% y_ho
  q_ho_ho    <- colSums(y_ho_ho^2)
  q_ho_ha    <- colSums(y_ho_ha^2)
  lr_ho      <- q_ho_ho / q_ho_ha               # reference LR draws

  ## critical values (1 %, 5 %, 10 %) ---------------------------
  sz_vec <- c(0.01, 0.05, 0.10)
  cv_vec <- quantile(lr_ho, probs = 1 - sz_vec, names = FALSE)

  ## LR test for each series in Y ------------------------------
  n_y   <- ncol(Y)
  LR    <- numeric(n_y)
  pval  <- numeric(n_y)

  for (i in seq_len(n_y)) {
    X <- Y[ , i] - mean(Y[ , i])               # demean series
    P <- t(R) %*% X
    LR[i]   <- ( t(P) %*% omi_ho %*% P ) /
      ( t(P) %*% omi_ha %*% P )
    pval[i] <- mean(lr_ho > LR[i])             # Monte‑Carlo p‑value
  }

  ## return -----------------------------------------------------
  list(
    LR      = as.vector(LR),
    pvalue  = as.vector(pval),
    cv_vec  = as.vector(cv_vec),
    ha_parm = ha_parm
  )
}
#---------------------------------------------------------------
