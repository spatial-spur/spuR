#' Spatial Persistence Core Routine
#'
#' Internal translation of SPUR Mata `spatial_persistence` and `c_ci`.
#'
#' @param Z Outcome matrix.
#' @param distmat Normalized distance matrix.
#' @param emat Monte Carlo draw matrix.
#' @param level Confidence level between 0 and 1.
#' @return List with lower/upper half-life confidence bounds.
#' @keywords internal
spur_persistence <- function(Z, distmat, emat, level) {
  n_hl <- 100
  hl_grid_ho <- c(
    seq(from = 0.001, to = 1, length.out = n_hl),
    seq(from = 1.01, to = 3, length.out = 30),
    100
  )
  n_hl_ha <- 50
  hl_grid_ha <- seq(from = 0.001, to = 1, length.out = n_hl_ha)

  c_grid_ho <- -log(0.5) / hl_grid_ho
  c_grid_ha <- -log(0.5) / hl_grid_ha

  pv_mat <- c_ci(
    Y = Z,
    distmat = distmat,
    emat = emat,
    c_grid_ho = c_grid_ho,
    c_grid_ha = c_grid_ha
  )

  ii <- as.vector(pv_mat > (1 - level))
  hl <- hl_grid_ho[ii]

  list(
    hl_ci_lower = min(hl),
    hl_ci_upper = max(hl)
  )
}

#' @keywords internal
c_ci <- function(Y, distmat, emat, c_grid_ho, c_grid_ha) {
  q <- nrow(emat)
  n_c_ho <- length(c_grid_ho)
  n_c_ha <- length(c_grid_ha)

  rho_bm <- 0.999
  c_bm <- get_cbar(rho_bm, distmat)
  sigdm_bm <- get_sigma_dm(distmat, c_bm)
  R <- get_r(sigdm_bm, q)$R

  ch_om_mat <- vector("list", n_c_ho)
  ch_omi_mat <- vector("list", n_c_ho)
  const_den_mat <- numeric(n_c_ho)
  for (i in seq_len(n_c_ho)) {
    sigdm <- get_sigma_dm(distmat, c_grid_ho[i])
    om <- t(R) %*% sigdm %*% R
    ch_om_mat[[i]] <- chol(om)
    omi <- solve(om)
    ch_omi_mat[[i]] <- chol(omi)
    const_den_mat[i] <- sqrt(det(omi)) * 0.5 * gamma(q / 2) / (pi^(q / 2))
  }

  ch_omi_mat_ha <- vector("list", n_c_ha)
  const_den_mat_ha <- numeric(n_c_ha)
  for (i in seq_len(n_c_ha)) {
    sigdm <- get_sigma_dm(distmat, c_grid_ha[i])
    om <- t(R) %*% sigdm %*% R
    omi <- solve(om)
    ch_omi_mat_ha[[i]] <- chol(omi)
    const_den_mat_ha[i] <- sqrt(det(omi)) * 0.5 * gamma(q / 2) / (pi^(q / 2))
  }

  pv_mat <- matrix(NA_real_, nrow = n_c_ho, ncol = 1)
  X <- t(R) %*% Y

  for (i in seq_len(n_c_ho)) {
    ch_null <- ch_om_mat[[i]]
    e <- t(ch_null) %*% emat

    ch_omi <- ch_omi_mat[[i]]
    const_den <- const_den_mat[i]
    xc <- ch_omi %*% e
    den_ho <- const_den * (colSums(xc^2)^(-q / 2))

    Xc <- ch_omi %*% X
    den_ho_X <- as.numeric(const_den * (colSums(Xc^2)^(-q / 2)))

    den_ha_mat <- matrix(NA_real_, nrow = ncol(emat), ncol = n_c_ha)
    den_ha_mat_X <- numeric(n_c_ha)

    for (j in seq_len(n_c_ha)) {
      ch_omi_ha <- ch_omi_mat_ha[[j]]
      const_den_ha <- const_den_mat_ha[j]

      xc <- ch_omi_ha %*% e
      den_ha_mat[, j] <- const_den_ha * (colSums(xc^2)^(-q / 2))

      Xc <- ch_omi_ha %*% X
      den_ha_mat_X[j] <- as.numeric(const_den_ha * (colSums(Xc^2)^(-q / 2)))
    }

    den_ha_avg <- rowMeans(den_ha_mat)
    lr <- den_ha_avg / den_ho

    den_ha_avg_X <- mean(den_ha_mat_X)
    lr_X <- den_ha_avg_X / den_ho_X
    pv_mat[i, 1] <- mean(lr > lr_X)
  }

  pv_mat
}
