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
  storage.mode(distmat) <- "double"
  storage.mode(emat) <- "double"

  rho_bm <- 0.999
  c_bm <- get_cbar(rho_bm, distmat)
  sigdm_bm <- get_sigma_dm(distmat, c_bm)
  R <- get_r(sigdm_bm, q)$R
  R_dm <- sweep(R, 2L, colMeans(R), "-")
  storage.mode(R_dm) <- "double"

  om_all <- .Call(cpp_om_grid, distmat, R_dm, as.numeric(c(c_grid_ho, c_grid_ha)))
  const_scale <- 0.5 * gamma(q / 2) / (pi^(q / 2))

  ch_om_ho_arr <- array(NA_real_, dim = c(q, q, n_c_ho))
  ch_omi_ho_arr <- array(NA_real_, dim = c(q, q, n_c_ho))
  const_den_mat <- numeric(n_c_ho)
  for (i in seq_len(n_c_ho)) {
    om <- om_all[, , i]
    ch_om_ho_arr[, , i] <- chol(om)
    omi <- solve(om)
    ch_omi_ho_arr[, , i] <- chol(omi)
    const_den_mat[i] <- sqrt(det(omi)) * const_scale
  }

  ch_omi_ha_arr <- array(NA_real_, dim = c(q, q, n_c_ha))
  omi_ha_arr <- array(NA_real_, dim = c(q, q, n_c_ha))
  const_den_mat_ha <- numeric(n_c_ha)
  for (i in seq_len(n_c_ha)) {
    om <- om_all[, , n_c_ho + i]
    omi <- solve(om)
    omi_ha_arr[, , i] <- omi
    ch_omi_ha_arr[, , i] <- chol(omi)
    const_den_mat_ha[i] <- sqrt(det(omi)) * const_scale
  }

  Y <- sweep(as.matrix(Y), 2L, colMeans(Y), "-")
  y_scale <- max(abs(Y))
  if (!is.finite(y_scale) || y_scale == 0) stop("Half-life inference requires a finite, nonconstant outcome.")
  X <- crossprod(R, Y / y_scale)
  storage.mode(X) <- "double"

  .Call(
    cpp_c_ci_mc,
    ch_om_ho_arr,
    ch_omi_ho_arr,
    const_den_mat,
    ch_omi_ha_arr,
    omi_ha_arr,
    const_den_mat_ha,
    emat,
    X
  )
}
