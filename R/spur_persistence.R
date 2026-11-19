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
  dims <- .validate_core_spur_inputs(distmat, emat, min_q = 2L)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("Invalid level.")
  }
  y_info <- .validate_spur_y(Z, name = "Z")
  if (ncol(y_info$Y_dm) != 1L || nrow(y_info$Y_dm) != dims$n) {
    stop("Half-life inference requires univariate `Z` of length nrow(distmat).")
  }

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
    c_grid_ha = c_grid_ha,
    y_info = y_info
  )

  # Compare the complement from the simulation count so equality at nominal
  # levels such as 90% stays excluded by the strict acceptance rule.
  non_exceedances <- dims$nrep - round(pv_mat * dims$nrep)
  ii <- as.vector(non_exceedances / dims$nrep < level)
  hl <- hl_grid_ho[ii]
  if (!length(hl)) {
    return(list(
      hl_ci_lower = NA_real_,
      hl_ci_upper = NA_real_,
      pv_mat = as.numeric(pv_mat),
      hl_grid_ho = hl_grid_ho
    ))
  }

  list(
    hl_ci_lower = min(hl),
    hl_ci_upper = max(hl),
    pv_mat = as.numeric(pv_mat),
    hl_grid_ho = hl_grid_ho
  )
}

#' @keywords internal
c_ci <- function(Y, distmat, emat, c_grid_ho, c_grid_ha, y_info = NULL) {
  dims <- .validate_core_spur_inputs(distmat, emat, min_q = 2L)
  q <- dims$q
  if (is.null(y_info)) {
    y_info <- .validate_spur_y(Y, name = "Y")
  }
  if (ncol(y_info$Y_dm) != 1L || nrow(y_info$Y_dm) != dims$n) {
    stop("Half-life inference requires univariate `Y` of length nrow(distmat).")
  }

  n_c_ho <- length(c_grid_ho)
  n_c_ha <- length(c_grid_ha)

  rho_bm <- 0.999
  c_bm <- get_cbar(rho_bm, distmat)
  sigdm_bm <- get_sigma_dm(distmat, c_bm)
  R <- .extract_spur_basis(sigdm_bm, q)
  P_mat <- .project_spur_y(y_info, R)
  X <- P_mat[, 1L, drop = FALSE]

  ch_om_mat <- vector("list", n_c_ho)
  ch_omi_mat <- vector("list", n_c_ho)
  log_det_mat <- numeric(n_c_ho)
  reference_shape <- NULL
  informative <- FALSE
  for (i in seq_len(n_c_ho)) {
    sigdm <- get_sigma_dm(distmat, c_grid_ho[i])
    om <- crossprod(R, sigdm %*% R)
    om <- 0.5 * (om + t(om))
    if (is.null(reference_shape)) reference_shape <- om
    informative <- informative || .spur_shapes_differ(om, reference_shape)
    ch_om_mat[[i]] <- chol(om)
    omi <- solve(om)
    omi <- 0.5 * (omi + t(omi))
    ch_omi <- chol(omi)
    ch_omi_mat[[i]] <- ch_omi
    log_det_mat[i] <- sum(log(diag(ch_omi)))
  }

  ch_omi_mat_ha <- vector("list", n_c_ha)
  log_det_mat_ha <- numeric(n_c_ha)
  log_ha_X <- numeric(n_c_ha)
  for (j in seq_len(n_c_ha)) {
    sigdm <- get_sigma_dm(distmat, c_grid_ha[j])
    om <- crossprod(R, sigdm %*% R)
    om <- 0.5 * (om + t(om))
    informative <- informative || .spur_shapes_differ(om, reference_shape)
    omi <- solve(om)
    omi <- 0.5 * (omi + t(omi))
    ch_omi_ha <- chol(omi)
    ch_omi_mat_ha[[j]] <- ch_omi_ha
    log_det_ha <- sum(log(diag(ch_omi_ha)))
    log_det_mat_ha[j] <- log_det_ha
    q_x <- sum((ch_omi_ha %*% X)^2)
    if (!is.finite(q_x) || q_x <= 0) {
      stop("Projected SPUR vector has zero quadratic form.")
    }
    log_ha_X[j] <- log_det_ha - 0.5 * q * log(q_x)
  }

  if (!informative) {
    stop("Half-life experiment is uninformative: covariance shapes are proportional.")
  }
  pv_mat <- matrix(NA_real_, nrow = n_c_ho, ncol = 1)
  nrep <- ncol(emat)

  for (i in seq_len(n_c_ho)) {
    ch_null <- ch_om_mat[[i]]
    e <- t(ch_null) %*% emat

    ch_omi <- ch_omi_mat[[i]]
    log_det_ho <- log_det_mat[i]
    log_ho <- log_det_ho - 0.5 * q * log(colSums((ch_omi %*% e)^2))

    q_x_ho <- sum((ch_omi %*% X)^2)
    if (!is.finite(q_x_ho) || q_x_ho <= 0) {
      stop("Projected SPUR vector has zero quadratic form.")
    }
    log_ho_X <- log_det_ho - 0.5 * q * log(q_x_ho)

    log_ha_mat <- matrix(NA_real_, nrow = nrep, ncol = n_c_ha)
    for (j in seq_len(n_c_ha)) {
      ch_omi_ha <- ch_omi_mat_ha[[j]]
      log_ha_mat[, j] <- log_det_mat_ha[j] - 0.5 * q * log(colSums((ch_omi_ha %*% e)^2))
    }

    lr <- rowMeans(exp(log_ha_mat - log_ho))
    lr_X <- mean(exp(log_ha_X - log_ho_X))
    pv_mat[i, 1] <- mean(lr > lr_X)
  }

  pv_mat
}
