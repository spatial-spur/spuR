.run_spur_test <- function(Y, distmat, emat, is_i0, X_in = NULL) {
  dims <- .validate_core_spur_inputs(distmat, emat)
  design <- .prepare_spur_diagnostics(Y, distmat, dims$q, X_in = X_in)
  if (is_i0) .spur_i0(design, distmat, emat) else .spur_i1(design, distmat, emat)
}

.prepare_spur_diagnostics <- function(Y, distmat, q, X_in = NULL,
                                      decay = .make_spur_decay_solver(distmat),
                                      bm = NULL, eig = NULL, q_x = NULL) {
  n <- nrow(distmat)
  if (NROW(Y) != n) {
    stop("`Y` row count must match `distmat`.")
  }

  if (is.null(X_in)) {
    y_info <- .validate_spur_y(Y)
    if (is.null(bm)) bm <- get_sigma_lbm_dm(distmat)
    basis <- .extract_spur_basis(bm, q, eig)
    projected <- sweep(basis, 2L, colMeans(basis), "-")
  } else {
    if (NROW(X_in) != n) {
      stop("`X_in` row count must match `distmat`.")
    }
    if (is.null(q_x)) q_x <- .orthonormal_design_basis(X_in)
    y_info <- .validate_spur_residual_y(Y, q_x)
    resid_dof <- n - ncol(q_x)
    if (q > resid_dof) {
      stop(sprintf(
        "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
        q, resid_dof
      ))
    }
    # Apply M = I - Q_x Q_x' on both sides of the kernel without building M.
    sigma <- exp(-decay(0.999) * distmat)
    left <- .residualize_spur(sigma, q_x)
    bm <- t(.residualize_spur(t(left), q_x))
    bm <- 0.5 * (bm + t(bm))
    basis <- .extract_spur_basis(bm, q)
    projected <- .residualize_spur(basis, q_x)
  }

  # Reuse M R (or demeaned R) so each new covariance only needs R' M Sigma M R.
  om_bm <- crossprod(basis, bm %*% basis)
  force(decay)
  list(
    response = .project_spur_y(y_info, basis),
    bm = 0.5 * (om_bm + t(om_bm)),
    basis = projected,
    decay = decay
  )
}

.spur_projected_covariance <- function(distmat, basis, rate) {
  om <- crossprod(basis, exp(-rate * distmat) %*% basis)
  0.5 * (om + t(om))
}

.spur_lr <- function(values, ch_omi_ho, ch_omi_ha, observed = FALSE) {
  numerator <- colSums((ch_omi_ho %*% values)^2)
  denominator <- colSums((ch_omi_ha %*% values)^2)
  if (observed && any(!is.finite(numerator) | !is.finite(denominator) | denominator <= 0)) {
    stop("SPUR LR denominator is non-positive.")
  }
  numerator / denominator
}

.spur_i0 <- function(design, distmat, emat) {
  om_ho <- .spur_projected_covariance(distmat, design$basis, design$decay(0.001))
  if (!.spur_shapes_differ(om_ho, design$bm)) {
    stop("SPUR experiment is uninformative: covariance shapes are proportional.")
  }
  ha_parm <- get_ha_param_i0(om_ho, om_ho, design$bm, emat)
  om_ha <- om_ho + ha_parm * design$bm
  ch_omi_ho <- chol(solve(om_ho))
  ch_omi_ha <- chol(solve(om_ha))
  lr <- .spur_lr(design$response, ch_omi_ho, ch_omi_ha, observed = TRUE)

  rho_grid <- seq(0.0001, 0.03, length.out = 30L)
  probs <- 1 - c(0.01, 0.05, 0.10)
  pvalues <- matrix(NA_real_, length(rho_grid), length(lr))
  cvalues <- matrix(NA_real_, length(rho_grid), length(probs))
  for (i in seq_along(rho_grid)) {
    omega <- .spur_projected_covariance(distmat, design$basis, design$decay(rho_grid[i]))
    lr_ho <- .spur_lr(t(chol(omega)) %*% emat, ch_omi_ho, ch_omi_ha)
    cvalues[i, ] <- stats::quantile(lr_ho, probs, names = FALSE)
    pvalues[i, ] <- vapply(lr, function(obs) mean(lr_ho > obs), 0.0)
  }

  list(
    LR = lr,
    pvalue = apply(pvalues, 2L, max),
    cvalue = apply(cvalues, 2L, max),
    ha_parm = ha_parm,
    cvalue_mat = cvalues,
    pvalue_mat = pvalues,
    rho_grid = rho_grid
  )
}

.spur_i1 <- function(design, distmat, emat) {
  covariance <- function(rate) .spur_projected_covariance(distmat, design$basis, rate)
  rates <- vapply(c(0.01, 0.5, 0.95), design$decay, 0.0)
  shapes <- lapply(rates, covariance)
  if (!.spur_shapes_differ(shapes[[1L]], shapes[[2L]]) &&
      !.spur_shapes_differ(shapes[[1L]], shapes[[3L]])) {
    stop("SPUR experiment is uninformative: covariance shapes are proportional.")
  }

  ha <- .fit_spur_i1(design$bm, emat, covariance, rates, shapes)
  ch_omi_ho <- chol(solve(design$bm))
  ch_omi_ha <- chol(solve(ha$covariance))
  lr_ho <- .spur_lr(t(chol(design$bm)) %*% emat, ch_omi_ho, ch_omi_ha)
  lr <- .spur_lr(design$response, ch_omi_ho, ch_omi_ha, observed = TRUE)

  list(
    LR = lr,
    pvalue = vapply(lr, function(obs) mean(lr_ho > obs), 0.0),
    cv_vec = stats::quantile(lr_ho, 1 - c(0.01, 0.05, 0.10), names = FALSE),
    ha_parm = ha$rate
  )
}

.fit_spur_i1 <- function(om_ho, emat, covariance, rates, shapes) {
  power <- .make_pow_qf_evaluator(om0 = om_ho, e = emat)

  lower <- rates[3L]
  pow <- power(shapes[[3L]])
  steps <- 1L
  while (pow > 0.5) {
    if (steps >= 200L) stop("SPUR I(1) step-1 bracketing did not converge.")
    lower <- lower / 2
    pow <- power(covariance(lower))
    steps <- steps + 1L
  }

  upper <- rates[1L]
  omega <- shapes[[1L]]
  pow <- power(omega)
  steps <- 1L
  while (pow < 0.5) {
    if (steps >= 200L) stop("SPUR I(1) step-2 bracketing did not converge.")
    upper <- 2 * upper
    omega <- covariance(upper)
    pow <- power(omega)
    steps <- steps + 1L
  }

  rate <- upper
  for (i in seq_len(20L)) {
    if (abs(pow - 0.5) <= 0.01) break
    rate <- (lower + upper) / 2
    omega <- covariance(rate)
    pow <- power(omega)
    if (pow > 0.5) {
      upper <- rate
    } else {
      lower <- rate
    }
  }

  list(rate = rate, covariance = omega)
}
