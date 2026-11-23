.run_spur_test <- function(Y, distmat, emat, is_i0, X_in = NULL) {
  dims <- .validate_core_spur_inputs(distmat, emat, blocks = 2L)
  design <- .prepare_spur_diagnostics(Y, distmat, dims$q, X_in = X_in)
  if (is_i0) .spur_i0(design, distmat, emat) else .spur_i1(design, distmat, emat)
}

.prepare_spur_diagnostics <- function(Y, distmat, q, X_in = NULL,
                                      decay = .make_spur_decay_solver(distmat),
                                      bm = NULL, eig = NULL) {
  n <- nrow(distmat)
  if (NROW(Y) != n) {
    stop("`Y` row count must match `distmat`.")
  }

  if (is.null(X_in)) {
    y_info <- .validate_spur_y(Y)
    if (is.null(bm)) bm <- get_sigma_lbm_dm(distmat)
    selection <- .extract_spur_basis(bm, q, distmat, eig)
    basis <- selection$R
    projected <- sweep(basis, 2L, colMeans(basis), "-")
  } else {
    if (NROW(X_in) != n) {
      stop("`X_in` row count must match `distmat`.")
    }
    q_x <- .orthonormal_design_basis(X_in)
    y_info <- .validate_spur_residual_y(Y, q_x, .spur_constant_residual(q_x, X_in))
    resid_dof <- n - ncol(q_x)
    if (q > resid_dof) {
      stop(sprintf(
        "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
        q, resid_dof
      ))
    }
    # The constant direction cancels the origin terms of the LBM covariance.
    left <- .residualize_spur(-0.5 * distmat, q_x)
    bm <- t(.residualize_spur(t(left), q_x))
    bm <- 0.5 * (bm + t(bm))
    selection <- .extract_spur_basis(bm, q, distmat)
    basis <- selection$R
    projected <- .residualize_spur(basis, q_x)
  }

  # Reuse M R (or demeaned R) so each new covariance only needs R' M Sigma M R.
  om_bm <- crossprod(basis, bm %*% basis)
  force(decay)
  list(
    response = .project_spur_y(y_info, basis),
    bm = 0.5 * (om_bm + t(om_bm)),
    basis = projected,
    tied = selection$tied, all_directions = selection$all_directions,
    decay = decay
  )
}

.spur_projected_covariance <- function(distmat, basis, rate) {
  kernel <- if (is.infinite(rate)) {
    distmat == 0
  } else if (rate < 0.01) {
    expm1(-rate * distmat)
  } else {
    exp(-rate * distmat)
  }
  om <- crossprod(basis, kernel %*% basis)
  0.5 * (om + t(om))
}

.spur_precision_factor <- function(lower) {
  forwardsolve(lower, diag(nrow(lower)))
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
  q <- ncol(design$basis)
  tuning <- emat[q + seq_len(q), , drop = FALSE]
  emat <- emat[seq_len(q), , drop = FALSE]
  floor_corr <- mean(lvech(distmat) == 0)
  if (floor_corr > 0.03) {
    stop("I(0) null requires a co-located-pair fraction at or below 0.03. ", .spur_site_summary(distmat))
  }
  reference_rate <- if (floor_corr < 0.001) design$decay(0.001) else Inf
  om_ho <- .spur_projected_covariance(distmat, design$basis, reference_rate)
  if (!.spur_shapes_differ(om_ho, design$bm)) {
    stop("SPUR experiment is uninformative: covariance shapes are proportional.")
  }
  ha <- .fit_spur_i0(om_ho, om_ho, design$bm, tuning)
  om_ha <- om_ho + ha$rate * design$bm
  ch_omi_ho <- .spur_precision_factor(t(chol(om_ho)))
  ch_omi_ha <- .spur_precision_factor(t(chol(om_ha)))
  lr <- .spur_lr(design$response, ch_omi_ho, ch_omi_ha, observed = TRUE)

  first <- if (floor_corr < 0.0001) 0.0001 else floor_corr + (0.03 - floor_corr) / 30
  finite_rhos <- seq(first, 0.03, length.out = 30L)
  rho_grid <- c(floor_corr, finite_rhos[finite_rhos > floor_corr])
  pvalues <- matrix(NA_real_, length(rho_grid), length(lr))
  cvalues <- matrix(NA_real_, length(rho_grid), 3L)
  for (i in seq_along(rho_grid)) {
    rate <- if (i == 1L) Inf else design$decay(rho_grid[i])
    omega <- .spur_projected_covariance(distmat, design$basis, rate)
    lr_ho <- .spur_lr(t(chol(omega)) %*% emat, ch_omi_ho, ch_omi_ha)
    cvalues[i, ] <- .spur_rank_cutoffs(lr_ho)
    pvalues[i, ] <- vapply(lr, function(obs) (1 + sum(lr_ho >= obs)) / (length(lr_ho) + 1), 0.0)
  }

  list(
    LR = lr,
    pvalue = apply(pvalues, 2L, max),
    cvalue = apply(cvalues, 2L, max),
    ha_parm = ha$rate,
    best_tuning_power = ha$power, tied = design$tied, all_directions = design$all_directions,
    cvalue_mat = cvalues,
    pvalue_mat = pvalues,
    rho_grid = rho_grid
  )
}

.spur_i1 <- function(design, distmat, emat) {
  q <- ncol(design$basis)
  tuning <- emat[q + seq_len(q), , drop = FALSE]
  emat <- emat[seq_len(q), , drop = FALSE]
  covariance <- function(rate) .spur_projected_covariance(distmat, design$basis, rate)
  floor_corr <- mean(lvech(distmat) == 0)
  rates <- vapply(floor_corr + (1 - floor_corr) * c(0.01, 0.5, 0.95), design$decay, 0.0)
  shapes <- lapply(rates, covariance)
  if (!any(vapply(shapes, .spur_shapes_differ, logical(1L), b = design$bm))) {
    stop("SPUR experiment is uninformative: covariance shapes are proportional.")
  }

  min_distance <- min(distmat[distmat > 0])
  ha <- .fit_spur_i1(design$bm, tuning, covariance, rates, shapes, min_distance)
  lower <- t(chol(design$bm))
  ch_omi_ho <- .spur_precision_factor(lower)
  ch_omi_ha <- .spur_precision_factor(t(chol(ha$covariance)))
  lr_ho <- .spur_lr(lower %*% emat, ch_omi_ho, ch_omi_ha)
  lr <- .spur_lr(design$response, ch_omi_ho, ch_omi_ha, observed = TRUE)

  list(
    LR = lr,
    pvalue = vapply(lr, function(obs) (1 + sum(lr_ho >= obs)) / (length(lr_ho) + 1), 0.0),
    cv_vec = .spur_rank_cutoffs(lr_ho),
    ha_parm = ha$rate,
    best_tuning_power = ha$power, tied = design$tied, all_directions = design$all_directions
  )
}

.fit_spur_i1 <- function(om_ho, emat, covariance, rates, shapes, min_distance) {
  power <- .make_pow_qf_evaluator(om0 = om_ho, e = emat)

  upper <- rates[1L]
  omega <- shapes[[1L]]
  lower <- rates[3L]
  low_cov <- shapes[[3L]]
  pow <- power(low_cov)
  steps <- 1L
  while (pow > 0.5) {
    if (steps >= 200L) stop("SPUR I(1) step-1 bracketing did not converge.")
    upper <- lower
    omega <- low_cov
    lower <- lower / 2
    low_cov <- covariance(lower)
    pow <- power(low_cov)
    steps <- steps + 1L
  }

  if (pow == 0.5) return(list(rate = lower, covariance = low_cov, power = NA_real_))
  best <- list(rate = lower, covariance = low_cov, power = pow)
  pow <- power(omega)
  steps <- 1L
  if (pow > best$power) best <- list(rate = upper, covariance = omega, power = pow)
  while (pow < 0.5) {
    if (steps >= 200L || exp(-upper * min_distance) == 0) return(best)
    upper <- 2 * upper
    omega <- covariance(upper)
    pow <- power(omega)
    if (pow > best$power) best <- list(rate = upper, covariance = omega, power = pow)
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

  list(rate = rate, covariance = omega, power = NA_real_)
}

.spur_rank_cutoffs <- function(draws) {
  ranks <- (length(draws) + 1L) %/% c(100L, 20L, 10L)
  active <- ranks > 0L
  cutoffs <- rep(Inf, 3L)
  if (any(active)) {
    indices <- length(draws) - ranks[active] + 1L
    cutoffs[active] <- sort.int(draws, partial = indices)[indices]
  }
  cutoffs
}
