.half_life_mixture <- function() {
  nodes <- c(-0.906179845938664, -0.5384693101056831, 0,
             0.5384693101056831, 0.906179845938664)
  weights <- c(0.2369268850561891, 0.4786286704993665, 128 / 225,
               0.4786286704993665, 0.2369268850561891)
  edges <- c(0, 2^(-9:0))
  width <- diff(edges) / 2
  list(
    points = as.vector(outer(1 + nodes, width) + rep(head(edges, -1L), each = 5L)),
    weights = as.vector(outer(weights, width))
  )
}

.invert_half_life <- function(grid, accepted, accepts) {
  components <- matrix(numeric(), 0L, 2L)
  start <- 0
  transitions <- which(head(accepted, -1L) != tail(accepted, -1L))
  for (i in transitions) {
    left <- grid[i]
    right <- grid[i + 1L]
    while (right - left > 1e-7) {
      mid <- (left + right) / 2
      if (accepts(mid) == accepted[i]) left <- mid else right <- mid
    }
    if (accepted[i]) {
      components <- rbind(components, c(start, right))
    } else {
      start <- left
    }
  }
  if (tail(accepted, 1L)) components <- rbind(components, c(start, 1))
  colnames(components) <- c("lower", "upper")
  components / (1 - components)
}

#' @keywords internal
spur_persistence <- function(Z, distmat, emat, level) {
  dims <- .validate_core_spur_inputs(distmat, emat)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("Invalid level.")
  }
  y_info <- .validate_spur_y(Z, name = "Z")
  if (ncol(y_info$Y_dm) != 1L || nrow(y_info$Y_dm) != dims$n) {
    stop("Half-life inference requires univariate Z of length nrow(distmat).")
  }

  q <- dims$q
  nrep <- dims$nrep
  basis <- .extract_spur_basis(get_sigma_lbm_dm(distmat), q)
  projected <- sweep(basis, 2L, colMeans(basis), "-")
  observed <- .project_spur_y(y_info, basis)
  covariance <- function(h) {
    kernel <- if (h == 0) {
      distmat == 0
    } else if (is.infinite(h)) {
      -0.5 * distmat
    } else {
      rate <- log(2) / h
      if (rate < 0.01) 0.5 * expm1(-rate * distmat) / rate else exp(-rate * distmat)
    }
    omega <- crossprod(projected, kernel %*% projected)
    0.5 * (omega + t(omega))
  }

  grid <- seq(0, 1, length.out = 129L)
  h_grid <- grid / (1 - grid)
  nulls <- lapply(h_grid, covariance)
  mixture <- .half_life_mixture()
  alternatives <- lapply(mixture$points, covariance)
  if (!any(vapply(c(nulls[-1L], alternatives), .spur_shapes_differ, logical(1L), b = nulls[[1L]]))) {
    stop("Half-life experiment is uninformative: covariance shapes are proportional.")
  }
  precision_factor <- function(omega) {
    inverse <- solve(omega)
    chol(0.5 * (inverse + t(inverse)))
  }
  factors <- lapply(alternatives, precision_factor)
  log_weights <- log(mixture$weights) + vapply(factors, function(f) sum(log(diag(f))), 0.0)
  log_density <- function(upper, values, log_weight) {
    log_weight - 0.5 * q * log(colSums((upper %*% values)^2))
  }
  observed_terms <- vapply(seq_along(factors), function(j) {
    log_density(factors[[j]], observed, log_weights[j])
  }, 0.0)
  peak <- max(observed_terms)
  observed_mix <- peak + log(sum(exp(observed_terms - peak)))
  draw_norm <- 0.5 * q * log(colSums(emat^2))

  exceedances <- function(omega) {
    simulated <- t(chol(omega)) %*% emat
    upper <- precision_factor(omega)
    observed_score <- observed_mix + 0.5 * q * log(sum((upper %*% observed)^2))
    log_terms <- vapply(seq_along(factors), function(j) {
      log_density(factors[[j]], simulated, log_weights[j])
    }, numeric(nrep))
    dim(log_terms) <- c(nrep, length(factors))
    peak <- apply(log_terms, 1L, max)
    scores <- peak + log(rowSums(exp(log_terms - peak))) + draw_norm
    sum(scores >= observed_score)
  }

  counts <- vapply(nulls, exceedances, 0.0)
  accepted <- (nrep - counts) / (nrep + 1) < level
  components <- .invert_half_life(grid, accepted, function(t) {
    (nrep - exceedances(covariance(t / (1 - t)))) / (nrep + 1) < level
  })
  list(
    hl_ci_lower = if (nrow(components)) unname(components[1L, 1L]) else NA_real_,
    hl_ci_upper = if (nrow(components)) unname(components[nrow(components), 2L]) else NA_real_,
    components = components,
    pv_mat = (1 + counts) / (nrep + 1),
    hl_grid_ho = h_grid
  )
}
