.half_life_mixture <- function(min_distance) {
  nodes <- c(-0.906179845938664, -0.5384693101056831, 0,
             0.5384693101056831, 0.906179845938664)
  weights <- c(0.2369268850561891, 0.4786286704993665, 128 / 225,
               0.4786286704993665, 0.2369268850561891)
  near_zero <- max(min_distance / 32, .Machine$double.xmin)
  exponent <- -9L
  while (2^exponent > near_zero) exponent <- exponent - 1L
  edges <- c(0, 2^(exponent:0))
  width <- diff(edges) / 2
  list(
    points = as.vector(outer(1 + nodes, width) + rep(utils::head(edges, -1L), each = 5L)),
    weights = as.vector(outer(weights, width))
  )
}

.half_life_search_grid <- function(min_distance) {
  short <- numeric()
  h <- 1 / 127
  near_zero <- max(min_distance / 32, .Machine$double.xmin)
  while (h > near_zero) {
    h <- h * 0.8408964152537145
    short <- c(short, h / (1 + h))
  }
  c(0, rev(short), seq_len(128L) / 128)
}

.invert_half_life <- function(grid, accepted, accepts) {
  components <- matrix(numeric(), 0L, 2L)
  start <- 0
  transitions <- which(utils::head(accepted, -1L) != utils::tail(accepted, -1L))
  for (i in transitions) {
    left <- grid[i]
    right <- grid[i + 1L]
    tolerance <- min(1e-7, 1e-5 * (right - left))
    while (right - left > tolerance) {
      mid <- (left + right) / 2
      if (accepts(mid) == accepted[i]) left <- mid else right <- mid
    }
    if (accepted[i]) {
      components <- rbind(components, c(start, right))
    } else {
      start <- left
    }
  }
  if (utils::tail(accepted, 1L)) components <- rbind(components, c(start, 1))
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
  if (q > 100L) stop(sprintf("q=%d exceeds the validated half-life quadrature limit of 100.", q))
  selection <- .extract_spur_basis(get_sigma_lbm_dm(distmat), q, distmat)
  basis <- selection$R
  projected <- sweep(basis, 2L, colMeans(basis), "-")
  observed <- .project_spur_y(y_info, basis)
  covariance <- function(h) {
    kernel <- if (h == 0) {
      distmat == 0
    } else if (is.infinite(h)) {
      -0.5 * distmat
    } else {
      rate <- log(2) / h
      if (rate < 0.01) 0.5 * expm1(-rate * distmat) / rate else exp(-log(2) * (distmat / h))
    }
    omega <- crossprod(projected, kernel %*% projected)
    0.5 * (omega + t(omega))
  }

  min_distance <- min(distmat[distmat > 0])
  grid <- .half_life_search_grid(min_distance)
  h_grid <- grid / (1 - grid)
  nulls <- lapply(h_grid, covariance)
  mixture <- .half_life_mixture(min_distance)
  alternatives <- lapply(mixture$points, covariance)
  if (!any(vapply(c(nulls[-1L], alternatives), .spur_shapes_differ, logical(1L), b = nulls[[1L]]))) {
    stop("Half-life experiment is uninformative: covariance shapes are proportional.")
  }
  factors <- lapply(alternatives, function(omega) .spur_precision_factor(t(chol(omega))))
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
    lower <- t(chol(omega))
    simulated <- lower %*% emat
    observed_score <- observed_mix + 0.5 * q * log(sum(forwardsolve(lower, observed)^2))
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
    tied = selection$tied, all_directions = selection$all_directions,
    pv_mat = (1 + counts) / (nrep + 1),
    hl_grid_ho = h_grid
  )
}
