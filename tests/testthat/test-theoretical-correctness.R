testthat::test_that("half-life output preserves Euclidean units and nuisance scale", {
  df <- spur_fixture_data(n = 60L)
  args <- list(formula = am ~ 1, data = df, coords_euclidean = c("s_1", "s_2"),
               q = 6L, nrep = 300L, seed = 42L)
  normalized <- do.call(spurhalflife, c(args, list(normdist = TRUE)))
  raw <- do.call(spurhalflife, args)
  diameter <- max(as.matrix(dist(df[, c("s_1", "s_2")])))
  testthat::expect_equal(raw$max_dist, diameter, tolerance = 1e-12)
  testthat::expect_equal(raw$ci_l, normalized$ci_l * diameter, tolerance = 1e-10)
  testthat::expect_equal(raw$ci_u, normalized$ci_u * diameter, tolerance = 1e-10)
  args$data$am <- args$data$am * 1e120
  scaled <- do.call(spurhalflife, args)
  testthat::expect_equal(c(scaled$ci_l, scaled$ci_u), c(raw$ci_l, raw$ci_u), tolerance = 1e-10)
})

testthat::test_that("correlation calibration solves the equation across coordinate units", {
  D <- as.matrix(dist(cbind(1:8, (1:8)^2)))
  for (scale in c(1e-12, 1, 1e12)) {
    c0 <- spuR:::get_cbar(0.03, D * scale)
    testthat::expect_equal(mean(exp(-c0 * D[lower.tri(D)] * scale)), 0.03, tolerance = 1e-9)
  }
  testthat::expect_error(spuR:::get_cbar(0.01, matrix(0, 4, 4)), "co-located")
})

testthat::test_that("I(1) calibration cannot hang on an uninformative experiment", {
  testthat::expect_error(spuR:::.calibrate_i1_power(function(c) 0.1, 1, 10), "Cannot bracket")
  testthat::expect_error(spuR:::get_r(diag(c(3, 2, 0)), 3L), "positive rank")
})

testthat::test_that("I(0) calibration terminates on uninformative covariance shapes", {
  set.seed(724)
  E <- matrix(rnorm(3L * 1000L), 3L)
  testthat::expect_error(get_ha_param_i0(diag(3), diag(3), diag(3), E),
    "Cannot bracket 50% power")
  alternative <- diag(c(1e-4, 1, 1e4))
  g <- get_ha_param_i0(diag(3), diag(3), alternative, E)
  power <- get_pow_qf(diag(3), diag(3) + g * alternative, E)
  testthat::expect_lte(abs(power - 0.5), 0.01)
})

testthat::test_that("great-circle distances remain finite at antipodes", {
  S <- rbind(c(35.123, 45.987), c(-35.123, -134.013))
  D <- spuR:::get_distmat_lat_lon(S)
  testthat::expect_true(all(is.finite(D)))
  testthat::expect_equal(D[1, 2], 0.5, tolerance = 1e-8)
})

testthat::test_that("residual tests are invariant to regressor units", {
  df <- spur_fixture_data(n = 80L)
  a <- spurtest_i0resid(am ~ fracblack + racseg, data = df, coords_euclidean = c("s_1", "s_2"), q = 6L, nrep = 300L)
  df$fracblack <- df$fracblack * 1e8
  df$racseg <- df$racseg * 1e-8
  b <- spurtest_i0resid(am ~ fracblack + racseg, data = df, coords_euclidean = c("s_1", "s_2"), q = 6L, nrep = 300L)
  testthat::expect_equal(a$teststat, b$teststat, tolerance = 1e-8)
  testthat::expect_equal(a$pvalue, b$pvalue, tolerance = 1e-8)
})


testthat::test_that("the full pipeline refuses nonlinear transformations applied in the wrong order", {
  df <- spur_fixture_data(n = 80L)
  testthat::expect_error(spur(am ~ I(fracblack^2), data = df,
    coords_euclidean = c("s_1", "s_2")), "create nonlinear terms")
  testthat::expect_error(spurtest_i0(am ~ 1, data = df,
    coords_euclidean = c("s_1", "s_2"), q = 1L), "at least two")
})


testthat::test_that("half-life Monte Carlo uses the angular Gaussian likelihood", {
  set.seed(710)
  q <- 3L
  E <- matrix(rnorm(q * 4000L), q)
  X <- matrix(c(0.4, -0.8, 0.3), q)
  covariance <- function() {
    A <- matrix(rnorm(q*q), q)
    crossprod(A) + diag(q)
  }
  nulls <- replicate(2, covariance(), simplify = FALSE)
  alternatives <- replicate(3, covariance(), simplify = FALSE)
  array_of <- function(values) array(unlist(values), c(q, q, length(values)))
  inverse <- function(values) lapply(values, solve)
  inv_null <- inverse(nulls)
  inv_alt <- inverse(alternatives)
  constants <- function(values) vapply(values, function(om) sqrt(det(om)), 0.0)
  # Angular Gaussian density after profiling out the unknown variance.
  log_density <- function(Y, Om) {
    -0.5 * as.numeric(determinant(Om, logarithm = TRUE)$modulus) -
      q / 2 * log(colSums((chol(solve(Om)) %*% Y)^2))
  }
  log_mixture <- function(Y) {
    values <- matrix(vapply(alternatives, function(om) log_density(Y, om), numeric(ncol(Y))), nrow = ncol(Y))
    apply(values, 1, function(row) {m <- max(row); m + log(mean(exp(row - m)))})
  }
  expected <- vapply(nulls, function(om) {
    Y <- t(chol(om)) %*% E
    mean(log_mixture(Y) - log_density(Y, om) > log_mixture(X) - log_density(X, om))
  }, 0.0)
  args <- list(array_of(lapply(nulls, chol)), array_of(lapply(inv_null, chol)),
    constants(inv_null), array_of(lapply(inv_alt, chol)), array_of(inv_alt),
    constants(inv_alt), E, X)
  actual <- do.call(spuR:::rcpp_c_ci_mc, c(args, list(nthreads = 2L)))
  testthat::expect_equal(as.numeric(actual), expected, tolerance = 1e-12)
})
