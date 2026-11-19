review_backends <- function() {
  c("fallback", if (is.loaded("wrap__rust_run_spur_i0", PACKAGE = "spuR")) "rust")
}

review_spur_data <- function(n = 40L) {
  set.seed(788)
  data.frame(y = rnorm(n), x = rnorm(n), sx = runif(n), sy = runif(n))
}

test_that("ordinary seeded SPUR calls agree and leave the same RNG state", {
  d <- review_spur_data(100L)
  if (!"rust" %in% review_backends()) skip("native SPUR is not installed")
  for (fn in list(spurtest_i0, spurtest_i1, spurtest_i0resid, spurtest_i1resid)) {
    form <- if (identical(fn, spurtest_i0resid) || identical(fn, spurtest_i1resid)) y ~ x else y ~ 1
    a <- fn(form, d, coords_euclidean = c("sx", "sy"), q = 8L, nrep = 200L,
            seed = 42L, backend = "fallback")
    state <- .Random.seed
    b <- fn(form, d, coords_euclidean = c("sx", "sy"), q = 8L, nrep = 200L,
            seed = 42L, backend = "rust", nthreads = 1L)
    expect_identical(.Random.seed, state)
    expect_identical(b$backend, "rust")
    expect_equal(b$teststat, a$teststat, tolerance = 1e-9)
    expect_equal(b$ha_param, a$ha_param, tolerance = 1e-8)
    expect_equal(b$pvalue, a$pvalue, tolerance = 1e-12)
    expect_equal(b$cv, a$cv, tolerance = 1e-9)
  }
  a <- spurhalflife(~y, d, coords_euclidean = c("sx", "sy"), q = 4L, nrep = 100L,
                    seed = 42L, backend = "fallback")
  state <- .Random.seed
  b <- spurhalflife(~y, d, coords_euclidean = c("sx", "sy"), q = 4L, nrep = 100L,
                    seed = 42L, backend = "rust", nthreads = 1L)
  expect_identical(.Random.seed, state)
  expect_equal(b$full$pv_mat, a$full$pv_mat, tolerance = 1e-12)
  expect_equal(b$ci_l, a$ci_l, tolerance = 1e-12)
  expect_equal(b$ci_u, a$ci_u, tolerance = 1e-12)
})

test_that("half-life uses Euclidean coordinate units and geographic metres", {
  d <- data.frame(sx = 0:29, sy = 0, y = sin(0:29))
  for (backend in review_backends()) {
    calc <- function(data, normdist, geographic = FALSE) {
      coords <- if (geographic) list(lon = "sx", lat = "sy") else list(coords_euclidean = c("sx", "sy"))
      do.call(spurhalflife, c(list(formula = ~y, data = data, q = 4L, nrep = 99L,
               seed = 42L, backend = backend, normdist = normdist), coords))
    }
    a <- calc(d, FALSE)
    z <- calc(d, TRUE)
    expect_equal(a$max_dist, 29, tolerance = 1e-12)
    expect_equal(a$ci_l, 0.029, tolerance = 1e-12)
    expect_equal(a$ci_l, z$ci_l * 29, tolerance = 1e-12)
    scaled <- d
    scaled$sx <- 7 * d$sx
    b <- calc(scaled, FALSE)
    expect_equal(b$max_dist, 7 * a$max_dist, tolerance = 1e-12)
    expect_equal(b$ci_l, 7 * a$ci_l, tolerance = 1e-12)
    geo <- calc(d, FALSE, TRUE)
    expect_equal(geo$max_dist, 29 / 360 * .SPUR_HALFLIFE_SCALE_M, tolerance = 1e-10)
  }
})

test_that("SPUR rejects invalid ranks and uninformative experiments", {
  d <- review_spur_data(12L)
  for (backend in review_backends()) {
    for (fn in list(spurtest_i0resid, spurtest_i1resid)) {
      expect_error(fn(y ~ x, d, coords_euclidean = c("sx", "sy"), q = 11L,
                      nrep = 50L, seed = 3L, backend = backend), "residual degrees")
    }
    for (fn in list(spurtest_i0, spurtest_i1, spurhalflife)) {
      expect_error(fn(~y, d, coords_euclidean = c("sx", "sy"), q = 1L,
                      nrep = 50L, seed = 3L, backend = backend), "at least 2")
    }
    constant <- d
    constant$y <- 1
    expect_error(spurtest_i0(~y, constant, coords_euclidean = c("sx", "sy"),
                             q = 4L, nrep = 50L, backend = backend), "variation|constant")
    dependent <- d
    dependent$x2 <- 2 * dependent$x
    expect_error(spurtest_i0resid(y ~ x + x2, dependent, coords_euclidean = c("sx", "sy"),
                                  q = 4L, nrep = 50L, backend = backend), "rank-deficient")
    eq <- data.frame(y = c(1, -2, 3), sx = c(0, 1, 0.5), sy = c(0, 0, sqrt(3) / 2))
    for (fn in list(spurtest_i0, spurtest_i1, spurhalflife)) {
      expect_error(fn(~y, eq, coords_euclidean = c("sx", "sy"), q = 2L,
                      nrep = 50L, seed = 3L, backend = backend), "uninformative")
    }
    line <- data.frame(y = c(1, -2, 3), sx = c(0, 1, 3), sy = 0)
    expect_s3_class(spurhalflife(~y, line, coords_euclidean = c("sx", "sy"),
                                  q = 2L, nrep = 50L, seed = 3L, backend = backend), "spur_halflife")
  }
})

test_that("SPUR residual projection is invariant to regressor units and response location", {
  d <- review_spur_data(100L)
  for (backend in review_backends()) {
    a <- spurtest_i0resid(y ~ x, d, coords_euclidean = c("sx", "sy"), q = 8L,
                          nrep = 200L, seed = 7L, backend = backend)
    for (scale in c(1e-12, 1e12)) {
      changed <- d
      changed$x <- scale * d$x
      b <- spurtest_i0resid(y ~ x, changed, coords_euclidean = c("sx", "sy"), q = 8L,
                            nrep = 200L, seed = 7L, backend = backend)
      expect_equal(b$teststat, a$teststat, tolerance = 1e-10)
      expect_equal(b$pvalue, a$pvalue, tolerance = 1e-12)
      expect_equal(b$ha_param, a$ha_param, tolerance = 1e-10)
    }
    # A large offset used to trip a check against the uncentered response scale.
    changed <- d
    changed$y <- d$y + 1e9
    b <- spurtest_i0resid(y ~ x, changed, coords_euclidean = c("sx", "sy"), q = 8L,
                          nrep = 200L, seed = 7L, backend = backend)
    expect_equal(b$teststat, a$teststat, tolerance = 1e-6)
    expect_equal(b$pvalue, a$pvalue, tolerance = 1e-12)
  }
})

test_that("partial fallback eigenvectors do not advance the host RNG", {
  set.seed(43)
  state <- .Random.seed
  sig <- diag(seq_len(500L))
  r <- get_r(sig, 3L)
  expect_identical(.Random.seed, state)
  expect_equal(r$DS, c(500, 499, 498), tolerance = 1e-8)
})

test_that("failed fallback power bracketing terminates explicitly", {
  calls <- 0L
  testthat::with_mocked_bindings(
    .make_pow_qf_evaluator = function(om0, e) function(om1) {
      calls <<- calls + 1L
      0
    },
    .package = "spuR",
    expect_error(get_ha_param_i0(diag(2), diag(2), diag(2), matrix(1, 2, 3)),
                 "step-2 bracketing did not converge")
  )
  expect_equal(calls, 201L)
})
testthat::test_that("half-life inversion excludes exact simulated equality at the nominal level", {
  y <- matrix(c(0, 1, 3, 2), ncol = 1L)
  distances <- as.matrix(stats::dist(0:3)) / 3
  draws <- matrix(seq_len(200), nrow = 2L)
  for (alpha in c(0.05, 0.10)) {
    local({
      pvalues <- c(alpha, alpha + 0.01, alpha - 0.01, rep(0, 128L))
      testthat::local_mocked_bindings(c_ci = function(...) pvalues, .package = "spuR")
      ci <- spur_persistence(y, distances, draws, level = 1 - alpha)
      testthat::expect_identical(ci$hl_ci_lower, ci$hl_grid_ho[2L])
      testthat::expect_identical(ci$hl_ci_upper, ci$hl_grid_ho[2L])
      testthat::expect_identical(ci$pv_mat, pvalues)
    })
  }
})
