review_backends <- function() {
  c("fallback", if (is.loaded("wrap__rust_run_spur_i0", PACKAGE = "spuR")) "rust")
}

review_spur_data <- function(n = 40L) {
  set.seed(788)
  data.frame(y = rnorm(n), x = rnorm(n), sx = runif(n), sy = runif(n))
}

test_that("SPUR calibration brackets large rates and rejects unattainable averages", {
  locations <- c((0:9) * 1e-8, 1)
  distances <- abs(outer(locations, locations, "-"))
  solver <- .make_spur_decay_solver(distances)
  for (target in c(0.03, 0.001, 0.0001)) {
    rate <- solver(target)
    expect_gt(rate, 10240)
    expect_equal(mean(exp(-rate * lvech(distances))), target, tolerance = 1e-11)
    expect_identical(solver(target), rate)
  }
  with_duplicates <- abs(outer(c(0, 0, 1), c(0, 0, 1), "-"))
  for (target in c(0.3, 1 / 3)) {
    expect_error(get_cbar(target, with_duplicates), "co-located-pair fraction")
  }
  endpoint <- matrix(c(0, 1e-4, 1e-4, 0), 2L)
  expect_identical(get_cbar(exp(-10240 * 1e-4), endpoint), 10240)
})

test_that("SPUR diagnostics accept attainable correlations on co-located sites", {
  sites <- c(0, 0, 0.04, 0.17, 0.19, 0.25, 0.42, 0.60, 0.78, 1)
  data <- data.frame(s = sites, y = sin(7 * sites))
  for (backend in review_backends()) {
    i0 <- spurtest_i0(~y, data, coords_euclidean = "s", q = 2L, nrep = 199L,
                      seed = 17L, backend = backend)
    i1 <- spurtest_i1(~y, data, coords_euclidean = "s", q = 2L, nrep = 199L,
                      seed = 17L, backend = backend)
    expect_true(is.finite(i0$teststat))
    expect_true(is.finite(i1$teststat))
  }
})

test_that("residual inference agrees for equivalent large-offset group designs", {
  data <- review_spur_data(36L)
  data$g <- factor(rep(letters[1:3], each = 12L))
  data$y <- 1e13 + round(4 * data$y) / 4
  for (backend in review_backends()) {
    for (fn in list(spurtest_i0resid, spurtest_i1resid)) {
      expected <- fn(y ~ g, data, coords_euclidean = "sx", q = 4L, nrep = 499L,
                     seed = 42L, backend = backend)
      actual <- fn(y ~ 0 + g, data, coords_euclidean = "sx", q = 4L, nrep = 499L,
                   seed = 42L, backend = backend)
      expect_equal(actual$teststat, expected$teststat, tolerance = 1e-12)
      expect_identical(actual$pvalue, expected$pvalue)
    }
  }
})

test_that("NN compares raw distances before averaging ties", {
  d <- data.frame(
    s = c(0, 0.66219144446564393, -0.66219144446564404, 86012.158427686139, 43006.07921384307),
    y = seq(0, 40, by = 10)
  )
  for (backend in review_backends()) {
    res <- spurtransform(~y, d, coords_euclidean = "s", transformation = "nn", backend = backend)
    expect_equal(res$h_y[1L], -10)
  }
})

test_that("geographic ISO includes a pair at the radius", {
  d <- data.frame(
    lat = c(-1.0036422005124446, 37.11177842377785),
    lon = c(10.076763636412846, -142.9944155287363),
    y = c(1, 4)
  )
  radius <- get_distmat_lat_lon(as.matrix(d[c("lat", "lon")]))[1, 2] * .SPUR_EARTH_CIRCUMFERENCE_M
  for (backend in review_backends()) {
    transform <- function(r) {
      spurtransform(~y, d, lat = "lat", lon = "lon", transformation = "iso",
                    radius = r, backend = backend, nthreads = 1L)$h_y
    }
    expect_equal(transform(radius), c(0, 0))
    expect_equal(transform(radius * (1 + 1e-12)), c(-3, 3))
    expect_equal(transform(radius * (1 - 1e-12)), c(0, 0))
  }
})

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
    expect_equal(a$ci_l, z$ci_l * 29, tolerance = 1e-12)
    scaled <- d
    scaled$sx <- 7 * d$sx
    b <- calc(scaled, FALSE)
    expect_equal(b$max_dist, 7 * a$max_dist, tolerance = 1e-12)
    expect_equal(b$ci_l, 7 * a$ci_l, tolerance = 1e-12)
    geo <- calc(d, FALSE, TRUE)
    expect_equal(geo$max_dist, 29 * pi / 180 * 6371000.009, tolerance = 1e-10)
    expect_output(print(z), "confidence set")
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
  # These increments are exactly representable before and after adding 1e13.
  d$y <- round(256 * d$y) / 256
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
    changed$y <- d$y + 1e13
    for (fn in list(spurtest_i0resid, spurtest_i1resid)) {
      base <- fn(y ~ x, d, coords_euclidean = c("sx", "sy"), q = 8L,
                 nrep = 200L, seed = 7L, backend = backend)
      shifted <- fn(y ~ x, changed, coords_euclidean = c("sx", "sy"), q = 8L,
                    nrep = 200L, seed = 7L, backend = backend)
      expect_equal(shifted$teststat, base$teststat, tolerance = 1e-12)
      expect_equal(shifted$pvalue, base$pvalue, tolerance = 1e-12)
    }
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
    expect_equal(get_ha_param_i0(diag(2), diag(2), diag(2), matrix(1, 2, 3)), 1)
  )
  expect_equal(calls, 2L)
})
testthat::test_that("half-life inversion refines disconnected components", {
  grid <- seq(0, 1, by = 0.25)
  accept <- function(t) t <= 0.1 | (t >= 0.4 & t <= 0.6) | t >= 0.9
  actual <- .invert_half_life(grid, accept(grid), accept)
  testthat::expect_equal(unname(actual), rbind(c(0, 1/9), c(2/3, 1.5), c(9, Inf)),
                         tolerance = 2e-6)
  testthat::expect_equal(.spur_rank_cutoffs(1:99), c(99, 95, 90))
})

test_that("half-life rank bounds, complete ties and co-located transformations are reported", {
  for (backend in review_backends()) {
    set.seed(207)
    d <- data.frame(s = runif(25), y = rnorm(25))
    hl <- function(B) spurhalflife(~y, d, coords_euclidean="s", q=5, nrep=B, seed=1L, backend=backend)
    expect_warning(out <- hl(18), "necessarily covers")
    expect_equal(unname(out$components), matrix(c(0, Inf), 1, 2))
    expect_warning(hl(19), "few spatial contrasts")
    expect_error(spurhalflife(~y, d, coords_euclidean="s", q=101, backend=backend), "quadrature limit")
    t <- 2*pi*(0:19)/20
    polygon <- data.frame(s=cos(t), z=sin(t), y=cos(t)+0.4*sin(2*t))
    run <- function(q) spurtest_i0(~y, polygon, coords_euclidean=c("s","z"), q=q, nrep=199, seed=31, backend=backend)
    expect_true(run(14)$tied)
    expect_error(run(15), "splits a tied eigenspace")
    d <- data.frame(s=c(0,0,1), y=c(2,4,10))
    for (radius in c(0.5,1)) expect_equal(spurtransform(~y, d, coords_euclidean="s", transformation="iso", radius=radius, backend=backend)$h_y, c(-2,2,0))
    d <- data.frame(s=c(0,1e-6,.3,.6,1), y=c(1,-1,0,0,0), x=c(1,-1,2,0,-1))
    whiten <- function(d) spurtransform(~y+x, d, coords_euclidean="s", backend=backend)[c("h_y","h_x")]
    expect_equal(whiten(d), whiten(transform(d, y=y+1e9, x=x+1e9)), tolerance=1e-12)
  }
})
