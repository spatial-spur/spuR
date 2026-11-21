testthat::test_that("transformation validation is consistent across available backends", {
  backends <- "fallback"
  if (is.loaded("wrap__rust_run_spur_transform", PACKAGE = "spuR")) {
    backends <- c(backends, "rust")
  }
  for (backend in backends) {
    d <- data.frame(y = 4, group = "a", h_y = 99)
    one <- spurtransform(~y, d, transformation = "cluster", clustvar = "group",
                         replace = TRUE, backend = backend)
    testthat::expect_identical(one$h_y, 0)

    d$y <- NA_real_
    testthat::expect_warning(
      none <- spurtransform(~y, d, transformation = "cluster", clustvar = "group",
                            replace = TRUE, backend = backend),
      "No valid observations"
    )
    testthat::expect_identical(none$h_y, NA_real_)
    testthat::expect_error(
      spurtransform(~y, d, transformation = "cluster", clustvar = "group",
                    prefix = NA_character_, backend = backend),
      "valid name prefix"
    )
    testthat::expect_error(
      spurtransform(~y, d, transformation = "iso", radius = Inf, backend = backend),
      "positive and finite"
    )
  }
})

testthat::test_that("canonical formula APIs run across test variants", {
  df <- spur_fixture_data(n = 180L)

  set.seed(123)
  out_i0 <- spurtest_i0(
    formula = am ~ 1,
    data = df,
    q = 8L,
    nrep = 500L,
    lon = "lon",
    lat = "lat"
  )

  set.seed(123)
  out_i1 <- spurtest_i1(
    formula = ~am,
    data = df,
    q = 8L,
    nrep = 500L,
    coords_euclidean = c("s_1", "s_2")
  )

  set.seed(123)
  out_i0r <- spurtest_i0resid(
    formula = am ~ fracblack + racseg,
    data = df,
    q = 8L,
    nrep = 500L,
    lon = "lon",
    lat = "lat"
  )

  set.seed(123)
  out_i1r <- spurtest_i1resid(
    formula = am ~ fracblack + racseg,
    data = df,
    q = 8L,
    nrep = 500L,
    coords_euclidean = c("s_1", "s_2")
  )

  for (out in list(out_i0, out_i1, out_i0r, out_i1r)) {
    testthat::expect_true(is.list(out))
    testthat::expect_true(all(c("pvalue", "teststat", "ha_param", "cv", "full") %in% names(out)))
    testthat::expect_true(is.finite(as.numeric(out$pvalue)))
    testthat::expect_true(is.finite(as.numeric(out$teststat)))
  }
})

testthat::test_that("canonical half-life and transform APIs run", {
  df <- spur_fixture_data(n = 180L)

  set.seed(123)
  out_hl <- spurhalflife(
    formula = am ~ 1,
    data = df,
    q = 8L,
    nrep = 500L,
    level = 90,
    normdist = FALSE,
    lon = "lon",
    lat = "lat"
  )

  testthat::expect_true(is.list(out_hl))
  testthat::expect_true(all(c("ci_l", "ci_u", "max_dist", "full") %in% names(out_hl)))
  testthat::expect_true(is.finite(as.numeric(out_hl$ci_l)))
  testthat::expect_true(is.finite(as.numeric(out_hl$max_dist)))

  out_tr <- spurtransform(
    formula = ~ am + fracblack,
    data = df,
    prefix = "h_",
    transformation = "lbmgls",
    coords_euclidean = c("s_1", "s_2"),
    separately = FALSE,
    replace = TRUE
  )

  testthat::expect_true(all(c("h_am", "h_fracblack") %in% names(out_tr)))
  testthat::expect_equal(nrow(out_tr), nrow(df))
})

testthat::test_that("coordinate alias and mode parsing behave as expected", {
  df <- spur_fixture_data(n = 180L)

  set.seed(123)
  out_alias <- spurtest_i0(
    am ~ 1,
    data = df,
    q = 8L,
    nrep = 500L,
    coords_euclidean = c("s_1", "s_2")
  )

  testthat::expect_error(
    spurtest_i0(am ~ 1, data = df, q = 5L, nrep = 100L),
    "Specify coordinates via `lon`/`lat` or `coords_euclidean`."
  )

  testthat::expect_error(
    spurtest_i0(
      am ~ 1,
      data = df,
      q = 5L,
      nrep = 100L,
      lon = "lon",
      lat = "lat",
      coords_euclidean = c("s_1", "s_2")
    ),
    "Specify either `lon`/`lat` or `coords_euclidean`, not both."
  )

  testthat::expect_error(
    spurtest_i0(am ~ 1, data = df, q = 5L, nrep = 100L, lon = "lon"),
    "For geodesic coordinates, provide both `lon` and `lat`."
  )
})

testthat::test_that("formula shape validation remains strict", {
  df <- spur_fixture_data(n = 180L)

  testthat::expect_error(
    spurtest_i0(am ~ fracblack, data = df, q = 5L, nrep = 100L, lon = "lon", lat = "lat"),
    "single-variable formula"
  )

  testthat::expect_error(
    spurtest_i0resid(~am, data = df, q = 5L, nrep = 100L, lon = "lon", lat = "lat"),
    "must be two-sided"
  )
})

testthat::test_that("input validation for q, nrep, data, and rank deficiency", {
  df <- spur_fixture_data(n = 180L)

  testthat::expect_error(
    spurtest_i0(am ~ 1, data = df, q = 0, nrep = 100L, lon = "lon", lat = "lat"),
    "`q` must be a positive integer."
  )
  testthat::expect_error(
    spurtest_i0(am ~ 1, data = df, q = 5L, nrep = 0, lon = "lon", lat = "lat"),
    "`nrep` must be a positive integer."
  )
  testthat::expect_error(
    spurtest_i0(am ~ 1, data = as.matrix(df), q = 5L, nrep = 100L, lon = "lon", lat = "lat"),
    "`data` must be a data frame."
  )

  bad <- df
  bad$xdup <- bad$fracblack
  testthat::expect_error(
    spurtest_i0resid(am ~ fracblack + xdup, data = bad, q = 5L, nrep = 100L, lon = "lon", lat = "lat"),
    "rank-deficient"
  )
})
testthat::test_that("spurtransform transforms evaluated response and interaction columns", {
  d <- data.frame(y = exp(1:6), x = c(2, 1, NA, 3, 2, 5), z = c(3, 1, 2, 1, 2, 2),
                  g = c("a", "a", "a", "b", "b", "b"))
  backends <- "fallback"
  if (is.loaded("wrap__rust_run_spur_transform", PACKAGE = "spuR")) backends <- c(backends, "rust")
  for (backend in backends) {
    out <- spurtransform(log(y) ~ x:z, d, transformation = "cluster", clustvar = "g",
                         backend = backend)
    testthat::expect_identical(names(out), c(names(d), "h_log.y.", "h_x.z"))
    testthat::expect_equal(out$h_log.y., c(-0.5, 0.5, NA, -1, 0, 1))
    testthat::expect_equal(out$h_x.z, c(2.5, -2.5, NA, 3 - 17 / 3, 4 - 17 / 3, 10 - 17 / 3))
    separate <- spurtransform(log(y) ~ x:z, d, transformation = "cluster", clustvar = "g",
                              backend = backend, separately = TRUE)
    testthat::expect_equal(separate$h_log.y., rep(c(-1, 0, 1), 2))
    testthat::expect_equal(separate$h_x.z, out$h_x.z)
  }
})

testthat::test_that("formula LHS expressions evaluate in model.frame and spur allows compound LHS", {
  df <- utils::read.csv(testthat::test_path("fixtures", "bench_200.csv"))
  df$pos_am <- exp(df$x1 + 0.3 * df$x2)
  df$log_pos_am <- log(df$pos_am)
  df$y1 <- 0.75 * df$y
  df$y2 <- 0.25 * df$y
  df$y_sum <- df$y1 + df$y2
  backends <- "fallback"
  if (is.loaded("wrap__rust_run_spur_i0", PACKAGE = "spuR")) {
    backends <- c(backends, "rust")
  }

  for (backend in backends) {
    set.seed(10)
    res_expr_i0 <- spurtest_i0(log(pos_am) ~ 1, data = df, q = 8L, nrep = 250L,
                               coords_euclidean = c("s_1", "s_2"), backend = backend)
    set.seed(10)
    res_col_i0 <- spurtest_i0(log_pos_am ~ 1, data = df, q = 8L, nrep = 250L,
                              coords_euclidean = c("s_1", "s_2"), backend = backend)
    testthat::expect_equal(res_expr_i0$teststat, res_col_i0$teststat, tolerance = 1e-12)
    testthat::expect_equal(res_expr_i0$pvalue, res_col_i0$pvalue, tolerance = 1e-12)

    set.seed(10)
    res_expr_i1 <- spurtest_i1(log(pos_am) ~ 1, data = df, q = 8L, nrep = 250L,
                               coords_euclidean = c("s_1", "s_2"), backend = backend)
    set.seed(10)
    res_col_i1 <- spurtest_i1(log_pos_am ~ 1, data = df, q = 8L, nrep = 250L,
                              coords_euclidean = c("s_1", "s_2"), backend = backend)
    testthat::expect_equal(res_expr_i1$teststat, res_col_i1$teststat, tolerance = 1e-12)
    testthat::expect_equal(res_expr_i1$pvalue, res_col_i1$pvalue, tolerance = 1e-12)
  }

  set.seed(123)
  hl_expr <- spurhalflife(log(pos_am) ~ 1, data = df, q = 6L, nrep = 100L,
                          coords_euclidean = c("s_1", "s_2"), backend = "fallback")
  set.seed(123)
  hl_col <- spurhalflife(log_pos_am ~ 1, data = df, q = 6L, nrep = 100L,
                         coords_euclidean = c("s_1", "s_2"), backend = "fallback")
  testthat::expect_equal(hl_expr$ci_l, hl_col$ci_l, tolerance = 1e-12)
  testthat::expect_equal(hl_expr$ci_u, hl_col$ci_u, tolerance = 1e-12)

  sample <- df[seq_len(120L), ]
  pipe_expr <- spur(I(y1 + y2) ~ x1 + x2, data = sample, q = 8L, nrep = 400L,
                    lon = "s_1", lat = "s_2", seed = 42L, backend = "fallback")
  pipe_col <- spur(y_sum ~ x1 + x2, data = sample, q = 8L, nrep = 400L,
                   lon = "s_1", lat = "s_2", seed = 42L, backend = "fallback")
  testthat::expect_equal(
    unname(stats::coef(pipe_expr$fits$transformed$model)),
    unname(stats::coef(pipe_col$fits$transformed$model)),
    tolerance = 1e-12
  )
  testthat::expect_equal(pipe_expr$tests$i0$teststat, pipe_col$tests$i0$teststat, tolerance = 1e-12)
})

testthat::test_that("spurhalflife avoids overflow at high q and matches fallback vs rust", {
  set.seed(2026)
  n <- 70L
  q <- 40L
  nrep <- 100L
  df <- data.frame(
    s_1 = stats::runif(n),
    s_2 = stats::runif(n)
  )
  df$y <- sin(3 * df$s_1) + 0.4 * cos(2 * df$s_2) + 0.1 * stats::rnorm(n)
  emat <- matrix(stats::rnorm(q * nrep), nrow = q, ncol = nrep)

  fb <- testthat::with_mocked_bindings(
    .draw_emat = function(q, nrep) emat,
    .package = "spuR",
    spurhalflife(y ~ 1, data = df, q = q, nrep = nrep, level = 95,
                 coords_euclidean = c("s_1", "s_2"), backend = "fallback")
  )
  testthat::expect_true(is.finite(fb$ci_l))
  testthat::expect_true(all(is.finite(fb$full$pv_mat)))

  if (is.loaded("wrap__rust_run_spur_halflife", PACKAGE = "spuR")) {
    ru <- testthat::with_mocked_bindings(
      .draw_emat = function(q, nrep) emat,
      .package = "spuR",
      spurhalflife(y ~ 1, data = df, q = q, nrep = nrep, level = 95,
                   coords_euclidean = c("s_1", "s_2"), backend = "rust", nthreads = 1L)
    )
    testthat::expect_equal(ru$ci_l, fb$ci_l, tolerance = 1e-12)
    testthat::expect_equal(ru$ci_u, fb$ci_u, tolerance = 1e-12)
    testthat::expect_equal(ru$full$pv_mat, fb$full$pv_mat, tolerance = 1e-12)
  }
})
