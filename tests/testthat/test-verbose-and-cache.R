testthat::test_that("verbose flag for spurtest_i0 only affects printing", {
  df <- spur_fixture_data(n = 140L)

  testthat::expect_output(
    {
      set.seed(42)
      silent_out <- spurtest_i0(
        formula = am ~ 1,
        data = df,
        q = 8L,
        nrep = 2000L,
        lon = "lon",
        lat = "lat"
      )
    },
    regexp = NA
  )

  testthat::expect_output(
    {
      set.seed(42)
      verbose_out <- spurtest_i0(
        formula = am ~ 1,
        data = df,
        q = 8L,
        nrep = 2000L,
        lon = "lon",
        lat = "lat",
        verbose = TRUE
      )
    },
    regexp = "Spatial I\\(0\\) Test Results"
  )

  testthat::expect_equal(silent_out$pvalue, verbose_out$pvalue, tolerance = 0)
  testthat::expect_equal(silent_out$teststat, verbose_out$teststat, tolerance = 0)
  testthat::expect_equal(silent_out$ha_param, verbose_out$ha_param, tolerance = 0)
  testthat::expect_equal(as.numeric(silent_out$cv), as.numeric(verbose_out$cv), tolerance = 0)
})

testthat::test_that("verbose flag for spurhalflife only affects printing", {
  df <- spur_fixture_data(n = 140L)

  testthat::expect_output(
    {
      set.seed(42)
      silent_out <- spurhalflife(
        formula = am ~ 1,
        data = df,
        q = 8L,
        nrep = 2000L,
        level = 90,
        normdist = TRUE,
        lon = "lon",
        lat = "lat"
      )
    },
    regexp = NA
  )

  testthat::expect_output(
    {
      set.seed(42)
      verbose_out <- spurhalflife(
        formula = am ~ 1,
        data = df,
        q = 8L,
        nrep = 2000L,
        level = 90,
        normdist = TRUE,
        lon = "lon",
        lat = "lat",
        verbose = TRUE
      )
    },
    regexp = "Spatial half life"
  )

  testthat::expect_equal(silent_out$ci_l, verbose_out$ci_l, tolerance = 0)
  testthat::expect_equal(silent_out$ci_u, verbose_out$ci_u, tolerance = 0)
  testthat::expect_equal(silent_out$max_dist, verbose_out$max_dist, tolerance = 0)
})

testthat::test_that("verbose flag for spurtransform only affects printing", {
  df <- spur_fixture_data(n = 180L)
  vars <- c("am", "fracblack")

  testthat::expect_output(
    {
      silent_out <- spurtransform(
        formula = ~ am + fracblack,
        data = df,
        prefix = "h_",
        transformation = "lbmgls",
        lon = "lon",
        lat = "lat",
        replace = TRUE,
        separately = FALSE
      )
    },
    regexp = NA
  )

  testthat::expect_output(
    {
      verbose_out <- spurtransform(
        formula = ~ am + fracblack,
        data = df,
        prefix = "h_",
        transformation = "lbmgls",
        lon = "lon",
        lat = "lat",
        replace = TRUE,
        separately = FALSE,
        verbose = TRUE
      )
    },
    regexp = "Found"
  )

  testthat::expect_equal(
    silent_out[order(silent_out$row_id), c("row_id", "h_am", "h_fracblack")],
    verbose_out[order(verbose_out$row_id), c("row_id", "h_am", "h_fracblack")],
    tolerance = 0
  )
})

testthat::test_that("cached power evaluator matches get_pow_qf", {
  set.seed(1)
  q <- 7L
  nrep <- 2000L

  A <- crossprod(matrix(stats::rnorm(q * q), q, q)) + diag(q)
  B <- crossprod(matrix(stats::rnorm(q * q), q, q)) + diag(q)
  e <- matrix(stats::rnorm(q * nrep), q, nrep)

  make_eval <- get(".make_pow_qf_evaluator", envir = asNamespace("spuR"))
  pow_qf <- get("get_pow_qf", envir = asNamespace("spuR"))
  eval_fn <- make_eval(om0 = A, e = e)

  pow_direct <- pow_qf(om0 = A, om1 = B, e = e)
  pow_cached <- eval_fn(B)

  testthat::expect_equal(pow_direct, pow_cached, tolerance = 0)
})
