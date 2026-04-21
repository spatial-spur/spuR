testthat::test_that("rewrite_formula_with_prefix prefixes a simple additive formula", {
  rewrite_formula <- getFromNamespace(".rewrite_formula_with_prefix", "spuR")
  out <- rewrite_formula(am ~ fracblack + racseg, "h_")

  testthat::expect_identical(
    format(out),
    format(h_am ~ h_fracblack + h_racseg)
  )
})

testthat::test_that("spur returns the same pipeline as manual composition", {
  df <- spur_fixture_data(n = 180L)

  out <- spur(
    am ~ fracblack + racseg,
    data = df,
    q = 8L,
    nrep = 500L,
    lon = "lon",
    lat = "lat",
    seed = 42L,
    avc = 0.03,
    uncond = FALSE,
    cvs = FALSE
  )

  i0 <- spurtest_i0(
    am ~ 1,
    data = df,
    q = 8L,
    nrep = 500L,
    lon = "lon",
    lat = "lat",
    seed = 42L
  )
  i1 <- spurtest_i1(
    am ~ 1,
    data = df,
    q = 8L,
    nrep = 500L,
    lon = "lon",
    lat = "lat",
    seed = 42L
  )
  i0r <- spurtest_i0resid(
    am ~ fracblack + racseg,
    data = df,
    q = 8L,
    nrep = 500L,
    lon = "lon",
    lat = "lat",
    seed = 42L
  )
  i1r <- spurtest_i1resid(
    am ~ fracblack + racseg,
    data = df,
    q = 8L,
    nrep = 500L,
    lon = "lon",
    lat = "lat",
    seed = 42L
  )

  levels_model <- stats::lm(am ~ fracblack + racseg, data = df)
  levels_scpc <- scpcR::scpc(
    levels_model,
    data = df,
    lon = "lon",
    lat = "lat",
    avc = 0.03,
    uncond = FALSE,
    cvs = FALSE
  )

  transformed_data <- spurtransform(
    am ~ fracblack + racseg,
    data = df,
    prefix = "h_",
    transformation = "lbmgls",
    lon = "lon",
    lat = "lat"
  )
  transformed_model <- stats::lm(h_am ~ h_fracblack + h_racseg, data = transformed_data)
  transformed_scpc <- scpcR::scpc(
    transformed_model,
    data = transformed_data,
    lon = "lon",
    lat = "lat",
    avc = 0.03,
    uncond = FALSE,
    cvs = FALSE
  )

  testthat::expect_s3_class(out, "spur_result")
  testthat::expect_true(all(c("tests", "fits") %in% names(out)))
  testthat::expect_true(all(c("i0", "i1", "i0resid", "i1resid") %in% names(out$tests)))
  testthat::expect_true(all(c("levels", "transformed") %in% names(out$fits)))

  testthat::expect_equal(out$tests$i0$teststat, i0$teststat)
  testthat::expect_equal(out$tests$i0$pvalue, i0$pvalue)
  testthat::expect_equal(out$tests$i1$teststat, i1$teststat)
  testthat::expect_equal(out$tests$i1$pvalue, i1$pvalue)
  testthat::expect_equal(out$tests$i0resid$teststat, i0r$teststat)
  testthat::expect_equal(out$tests$i0resid$pvalue, i0r$pvalue)
  testthat::expect_equal(out$tests$i1resid$teststat, i1r$teststat)
  testthat::expect_equal(out$tests$i1resid$pvalue, i1r$pvalue)

  testthat::expect_equal(
    unname(out$fits$levels$scpc$scpcstats),
    unname(levels_scpc$scpcstats),
    tolerance = 1e-8
  )
  testthat::expect_equal(out$fits$levels$scpc$q, levels_scpc$q)
  testthat::expect_equal(out$fits$levels$scpc$cv, levels_scpc$cv, tolerance = 1e-8)

  testthat::expect_equal(
    unname(out$fits$transformed$scpc$scpcstats),
    unname(transformed_scpc$scpcstats),
    tolerance = 1e-8
  )
  testthat::expect_equal(out$fits$transformed$scpc$q, transformed_scpc$q)
  testthat::expect_equal(out$fits$transformed$scpc$cv, transformed_scpc$cv, tolerance = 1e-8)
})

testthat::test_that("summary.spur_result prints the full table", {
  df <- spur_fixture_data(n = 180L)

  out <- spur(
    am ~ fracblack + racseg,
    data = df,
    q = 8L,
    nrep = 500L,
    lon = "lon",
    lat = "lat",
    seed = 42L
  )

  txt <- paste(capture.output(summary(out)), collapse = "\n")

  testthat::expect_match(txt, "SPUR Diagnostics")
  testthat::expect_match(txt, "Regression results")
  testthat::expect_match(txt, "Levels")
  testthat::expect_match(txt, "Transformed")
  testthat::expect_match(txt, "R-squared")
  testthat::expect_match(txt, "Adj\\. R-squared")
  testthat::expect_match(txt, "SCPC q")
  testthat::expect_match(txt, "SCPC cv")
  testthat::expect_match(txt, "SCPC avc")
  testthat::expect_false(grepl("ha_param", txt, fixed = TRUE))
  testthat::expect_false(grepl("h_fracblack", txt, fixed = TRUE))
  testthat::expect_false(grepl("h_racseg", txt, fixed = TRUE))
})

testthat::test_that("print.spur_result matches summary.spur_result", {
  df <- spur_fixture_data(n = 180L)

  out <- spur(
    am ~ fracblack + racseg,
    data = df,
    q = 8L,
    nrep = 500L,
    coords_euclidean = c("s_1", "s_2"),
    seed = 42L
  )

  print_txt <- paste(capture.output(print(out)), collapse = "\n")
  summary_txt <- paste(capture.output(summary(out)), collapse = "\n")

  testthat::expect_identical(print_txt, summary_txt)
})
