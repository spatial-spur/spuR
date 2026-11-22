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
testthat::test_that("spur transforms the evaluated model on one sample and consumes one draw matrix", {
  df <- spur_fixture_data(n = 75L)
  df$positive <- df$am - min(df$am) + 1
  df$positive[c(1L, 4L)] <- NA_real_
  df$fracblack[9L] <- NA_real_
  df$lon[c(1L, 9L)] <- NA_real_
  f <- log(positive) ~ fracblack * racseg
  q <- 8L
  nrep <- 500L
  set.seed(42)
  .draw_emat(q, nrep)
  expected_rng <- .Random.seed

  out <- spur(f, df, q, nrep, lon = "lon", lat = "lat", seed = 42,
              avc = 0.2, uncond = TRUE, backend = "fallback")
  testthat::expect_identical(.Random.seed, expected_rng)
  original <- stats::lm(f, df)
  sample <- match(rownames(model.frame(original)), rownames(df))
  transformed <- out$fits$transformed$model
  testthat::expect_identical(rownames(model.frame(transformed)), rownames(model.frame(original)))
  complete <- df[sample, , drop = FALSE]
  complete$testing_response <- as.numeric(model.response(model.frame(original)))
  expected <- spurtransform(~ testing_response, complete, transformation = "lbmgls",
                            lon = "lon", lat = "lat", backend = "fallback")
  testthat::expect_equal(
    as.numeric(model.response(model.frame(transformed))),
    expected$h_testing_response, tolerance = 1e-10
  )
  diagnostics <- list(i0 = spurtest_i0, i1 = spurtest_i1,
                      i0resid = spurtest_i0resid, i1resid = spurtest_i1resid)
  for (kind in names(diagnostics)) {
    formula <- if (endsWith(kind, "resid")) testing_response ~ fracblack * racseg else ~ testing_response
    individual <- diagnostics[[kind]](formula, complete, q, nrep, lon = "lon", lat = "lat",
                                       seed = 42L, backend = "fallback")
    testthat::expect_equal(out$tests[[kind]]$full, individual$full, tolerance = 1e-12)
  }
  testthat::expect_true("fracblack:racseg" %in% out$fits$transformed$term_names)
  testthat::expect_identical(unname(out$stage_backends), rep("fallback", 7L))
  if (is.loaded("wrap__rust_run_spur_pipeline", PACKAGE = "spuR")) {
    native <- spur(f, df, q, nrep, lon = "lon", lat = "lat", seed = 42,
                   avc = 0.2, uncond = TRUE, backend = "rust", nthreads = 1L)
    testthat::expect_identical(.Random.seed, expected_rng)
    testthat::expect_identical(unname(native$stage_backends), rep("rust", 7L))
    testthat::expect_equal(
      coef(native$fits$transformed$model), coef(transformed), tolerance = 1e-9
    )
  }
})
