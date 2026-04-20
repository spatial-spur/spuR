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
