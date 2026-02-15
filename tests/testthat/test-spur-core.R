testthat::test_that("spurtransform returns transformed columns", {
  set.seed(123)
  n <- 40L
  df <- data.frame(
    y = as.numeric(stats::rnorm(n)),
    x = as.numeric(stats::rnorm(n)),
    lat = seq(30, 45, length.out = n),
    lon = seq(-120, -80, length.out = n)
  )

  out <- spurtransform(
    formula = ~ y + x,
    data = df,
    prefix = "h_",
    transformation = "lbmgls",
    lon = "lon",
    lat = "lat"
  )

  testthat::expect_true(all(c("h_y", "h_x") %in% names(out)))
  testthat::expect_equal(nrow(out), n)
  testthat::expect_true(all(is.finite(out$h_y)))
  testthat::expect_true(all(is.finite(out$h_x)))
})

testthat::test_that("spurtest_i0 and spurtest_i1 run on a simple fixture", {
  n <- 100L
  idx <- seq_len(n)
  lat <- seq(from = 25, to = 49, length.out = n) + sin(idx / 3) * 0.1
  lon <- seq(from = -124, to = -67, length.out = n) + cos(idx / 4) * 0.1
  x <- as.numeric(scale(sin(idx / 5) + lat * 0.02 - lon * 0.01))
  y <- 1.0 * x + 0.5 * as.numeric(scale(lat))
  df <- data.frame(y = y, x = x, lat = lat, lon = lon)

  set.seed(123)
  res_i0 <- spurtest_i0(y ~ 1, data = df, q = 8L, nrep = 500L, lon = "lon", lat = "lat")
  set.seed(123)
  res_i1 <- spurtest_i1(y ~ 1, data = df, q = 8L, nrep = 500L, lon = "lon", lat = "lat")

  for (res in list(res_i0, res_i1)) {
    testthat::expect_true(is.list(res))
    testthat::expect_true(all(c("pvalue", "teststat", "ha_param", "cv", "full") %in% names(res)))
    testthat::expect_true(is.finite(as.numeric(res$pvalue)))
    testthat::expect_true(is.finite(as.numeric(res$teststat)))
  }
})
