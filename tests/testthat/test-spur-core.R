testthat::test_that("spur_transform returns transformed columns", {
  set.seed(123)
  n <- 40L
  df <- data.frame(
    y = as.numeric(stats::rnorm(n)),
    x = as.numeric(stats::rnorm(n)),
    lat = seq(30, 45, length.out = n),
    lon = seq(-120, -80, length.out = n)
  )

  out <- spur_transform(
    data = df,
    vars = c("y", "x"),
    prefix = "h_",
    transformation = "lbmgls",
    latlong = TRUE
  )

  testthat::expect_true(all(c("h_y", "h_x") %in% names(out)))
  testthat::expect_equal(nrow(out), n)
  testthat::expect_true(all(is.finite(out$h_y)))
  testthat::expect_true(all(is.finite(out$h_x)))
})

testthat::test_that("spur_i0 and spur_i1 run on a simple fixture", {
  n <- 100L
  idx <- seq_len(n)
  lat <- seq(from = 25, to = 49, length.out = n) + sin(idx / 3) * 0.1
  lon <- seq(from = -124, to = -67, length.out = n) + cos(idx / 4) * 0.1
  x <- as.numeric(scale(sin(idx / 5) + lat * 0.02 - lon * 0.01))
  y <- 1.0 * x + 0.5 * as.numeric(scale(lat))
  df <- data.frame(y = y, x = x, lat = lat, lon = lon)

  set.seed(123)
  res_i0 <- spur_i0(var = "y", q = 8L, nrep = 500L, latlong = TRUE, data = df)
  set.seed(123)
  res_i1 <- spur_i1(var = "y", q = 8L, nrep = 500L, latlong = TRUE, data = df)

  for (res in list(res_i0, res_i1)) {
    testthat::expect_true(is.list(res))
    testthat::expect_true(all(c("pvalue", "teststat", "ha_param", "cv", "full") %in% names(res)))
    testthat::expect_true(is.finite(as.numeric(res$pvalue)))
    testthat::expect_true(is.finite(as.numeric(res$teststat)))
  }
})
