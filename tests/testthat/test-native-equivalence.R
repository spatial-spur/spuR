testthat::test_that("native distance and covariance helpers match reference R definitions", {
  df <- spur_fixture_data(n = 80L)
  s_ll <- as.matrix(df[, c("lat", "lon")])
  s_eu <- as.matrix(df[, c("s_1", "s_2")])

  # Reference Euclidean distance
  d_eu_ref <- unname(as.matrix(stats::dist(s_eu, method = "euclidean")))
  d_eu_cpp <- spuR:::get_distmat_euclidean(s_eu)
  testthat::expect_equal(d_eu_cpp, d_eu_ref, tolerance = 1e-14)

  # Reference great-circle distance
  n <- nrow(s_ll)
  d_ll_ref <- matrix(0, nrow = n, ncol = n)
  c_rad <- pi / 180
  for (i in seq_len(n)) {
    d <- 0.5 * c_rad * sweep(s_ll, 2, s_ll[i, ], "-")
    d_ll_ref[, i] <- asin(sqrt(sin(d[, 1])^2 +
      cos(c_rad * s_ll[i, 1]) * (cos(c_rad * s_ll[, 1]) * sin(d[, 2])^2))) / pi
  }
  d_ll_cpp <- spuR:::get_distmat_lat_lon(s_ll)
  testthat::expect_equal(d_ll_cpp, d_ll_ref, tolerance = 1e-14)

  distmat <- d_ll_cpp / max(d_ll_cpp)

  # Reference demean_sigma and get_sigma_dm
  c_val <- 2.5
  sig_raw <- exp(-c_val * distmat)
  sig_dm_ref <- sig_raw - matrix(colMeans(sig_raw), n, n, byrow = TRUE)
  sig_dm_ref <- sig_dm_ref - matrix(rowMeans(sig_dm_ref), n, n)
  testthat::expect_equal(spuR:::demean_sigma(sig_raw), sig_dm_ref, tolerance = 1e-14)
  testthat::expect_equal(spuR:::get_sigma_dm(distmat, c_val), sig_dm_ref, tolerance = 1e-14)

  # Reference get_sigma_lbm and get_sigma_lbm_dm
  first_col <- matrix(distmat[, 1], n, n)
  first_row <- matrix(distmat[1, ], n, n, byrow = TRUE)
  lbm_ref <- 0.5 * (first_col + first_row - distmat)
  lbm_dm_ref <- lbm_ref - matrix(colMeans(lbm_ref), n, n, byrow = TRUE)
  lbm_dm_ref <- lbm_dm_ref - matrix(rowMeans(lbm_dm_ref), n, n)
  testthat::expect_equal(spuR:::get_sigma_lbm(distmat), lbm_ref, tolerance = 1e-14)
  testthat::expect_equal(spuR:::get_sigma_lbm_dm(distmat), lbm_dm_ref, tolerance = 1e-14)
})

testthat::test_that("native get_cbar and transformation matrices match reference R definitions", {
  df <- spur_fixture_data(n = 75L)
  s_ll <- as.matrix(df[, c("lat", "lon")])
  distmat <- spuR:::get_distmat_lat_lon(s_ll)
  distmat_norm <- distmat / max(distmat)
  vd <- distmat_norm[lower.tri(distmat_norm)]

  ref_get_cbar <- function(rhobar) {
    c0 <- 10
    c1 <- 10
    while (mean(exp(-c0 * vd)) <= rhobar) {
      c1 <- c0
      c0 <- c0 / 2
    }
    while (mean(exp(-c1 * vd)) >= rhobar && c1 <= 10000) {
      c0 <- c1
      c1 <- 2 * c1
    }
    while ((c1 - c0) > 0.001) {
      cm <- sqrt(c0 * c1)
      if (mean(exp(-cm * vd)) < rhobar) c1 <- cm else c0 <- cm
    }
    sqrt(c0 * c1)
  }

  rhos <- c(0.0001, 0.001, 0.01, 0.03, 0.95, 0.999)
  cb_ref <- vapply(rhos, ref_get_cbar, numeric(1))
  cb_cpp <- spuR:::get_cbar(rhos, distmat_norm)
  testthat::expect_equal(cb_cpp, cb_ref, tolerance = 1e-14)

  # Transformation matrices
  n <- nrow(s_ll)
  d_noz <- distmat_norm + (distmat_norm == 0) * 1e10
  rmins <- apply(d_noz, 1L, min)
  nn_ref <- d_noz == matrix(rmins, n, n)
  nn_ref <- diag(n) - nn_ref / rowSums(nn_ref)
  testthat::expect_equal(spuR:::nn_matrix(s_ll, latlong = TRUE), nn_ref, tolerance = 1e-14)

  radius <- 150000
  d_m <- distmat * 3.14159265359 * 6371000.009 * 2
  d_m_noz <- d_m + (d_m == 0) * 1e10
  below <- d_m_noz <= radius
  below <- below / rowSums(below)
  iso_ref <- diag(n) - below
  iso_ref[is.na(iso_ref)] <- 0
  testthat::expect_equal(spuR:::iso_matrix(s_ll, radius = radius, latlong = TRUE), iso_ref, tolerance = 1e-14)

  cl <- as.numeric(factor(df$state))
  cl_ref <- outer(cl, cl, "==")
  cl_ref <- diag(n) - cl_ref / rowSums(cl_ref)
  testthat::expect_equal(spuR:::cluster_matrix(cl), cl_ref, tolerance = 1e-14)
})
