make_det_emat_r <- function(q, nrep) {
  out <- numeric(q * nrep)
  for (i in seq_len(length(out))) {
    out[i] <- sin(i * 1.7320508075688772) + 0.5 * cos(i * 0.6180339887498949)
  }
  matrix(out, nrow = q, ncol = nrep, byrow = TRUE)
}

read_bench_200 <- function() {
  utils::read.csv(testthat::test_path("fixtures", "bench_200.csv"))
}

testthat::test_that("spuR diagnostics match fallback with frozen emat (geodesic and euclidean)", {
  testthat::skip_if_not(is.loaded("wrap__rust_run_spur_i0", PACKAGE = "spuR"), "spuR native Rust library not loaded")
  df <- read_bench_200()
  q <- 8L
  nrep <- 1000L
  set.seed(42)
  base_emat <- matrix(stats::rnorm(q * nrep), nrow = q, ncol = nrep)
  x_in <- cbind(1, df$x1, df$x2)

  configs <- list(
    list(latlong = TRUE, lon = "s_1", lat = "s_2", euc = NULL, coords = as.matrix(df[, c("s_2", "s_1")])),
    list(latlong = FALSE, lon = NULL, lat = NULL, euc = c("s_1", "s_2"), coords = as.matrix(df[, c("s_1", "s_2")]))
  )

  for (cfg in configs) {
    e_lbm <- list(fb = base_emat, ru = base_emat)
    e_res <- list(fb = base_emat, ru = base_emat)

    for (pair in list(
      list(name = "i0", fn = spuR::spurtest_i0),
      list(name = "i1", fn = spuR::spurtest_i1)
    )) {
      fb <- testthat::with_mocked_bindings(
        .draw_emat = function(q, nrep) e_lbm$fb,
        .package = "spuR",
        pair$fn(
          y ~ 1,
          data = df,
          q = q,
          nrep = nrep,
          lon = cfg$lon,
          lat = cfg$lat,
          coords_euclidean = cfg$euc,
          seed = 7L,
          backend = "fallback"
        )
      )
      ru <- testthat::with_mocked_bindings(
        .draw_emat = function(q, nrep) e_lbm$ru,
        .package = "spuR",
        pair$fn(
          y ~ 1,
          data = df,
          q = q,
          nrep = nrep,
          lon = cfg$lon,
          lat = cfg$lat,
          coords_euclidean = cfg$euc,
          seed = 7L,
          backend = "rust",
          nthreads = 2L
        )
      )
      testthat::expect_identical(fb$backend, "fallback")
      testthat::expect_identical(ru$backend, "rust")
      testthat::expect_false(ru$tied)
      testthat::expect_equal(ru$teststat, fb$teststat, tolerance = 1e-10)
      testthat::expect_equal(ru$ha_param, fb$ha_param, tolerance = 1e-9)
      testthat::expect_equal(ru$cv, fb$cv, tolerance = 1e-9)
      testthat::expect_equal(ru$pvalue, fb$pvalue, tolerance = 1e-12)
      if (identical(pair$name, "i0")) {
        testthat::expect_length(ru$full$rho_grid, 30L)
        testthat::expect_equal(ru$full$pvalue_mat, fb$full$pvalue_mat, tolerance = 1e-12)
        testthat::expect_equal(ru$full$cvalue_mat, fb$full$cvalue_mat, tolerance = 1e-9)
        testthat::expect_equal(ru$pvalue, max(ru$full$pvalue_mat), tolerance = 1e-15)
        testthat::expect_equal(ru$cv, apply(ru$full$cvalue_mat, 2, max), tolerance = 1e-15)
      }
    }

    for (pair in list(
      list(name = "i0resid", fn = spuR::spurtest_i0resid),
      list(name = "i1resid", fn = spuR::spurtest_i1resid)
    )) {
      fb <- testthat::with_mocked_bindings(
        .draw_emat = function(q, nrep) e_res$fb,
        .package = "spuR",
        pair$fn(
          y ~ x1 + x2,
          data = df,
          q = q,
          nrep = nrep,
          lon = cfg$lon,
          lat = cfg$lat,
          coords_euclidean = cfg$euc,
          seed = 7L,
          backend = "fallback"
        )
      )
      ru <- testthat::with_mocked_bindings(
        .draw_emat = function(q, nrep) e_res$ru,
        .package = "spuR",
        pair$fn(
          y ~ x1 + x2,
          data = df,
          q = q,
          nrep = nrep,
          lon = cfg$lon,
          lat = cfg$lat,
          coords_euclidean = cfg$euc,
          seed = 7L,
          backend = "rust",
          nthreads = 2L
        )
      )
      testthat::expect_identical(ru$backend, "rust")
      testthat::expect_false(ru$tied)
      testthat::expect_equal(ru$teststat, fb$teststat, tolerance = 1e-9)
      testthat::expect_equal(ru$ha_param, fb$ha_param, tolerance = 1e-8)
      testthat::expect_equal(ru$cv, fb$cv, tolerance = 1e-9)
      testthat::expect_equal(ru$pvalue, fb$pvalue, tolerance = 1e-12)
    }
  }
})

testthat::test_that("spurhalflife matches fallback with frozen emat and preserves R scaling", {
  testthat::skip_if_not(is.loaded("wrap__rust_run_spur_halflife", PACKAGE = "spuR"), "spuR native Rust library not loaded")
  df <- read_bench_200()[seq_len(120L), ]
  q <- 6L
  nrep <- 200L
  set.seed(77)
  base_emat <- matrix(stats::rnorm(q * nrep), nrow = q, ncol = nrep)

  configs <- list(
    list(latlong = TRUE, lon = "s_1", lat = "s_2", euc = NULL, coords = as.matrix(df[, c("s_2", "s_1")])),
    list(latlong = FALSE, lon = NULL, lat = NULL, euc = c("s_1", "s_2"), coords = as.matrix(df[, c("s_1", "s_2")]))
  )
  for (cfg in configs) {
    e_hl <- list(fb = base_emat, ru = base_emat)
    for (normdist in c(FALSE, TRUE)) {
      fb <- testthat::with_mocked_bindings(
        .draw_emat = function(q, nrep) e_hl$fb,
        .package = "spuR",
        spuR::spurhalflife(
          y ~ 1,
          data = df,
          q = q,
          nrep = nrep,
          level = 95,
          normdist = normdist,
          lon = cfg$lon,
          lat = cfg$lat,
          coords_euclidean = cfg$euc,
          backend = "fallback"
        )
      )
      ru <- testthat::with_mocked_bindings(
        .draw_emat = function(q, nrep) e_hl$ru,
        .package = "spuR",
        spuR::spurhalflife(
          y ~ 1,
          data = df,
          q = q,
          nrep = nrep,
          level = 95,
          normdist = normdist,
          lon = cfg$lon,
          lat = cfg$lat,
          coords_euclidean = cfg$euc,
          backend = "rust",
          nthreads = 2L
        )
      )
      testthat::expect_identical(ru$backend, "rust")
      testthat::expect_length(ru$full$pv_mat, 131L)
      testthat::expect_equal(ru$max_dist, fb$max_dist, tolerance = 1e-10)
      testthat::expect_equal(ru$ci_l, fb$ci_l, tolerance = 1e-10)
      if (is.na(fb$ci_u)) {
        testthat::expect_true(is.na(ru$ci_u))
      } else {
        testthat::expect_equal(ru$ci_u, fb$ci_u, tolerance = 1e-10)
      }
    }
  }
})

testthat::test_that("spurtransform parity holds across all 4 transforms and options", {
  testthat::skip_if_not(is.loaded("wrap__rust_run_spur_transform", PACKAGE = "spuR"), "spuR native Rust library not loaded")
  df <- read_bench_200()[seq_len(60L), ]
  df$s_1[2L] <- df$s_1[1L]
  df$s_2[2L] <- df$s_2[1L]
  df$grp <- rep(c("A", "B", "C", "D"), length.out = nrow(df))

  for (cfg in list(
    list(lon = "s_1", lat = "s_2", euc = NULL, radius = 250000),
    list(lon = NULL, lat = NULL, euc = c("s_1", "s_2"), radius = 1.5)
  )) {
    for (tr in c("lbmgls", "nn", "iso")) {
      rad <- if (identical(tr, "iso")) cfg$radius else NULL
      fb <- spuR::spurtransform(
        y ~ x1, df, transformation = tr, radius = rad,
        lon = cfg$lon, lat = cfg$lat, coords_euclidean = cfg$euc, backend = "fallback"
      )
      ru <- spuR::spurtransform(
        y ~ x1, df, transformation = tr, radius = rad,
        lon = cfg$lon, lat = cfg$lat, coords_euclidean = cfg$euc, backend = "rust", nthreads = 2L
      )
      testthat::expect_equal(ru$h_y, fb$h_y, tolerance = 1e-9)
      testthat::expect_equal(ru$h_x1, fb$h_x1, tolerance = 1e-9)
      if (identical(tr, "lbmgls")) {
        testthat::expect_equal(mean(ru$h_y), 0, tolerance = 1e-12)
        testthat::expect_equal(ru$h_y[1L], ru$h_y[2L], tolerance = 1e-9)
      }
    }
  }

  # Cluster transform + separately + missing values
  df$y[5L] <- NA_real_
  fb_cl <- spuR::spurtransform(
    y ~ x1, df, transformation = "cluster", clustvar = "grp",
    lon = "s_1", lat = "s_2", separately = TRUE, backend = "fallback"
  )
  ru_cl <- spuR::spurtransform(
    y ~ x1, df, transformation = "cluster", clustvar = "grp",
    lon = "s_1", lat = "s_2", separately = TRUE, backend = "rust"
  )
  testthat::expect_equal(ru_cl$h_x1, fb_cl$h_x1, tolerance = 1e-12)
  testthat::expect_equal(ru_cl$h_y, fb_cl$h_y, tolerance = 1e-12)
})

testthat::test_that("spur pipeline is thread-invariant and matches individual rust calls", {
  testthat::skip_if_not(is.loaded("wrap__rust_run_spur_pipeline", PACKAGE = "spuR"), "spuR native Rust library not loaded")
  df <- read_bench_200()[seq_len(120L), ]
  p1 <- spuR::spur(
    y ~ x1 + x2, data = df, q = 8L, nrep = 400L,
    lon = "s_1", lat = "s_2", seed = 42L, backend = "rust", nthreads = 1L
  )
  p2 <- spuR::spur(
    y ~ x1 + x2, data = df, q = 8L, nrep = 400L,
    lon = "s_1", lat = "s_2", seed = 42L, backend = "rust", nthreads = 2L
  )
  testthat::expect_identical(p1$backend, "rust")
  testthat::expect_identical(p1$fits$levels$scpc$backend, "rust")
  testthat::expect_identical(p1$fits$transformed$scpc$backend, "rust")
  for (nm in c("i0", "i1", "i0resid", "i1resid")) {
    testthat::expect_equal(p1$tests[[nm]]$teststat, p2$tests[[nm]]$teststat, tolerance = 1e-13)
    testthat::expect_equal(p1$tests[[nm]]$pvalue, p2$tests[[nm]]$pvalue, tolerance = 1e-15)
  }
})

testthat::test_that("tied eigenspace triggers auto fallback and backend validation works", {
  k <- 0:15
  df <- data.frame(
    x_c = cos(2 * pi * k / 16),
    y_c = sin(2 * pi * k / 16),
    val = sin(k * 1.3) + cos(k * 0.7)
  )
  if (is.loaded("wrap__rust_run_spur_i0", PACKAGE = "spuR")) {
    ru_i0 <- spuR::spurtest_i0(
      val ~ 1, data = df, coords_euclidean = c("x_c", "y_c"),
      q = 8L, nrep = 2000L, seed = 42L, backend = "rust"
    )
    auto_i0 <- spuR::spurtest_i0(
      val ~ 1, data = df, coords_euclidean = c("x_c", "y_c"),
      q = 8L, nrep = 2000L, seed = 42L, backend = "auto"
    )
    testthat::expect_identical(ru_i0$backend, "rust")
    testthat::expect_true(ru_i0$tied)
    testthat::expect_identical(auto_i0$backend, "fallback")
    testthat::expect_true(auto_i0$tied)
  }

  testthat::expect_error(
    spuR::spurtest_i0(val ~ 1, data = df, coords_euclidean = c("x_c", "y_c"), q = 8L, nrep = 2000L, backend = "bogus"),
    "auto"
  )
  testthat::expect_error(
    spuR::spurtest_i0(val ~ 1, data = df, coords_euclidean = c("x_c", "y_c"), q = 8L, nrep = 2000L, nthreads = 0L),
    "positive integer"
  )
  testthat::with_mocked_bindings(
    .is_spur_rust_loaded = function(sym = "wrap__rust_run_spur_i0") FALSE,
    .package = "spuR",
    {
      auto_fb <- spuR::spurtest_i0(
        val ~ 1, data = df, coords_euclidean = c("x_c", "y_c"),
        q = 8L, nrep = 2000L, seed = 42L, backend = "auto"
      )
      testthat::expect_identical(auto_fb$backend, "fallback")
      testthat::expect_error(
        spuR::spurtest_i0(
          val ~ 1, data = df, coords_euclidean = c("x_c", "y_c"),
          q = 8L, nrep = 2000L, backend = "rust"
        ),
        "not loaded"
      )
    }
  )
})
