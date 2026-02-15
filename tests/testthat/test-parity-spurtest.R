testthat::test_that("SPUR test parity: i0 and i1 (latlong and euclidean)", {
  spur_skip_if_no_stata_or_spur()

  df <- spur_fixture_data(n = 180L)
  q <- 15L
  nrep <- 100000L

  cases <- list(
    list(test_type = "i0", latlong = TRUE),
    list(test_type = "i1", latlong = TRUE),
    list(test_type = "i0", latlong = FALSE),
    list(test_type = "i1", latlong = FALSE)
  )

  for (case in cases) {
    stata <- spur_run_stata_spurtest(
      df = df,
      test_type = case$test_type,
      varlist = "am",
      q = q,
      nrep = nrep,
      latlong = case$latlong
    )

    set.seed(123)
    coord_args <- if (isTRUE(case$latlong)) {
      list(lon = "lon", lat = "lat")
    } else {
      list(coords_euclidean = c("s_1", "s_2"))
    }
    r_out <- if (case$test_type == "i0") {
      do.call(
        spurtest_i0,
        c(
          list(
            formula = am ~ 1,
            data = df,
            q = q,
            nrep = nrep
          ),
          coord_args
        )
      )
    } else {
      do.call(
        spurtest_i1,
        c(
          list(
            formula = am ~ 1,
            data = df,
            q = q,
            nrep = nrep
          ),
          coord_args
        )
      )
    }

    teststat_rel <- abs(as.numeric(r_out$teststat) - as.numeric(stata$teststat)) / max(abs(as.numeric(stata$teststat)), 1e-8)
    ha_rel <- abs(as.numeric(r_out$ha_param) - as.numeric(stata$ha_param)) / max(abs(as.numeric(stata$ha_param)), 1e-8)
    cv_rel <- abs(unname(as.numeric(r_out$cv)) - as.numeric(stata[1, c("cv1", "cv2", "cv3")])) /
      pmax(abs(as.numeric(stata[1, c("cv1", "cv2", "cv3")])), 1e-8)

    testthat::expect_lte(teststat_rel, 0.10)
    testthat::expect_lte(abs(as.numeric(r_out$pvalue) - as.numeric(stata$pvalue)), 0.005)
    testthat::expect_lte(ha_rel, 0.20)
    testthat::expect_lte(
      max(cv_rel),
      0.10
    )
  }
})

testthat::test_that("SPUR test parity: i0resid and i1resid (latlong and euclidean)", {
  spur_skip_if_no_stata_or_spur()

  df <- spur_fixture_data(n = 180L)
  q <- 15L
  nrep <- 100000L
  indep <- c("fracblack", "racseg")

  cases <- list(
    list(test_type = "i0resid", latlong = TRUE),
    list(test_type = "i1resid", latlong = TRUE),
    list(test_type = "i0resid", latlong = FALSE),
    list(test_type = "i1resid", latlong = FALSE)
  )

  for (case in cases) {
    stata <- spur_run_stata_spurtest(
      df = df,
      test_type = case$test_type,
      varlist = paste(c("am", indep), collapse = " "),
      q = q,
      nrep = nrep,
      latlong = case$latlong
    )

    set.seed(123)
    coord_args <- if (isTRUE(case$latlong)) {
      list(lon = "lon", lat = "lat")
    } else {
      list(coords_euclidean = c("s_1", "s_2"))
    }
    r_out <- if (case$test_type == "i0resid") {
      do.call(
        spurtest_i0resid,
        c(
          list(
            formula = am ~ fracblack + racseg,
            data = df,
            q = q,
            nrep = nrep
          ),
          coord_args
        )
      )
    } else {
      do.call(
        spurtest_i1resid,
        c(
          list(
            formula = am ~ fracblack + racseg,
            data = df,
            q = q,
            nrep = nrep
          ),
          coord_args
        )
      )
    }

    teststat_rel <- abs(as.numeric(r_out$teststat) - as.numeric(stata$teststat)) / max(abs(as.numeric(stata$teststat)), 1e-8)
    ha_rel <- abs(as.numeric(r_out$ha_param) - as.numeric(stata$ha_param)) / max(abs(as.numeric(stata$ha_param)), 1e-8)
    cv_rel <- abs(unname(as.numeric(r_out$cv)) - as.numeric(stata[1, c("cv1", "cv2", "cv3")])) /
      pmax(abs(as.numeric(stata[1, c("cv1", "cv2", "cv3")])), 1e-8)

    testthat::expect_lte(teststat_rel, 0.10)
    testthat::expect_lte(abs(as.numeric(r_out$pvalue) - as.numeric(stata$pvalue)), 0.005)
    testthat::expect_lte(ha_rel, 0.20)
    testthat::expect_lte(
      max(cv_rel),
      0.10
    )
  }
})
