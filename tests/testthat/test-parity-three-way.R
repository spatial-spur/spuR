spur_three_way_rel_spread <- function(values, floor = 1e-8) {
  (max(values) - min(values)) / max(max(abs(values)), floor)
}

testthat::test_that("SPUR 3-way parity for spurtest core outputs", {
  spur_skip_if_no_stata_or_spur()
  spur_skip_if_no_matlab_or_mw()

  df <- spur_fixture_data(n = 550L)
  q <- 15L
  nrep <- 100000L
  indep <- c("fracblack", "racseg")

  cases <- list(
    list(name = "i0_latlong", test_type = "i0", latlong = TRUE, dep_var = "am", indep_vars = character(0)),
    list(name = "i1_latlong", test_type = "i1", latlong = TRUE, dep_var = "am", indep_vars = character(0)),
    list(name = "i0resid_latlong", test_type = "i0resid", latlong = TRUE, dep_var = "am", indep_vars = indep),
    list(name = "i1resid_latlong", test_type = "i1resid", latlong = TRUE, dep_var = "am", indep_vars = indep),
    list(name = "i0_euclidean", test_type = "i0", latlong = FALSE, dep_var = "am", indep_vars = character(0)),
    list(name = "i1_euclidean", test_type = "i1", latlong = FALSE, dep_var = "am", indep_vars = character(0))
  )

  for (case in cases) {
    stata <- spur_run_stata_spurtest(
      df = df,
      test_type = case$test_type,
      varlist = paste(c(case$dep_var, case$indep_vars), collapse = " "),
      q = q,
      nrep = nrep,
      latlong = case$latlong
    )

    matlab <- spur_run_matlab_spurtest(
      df = df,
      test_type = case$test_type,
      dep_var = case$dep_var,
      indep_vars = case$indep_vars,
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
    } else if (case$test_type == "i1") {
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
    } else if (case$test_type == "i0resid") {
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

    teststat_vals <- c(
      as.numeric(r_out$teststat),
      as.numeric(stata$teststat),
      as.numeric(matlab$teststat)
    )
    pvalue_vals <- c(
      as.numeric(r_out$pvalue),
      as.numeric(stata$pvalue),
      as.numeric(matlab$pvalue)
    )
    ha_vals <- c(
      as.numeric(r_out$ha_param),
      as.numeric(stata$ha_param),
      as.numeric(matlab$ha_param)
    )

    cv_r <- unname(as.numeric(r_out$cv))
    cv_s <- as.numeric(stata[1, c("cv1", "cv2", "cv3")])
    cv_m <- as.numeric(matlab[1, c("cv1", "cv2", "cv3")])

    testthat::expect_lte(
      spur_three_way_rel_spread(teststat_vals),
      0.08
    )
    testthat::expect_lte(
      max(pvalue_vals) - min(pvalue_vals),
      0.015
    )
    testthat::expect_lte(
      spur_three_way_rel_spread(ha_vals),
      0.20
    )

    for (i in seq_len(3)) {
      cv_vals <- c(cv_r[i], cv_s[i], cv_m[i])
      testthat::expect_lte(
        spur_three_way_rel_spread(cv_vals),
        0.08
      )
    }
  }
})

testthat::test_that("SPUR 3-way parity for half-life outputs", {
  spur_skip_if_no_stata_or_spur()
  spur_skip_if_no_matlab_or_mw()

  df <- spur_fixture_data(n = 550L)
  q <- 15L
  nrep <- 100000L
  level <- 90

  cases <- list(
    list(name = "latlong_normdist", latlong = TRUE, normdist = TRUE),
    list(name = "latlong_raw", latlong = TRUE, normdist = FALSE),
    list(name = "euclidean_normdist", latlong = FALSE, normdist = TRUE),
    list(name = "euclidean_raw", latlong = FALSE, normdist = FALSE)
  )

  for (case in cases) {
    stata <- spur_run_stata_spurhalflife(
      df = df,
      var = "am",
      q = q,
      nrep = nrep,
      level = level,
      latlong = case$latlong,
      normdist = case$normdist
    )

    matlab <- spur_run_matlab_halflife(
      df = df,
      var = "am",
      q = q,
      nrep = nrep,
      level = level,
      latlong = case$latlong,
      normdist = case$normdist
    )

    set.seed(123)
    coord_args <- if (isTRUE(case$latlong)) {
      list(lon = "lon", lat = "lat")
    } else {
      list(coords_euclidean = c("s_1", "s_2"))
    }
    r_out <- do.call(
      spurhalflife,
      c(
        list(
          formula = am ~ 1,
          data = df,
          q = q,
          nrep = nrep,
          level = level,
          normdist = case$normdist
        ),
        coord_args
      )
    )

    max_dist_vals <- c(
      as.numeric(r_out$max_dist),
      as.numeric(stata$max_dist),
      as.numeric(matlab$max_dist)
    )
    testthat::expect_lte(
      spur_three_way_rel_spread(max_dist_vals),
      5e-5
    )

    ci_l_vals <- c(
      as.numeric(if (case$normdist) r_out$ci_l else r_out$ci_l / r_out$max_dist),
      as.numeric(if (case$normdist) stata$ci_l else stata$ci_l / stata$max_dist),
      as.numeric(if (case$normdist) matlab$ci_l else matlab$ci_l / matlab$max_dist)
    )
    testthat::expect_lte(
      max(ci_l_vals) - min(ci_l_vals),
      0.05
    )

    ci_u_vals <- c(
      if (case$normdist) as.numeric(r_out$ci_u) else as.numeric(r_out$ci_u / r_out$max_dist),
      if (case$normdist) as.numeric(stata$ci_u) else as.numeric(stata$ci_u / stata$max_dist),
      if (case$normdist) as.numeric(matlab$ci_u) else as.numeric(matlab$ci_u / matlab$max_dist)
    )
    if (any(is.na(ci_u_vals))) {
      testthat::expect_true(all(is.na(ci_u_vals)))
    } else {
      testthat::expect_lte(
        max(ci_u_vals) - min(ci_u_vals),
        0.08
      )
    }
  }
})
