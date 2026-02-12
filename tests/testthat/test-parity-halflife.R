testthat::test_that("SPUR half-life parity across options", {
  spur_skip_if_no_stata_or_spur()

  df <- spur_fixture_data(n = 180L)
  q <- 15L
  nrep <- 20000L
  level <- 90

  cases <- list(
    list(latlong = TRUE, normdist = TRUE),
    list(latlong = TRUE, normdist = FALSE),
    list(latlong = FALSE, normdist = TRUE),
    list(latlong = FALSE, normdist = FALSE)
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

    set.seed(123)
    r_out <- spur_halflife(
      var = "am",
      q = q,
      nrep = nrep,
      level = level,
      latlong = case$latlong,
      normdist = case$normdist,
      data = df,
      coord_cols = c("s_1", "s_2")
    )

    testthat::expect_equal(
      as.numeric(r_out$max_dist),
      as.numeric(stata$max_dist),
      tolerance = 1e-6
    )

    r_ci_l <- if (case$normdist) r_out$ci_l else r_out$ci_l / r_out$max_dist
    s_ci_l <- if (case$normdist) as.numeric(stata$ci_l) else as.numeric(stata$ci_l / stata$max_dist)
    testthat::expect_lte(abs(as.numeric(r_ci_l) - as.numeric(s_ci_l)), 0.05)

    if (is.na(r_out$ci_u) || is.na(stata$ci_u)) {
      testthat::expect_true(is.na(r_out$ci_u) && is.na(stata$ci_u))
    } else {
      r_ci_u <- if (case$normdist) r_out$ci_u else r_out$ci_u / r_out$max_dist
      s_ci_u <- if (case$normdist) as.numeric(stata$ci_u) else as.numeric(stata$ci_u / stata$max_dist)
      testthat::expect_lte(abs(as.numeric(r_ci_u) - as.numeric(s_ci_u)), 0.10)
    }
  }
})

testthat::test_that("SPUR half-life level validation parity", {
  spur_skip_if_no_stata_or_spur()

  df <- spur_fixture_data(n = 180L)
  in_csv <- tempfile(fileext = ".csv")
  utils::write.csv(df, in_csv, row.names = FALSE)
  run <- spur_run_stata_do(c(
    "clear all",
    sprintf("import delimited using \"%s\", varnames(1) clear", spur_stata_path(in_csv)),
    "spurhalflife am, level(101) latlong",
    "exit, clear"
  ))
  testthat::expect_true(any(grepl("Invalid level\\.", run$log_lines, fixed = FALSE)))

  testthat::expect_error(
    spur_halflife(
      var = "am",
      level = 101,
      latlong = TRUE,
      data = df,
      coord_cols = c("s_1", "s_2")
    )
  )
})
