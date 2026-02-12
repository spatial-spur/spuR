testthat::test_that("SPUR transform parity across transformation types", {
  spur_skip_if_no_stata_or_spur()

  df <- spur_fixture_data(n = 180L)
  vars <- c("am", "fracblack")

  cases <- list(
    list(transformation = "lbmgls", latlong = TRUE, radius = NULL, clustvar = NULL),
    list(transformation = "nn", latlong = TRUE, radius = NULL, clustvar = NULL),
    list(transformation = "iso", latlong = TRUE, radius = 150000, clustvar = NULL),
    list(transformation = "cluster", latlong = TRUE, radius = NULL, clustvar = "state"),
    list(transformation = "lbmgls", latlong = FALSE, radius = NULL, clustvar = NULL),
    list(transformation = "nn", latlong = FALSE, radius = NULL, clustvar = NULL),
    list(transformation = "iso", latlong = FALSE, radius = 0.1, clustvar = NULL)
  )

  for (case in cases) {
    stata <- spur_run_stata_spurtransform(
      df = df,
      vars = vars,
      prefix = "h_",
      transformation = case$transformation,
      latlong = case$latlong,
      radius = case$radius,
      clustvar = case$clustvar,
      separately = FALSE,
      replace = TRUE
    )

    r_out <- spur_transform(
      data = df,
      vars = vars,
      prefix = "h_",
      transformation = case$transformation,
      radius = case$radius,
      clustvar = case$clustvar,
      latlong = case$latlong,
      coord_cols = c("s_1", "s_2"),
      replace = TRUE,
      separately = FALSE
    )

    r_cmp <- r_out[order(r_out$row_id), c("row_id", "h_am", "h_fracblack")]
    s_cmp <- stata[order(stata$row_id), c("row_id", "h_am", "h_fracblack")]

    testthat::expect_equal(r_cmp$row_id, s_cmp$row_id)
    testthat::expect_equal(r_cmp$h_am, s_cmp$h_am, tolerance = 1e-4)
    testthat::expect_equal(r_cmp$h_fracblack, s_cmp$h_fracblack, tolerance = 1e-4)
  }
})

testthat::test_that("SPUR transform parity for separately option with missing values", {
  spur_skip_if_no_stata_or_spur()

  df <- spur_fixture_data(n = 180L)
  df$fracblack[c(3, 7, 11, 25)] <- NA_real_
  df$racseg[c(5, 9, 13)] <- NA_real_
  vars <- c("fracblack", "racseg")

  for (separately in c(FALSE, TRUE)) {
    stata <- spur_run_stata_spurtransform(
      df = df,
      vars = vars,
      prefix = "h_",
      transformation = "lbmgls",
      latlong = TRUE,
      separately = separately,
      replace = TRUE
    )

    r_out <- spur_transform(
      data = df,
      vars = vars,
      prefix = "h_",
      transformation = "lbmgls",
      latlong = TRUE,
      coord_cols = c("s_1", "s_2"),
      separately = separately,
      replace = TRUE
    )

    r_cmp <- r_out[order(r_out$row_id), c("row_id", "h_fracblack", "h_racseg")]
    s_cmp <- stata[order(stata$row_id), c("row_id", "h_fracblack", "h_racseg")]

    testthat::expect_equal(r_cmp$row_id, s_cmp$row_id)
    testthat::expect_equal(r_cmp$h_fracblack, s_cmp$h_fracblack, tolerance = 1e-4)
    testthat::expect_equal(r_cmp$h_racseg, s_cmp$h_racseg, tolerance = 1e-4)
  }
})

testthat::test_that("SPUR transform replace and option-validation parity", {
  spur_skip_if_no_stata_or_spur()

  df <- spur_fixture_data(n = 180L)
  df$h_am <- 0

  stata_rc_no_replace <- spur_run_stata_command_rc(
    df = df,
    command_line = "spurtransform am, prefix(\"h_\") transformation(lbmgls) latlong"
  )
  testthat::expect_true(stata_rc_no_replace != 0L)

  testthat::expect_error(
    spur_transform(
      data = df,
      vars = "am",
      prefix = "h_",
      transformation = "lbmgls",
      latlong = TRUE,
      coord_cols = c("s_1", "s_2"),
      replace = FALSE
    )
  )

  stata_rc_replace <- spur_run_stata_command_rc(
    df = df,
    command_line = "spurtransform am, prefix(\"h_\") transformation(lbmgls) latlong replace"
  )
  testthat::expect_equal(stata_rc_replace, 0L)

  testthat::expect_no_error(
    spur_transform(
      data = df,
      vars = "am",
      prefix = "h_",
      transformation = "lbmgls",
      latlong = TRUE,
      coord_cols = c("s_1", "s_2"),
      replace = TRUE
    )
  )

  testthat::expect_error(
    spur_transform(
      data = df,
      vars = "am",
      prefix = "h_",
      transformation = "lbmgls",
      radius = 1,
      latlong = TRUE,
      coord_cols = c("s_1", "s_2")
    )
  )
  testthat::expect_error(
    spur_transform(
      data = df,
      vars = "am",
      prefix = "h_",
      transformation = "iso",
      clustvar = "state",
      radius = 1,
      latlong = TRUE,
      coord_cols = c("s_1", "s_2")
    )
  )
})
