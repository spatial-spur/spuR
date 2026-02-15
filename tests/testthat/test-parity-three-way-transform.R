spur_three_way_max_abs_diff <- function(a, b, c) {
  diffs <- c(abs(a - b), abs(a - c), abs(b - c))
  diffs <- diffs[is.finite(diffs)]
  if (!length(diffs)) {
    return(0)
  }
  max(diffs)
}

testthat::test_that("SPUR 3-way transform parity: lbmgls across distance types", {
  spur_skip_if_no_stata_or_spur()
  spur_skip_if_no_matlab_or_mw()

  df <- spur_fixture_data(n = 550L)
  vars <- c("am", "fracblack")

  for (latlong in c(TRUE, FALSE)) {
    stata <- spur_run_stata_spurtransform(
      df = df,
      vars = vars,
      prefix = "h_",
      transformation = "lbmgls",
      latlong = latlong,
      separately = FALSE,
      replace = TRUE
    )

    matlab <- spur_run_matlab_spurtransform(
      df = df,
      vars = vars,
      prefix = "h_",
      transformation = "lbmgls",
      latlong = latlong,
      separately = FALSE,
      replace = TRUE
    )

    coord_args <- if (isTRUE(latlong)) {
      list(lon = "lon", lat = "lat")
    } else {
      list(coords_euclidean = c("s_1", "s_2"))
    }
    r_out <- do.call(
      spurtransform,
      c(
        list(
          formula = ~ am + fracblack,
          data = df,
          prefix = "h_",
          transformation = "lbmgls",
          replace = TRUE,
          separately = FALSE
        ),
        coord_args
      )
    )

    r_cmp <- r_out[order(r_out$row_id), c("row_id", "h_am", "h_fracblack")]
    s_cmp <- stata[order(stata$row_id), c("row_id", "h_am", "h_fracblack")]
    m_cmp <- matlab[order(matlab$row_id), c("row_id", "h_am", "h_fracblack")]

    testthat::expect_equal(r_cmp$row_id, s_cmp$row_id)
    testthat::expect_equal(r_cmp$row_id, m_cmp$row_id)

    testthat::expect_lte(
      spur_three_way_max_abs_diff(r_cmp$h_am, s_cmp$h_am, m_cmp$h_am),
      2.5e-4
    )
    testthat::expect_lte(
      spur_three_way_max_abs_diff(r_cmp$h_fracblack, s_cmp$h_fracblack, m_cmp$h_fracblack),
      2.5e-4
    )
  }
})

testthat::test_that("SPUR 3-way transform parity: lbmgls with separately and missing values", {
  spur_skip_if_no_stata_or_spur()
  spur_skip_if_no_matlab_or_mw()

  df <- spur_fixture_data(n = 550L)
  df$fracblack[c(3, 7, 11, 25, 41, 59)] <- NA_real_
  df$racseg[c(5, 9, 13, 31, 47)] <- NA_real_
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

    matlab <- spur_run_matlab_spurtransform(
      df = df,
      vars = vars,
      prefix = "h_",
      transformation = "lbmgls",
      latlong = TRUE,
      separately = separately,
      replace = TRUE
    )

    r_out <- spurtransform(
      formula = ~ fracblack + racseg,
      data = df,
      prefix = "h_",
      transformation = "lbmgls",
      lon = "lon",
      lat = "lat",
      separately = separately,
      replace = TRUE
    )

    r_cmp <- r_out[order(r_out$row_id), c("row_id", "h_fracblack", "h_racseg")]
    s_cmp <- stata[order(stata$row_id), c("row_id", "h_fracblack", "h_racseg")]
    m_cmp <- matlab[order(matlab$row_id), c("row_id", "h_fracblack", "h_racseg")]

    testthat::expect_equal(r_cmp$row_id, s_cmp$row_id)
    testthat::expect_equal(r_cmp$row_id, m_cmp$row_id)

    testthat::expect_lte(
      spur_three_way_max_abs_diff(r_cmp$h_fracblack, s_cmp$h_fracblack, m_cmp$h_fracblack),
      2.5e-4
    )
    testthat::expect_lte(
      spur_three_way_max_abs_diff(r_cmp$h_racseg, s_cmp$h_racseg, m_cmp$h_racseg),
      2.5e-4
    )
  }
})
