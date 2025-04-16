library(testthat)
library(dplyr)

rm(list = c("getcbar", "lvech"))
devtools::load_all()

test_that("spur_i0_test returns expected values", {
  # Load your dataset from the data/ directory
  load("data/chetty.Rda")  # This should load an object called chetty

  # Rename coordinate columns as expected by the function (e.g., "Lat" -> "lat", "Lon" -> "lon")
  chetty <- chetty %>%
    rename(lat = Lat,
           lon = Lon)

  # Call the spur_i0 wrapper function.
  result <- spur_i0(var = "AM", q = 15, nrep = 100000, latlong = TRUE, data = chetty)

  # For spur_i0_test, you expect:
  #   Test Statistic (LFST) :    3.1451
  #   P-value               :    0.0013
  # Allow a small tolerance for numerical differences (e.g., tolerance = 1e-4)
  cat("Test Statistic (LFST):", result$teststat, "\n")
  cat("P-value:", result$pvalue, "\n")

  expect_equal(round(result$teststat, 4), 3.1451, tolerance = 1e-4)
  expect_equal(round(result$pvalue, 4), 0.0013, tolerance = 1e-4)
})
