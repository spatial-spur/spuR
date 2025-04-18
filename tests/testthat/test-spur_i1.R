library(testthat)
library(dplyr)
library(here)

#rm(list = c("getcbar", "lvech"))
rm(list = ls())
devtools::load_all()

test_that("spur_i1_test returns expected values", {
  # Load your dataset from the data/ directory
  load(here("data/chetty.Rda"))  # This should load an object called chetty

  # Rename coordinate columns as expected by the function (e.g., "Lat" -> "lat", "Lon" -> "lon")
  chetty <- chetty %>%
    rename(lat = Lat,
           lon = Lon)

  # drop if state Hawaii or Alaska
  chetty <- chetty %>%
    filter(!(State == "HI" | State == "AK"))

  # Call the spur_i0 wrapper function.
  result <- spur_i1(var = "AM", q = 15, nrep = 100000, latlong = TRUE, data = chetty)

  # For spur_i1_test, you expect:
  #   Test Statistic (LFST) :    5.8605
  #   P-value               :    0.3788
  # Allow a small tolerance for numerical differences (e.g., tolerance = 1e-4)
  cat("Test Statistic (LFST):", result$teststat, "\n")
  cat("P-value:", result$pvalue, "\n")

  expect_equal(round(result$teststat, 4), 5.8605, tolerance = 1e-4)
  expect_equal(round(result$pvalue, 4), 0.3788, tolerance = 1e-3)
})
