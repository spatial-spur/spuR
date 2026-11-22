#' @export
print.spur_test <- function(x, ...) {
  type <- attr(x, "test_type")
  if (is.null(type)) type <- "SPUR"
  cat(sprintf("Spatial %s Test\n", type))
  cat("---------------------------------------\n")
  cat(sprintf("Test statistic : %9.4f\n", x$teststat))
  cat(sprintf("P-value        : %9.4f\n", x$pvalue))
  cat("---------------------------------------\n")
  invisible(x)
}

#' @export
print.spur_halflife <- function(x, ...) {
  level <- attr(x, "level")
  if (is.null(level)) level <- 95
  cat(sprintf("Spatial half-life %s%% confidence grid hull (Table 1)\n", level))
  cat("---------------------------------------\n")
  .print_halflife_bounds(x$ci_l, x$ci_u)
  invisible(x)
}

.print_halflife_bounds <- function(lower, upper) {
  if (is.na(lower)) {
    cat("No tested half-life was accepted.\n")
  } else {
    cat(sprintf("Lower bound : %.6g\n", lower))
    cat(sprintf("Upper bound : %s\n", if (is.na(upper)) {
      "Inf (100 sample diameters accepted on the grid)"
    } else {
      sprintf("%.6g", upper)
    }))
  }
}
