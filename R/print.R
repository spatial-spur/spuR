#' @export
print.spur_test <- function(x, ...) {
  type <- attr(x, "test_type") %||% "SPUR"
  cat(sprintf("Spatial %s Test\n", type))
  cat("---------------------------------------\n")
  cat(sprintf("Test statistic : %9.4f\n", x$teststat))
  cat(sprintf("P-value        : %9.4f\n", x$pvalue))
  cat("---------------------------------------\n")
  invisible(x)
}

#' @export
print.spur_halflife <- function(x, ...) {
  level <- attr(x, "level") %||% 95
  cat(sprintf("Spatial Half-Life %s%% Confidence Interval\n", level))
  cat("---------------------------------------\n")
  .print_halflife_bounds(x$ci_l, x$ci_u)
  invisible(x)
}

.print_halflife_bounds <- function(lower, upper) {
  if (is.na(lower)) {
    cat("Confidence set is empty.\n")
  } else {
    cat(sprintf("Lower bound : %.6g\n", lower))
    cat(sprintf("Upper bound : %s\n", if (is.na(upper)) "Inf" else sprintf("%.6g", upper)))
  }
}
