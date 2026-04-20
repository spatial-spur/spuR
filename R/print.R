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
  cat(sprintf("Lower bound : %9.4f\n", x$ci_l))
  if (is.na(x$ci_u)) {
    cat("Upper bound :       Inf\n")
  } else {
    cat(sprintf("Upper bound : %9.4f\n", x$ci_u))
  }
  cat("---------------------------------------\n")
  invisible(x)
}
