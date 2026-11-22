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
  cat(sprintf("Spatial half-life %s%% confidence set\n", level))
  cat("---------------------------------------\n")
  .print_halflife_components(x$components)
  invisible(x)
}

.print_halflife_components <- function(components) {
  if (nrow(components) == 0L) {
    cat("Empty confidence set.\n")
  } else {
    cat(sprintf("[%.6g, %.6g]\n", components[, 1L], components[, 2L]), sep = "")
  }
}
