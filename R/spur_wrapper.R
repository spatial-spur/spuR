#' @keywords internal
.display_term_name <- function(name) {
  sub("^h_", "", name)
}

#' @keywords internal
.format_decimal <- function(x) {
  sprintf("%.4f", as.numeric(x))
}

#' @keywords internal
.format_count <- function(x) {
  format(round(as.numeric(x)), trim = TRUE, scientific = FALSE)
}

#' @keywords internal
.center_text <- function(text, width) {
  text <- as.character(text)
  if (nchar(text) >= width) {
    return(text)
  }

  pad <- width - nchar(text)
  left <- floor(pad / 2)
  right <- pad - left
  paste0(strrep(" ", left), text, strrep(" ", right))
}

#' @keywords internal
.render_row <- function(label,
                        left,
                        right,
                        label_width,
                        value_width,
                        gap = 4L,
                        header = FALSE) {
  pad <- strrep(" ", gap)

  if (isTRUE(header)) {
    sprintf(
      "%-*s%s%-*s%s%-*s",
      label_width,
      label,
      pad,
      value_width,
      left,
      pad,
      value_width,
      right
    )
  } else {
    sprintf(
      "%-*s%s%*s%s%*s",
      label_width,
      label,
      pad,
      value_width,
      left,
      pad,
      value_width,
      right
    )
  }
}

#' @keywords internal
.collect_coefficient_rows <- function(x) {
  levels_stats <- x$fits$levels$scpc$scpcstats
  transformed_stats <- x$fits$transformed$scpc$scpcstats

  levels_names <- rownames(levels_stats)
  transformed_names <- rownames(transformed_stats)
  term_names <- x$fits$transformed$term_names
  if (is.null(term_names)) {
    transformed_names <- .display_term_name(transformed_names)
  } else {
    transformed_names <- unname(term_names[transformed_names])
  }
  row_order <- unique(c(levels_names, transformed_names))

  out <- vector("list", length = 2L * length(row_order))
  k <- 1L

  for (name in row_order) {
    i_level <- match(name, levels_names)
    i_trans <- match(name, transformed_names)

    out[[k]] <- c(
      name,
      if (!is.na(i_level)) .format_decimal(levels_stats[i_level, "Coef"]) else "",
      if (!is.na(i_trans)) .format_decimal(transformed_stats[i_trans, "Coef"]) else ""
    )
    k <- k + 1L

    out[[k]] <- c(
      "",
      if (!is.na(i_level)) sprintf("(%s)", .format_decimal(levels_stats[i_level, "Std_Err"])) else "",
      if (!is.na(i_trans)) sprintf("(%s)", .format_decimal(transformed_stats[i_trans, "Std_Err"])) else ""
    )
    k <- k + 1L
  }

  out
}

#' @keywords internal
.collect_diagnostic_rows <- function(x) {
  list(
    c("i0", .format_decimal(x$tests$i0$teststat), .format_decimal(x$tests$i0$pvalue)),
    c("i1", .format_decimal(x$tests$i1$teststat), .format_decimal(x$tests$i1$pvalue)),
    c("i0resid", .format_decimal(x$tests$i0resid$teststat), .format_decimal(x$tests$i0resid$pvalue)),
    c("i1resid", .format_decimal(x$tests$i1resid$teststat), .format_decimal(x$tests$i1resid$pvalue))
  )
}

#' @keywords internal
.collect_model_rows <- function(x) {
  levels_sum <- summary(x$fits$levels$model)
  transformed_sum <- summary(x$fits$transformed$model)

  list(
    c("N", .format_count(stats::nobs(x$fits$levels$model)), .format_count(stats::nobs(x$fits$transformed$model))),
    c("R-squared", .format_decimal(levels_sum$r.squared), .format_decimal(transformed_sum$r.squared)),
    c("Adj. R-squared", .format_decimal(levels_sum$adj.r.squared), .format_decimal(transformed_sum$adj.r.squared)),
    c("SCPC q", .format_count(x$fits$levels$scpc$q), .format_count(x$fits$transformed$scpc$q)),
    c("SCPC cv", .format_decimal(x$fits$levels$scpc$cv), .format_decimal(x$fits$transformed$scpc$cv)),
    c("SCPC avc", .format_decimal(x$fits$levels$scpc$avc), .format_decimal(x$fits$transformed$scpc$avc))
  )
}

#' @keywords internal
.render_spur_summary <- function(x) {
  coefficient_rows <- .collect_coefficient_rows(x)
  diagnostic_rows <- .collect_diagnostic_rows(x)
  model_rows <- .collect_model_rows(x)

  label_width <- max(
    nchar("Coefficient"),
    max(nchar(c(
      vapply(coefficient_rows, `[[`, character(1), 1L),
      vapply(diagnostic_rows, `[[`, character(1), 1L),
      vapply(model_rows, `[[`, character(1), 1L)
    )))
  )
  value_width <- max(
    nchar("Transformed"),
    max(nchar(c(
      unlist(lapply(coefficient_rows, `[`, 2:3), use.names = FALSE),
      unlist(lapply(diagnostic_rows, `[`, 2:3), use.names = FALSE),
      unlist(lapply(model_rows, `[`, 2:3), use.names = FALSE)
    )))
  )

  gap <- 4L
  table_width <- label_width + gap + value_width + gap + value_width
  rule <- strrep("-", table_width)

  depvar <- all.vars(stats::formula(x$fits$levels$model))[1]
  span_width <- value_width + gap + value_width
  span_indent <- label_width + gap

  lines <- c(
    rule,
    rule,
    .center_text("SPUR Diagnostics", table_width),
    rule,
    .render_row("Test", "LR", "p-value", label_width, value_width, gap, header = TRUE),
    vapply(
      diagnostic_rows,
      function(row) .render_row(row[1], row[2], row[3], label_width, value_width, gap),
      character(1)
    ),
    rule,
    rule,
    "",
    .center_text("Regression results", table_width),
    rule,
    rule,
    paste0(strrep(" ", span_indent), .center_text(depvar, span_width)),
    paste0(strrep(" ", span_indent), strrep("-", span_width)),
    .render_row("Coefficient", "Levels", "Transformed", label_width, value_width, gap, header = TRUE),
    rule,
    vapply(
      coefficient_rows,
      function(row) .render_row(row[1], row[2], row[3], label_width, value_width, gap),
      character(1)
    ),
    rule,
    vapply(
      model_rows,
      function(row) .render_row(row[1], row[2], row[3], label_width, value_width, gap),
      character(1)
    ),
    rule,
    rule
  )

  paste(lines, collapse = "\n")
}

#' Full SPUR pipeline
#'
#' Run the SPUR diagnostics, fit the levels and transformed regressions,
#' and compute SCPC inference for both branches.
#'
#' @param formula Two-sided regression formula.
#' @param data Data frame containing model and coordinate variables.
#' @param q Number of low-frequency weighted averages (at least 2).
#' @param nrep Number of Monte Carlo draws.
#' @param lon Optional longitude column name.
#' @param lat Optional latitude column name.
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param seed Optional integer seed for the draws shared by all four SPUR tests.
#' @param avc SCPC tuning parameter passed to \code{scpcR::scpc()}.
#' @param uncond Passed through to \code{scpcR::scpc()}.
#' @param cvs Passed through to \code{scpcR::scpc()}.
#' @param backend Backend selector: \code{"auto"} (default), \code{"rust"}, or
#'   \code{"fallback"}.
#' @param nthreads Optional positive integer thread count for parallel Rust
#'   execution.
#' @details Both regressions and the diagnostics use the levels model's estimation
#'   sample. The spatial transformation applies to the evaluated response and
#'   design columns. For example, a log response is logged before transformation.
#' @return A list of class \code{"spur_result"} with \code{tests}, \code{fits},
#'   \code{stage_backends}, and \code{fallback_reasons}.
#' @useDynLib spuR, .registration = TRUE
#' @export
spur <- function(formula,
                 data,
                 q = 15L,
                 nrep = 100000L,
                 lon = NULL,
                 lat = NULL,
                 coords_euclidean = NULL,
                 seed = 42L,
                 avc = 0.03,
                 uncond = FALSE,
                 cvs = FALSE,
                 backend = c("auto", "rust", "fallback"),
                 nthreads = NULL) {
  backend_resolved <- .resolve_spur_backend(backend, missing(backend))
  nthreads_resolved <- .resolve_spur_nthreads(nthreads)
  .validate_common_args(data = data, q = q, nrep = nrep)

  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` for spur() must be two-sided, e.g. `y ~ x1 + x2`.")
  }

  depvar <- all.vars(formula[[2L]])
  if (length(depvar) != 1L) {
    stop("`formula` for spur() must have exactly one dependent variable.")
  }

  # Both diagnostics and regressions use the levels fit's estimation sample.
  # Transform the evaluated response and design, including formula expressions.
  levels_model <- stats::lm(formula = formula, data = data, na.action = stats::na.omit)
  mf <- stats::model.frame(levels_model)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || (!is.null(dim(y)) && ncol(y) != 1L) || any(!is.finite(y))) {
    stop("The fitted response for spur() must be univariate, numeric and finite.")
  }
  if (!is.null(stats::model.offset(mf))) stop("Offsets are not supported by spur().")
  design <- stats::model.matrix(levels_model)
  has_intercept <- "(Intercept)" %in% colnames(design)
  x_in <- if (has_intercept) design else cbind("(Intercept)" = 1, design)
  resid_dof <- nrow(x_in) - ncol(x_in)
  if (q > resid_dof) {
    stop(sprintf(
      "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
      q, resid_dof
    ))
  }

  sample <- match(rownames(mf), rownames(data))
  if (anyNA(sample)) stop("Could not align the fitted SPUR sample with the data.")
  sample_data <- data[sample, , drop = FALSE]
  coord_info <- .resolve_spur_coords(
    sample_data, rep(TRUE, length(sample)), lon, lat, coords_euclidean
  )
  regressors <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  variables <- cbind(as.numeric(y), regressors)
  storage.mode(variables) <- "double"
  storage.mode(x_in) <- "double"

  rust_available <- !identical(backend_resolved, "fallback") &&
    .require_spur_rust_if_requested(backend_resolved, "wrap__rust_run_spur_pipeline")
  if (!is.null(seed)) set.seed(seed)
  emat <- .draw_emat(q = q, nrep = nrep)

  kinds <- c("i0", "i1", "i0resid", "i1resid")
  tied <- stats::setNames(rep(FALSE, length(kinds)), kinds)
  reason <- if (identical(backend_resolved, "fallback")) "requested" else "native_unavailable"
  diagnostics <- NULL
  if (rust_available) {
    pipe <- .call_spur_rust(
      "wrap__rust_run_spur_pipeline",
      .swap_coords_for_rust(coord_info$coords, coord_info$latlong),
      isTRUE(coord_info$latlong), as.numeric(y), x_in, emat, variables, TRUE,
      nthreads_resolved
    )
    tied[] <- vapply(pipe[kinds], function(test) isTRUE(test$tied), logical(1L))
    if (any(tied) && identical(backend_resolved, "auto")) {
      reason <- "tied_eigenspace"
    } else {
      diagnostics <- pipe[kinds]
      transformed <- pipe$transformed_vars
      reason <- NULL
    }
  }
  if (is.null(diagnostics)) {
    d <- .normalized_distmat(coord_info$coords, coord_info$latlong)
    y_mat <- matrix(y, ncol = 1L)
    diagnostics <- list(
      i0 = spur_i0_test(y_mat, d, emat),
      i1 = spur_i1_test(y_mat, d, emat),
      i0resid = spur_i0_resid_test(y_mat, x_in, d, emat),
      i1resid = spur_i1_resid_test(y_mat, x_in, d, emat)
    )
    transformed <- .apply_lbm_transform(d, variables)
  }

  operation_backend <- if (is.null(reason)) "rust" else "fallback"
  types <- c("I(0)", "I(1)", "I(0) Residual", "I(1) Residual")
  titles <- c(
    "Spatial I(0) Test Results", "Spatial I(1) Test Results",
    "Spatial I(0) Test Results for Residuals", "Spatial I(1) Test Results for Residuals"
  )
  tests <- stats::setNames(lapply(seq_along(kinds), function(i) {
    res <- diagnostics[[i]]
    is_i0 <- i %in% c(1L, 3L)
    .format_spur_test_from_full(
      res, test_type = types[i], title = titles[i],
      stat_label = if (is_i0) "LFST" else "LFUR",
      cv = if (is_i0) res$cvalue else res$cv_vec,
      verbose = FALSE, backend = operation_backend, tied = tied[[i]]
    )
  }), kinds)

  # Use valid distinct column names, and retain original labels for the summary.
  preferred <- paste0("h_", make.names(c(depvar, colnames(regressors))))
  transformed_names <- make.unique(c(names(sample_data), preferred))
  transformed_names <- transformed_names[ncol(sample_data) + seq_along(preferred)]
  transformed_data <- sample_data
  transformed_data[transformed_names] <- as.data.frame(transformed)
  rhs <- if (ncol(regressors)) transformed_names[-1L] else "1"
  transformed_formula <- stats::reformulate(rhs, transformed_names[1L], intercept = has_intercept)
  transformed_model <- stats::lm(transformed_formula, data = transformed_data, na.action = stats::na.fail)
  infer <- function(model, fitted_data) {
    scpcR::scpc(model, data = fitted_data, lon = lon, lat = lat,
                coords_euclidean = coords_euclidean, avc = avc, uncond = uncond,
                cvs = cvs, backend = backend_resolved, nthreads = nthreads_resolved)
  }
  levels_scpc <- infer(levels_model, data)
  transformed_scpc <- infer(transformed_model, transformed_data)
  stage_backends <- c(
    stats::setNames(rep(operation_backend, 5L), c(kinds, "transform")),
    levels_scpc = levels_scpc$backend, transformed_scpc = transformed_scpc$backend
  )
  term_names <- stats::setNames(
    c("(Intercept)", colnames(regressors)), c("(Intercept)", transformed_names[-1L])
  )
  structure(list(
    tests = tests,
    fits = list(
      levels = list(model = levels_model, scpc = levels_scpc),
      transformed = list(model = transformed_model, scpc = transformed_scpc, term_names = term_names)
    ),
    backend = if (all(stage_backends == "rust")) "rust" else "fallback",
    stage_backends = stage_backends,
    fallback_reasons = list(
      diagnostics_and_transform = reason, levels_scpc = levels_scpc$fallback_reason,
      transformed_scpc = transformed_scpc$fallback_reason
    ),
    tied = any(tied) || isTRUE(levels_scpc$tied) || isTRUE(transformed_scpc$tied)
  ), class = "spur_result")
}

#' @export
print.spur_result <- function(x, ...) {
  cat(.render_spur_summary(x), "\n", sep = "")
  invisible(x)
}

#' @export
summary.spur_result <- function(object, ...) {
  cat(.render_spur_summary(object), "\n", sep = "")
  invisible(object)
}
