#' @keywords internal
.prefix_formula_expr <- function(expr, prefix) {
  if (is.symbol(expr)) {
    name <- as.character(expr)
    if (name %in% c("0", "1")) {
      return(expr)
    }
    return(as.name(paste0(prefix, name)))
  }

  if (is.call(expr)) {
    # keep operators and function names as they are, and only rename variables.
    head <- expr[[1L]]
    args <- as.list(expr)[-1L]
    new_args <- lapply(args, .prefix_formula_expr, prefix = prefix)
    return(as.call(c(list(head), new_args)))
  }

  expr
}

#' @keywords internal
.rewrite_formula_with_prefix <- function(formula, prefix) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` for spur() must be two-sided, e.g. `y ~ x1 + x2`.")
  }

  lhs <- .prefix_formula_expr(formula[[2L]], prefix)
  rhs <- .prefix_formula_expr(formula[[3L]], prefix)

  stats::as.formula(call("~", lhs, rhs), env = environment(formula))
}

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
  transformed_names <- .display_term_name(rownames(transformed_stats))
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

  depvar <- all.vars(formula(x$fits$levels$model))[1]
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
#' @param q Number of low-frequency weighted averages.
#' @param nrep Number of Monte Carlo draws.
#' @param lon Optional longitude column name.
#' @param lat Optional latitude column name.
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param seed Integer seed passed to the SPUR test steps.
#' @param avc SCPC tuning parameter passed to \code{scpcR::scpc()}.
#' @param uncond Passed through to \code{scpcR::scpc()}.
#' @param cvs Passed through to \code{scpcR::scpc()}.
#' @return A list of class \code{"spur_result"} with \code{tests} and \code{fits}.
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
                 cvs = FALSE) {
  .validate_common_args(data = data, q = q, nrep = nrep)

  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` for spur() must be two-sided, e.g. `y ~ x1 + x2`.")
  }

  depvar <- all.vars(formula[[2L]])
  if (length(depvar) != 1L) {
    stop("`formula` for spur() must have exactly one dependent variable.")
  }

  depvar_formula <- stats::reformulate("1", response = depvar)

  i0 <- spurtest_i0(
    depvar_formula,
    data = data,
    q = q,
    nrep = nrep,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    seed = seed
  )
  i1 <- spurtest_i1(
    depvar_formula,
    data = data,
    q = q,
    nrep = nrep,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    seed = seed
  )
  i0resid <- spurtest_i0resid(
    formula,
    data = data,
    q = q,
    nrep = nrep,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    seed = seed
  )
  i1resid <- spurtest_i1resid(
    formula,
    data = data,
    q = q,
    nrep = nrep,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    seed = seed
  )

  levels_model <- stats::lm(formula = formula, data = data)
  levels_scpc <- scpcR::scpc(
    levels_model,
    data = data,
    lon = lon,
    lat = lat,
    coord_euclidean = coords_euclidean,
    avc = avc,
    uncond = uncond,
    cvs = cvs
  )

  transformed_data <- spurtransform(
    formula,
    data = data,
    prefix = "h_",
    transformation = "lbmgls",
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean
  )
  transformed_formula <- .rewrite_formula_with_prefix(formula, "h_")
  transformed_model <- stats::lm(formula = transformed_formula, data = transformed_data)
  transformed_scpc <- scpcR::scpc(
    transformed_model,
    data = transformed_data,
    lon = lon,
    lat = lat,
    coord_euclidean = coords_euclidean,
    avc = avc,
    uncond = uncond,
    cvs = cvs
  )

  structure(
    list(
      tests = list(
        i0 = i0,
        i1 = i1,
        i0resid = i0resid,
        i1resid = i1resid
      ),
      fits = list(
        levels = list(
          model = levels_model,
          scpc = levels_scpc
        ),
        transformed = list(
          model = transformed_model,
          scpc = transformed_scpc
        )
      )
    ),
    class = "spur_result"
  )
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
