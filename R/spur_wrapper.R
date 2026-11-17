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
#' @param q Number of low-frequency weighted averages (at least 2).
#' @param nrep Number of Monte Carlo draws.
#' @param lon Optional longitude column name.
#' @param lat Optional latitude column name.
#' @param coords_euclidean Optional Euclidean coordinate column names.
#' @param seed Integer seed passed to the SPUR test steps.
#' @param avc SCPC tuning parameter passed to \code{scpcR::scpc()}.
#' @param uncond Passed through to \code{scpcR::scpc()}.
#' @param cvs Passed through to \code{scpcR::scpc()}.
#' @param backend Backend selector: \code{"auto"} (default), \code{"rust"}, or
#'   \code{"fallback"}.
#' @param nthreads Optional positive integer thread count for parallel Rust
#'   execution.
#' @return A list of class \code{"spur_result"} with \code{tests} and \code{fits}.
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

  depvar_formula <- stats::reformulate("1", response = depvar)
  parsed <- .parse_residual_formula(formula = formula, data = data, fn_name = "spur()")
  .validate_spur_y(parsed$Y, name = depvar)
  q_x <- .orthonormal_design_basis(parsed$X_in)
  resid_dof <- nrow(parsed$Y) - ncol(q_x)
  if (q > resid_dof) {
    stop(sprintf(
      "`q` (%d) cannot exceed residual degrees of freedom (n - ncol(X) = %d).",
      q, resid_dof
    ))
  }
  rust_loaded <- .require_spur_rust_if_requested(backend_resolved, "wrap__rust_run_spur_pipeline")

  if (!identical(backend_resolved, "fallback") && rust_loaded) {
    tr_vars <- .parse_transform_formula(formula = formula, data = data)
    all_tr_complete <- length(tr_vars) > 0L && all(stats::complete.cases(data[tr_vars]))
    if (isTRUE(all(parsed$use_rows)) && isTRUE(all_tr_complete) && nrow(data) >= 2L) {
      coord_info <- .resolve_spur_coords(
        data = data,
        use_rows = parsed$use_rows,
        lon = lon,
        lat = lat,
        coords_euclidean = coords_euclidean,
        verbose = FALSE
      )
      if (!is.null(seed)) set.seed(seed)
      emat <- .draw_emat(q = q, nrep = nrep)
      e_m <- as.matrix(emat)
      storage.mode(e_m) <- "double"
      x_m <- as.matrix(parsed$X_in)
      storage.mode(x_m) <- "double"
      tr_mat_in <- as.matrix(data[, tr_vars, drop = FALSE])
      storage.mode(tr_mat_in) <- "double"

      pipe <- .call_spur_rust(
        "wrap__rust_run_spur_pipeline",
        .swap_coords_for_rust(coord_info$coords, coord_info$latlong),
        isTRUE(coord_info$latlong),
        as.double(parsed$Y),
        x_m,
        e_m,
        tr_mat_in,
        TRUE,
        nthreads_resolved
      )
      any_tied <- isTRUE(pipe$i0$tied) || isTRUE(pipe$i0resid$tied)
      if (!any_tied || identical(backend_resolved, "rust")) {
        i0 <- .format_spur_test_from_full(
          pipe$i0,
          test_type = "I(0)",
          title = "Spatial I(0) Test Results",
          stat_label = "LFST",
          cv = as.numeric(pipe$i0$cvalue),
          verbose = FALSE,
          backend = "rust",
          tied = isTRUE(pipe$i0$tied)
        )
        i1 <- .format_spur_test_from_full(
          pipe$i1,
          test_type = "I(1)",
          title = "Spatial I(1) Test Results",
          stat_label = "LFUR",
          cv = as.numeric(pipe$i1$cv_vec),
          verbose = FALSE,
          backend = "rust",
          tied = isTRUE(pipe$i1$tied)
        )
        i0resid <- .format_spur_test_from_full(
          pipe$i0resid,
          test_type = "I(0) Residual",
          title = "Spatial I(0) Test Results for Residuals",
          stat_label = "LFST",
          cv = as.numeric(pipe$i0resid$cvalue),
          verbose = FALSE,
          backend = "rust",
          tied = isTRUE(pipe$i0resid$tied)
        )
        i1resid <- .format_spur_test_from_full(
          pipe$i1resid,
          test_type = "I(1) Residual",
          title = "Spatial I(1) Test Results for Residuals",
          stat_label = "LFUR",
          cv = as.numeric(pipe$i1resid$cv_vec),
          verbose = FALSE,
          backend = "rust",
          tied = isTRUE(pipe$i1resid$tied)
        )

        transformed_data <- data
        tr_out <- pipe$transformed_vars
        for (j in seq_along(tr_vars)) {
          transformed_data[[paste0("h_", tr_vars[[j]])]] <- as.numeric(tr_out[, j])
        }

        levels_model <- stats::lm(formula = formula, data = data)
        levels_scpc <- scpcR::scpc(
          levels_model,
          data = data,
          lon = lon,
          lat = lat,
          coords_euclidean = coords_euclidean,
          avc = avc,
          uncond = uncond,
          cvs = cvs,
          backend = backend_resolved,
          nthreads = nthreads_resolved
        )

        transformed_formula <- .rewrite_formula_with_prefix(formula, "h_")
        transformed_model <- stats::lm(formula = transformed_formula, data = transformed_data)
        transformed_scpc <- scpcR::scpc(
          transformed_model,
          data = transformed_data,
          lon = lon,
          lat = lat,
          coords_euclidean = coords_euclidean,
          avc = avc,
          uncond = uncond,
          cvs = cvs,
          backend = backend_resolved,
          nthreads = nthreads_resolved
        )

        return(structure(
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
            ),
            backend = "rust",
            tied = any_tied
          ),
          class = "spur_result"
        ))
      }
    }
  }

  i0 <- spurtest_i0(
    depvar_formula,
    data = data,
    q = q,
    nrep = nrep,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    seed = seed,
    backend = backend_resolved,
    nthreads = nthreads_resolved
  )
  i1 <- spurtest_i1(
    depvar_formula,
    data = data,
    q = q,
    nrep = nrep,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    seed = seed,
    backend = backend_resolved,
    nthreads = nthreads_resolved
  )
  i0resid <- spurtest_i0resid(
    formula,
    data = data,
    q = q,
    nrep = nrep,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    seed = seed,
    backend = backend_resolved,
    nthreads = nthreads_resolved
  )
  i1resid <- spurtest_i1resid(
    formula,
    data = data,
    q = q,
    nrep = nrep,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    seed = seed,
    backend = backend_resolved,
    nthreads = nthreads_resolved
  )

  levels_model <- stats::lm(formula = formula, data = data)
  levels_scpc <- scpcR::scpc(
    levels_model,
    data = data,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    avc = avc,
    uncond = uncond,
    cvs = cvs,
    backend = backend_resolved,
    nthreads = nthreads_resolved
  )

  transformed_data <- spurtransform(
    formula,
    data = data,
    prefix = "h_",
    transformation = "lbmgls",
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    backend = backend_resolved,
    nthreads = nthreads_resolved
  )
  transformed_formula <- .rewrite_formula_with_prefix(formula, "h_")
  transformed_model <- stats::lm(formula = transformed_formula, data = transformed_data)
  transformed_scpc <- scpcR::scpc(
    transformed_model,
    data = transformed_data,
    lon = lon,
    lat = lat,
    coords_euclidean = coords_euclidean,
    avc = avc,
    uncond = uncond,
    cvs = cvs,
    backend = backend_resolved,
    nthreads = nthreads_resolved
  )

  all_rust <- identical(i0$backend, "rust") &&
    identical(i1$backend, "rust") &&
    identical(i0resid$backend, "rust") &&
    identical(i1resid$backend, "rust")
  any_tied <- isTRUE(i0$tied) || isTRUE(i1$tied) || isTRUE(i0resid$tied) || isTRUE(i1resid$tied)

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
      ),
      backend = if (all_rust) "rust" else "fallback",
      tied = any_tied
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
