#' Transform Variables Using SPUR Spatial Transformations
#'
#' R equivalent of Stata `spurtransform`.
#'
#' @param data A data frame containing variables and coordinates.
#' @param vars Character vector of numeric variable names to transform.
#' @param prefix Prefix for transformed variable names.
#' @param transformation One of `"lbmgls"`, `"nn"`, `"iso"`, or `"cluster"`.
#' @param radius Radius used by `transformation = "iso"`.
#' @param clustvar Cluster variable used by `transformation = "cluster"`.
#' @param latlong Logical; if `TRUE`, coordinates are latitude/longitude.
#' @param lat Name of latitude column when `coord_cols` is `NULL`.
#' @param lon Name of longitude column when `coord_cols` is `NULL`.
#' @param coord_cols Optional explicit coordinate column names.
#' @param replace Logical; allow overwriting existing output columns.
#' @param separately Logical; handle missingness separately by variable.
#' @return A data frame with transformed variables appended.
#' @export
spur_transform <- function(data,
                           vars,
                           prefix = "tr",
                           transformation = "lbmgls",
                           radius = NULL,
                           clustvar = NULL,
                           latlong = TRUE,
                           lat = "lat",
                           lon = "lon",
                           coord_cols = NULL,
                           replace = FALSE,
                           separately = FALSE) {
  if (!length(vars)) {
    stop("vars must contain at least one variable name.")
  }
  if (!all(vars %in% names(data))) {
    stop("At least one variable in vars was not found in data.")
  }
  if (!all(vapply(data[vars], is.numeric, logical(1)))) {
    stop("All variables in vars must be numeric.")
  }

  allowed <- c("lbmgls", "nn", "iso", "cluster")
  if (!transformation %in% allowed) {
    stop("Invalid transformation.")
  }

  if (transformation != "iso" && !is.null(radius)) {
    stop("Option radius only allowed with transformation(iso).")
  }
  if (transformation == "iso" && is.null(radius)) {
    stop("Radius missing.")
  }
  if (transformation == "iso" && radius <= 0) {
    stop("Radius must be positive.")
  }

  if (transformation != "cluster" && !is.null(clustvar)) {
    stop("Option clustvar only allowed with transformation(cluster).")
  }
  if (transformation == "cluster" && is.null(clustvar)) {
    stop("Clustvar missing.")
  }
  if (!is.null(clustvar) && !clustvar %in% names(data)) {
    stop("clustvar was not found in data.")
  }

  if (!is.character(prefix) || length(prefix) != 1L || !nzchar(prefix) || make.names(prefix) != prefix) {
    stop("prefix must be a valid name prefix.")
  }

  coord_cols <- .resolve_coord_cols(
    data = data,
    latlong = latlong,
    coord_cols = coord_cols,
    lat = lat,
    lon = lon
  )

  base_sample <- if (isTRUE(separately)) {
    rep(TRUE, nrow(data))
  } else {
    stats::complete.cases(data[vars])
  }

  cluster_group <- if (transformation == "cluster") {
    as.numeric(factor(data[[clustvar]]))
  } else {
    NULL
  }

  result <- data

  for (var in vars) {
    out_var <- paste0(prefix, var)
    if (out_var %in% names(result) && !isTRUE(replace)) {
      stop(sprintf("%s already defined or invalid name", out_var))
    }

    use_rows <- base_sample
    if (isTRUE(separately)) {
      use_rows <- use_rows & !is.na(data[[var]])
    }
    if (transformation == "cluster") {
      use_rows <- use_rows & !is.na(data[[clustvar]])
    }

    if (!any(use_rows)) {
      warning(sprintf("No valid observations for variable: %s", var))
      next
    }

    s <- .extract_coords(
      data = data,
      use_rows = use_rows,
      coord_cols = coord_cols,
      latlong = latlong
    )

    cluster_vals <- if (transformation == "cluster") {
      get_cluster_matrix(cluster_values = cluster_group, use_rows = use_rows)
    } else {
      NULL
    }

    H <- make_transform(
      s = s,
      transformation = transformation,
      radius = radius,
      cluster = cluster_vals,
      latlong = latlong
    )

    transformed <- spur_transform_vector(
      data = data,
      varname = var,
      H = H,
      transformation = transformation,
      use_rows = use_rows
    )

    result_vec <- rep(NA_real_, nrow(data))
    result_vec[use_rows] <- transformed
    result[[out_var]] <- result_vec
  }

  result
}
