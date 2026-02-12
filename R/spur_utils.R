#' @keywords internal
.detect_s_coord_columns <- function(data) {
  s_cols <- grep("^s_[0-9]+$", names(data), value = TRUE)
  if (!length(s_cols)) {
    return(character(0))
  }

  s_nums <- as.integer(sub("^s_", "", s_cols))
  ord <- order(s_nums)
  s_cols <- s_cols[ord]
  s_nums <- s_nums[ord]

  expected <- seq_len(length(s_nums))
  if (!identical(s_nums, expected)) {
    stop("s_* coordinate variables must be continuously numbered starting at s_1.")
  }

  s_cols
}

#' @keywords internal
.resolve_coord_cols <- function(data, latlong, coord_cols = NULL, lat = "lat", lon = "lon") {
  if (!is.null(coord_cols)) {
    if (!all(coord_cols %in% names(data))) {
      stop("Some coordinate columns were not found in data.")
    }
    cols <- coord_cols
  } else if (all(c(lat, lon) %in% names(data))) {
    cols <- c(lat, lon)
  } else {
    cols <- .detect_s_coord_columns(data)
    if (!length(cols)) {
      stop("Could not find coordinate columns. Supply coord_cols, or provide lat/lon or s_* columns.")
    }
  }

  if (isTRUE(latlong) && length(cols) != 2L) {
    stop("With latlong=TRUE there must be exactly 2 coordinate columns.")
  }

  cols
}

#' @keywords internal
.normalized_distmat <- function(coords, latlong) {
  distmat <- if (isTRUE(latlong)) {
    get_distmat_lat_lon(coords)
  } else {
    get_distmat_euclidean(coords)
  }
  distmat / max(distmat)
}

#' @keywords internal
.extract_coords <- function(data, use_rows, coord_cols, latlong) {
  coords <- as.matrix(data[use_rows, coord_cols, drop = FALSE])
  get_s_matrix(coords = coords, latlong = latlong)
}

#' @keywords internal
.format_spur_test_output <- function(title, stat_label, teststat, pvalue) {
  cat(title, "\n", sep = "")
  cat("---------------------------------------\n")
  cat(sprintf("Test Statistic (%s) : %9.4f\n", stat_label, teststat))
  cat(sprintf("P-value               : %9.4f\n", pvalue))
  cat("---------------------------------------\n")
}

#' @keywords internal
.draw_emat <- function(q, nrep) {
  matrix(stats::rnorm(q * nrep), nrow = q, ncol = nrep)
}
