#' Get Spatial Coordinates Matrix
#'
#' Validates coordinate data and reports the distance calculation mode.
#' R translation of the Stata Mata 'get_s_matrix' function.
#'
#' @param coords A matrix or data frame containing coordinates.
#' @param latlong Logical; if TRUE, coordinates are latitude/longitude.
#' @return The validated coordinates matrix.
get_s_matrix <- function(coords, latlong = FALSE) {
  # Check for missing values (equivalent to: if(sum(rowmissing(s) :== 0) < ns))
  if (any(is.na(coords))) {
    stop("Missing value(s) in coordinate variables; aborting")
  }
  
  # Check for minimum number of observations (equivalent to: if(ns<5))
  ns <- nrow(coords)
  if (ns < 5) {
    stop("Too few locations found; aborting")
  }
  
  # Print diagnostic message (equivalent to the stata("disp...") calls)
  cat("Found", ns, "observations and", ncol(coords), "-dimensional locations\n")
  
  # Print additional diagnostics
  if (!latlong) {
    cat("Using Euclidean norm to compute distance between locations\n")
  } else {
    if (ncol(coords) != 2) {
      stop("With latlong option, there must only be 2 coordinate columns")
    } else {
      cat("Computing distances on surface of sphere treating coordinates as latitude and longitude\n")
    }
  }
  
  return(coords)
}
