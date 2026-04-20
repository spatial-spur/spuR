
#' Transform One Variable with a Spatial Matrix
#'
#' Internal helper that applies a spatial transformation matrix to one variable.
#'
#' @param data Data frame containing the variable
#' @param varname Name of the variable to transform
#' @param H Transformation matrix
#' @param transformation Type of transformation
#' @param use_rows Logical vector indicating which rows to use (equivalent to touse)
#'
#' @return A numeric vector of transformed values
#' @keywords internal
spur_transform_vector <- function(data, varname, H, transformation, use_rows = NULL) {
  # Extract the variable
  y <- data[[varname]]
  
  # Apply subset if provided 
  if (!is.null(use_rows)) {
    y <- y[use_rows]
  }
  
  # Apply transformation
  hy <- as.vector(H %*% y)
  
  # For lbmgls transformation, demean the result
  if (transformation == "lbmgls") {
    hy <- hy - mean(hy)
  }
  
  return(hy)
}
