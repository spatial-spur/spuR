
#' Transform Variable Using Spatial Transformation Matrix
#'
#' This function transforms a variable using a spatial transformation matrix.
#' Direct R equivalent of the Stata transform() function.
#'
#' @param data Data frame containing the variable
#' @param varname Name of the variable to transform
#' @param H Transformation matrix
#' @param transformation Type of transformation
#' @param use_rows Logical vector indicating which rows to use (equivalent to touse)
#'
#' @return A numeric vector of transformed values
#' @export
transform <- function(data, varname, H, transformation, use_rows = NULL) {
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
