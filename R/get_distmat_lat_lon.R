#' Calculate Distance Matrix for Latitude/Longitude Points
#'
#' Computes the matrix of great-circle distances between all pairs of locations
#' on a sphere using latitude and longitude coordinates.
#'
#' @param s A numeric matrix where rows represent locations and columns are
#'   `lat` and `lon`.
#' @return A numeric matrix of pairwise great-circle distances in normalized
#'   units, matching the SPUR Mata `getdistmat()` implementation.
get_distmat_lat_lon <- function(s) {
  n <- nrow(s)
  distmat <- matrix(0, nrow = n, ncol = n)
  
  # Use exact same pi constant as Stata code
  c <- 3.14159265359/180
  
  # Ensure s is a matrix, not a data.frame
  s <- as.matrix(s)
  
  # Loop through each point exactly as in Stata
  for (i in 1:n) {
    # Calculate (s - s[i,])*0.5*c
    d <- 0.5 * c * (sweep(s, 2, s[i,], "-"))
    
    # Exact Stata haversine formula implementation
    distmat[, i] <- asin(sqrt(sin(d[, 1])^2 + 
                     cos(c * s[i, 1]) * 
                     (cos(c * s[, 1]) * sin(d[, 2])^2))) / 3.14159265359
  }
  
  return(distmat)
}

#' @rdname get_distmat_lat_lon
getdistmat_lat_lon <- function(s) {
  get_distmat_lat_lon(s = s)
}
