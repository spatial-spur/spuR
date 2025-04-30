#' Calculate Distance Matrix for Latitude/Longitude Points
#'
#' Computes the matrix of great-circle distances between all pairs of locations
#' on a sphere using latitude and longitude coordinates.
#'
#' @param s A numeric matrix where rows represent locations and columns are [lat, lon].
#' @return A numeric matrix of pairwise distances in meters between locations.
#' @export
getdistmat_lat_lon <- function(s) {
  # Check if required package is available
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is needed for this function. Please install it with install.packages('sf')")
  }
  
  # Convert matrix to sf points
  # Note: sf expects coordinates as (longitude, latitude)
  #points_sf <- sf::st_as_sf(
  #  as.data.frame(s), 
  #  coords = c(2, 1),  # lon, lat order for sf 
  #  crs = 4326         # WGS84 - standard for lat/lon
  #)
    # Calculate distance matrix
  #distmat <- sf::st_distance(points_sf)

  #distmat <- as.matrix(geodist::geodist(s[, c("lat", "lon")], measure = "haversine"))
  #distmat <- distmat / max(distmat)  # Normalize the distance matrix
  
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
  
  # Convert from units object to numeric matrix
  #distmat <- as.numeric(distmat)
  #dim(distmat) <- c(nrow(s), nrow(s))
  
  return(distmat)
}