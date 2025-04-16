#' SPUR I(0) Test Wrapper
#'
#' A user-friendly wrapper for performing the spatial I(0) test. The user only needs to supply
#' a data frame and the name of the outcome variable. The function assumes that the data frame also
#' contains spatial coordinates in columns `"lat"` and `"lon"` (if \code{latlong} is \code{TRUE}).
#'
#' @param data A data frame containing the data.
#' @param var A character string specifying the name of the numeric outcome variable to test.
#' @param q An integer specifying the number of eigenvectors to use (default is 15). # nolint: line_length_linter.
#' @param nrep An integer specifying the number of simulation replications for generating \code{emat} (default is 100000).
#' @param latlong Logical; if \code{TRUE} the data frame must have columns \code{"lat"} and \code{"lon"} to build the distance matrix.
#'
#' @return A list with elements:
#' \describe{
#'   \item{pvalue}{Scalar p-value of the test.}
#'   \item{teststat}{Test statistic (likelihood ratio) for the outcome variable.}
#'   \item{ha_param}{The Ha parameter that yields roughly 50\% power.}
#'   \item{cv}{A vector of critical values.}
#'   \item{full}{The complete output from the core SPUR I(0) test function.}
#' }
#'
#' @export
spur_i0 <- function(var, q = 15, nrep = 100000, latlong = TRUE, data) {
  # Check that the outcome variable exists in the data frame
  if (!var %in% names(data)) {
    stop("The variable specified does not exist in the data frame.")
  }

  # Extract the outcome variable as a matrix (n x 1)
  Y <- as.matrix(data[[var]])

  # Build the distance matrix from coordinates.
  # If latlong is TRUE, we expect the data frame to have "lat" and "lon" columns.
  if (latlong) {
    if (!all(c("lat", "lon") %in% names(data))) {
      stop("When latlong=TRUE, the data frame must contain 'lat' and 'lon' columns.")
    }
    coords <- as.matrix(data[, c("lat", "lon")])
  } else {
    stop("Currently only latlong=TRUE is supported.")
  }

  # Construct the distance matrix (using Euclidean distances)
  distmat <- as.matrix(dist(coords))

  # Generate the simulation matrix 'emat': dimensions q x nrep.
  # This mimics the Stata rnormal(q, nrep, 0, 1) call.
  set.seed(123)  # for reproducibility (optional)
  emat <- matrix(rnorm(q * nrep), nrow = q, ncol = nrep)

  # Call the core SPUR I(0) test function.
  # Note: spur_i0_test() implements the test exactly as in our previous translation.
  res <- spur_i0_test(Y, distmat, emat, q = q)

  # Extract key outputs.
  # For a single outcome variable, LR should be of length 1.
  teststat <- res$LR[1]
  # Here we assume pvalue is a vector (e.g., for different significance levels).
  pvalue <- res$pvalue[1]
  ha_param <- res$ha_parm
  cv <- res$cvalue

  # Display results (similar to the Stata display output)
  cat("Spatial I(0) Test Results\n")
  cat("---------------------------------------\n")
  cat(sprintf("Test Statistic (LFST) : %9.4f\n", teststat))
  cat(sprintf("P-value               : %9.4f\n", pvalue))
  cat("---------------------------------------\n")

  # Return the results in a list.
  return(list(
    pvalue = pvalue,
    teststat = teststat,
    ha_param = ha_param,
    cv = cv,
    full = res
  ))
}
