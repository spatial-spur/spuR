#' Extract Lower Triangle (Excluding Diagonal)
#'
#' Extracts the elements of a matrix that are strictly below the main diagonal
#' (excluding the diagonal) and returns them as a vector.
#'
#' @param S A numeric matrix.
#' @return A numeric vector containing the elements below the main diagonal.
lvech <- function(S) {
  # Extract elements strictly below the main diagonal (diag = FALSE)
  return(S[lower.tri(S, diag = FALSE)])
}