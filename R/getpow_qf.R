#' Compute Power Using Quadratic Forms
#'
#' Computes the power of the test by comparing likelihood ratios derived from two sets
#' of covariance matrices. It computes:
#' \deqn{
#'   \text{lr\_o} = \frac{qe}{qa\_o} \quad \text{and} \quad \text{lr\_a} = \frac{qo\_a}{qe},
#' }
#' where \eqn{qe} is the column sum of \eqn{e^2}. The critical value is the 95th percentile of
#' \eqn{lr\_o}, and power is computed as the proportion of \eqn{lr\_a} exceeding that critical value.
#'
#' @param om0 A numeric matrix representing the null covariance matrix.
#' @param om1 A numeric matrix representing the alternative covariance matrix.
#' @param e A numeric matrix of test statistics or error terms.
#' @return A numeric scalar representing the computed power.
#' @export
getpow_qf <- function(om0, om1, e) {
  # Compute inverses of om0 and om1
  om0i <- solve(om0)
  om1i <- solve(om1)
  
  # Cholesky decompositions
  ch_om0 <- chol(om0)
  ch_om1 <- chol(om1)
  ch_om0i <- chol(om0i)
  ch_om1i <- chol(om1i)
  
  # Compute transformation matrices
  ho <- ch_om1i %*% t(ch_om0)
  ha <- ch_om0i %*% t(ch_om1)
  
  # Compute quadratic forms
  qe <- colSums(e^2)
  ya_o <- ho %*% e
  yo_a <- ha %*% e
  qa_o <- colSums(ya_o^2)
  qo_a <- colSums(yo_a^2)
  
  # Likelihood ratios
  lr_o <- qe / qa_o
  lr_a <- qo_a / qe
  
  # Determine the 95th percentile critical value of lr_o
  cv <- as.numeric(quantile(lr_o, probs = 0.95))
  
  # Compute power as the proportion of lr_a values exceeding the critical value
  pow <- mean(lr_a > cv)
  
  return(pow)
}
