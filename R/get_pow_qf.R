#' @keywords internal
#' @noRd
.make_pow_qf_evaluator <- function(om0, e) {
  om0i <- solve(om0)
  ch_om0 <- chol(om0)
  ch_om0i <- chol(om0i)
  qe <- colSums(e^2)

  function(om1) {
    om1i <- solve(om1)

    # Mata::cholesky() returns a lower-triangular factor L with A = L %*% t(L),
    # while R::chol() returns upper-triangular U with A = t(U) %*% U.
    # In Mata code, factors are stored as cholesky(A)' (i.e., U), so use chol(A).
    ch_om1 <- chol(om1)
    ch_om1i <- chol(om1i)

    # Compute transformation matrices
    ho <- ch_om1i %*% t(ch_om0)
    ha <- ch_om0i %*% t(ch_om1)

    # Compute quadratic forms
    ya_o <- ho %*% e
    yo_a <- ha %*% e
    qa_o <- colSums(ya_o^2)
    qo_a <- colSums(yo_a^2)

    # Likelihood ratios
    lr_o <- qe / qa_o
    lr_a <- qo_a / qe

    # Determine the 95th percentile critical value of lr_o
    cv <- as.numeric(stats::quantile(lr_o, probs = 0.95))

    # Compute power as the proportion of lr_a values exceeding the critical value
    mean(lr_a > cv)
  }
}
