#' @keywords internal
#' @noRd
.make_pow_qf_evaluator <- function(om0, e) {
  lower0 <- t(chol(om0))
  qe <- colSums(e^2)

  function(om1) {
    lower1 <- t(chol(om1))
    ho <- forwardsolve(lower1, lower0)
    ha <- forwardsolve(lower0, lower1)

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
