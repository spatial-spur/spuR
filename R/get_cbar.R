#' Compute cbar Value
#'
#' Determines the value of \code{cbar} such that \eqn{mean(exp(-cbar \times vd))} is approximately equal to
#' a given \code{rhobar}, where \code{vd} is the vectorized lower-triangular portion (excluding the diagonal)
#' of the distance matrix \code{distmat}. This function uses iterative adjustments and bisection to find
#' an appropriate value.
#'
#' @param rhobar A numeric scalar representing the target value.
#' @param distmat A numeric matrix of distances.
#' @return A numeric scalar \code{cbar} satisfying the desired property.
#' @keywords internal
get_cbar <- function(rhobar, distmat) {
  .make_spur_decay_solver(distmat)(rhobar)
}

.make_spur_decay_solver <- function(distmat) {
  distances <- lvech(distmat)
  averages <- new.env(parent = emptyenv())
  roots <- new.env(parent = emptyenv())
  average <- function(c) {
    key <- sprintf("%a", c)
    out <- averages[[key]]
    if (is.null(out)) {
      out <- mean(exp(-c * distances))
      averages[[key]] <- out
    }
    out
  }

  function(rho) {
    key <- sprintf("%a", rho)
    out <- roots[[key]]
    if (is.null(out)) {
      out <- .solve_spur_decay(rho, average)
      roots[[key]] <- out
    }
    out
  }
}

.solve_spur_decay <- function(rhobar, average) {
  
  # Initialize c0 and c1
  c0 <- 10
  c1 <- 10
  
  # Adjust c0 so that mean(exp(-c0 * vd)) > rhobar
  i1 <- FALSE
  jj <- 0
  while (!i1) {
    v <- average(c0)
    i1 <- (v > rhobar)
    if (!i1) {
      c1 <- c0
      c0 <- c0 / 2
      jj <- jj + 1
    }
    if (jj > 500) {
      stop("rhobar too large")
    }
  }
  
  # Adjust c1 so that mean(exp(-c1 * vd)) < rhobar
  i1 <- FALSE
  jj <- 0
  while (!i1) {
    v <- average(c1)
    i1 <- (v < rhobar)
    if (!i1) {
      c0 <- c1
      c1 <- 2 * c1
      jj <- jj + 1
    }
    if (c1 > 10000) {  # safeguard to break if c1 becomes unreasonably large
      i1 <- TRUE
    }
    if (jj > 500) {
      stop("rhobar too small")
    }
  }
  
  # Bisection method to determine cbar
  while ((c1 - c0) > 0.001) {
    cm <- sqrt(c0 * c1)
    v <- average(cm)
    if (v < rhobar) {
      c1 <- cm
    } else {
      c0 <- cm
    }
  }
  
  cbar <- sqrt(c0 * c1)
  return(cbar)
}
