#' Set and get thread count for spuR
#'
#' Configure the number of threads used for parallel computation in \pkg{spuR}.
#'
#' @param num_threads Integer >= 1 specifying the number of threads, or \code{"auto"}
#'   / \code{NULL} to restore automatic thread detection.
#'
#' @return \code{spur_set_num_threads} invisibly returns the previous package option (NULL for automatic detection).
#'   \code{spur_get_num_threads} returns the currently configured or active thread count.
#'
#' @details
#' \pkg{spuR} supports multithreading via \pkg{RcppParallel} for grid evaluations,
#' Monte Carlo persistence intervals, and distance calculations.
#'
#' Thread resolution follows this priority:
#' \enumerate{
#'   \item Explicit argument \code{nthreads} in \code{\link{spur}} and test functions
#'   \item Package option \code{options(spuR.threads = N)} (set directly or via \code{spur_set_num_threads(N)})
#'   \item Environment variable \code{SPUR_NUM_THREADS}
#'   \item General R parallel option \code{options(mc.cores = N)}
#'   \item \pkg{RcppParallel} environment variable \code{RCPP_PARALLEL_NUM_THREADS}
#'   \item Automatic hardware detection (all available logical cores)
#' }
#'
#' The CRAN \code{_R_CHECK_LIMIT_CORES_} limit caps every resolved setting at two threads.
#'
#' Setting \code{num_threads = 1} (or \code{options(spuR.threads = 1)}) completely bypasses
#' background task submission, ensuring pure single-threaded execution. This is recommended
#' when parallelising outer loops with \code{parallel::mclapply}, \pkg{future}, or when running
#' on shared cluster nodes.
#'
#' @examples
#' \dontrun{
#' # Check currently active threads
#' spur_get_num_threads()
#'
#' # Force single-threaded execution (e.g. inside mclapply or cluster jobs)
#' old_th <- spur_set_num_threads(1)
#'
#' # Restore previous setting
#' spur_set_num_threads(old_th)
#'
#' # Or restore default automatic threading
#' spur_set_num_threads("auto")
#' }
#'
#' @export
spur_set_num_threads <- function(num_threads = "auto") {
  old <- getOption("spuR.threads")
  if (is.null(num_threads) || identical(num_threads, "auto")) {
    options(spuR.threads = NULL)
  } else {
    if (!is.numeric(num_threads) || length(num_threads) != 1L ||
        !is.finite(num_threads) || num_threads < 1 ||
        num_threads > .Machine$integer.max || num_threads != floor(num_threads)) {
      stop("`num_threads` must be a single positive integer, NULL, or \"auto\".")
    }
    options(spuR.threads = as.integer(num_threads))
  }
  invisible(old)
}

#' @rdname spur_set_num_threads
#' @export
spur_get_num_threads <- function() {
  .get_spur_threads()
}

.get_spur_threads <- function(nthreads = NULL) {
  valid_count <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value >= 1 && value <= .Machine$integer.max && value == floor(value)
  }
  if (!is.null(nthreads) && !valid_count(nthreads)) {
    stop("`nthreads` must be a single positive integer or NULL.")
  }
  count <- nthreads
  if (is.null(count)) count <- getOption("spuR.threads")
  if (!is.null(count) && !valid_count(count)) {
    stop("`options(spuR.threads)` must be a single positive integer.")
  }
  if (is.null(count)) {
    env <- Sys.getenv("SPUR_NUM_THREADS", "")
    if (nzchar(env)) {
      count <- suppressWarnings(as.numeric(env))
      if (!valid_count(count)) stop("`SPUR_NUM_THREADS` must be a positive integer.")
    }
  }
  if (is.null(count)) count <- getOption("mc.cores")
  if (!is.null(count) && !valid_count(count)) {
    stop("`options(mc.cores)` must be a single positive integer.")
  }
  if (is.null(count)) {
    env <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "")
    if (nzchar(env)) {
      count <- suppressWarnings(as.numeric(env))
      if (!valid_count(count)) stop("`RCPP_PARALLEL_NUM_THREADS` must be a positive integer.")
    }
  }
  if (is.null(count)) count <- RcppParallel::defaultNumThreads()
  cran_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(cran_limit) && !cran_limit %in% c("false", "no", "0")) {
    count <- min(count, 2L)
  }
  as.integer(count)
}
