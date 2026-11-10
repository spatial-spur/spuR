\name{spur_set_num_threads}
\alias{spur_set_num_threads}
\alias{spur_get_num_threads}
\title{Configure parallel worker threads}
\usage{
spur_set_num_threads(num_threads = "auto")
spur_get_num_threads()
}
\arguments{
\item{num_threads}{A positive integer, \code{NULL}, or \code{"auto"}.}
}
\value{The setter invisibly returns the previous package option, including \code{NULL} for automatic detection. The getter returns the resolved thread count.}
\description{Controls parallel worker threads without changing process-wide environment variables.}
\details{Resolution uses the explicit argument, \code{spuR.threads} option, \code{SPUR_NUM_THREADS} environment variable, \code{mc.cores} option, \code{RCPP_PARALLEL_NUM_THREADS}, then hardware detection. A true \code{_R_CHECK_LIMIT_CORES_} flag caps every result at two threads. BLAS thread settings are independent of worker threads.}
