<p align="center">
  <img src="assets/logo.png" alt="SPUR logo">
</p>

# spuR

`spuR` implements diagnostics and correction methods for spatial unit roots developed by Müller & Watson 2024 in R.

This package implements low-frequency spatial unit root tests and spatial
transformations for cross-sectional data as developed by
[Müller & Watson 2024](https://doi.org/10.3982/ECTA21654). A practical
guide to these methods can be found in
[Becker, Boll and Voth 2026](https://journals.sagepub.com/doi/10.1177/1536867X261449932).

## Installation

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("spatial-spur/scpcR@v0.1.3")
remotes::install_github("spatial-spur/spuR@v0.1.2")
```

GitHub installation does not guarantee the declared `scpcR` version is
updated, so make sure both versions are up to date by installing the tagged
version explicitly.

## Usage

```r
library(spuR)
data(spur_example)

# Spatial I(0) test
spurtest_i0(am ~ 1, data = spur_example, lon = "lon", lat = "lat", seed = 42)

# Spatial I(1) test
spurtest_i1(am ~ 1, data = spur_example, lon = "lon", lat = "lat", seed = 42)

# Residual-based I(0) test
spurtest_i0resid(am ~ gini, data = spur_example, lon = "lon", lat = "lat", seed = 42)

# Half-life confidence interval
spurhalflife(am ~ 1, data = spur_example, lon = "lon", lat = "lat", seed = 42)

# Spatial transformation — pass your analysis formula to transform all variables
out <- spurtransform(am ~ gini + fracblack, data = spur_example,
                     lon = "lon", lat = "lat", transformation = "lbmgls")
head(out[, c("h_am", "h_gini", "h_fracblack")])
```

## Multithreading and Performance

`spuR` uses `RcppParallel` for multi-core acceleration in computationally intensive steps (Omega grid evaluation, Monte Carlo persistence intervals, and distance calculations). By default, it automatically uses available cores.

You can configure or restrict multithreading using:

```r
# Restrict to a specific thread count via helper function:
spur_set_num_threads(2)

# Or set via standard R package option:
options(spuR.threads = 1)

# Or set via environment variable:
Sys.setenv(SPUR_NUM_THREADS = "1")

# Check currently active thread count:
spur_get_num_threads()

# Restore automatic multi-core threading:
spur_set_num_threads("auto")
```

> **Tip for cluster & batch workflows:** If you are running `spur()` or `spurtest_*()` inside parallel loops (e.g. `parallel::mclapply`, `future`, or a shared cluster node), set `spur_set_num_threads(1)` or `options(spuR.threads = 1)` to prevent CPU oversubscription. Setting thread count to `1` completely bypasses background task submission, ensuring pure single-threaded execution without overhead or thread creation issues.

## Documentation

Please refer to [the package documentation](https://spatial-spur.github.io/spuR/) for detailed information and other (R, Python, Stata) packages.

## Citation

If you use this package, please cite Becker, Boll and Voth 2026 and Müller & Watson 2024. See [CITATION.bib](CITATION.bib) for the BibTeX entries.

## References

Mueller, U. K. and Watson, M. W. (2024). Spatial Unit Roots and Spurious Regression.
*Econometrica*, 92(5), 1661-1695. doi:
[10.3982/ECTA21654](https://doi.org/10.3982/ECTA21654)
