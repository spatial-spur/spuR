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

For this development checkout, install from the workspace root with Cargo on your PATH.
Set `SPATIAL_RUST=0` for a pure R build; `backend = "fallback"` selects pure R at
runtime in either build.

```sh
R CMD INSTALL scpcR
R CMD INSTALL spuR
```

For a standalone source archive that can build Rust without network access, run
`bash spatial-rust/scripts/build_r_source.sh scpcR ../r-source` and
`bash spatial-rust/scripts/build_r_source.sh spuR ../r-source`. Install the
scpcR archive before the spuR archive.

On Windows install the Rust target matching Rtools and your R installation.
The build uses `x86_64-pc-windows-gnu` for x86-64 R or
`aarch64-pc-windows-gnullvm` for ARM64 R. Set `WIN_TARGET` to override
the selection. Windows native execution still needs validation.

For the published version:

```r
remotes::install_github("spatial-spur/scpcR@v0.1.3")
remotes::install_github("spatial-spur/spuR@v0.1.2")
```

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

# Half-life confidence grid hull
spurhalflife(am ~ 1, data = spur_example, lon = "lon", lat = "lat", seed = 42)

# Spatial transformation — pass your analysis formula to transform all variables
out <- spurtransform(am ~ gini + fracblack, data = spur_example,
                     lon = "lon", lat = "lat", transformation = "lbmgls")
head(out[, c("h_am", "h_gini", "h_fracblack")])
```

## Shared Rust backend conventions

Use `backend = "auto"`, `"rust"`, or `"fallback"` independently of `nthreads`.
Setting `nthreads = 1` keeps a single-threaded Rust computation when Rust is selected.
The fallback remains available when Rust is not installed.

SPUR inference requires at least two retained directions, positive retained covariance
rank, and a nonzero projected response. Residual tests require
`q <= n - rank(X)` and a fitted design spanning a constant. For example,
all indicators of a factor span a constant. The combined `spur()` pipeline
requires a formula intercept, which R formulas include by default, so both fits have
a common formula convention. Uninformative covariance shapes and decay
targets at or below the fraction of co-located pairs produce explicit errors.
The I(0) calibration grid starts at `0.0001`, so the I(0) tests and `spur()`
require the co-located share of observation pairs to lie below that value.

Euclidean absolute half-lives use the original coordinate units. Geographic absolute
half-lives use metres, and `normdist = TRUE` reports fractions of the sample diameter.
The bounds summarize the accepted candidates on the finite Table 1 grid of
Müller & Watson (2024). The grid jumps from 3 to 100 sample diameters.
Acceptance at 100 yields `ci_u = NA` and the replication display `Inf`.
Inspect `full$hl_grid_ho` and `full$pv_mat` for the tested candidates
and their Monte Carlo p-values.

A fixed seed reproduces the host draws. Finite simulations can still differ
with the basis chosen by the backend, especially in tied or nearly tied
subspaces.

A spatial transformation can also change IV exclusion and first-stage
relevance. For a manual IV fit on transformed variables, verify both
conditions for the transformed moments. Zero same-site instrument–error
covariances alone are insufficient after mixing observations across sites.

## Documentation

Please refer to [the package documentation](https://spatial-spur.github.io/spuR/) for detailed information and other (R, Python, Stata) packages.

## Citation

If you use this package, please cite Becker, Boll and Voth 2026 and Müller & Watson 2024. See [CITATION.bib](CITATION.bib) for the BibTeX entries.

## References

Mueller, U. K. and Watson, M. W. (2024). Spatial Unit Roots and Spurious Regression.
*Econometrica*, 92(5), 1661-1695. doi:
[10.3982/ECTA21654](https://doi.org/10.3982/ECTA21654)
