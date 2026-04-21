# spuR

Spatial Unit Root Tests and Transformations for R.

This package implements low-frequency spatial unit root tests and spatial
transformations for cross-sectional data as developed by
[Mueller and Watson (2024)](https://doi.org/10.3982/ECTA20682). A practical
guide to these methods can be found in
[Becker, Boll and Voth (2025)](https://warwick.ac.uk/fac/soc/economics/research/workingpapers/2025/twerp_1541-_becker.pdf).

## Installation

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("spatial-spur/spuR")
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

# Half-life confidence interval
spurhalflife(am ~ 1, data = spur_example, lon = "lon", lat = "lat", seed = 42)

# Spatial transformation — pass your analysis formula to transform all variables
out <- spurtransform(am ~ gini + fracblack, data = spur_example,
                     lon = "lon", lat = "lat", transformation = "lbmgls")
head(out[, c("h_am", "h_gini", "h_fracblack")])
```

## Citation

When using this package, please cite:

```bibtex
@TechReport{becker2025,
  author    = {Becker, Sascha O. and Boll, P. David and Voth, Hans-Joachim},
  title     = {Spatial Unit Roots in Regressions: A Practitioner's Guide and a Stata Package},
  year      = {2025},
  institution = {University of Warwick, Department of Economics},
  type      = {The Warwick Economics Research Paper Series (TWERPS)},
  number    = {1541},
  url       = {https://warwick.ac.uk/fac/soc/economics/research/workingpapers/2025/twerp_1541-_becker.pdf}
}
```

## References

Mueller, U. K. and Watson, M. W. (2024). Spatial Unit-Root Testing.
*Econometrica*, 92(6), 2121-2152. doi:
[10.3982/ECTA20682](https://doi.org/10.3982/ECTA20682)
