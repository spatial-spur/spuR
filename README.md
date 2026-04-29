<p align="center">
  <img src="assets/logo.png" alt="SPUR logo">
</p>

# spuR

`spuR` implements diagnostics and correction methods for spatial unit roots developed by Müller and Watson (2024) in R.

**When using this code, please cite [Becker, Boll and Voth (2026)](https://pauldavidboll.com/SPUR_Stata_Journal_website.pdf):**

```bibtex
@Article{becker2026,
  author    = {Becker, Sascha O. and Boll, P. David and Voth, Hans-Joachim},
  title     = {Testing and Correcting for Spatial Unit Roots in Regression Analysis},
  journal   = {Stata Journal},
  year      = {forthcoming},
  note      = {Forthcoming}
}
```

If you encounter any issues or have any questions, please open an issue on GitHub or contact the authors.

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

## Documentation

Please refer to [the package documentation](https://spatial-spur.github.io/scpcR/) for detailed information and other (R, Python, Stata) packages.

## References

Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Testing and Correcting for Spatial Unit Roots in Regression Analysis", Forthcoming at the Stata Journal.

Müller, Ulrich K. and Mark W. Watson "Spatial Unit Roots and Spurious Regression", Econometrica 92(5) (2024), 1661–1695. https://www.princeton.edu/~umueller/SPUR.pdf.
