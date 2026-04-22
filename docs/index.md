# spuR

`spuR` implements the SPUR workflow for diagnosing and correcting spatial unit
roots in cross-sectional regressions. It covers the diagnostic and
transformation stage of the workflow. For standalone SCPC inference on fitted
models, see `scpcR`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("spatial-spur/spuR@v0.1.2")
```

## Example: Chetty Dataset

This walkthrough follows the practitioner guide in Becker, Boll, and Voth
(2026). The branch decision is based on the dependent-variable `I(0)` and
`I(1)` tests, using a 10% significance level.

## 1. Prepare the sample

Construct the estimation sample and retain the variables used in the
diagnostics and regression.

```r
library(spuR)

data(spur_example)

df <- subset(
  spur_example,
  !state %in% c("AK", "HI"),
  select = c(am, gini, fracblack, lat, lon, state)
)

df <- stats::na.omit(df)
```

## 2. Test the dependent variable against the I(0) alternative

The first diagnostic tests the dependent variable under the spatial `I(0)`
null.

```r
i0 <- spurtest_i0(am ~ 1, df, lon = "lon", lat = "lat")

print(i0)
```

## 3. Test the dependent variable against the I(1) alternative

The second diagnostic tests the same variable under the spatial `I(1)` null.

```r
i1 <- spurtest_i1(am ~ 1, df, lon = "lon", lat = "lat")

print(i1)
```

## 4. Apply the decision rule

Using a 10% significance threshold:

- If you do **not** reject `I(0)` and you **do** reject `I(1)`, proceed in
  levels.
- In every other case, treat the specification as requiring spatial
  differencing and transform the dependent and independent variables together.

## 5. Levels branch

Use this branch only when the decision rule implies that the dependent variable
is consistent with `I(0)` and inconsistent with `I(1)`.

```r
fit_levels <- stats::lm(am ~ gini + fracblack, data = df)
scpc_levels <- scpcR::scpc(fit_levels, df, lon = "lon", lat = "lat")

summary(scpc_levels)
```

In this branch there is no SPUR transformation step; the regression is
estimated in levels and SCPC is used for inference.

## 6. Transformed branch

In every other case, transform the dependent and independent variables
together, re-estimate the regression on the transformed data, and use SCPC for
inference there.

```r
transformed <- spurtransform(
  am ~ gini + fracblack,
  df,
  lon = "lon",
  lat = "lat",
  transformation = "lbmgls"
)

fit_transformed <- stats::lm(
  h_am ~ h_gini + h_fracblack,
  data = transformed
)

scpc_transformed <- scpcR::scpc(
  fit_transformed,
  transformed,
  lon = "lon",
  lat = "lat"
)

summary(scpc_transformed)
```

The default empirical branch is the `lbmgls` transformation.

## 7. Residual diagnostics

`spuR` also provides residual-based `I(0)` and `I(1)` tests:

```r
i0resid <- spurtest_i0resid(
  am ~ gini + fracblack,
  df,
  lon = "lon",
  lat = "lat"
)

i1resid <- spurtest_i1resid(
  am ~ gini + fracblack,
  df,
  lon = "lon",
  lat = "lat"
)
```

## 8. Packaged shortcut

If you want the package’s default end-to-end implementation rather than the
manual branch logic, use `spur()`.

```r
pipeline <- spur(
  am ~ gini + fracblack,
  df,
  lon = "lon",
  lat = "lat"
)

summary(pipeline)
```

`spur()` returns all four SPUR diagnostics together with the levels and
transformed fits.

## Next Step

See [Reference](reference.md) for the full public API, parameter meanings, and
return objects.
