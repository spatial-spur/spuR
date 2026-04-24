# spuR

`spuR` implements the SPUR workflow for diagnosing and correcting spatial unit
roots in cross-sectional regressions. It covers the diagnostic and
transformation stage of the workflow. For standalone SCPC inference on fitted
models, see `scpcR`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("spatial-spur/scpcR@v0.1.3")
remotes::install_github("spatial-spur/spuR@v0.1.2")
```

We recommend installing the latest tagged version of both packages by 
pointing to the latest tagged version of each.

## Example: Chetty Dataset

In this example, we walk you through the workflow we recommend with the
packages step-by-step. We also provide a one-stop [pipeline wrapper](#pipeline-wrapper)
implementing the entire workflow in one step.

### Data preparation

For illustration, we load the Chetty dataset we ship as part of the package. Of
course, the analysis in principle follows the same logic on any other dataset.
In this specific case, we first omit the non-contiguous US states. We also drop
rows with missing values.

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

### Testing for a spatial unit root

Based on MW 2024, we suggest first testing for a spatial unit root setting
using the `I(0)` and `I(1)` tests on the dependent variable.

One way to do this is to use the `spurtest_i0()` and `spurtest_i1()` functions
directly:

```r
# am is the dependent variable
i0 <- spurtest_i0(am ~ 1, df, lon = "lon", lat = "lat")
i1 <- spurtest_i1(am ~ 1, df, lon = "lon", lat = "lat")

print(i0)
print(i1)
```

### Interpreting the test statistics

Using a 10% significance threshold, we suggest interpreting the results with the following heuristic:

- If you do **not** reject `I(0)` and you **do** reject `I(1)`, there is **likely no spatial unit root** and you can proceed in levels
- every other case implies a **possible spatial unit root** - in that case, we suggest transforming all dependent and independent variables before running regressions

We suggest always applying SCPC inference.

### Case 1: likely no spatial unit root

If the heuristic implies your scenario is unlikely to be a spatial unit root, we suggest proceeding in levels but applying SCPC inference:

```r
fit_levels <- stats::lm(am ~ gini + fracblack, data = df)
scpc_levels <- scpcR::scpc(fit_levels, df, lon = "lon", lat = "lat")

summary(scpc_levels)
```

### Case 2: likely spatial unit root

If you do have a likely spatial unit root according to the heuristic above, we suggest applying the transformation and running the regression on transformed variables with SCPC inference:

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

### Sanity check

As a sanity check, we recommend validating that your regression residuals do not
have a spatial unit root. You can do that using the `I(0) residual` and `I(1)
residual` tests:

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

### Pipeline wrapper

As a shortcut to implementing all of those steps individually, we also provide a
`spur()` wrapper that implements the entire pipeline in one step. It simply runs
all tests and returns all results.

```r
pipeline <- spur(
  am ~ gini + fracblack,
  df,
  lon = "lon",
  lat = "lat"
)

summary(pipeline)
```

## Next Step

See [Reference](reference.md) for the full public API, parameter meanings, and
return objects.
