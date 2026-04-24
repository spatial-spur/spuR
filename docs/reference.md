# Reference

This page documents the public `spuR` API.

## Overview

| Function | Description |
|---|---|
| `spur()` | Run the full SPUR workflow and return diagnostics plus levels and transformed branches |
| `spurtest_i0()` | Variable-level `I(0)` test |
| `spurtest_i1()` | Variable-level `I(1)` test |
| `spurtest_i0resid()` | Residual-based `I(0)` test for a regression formula |
| `spurtest_i1resid()` | Residual-based `I(1)` test for a regression formula |
| `spurtransform()` | Transform all variables referenced in a formula |
| `spurhalflife()` | Estimate a spatial half-life confidence interval |
| `spur_example` | Bundled Chetty-based example data |

## Conventions

### Coordinates

Functions that require spatial coordinates accept exactly one of:

- `lon` and `lat` for geographic coordinates
- `coords_euclidean` for planar coordinates

Do not pass both specifications in the same call.

### Simulation Controls

The SPUR tests and half-life interval are simulation-based.

- `q`: number of low-frequency weighted averages; default `15`
- `nrep`: number of Monte Carlo draws; default `100000`
- `seed`: integer seed for reproducibility; default `42L`

### Formula Inputs

- Variable-level diagnostics and half-life use a single-variable formula such
  as `am ~ 1` or `~ am`
- Residual diagnostics and the full workflow use a two-sided regression formula
  such as `am ~ gini + fracblack`
- `spurtransform()` accepts either a one-sided formula listing variables to
  transform or a two-sided analysis formula

## Full Workflow

### `spur()`

Run the full SPUR workflow.

**Signature**

```r
spur(
  formula,
  data,
  q = 15L,
  nrep = 100000L,
  lon = NULL,
  lat = NULL,
  coords_euclidean = NULL,
  seed = 42L,
  avc = 0.03,
  uncond = FALSE,
  cvs = FALSE
)
```

**Parameters**

- `formula` (`formula`): two-sided regression formula. The left-hand side is
  the dependent variable. The right-hand side defines the regressors used in
  both the levels and transformed branches.
- `data` (`data.frame`): input data containing the variables referenced in
  `formula` and the coordinate columns.
- `q`, `nrep`, `seed`: SPUR diagnostic simulation controls.
- `lon`, `lat` (`character | NULL`): longitude and latitude column names for
  geographic distance calculations.
- `coords_euclidean` (`character | NULL`): Euclidean coordinate column names.
  Use instead of `lon` and `lat`.
- `avc` (`numeric`): upper bound on average pairwise correlation passed to
  `scpcR::scpc()`.
- `uncond` (`logical`): passed through to `scpcR::scpc()`. If `TRUE`, use the
  unconditional SCPC branch.
- `cvs` (`logical`): passed through to `scpcR::scpc()`. If `TRUE`, store
  additional critical values in the SCPC result.

**Returns**

- object of class `"spur_result"`
  - `tests`: list containing `i0`, `i1`, `i0resid`, and `i1resid`
  - `fits`: list containing `levels` and `transformed`
  - each fit contains `model` and `scpc`

## Diagnostics

### `spurtest_i0()`

Variable-level `I(0)` test.

**Signature**

```r
spurtest_i0(
  formula,
  data,
  q = 15,
  nrep = 100000,
  lon = NULL,
  lat = NULL,
  coords_euclidean = NULL,
  seed = 42L,
  verbose = FALSE
)
```

**Parameters**

- `formula` (`formula`): single-variable formula such as `am ~ 1` or `~ am`.
- `data` (`data.frame`): input data containing the tested variable and the
  coordinate columns.
- `q`, `nrep`, `seed`: as above.
- `lon`, `lat`, `coords_euclidean`: coordinate specification.
- `verbose` (`logical`): if `TRUE`, print the test summary when the function
  runs.

**Returns**

- object of class `"spur_test"`
  - `pvalue`
  - `teststat`
  - `ha_param`
  - `cv`
  - `full`
  - attribute `test_type = "I(0)"`

### `spurtest_i1()`

Variable-level `I(1)` test.

**Signature**

```r
spurtest_i1(
  formula,
  data,
  q = 15,
  nrep = 100000,
  lon = NULL,
  lat = NULL,
  coords_euclidean = NULL,
  seed = 42L,
  verbose = FALSE
)
```

**Parameters**

- `formula` (`formula`): single-variable formula such as `am ~ 1` or `~ am`.
- `data` (`data.frame`): input data containing the tested variable and the
  coordinate columns.
- `q`, `nrep`, `seed`, `lon`, `lat`, `coords_euclidean`, `verbose`: as above.

**Returns**

- object of class `"spur_test"` with the same structure as `spurtest_i0()`
  and attribute `test_type = "I(1)"`

### `spurtest_i0resid()`

Residual-based `I(0)` test for a regression formula.

**Signature**

```r
spurtest_i0resid(
  formula,
  data,
  q = 15,
  nrep = 100000,
  lon = NULL,
  lat = NULL,
  coords_euclidean = NULL,
  seed = 42L,
  verbose = FALSE
)
```

**Parameters**

- `formula` (`formula`): two-sided regression formula defining the residual
  process to be tested.
- `data` (`data.frame`): input data.
- `q`, `nrep`, `seed`, `lon`, `lat`, `coords_euclidean`, `verbose`: as above.

**Returns**

- object of class `"spur_test"` with attribute `test_type = "I(0) Residual"`

### `spurtest_i1resid()`

Residual-based `I(1)` test for a regression formula.

**Signature**

```r
spurtest_i1resid(
  formula,
  data,
  q = 15,
  nrep = 100000,
  lon = NULL,
  lat = NULL,
  coords_euclidean = NULL,
  seed = 42L,
  verbose = FALSE
)
```

**Parameters**

- `formula` (`formula`): two-sided regression formula defining the residual
  process to be tested.
- `data` (`data.frame`): input data.
- `q`, `nrep`, `seed`, `lon`, `lat`, `coords_euclidean`, `verbose`: as above.

**Returns**

- object of class `"spur_test"` with attribute `test_type = "I(1) Residual"`

## Transformation

### `spurtransform()`

Transform all variables referenced in a formula and append the transformed
columns to a copy of the input data.

**Signature**

```r
spurtransform(
  formula,
  data,
  prefix = "h_",
  transformation = "lbmgls",
  radius = NULL,
  clustvar = NULL,
  lon = NULL,
  lat = NULL,
  coords_euclidean = NULL,
  separately = FALSE,
  replace = FALSE,
  verbose = FALSE
)
```

**Parameters**

- `formula` (`formula`): formula identifying the variables to transform. Every
  variable referenced in the formula is transformed.
- `data` (`data.frame`): input data.
- `prefix` (`character`): prefix applied to transformed variable names.
  Default: `h_`.
- `transformation` (`character`): transformation mode.
  - `"lbmgls"`: GLS-style transformation; default empirical branch
  - `"nn"`: nearest-neighbor differencing
  - `"iso"`: radius-based differencing
  - `"cluster"`: within-cluster demeaning
- `radius` (`numeric | NULL`): radius used by the isotropic transformation;
  required when `transformation = "iso"`.
- `clustvar` (`character | NULL`): cluster label column; required when
  `transformation = "cluster"`.
- `lon`, `lat`, `coords_euclidean`: coordinate specification for
  distance-based transformations.
- `separately` (`logical`): if `TRUE`, handle missingness separately by
  variable; otherwise use a common estimation sample across transformed
  variables.
- `replace` (`logical`): if `TRUE`, allow overwriting existing output columns.
- `verbose` (`logical`): if `TRUE`, print coordinate diagnostics.

**Returns**

- `data.frame`: copy of the input data with transformed columns appended

## Persistence

### `spurhalflife()`

Estimate a confidence interval for the spatial half-life of a variable.

**Signature**

```r
spurhalflife(
  formula,
  data,
  q = 15,
  nrep = 100000,
  level = 95,
  normdist = FALSE,
  lon = NULL,
  lat = NULL,
  coords_euclidean = NULL,
  seed = 42L,
  verbose = FALSE
)
```

**Parameters**

- `formula` (`formula`): single-variable formula such as `am ~ 1` or `~ am`.
- `data` (`data.frame`): input data containing the tested variable and the
  coordinate columns.
- `q`, `nrep`, `seed`, `lon`, `lat`, `coords_euclidean`, `verbose`: as above.
- `level` (`numeric`): confidence level in percent.
- `normdist` (`logical`): if `TRUE`, report the interval as fractions of the
  maximum pairwise distance; otherwise report it in distance units.

**Returns**

- object of class `"spur_halflife"`
  - `ci_l`
  - `ci_u`
  - `max_dist`
  - `full`
  - attribute `level`

## Example Data

### `spur_example`

Bundled example data derived from Chetty et al. (2014).

**Contents**

- `cz`, `czname`, `state`
- `lat`, `lon`
- `am`, `rm`, `gini`, `fracblack`

Use `data(spur_example)` to load it into the workspace.

## Return Objects

### `"spur_test"`

Returned by `spurtest_i0()`, `spurtest_i1()`, `spurtest_i0resid()`, and
`spurtest_i1resid()`.

**Fields**

- `pvalue`
- `teststat`
- `ha_param`
- `cv`
- `full`

**Methods**

- `print()`

The printed summary reports the test statistic and p-value.

### `"spur_halflife"`

Returned by `spurhalflife()`.

**Fields**

- `ci_l`
- `ci_u`
- `max_dist`
- `full`

**Methods**

- `print()`

The printed summary reports the lower and upper half-life bounds.

### `"spur_result"`

Returned by `spur()`.

**Fields**

- `tests`
- `fits`

`tests` contains:

- `i0`
- `i1`
- `i0resid`
- `i1resid`

`fits` contains:

- `levels`
- `transformed`

Each fit contains:

- `model`
- `scpc`

**Methods**

- `print()`
- `summary()`

Both methods render the combined diagnostic table and the side-by-side levels
and transformed regression summary.
