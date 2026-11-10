# Development changes

- Return Euclidean half-life bounds in the original coordinate units; apply Earth-radius conversion only to geographic coordinates.
- Solve the average-correlation equation to relative accuracy, rejecting unattainable targets instead of silently using a parameter cap.
- Use a QR basis for residual projection, avoiding normal-equation conditioning losses and dense residual-maker matrices.
- Reject insufficient covariance rank and scale-uninformative one-direction tests; bound I(1) calibration when 50% power cannot be bracketed.
- Normalize the observed half-life signal before angular-likelihood calculations and clamp antipodal distance roundoff.
- Reject nonlinear terms in the full `spur()` transformation pipeline: construct those numeric columns in the data first so that `H f(x)` is used rather than `f(H x)`. Standalone residual tests still accept general model-matrix formulas.

Upstream comparison suites are opt-in diagnostics (`SPUR_RUN_UPSTREAM_PARITY=true`). Default regression checks use analytic identities, defining matrix formulas, invariance, and Monte Carlo uncertainty.
