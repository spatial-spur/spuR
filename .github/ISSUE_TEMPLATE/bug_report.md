---
name: Bug report
about: Report a reproducible problem with spuR
title: "[bug] "
labels: bug
---

## Summary

Describe the problem in 1-3 sentences.

## Minimal reproduction

Provide a copy-pasteable example if possible.

```r

```

## Expected behavior

What did you expect to happen?

## Actual behavior

What happened instead? Include the full traceback, warning, or numerical mismatch if relevant.

```text

```

## Environment

- OS:
- R version:
- install method: `install.packages` / `remotes::install_github` / `devtools::load_all` / local source install
- spuR version or commit:
- affected function or workflow: `spurtest_i0` / `spurtest_i1` / `spurtest_i0resid` / `spurtest_i1resid` / `spurhalflife` / `spurtransform` / formula API / docs / CI
- parity tools available, if relevant: Stata / Matlab + MW replication package / neither

## Additional context

Anything else that might help reproduce or explain the issue.

- Does this reproduce on synthetic data, the Chetty fixture, or only your own data?
- Does this affect `lon`/`lat`, `coords_euclidean`, or both?
- If relevant, which mode or option is affected: `i0`, `i1`, `i0resid`, `i1resid`, `lbmgls`, `nn`, `iso`, `cluster`, `level`, `normdist`, or the formula API?
- Is the issue a crash, wrong result, parity mismatch, docs mismatch, install problem, or performance regression?
