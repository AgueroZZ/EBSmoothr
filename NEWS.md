# EBSmoothr 0.2.1

## Improvements

- Added a selected-inverse posterior variance path for large sparse Matérn
  precision matrices when the required sparsity pattern is covered.
- Made Matérn Laplace public backends TMB-or-error: `backend = "laplace"` and
  `backend = "laplace_fisher"` no longer silently fall back to the internal R
  reference implementation.
- Made Matérn log-link Step A multistart construction strictly respect
  `matern_n_starts`, including known-noise fits where `matern_n_starts = 1`
  now means exactly one start.

## Validation

- Added targeted tests for Matérn TMB-or-error backend semantics, internal
  `laplace_r` validation access, sparse selected-inverse posterior variance,
  and strict `matern_n_starts` behavior.

# EBSmoothr 0.2.0

## New Features

- Added softplus response-link support for L-GP and Matern smoothers on the
  Laplace backends.
- Added deterministic softplus posterior response moments under the marginal
  Gaussian Laplace approximation. These moments use fixed Gauss-Hermite
  quadrature and are designed to agree with `posterior_sampler()` in
  expectation without storing Monte Carlo noise.
- Added point-mass reference families to `eb_smoother()` and spatial scoring
  workflows, including `point_exponential`, `point_normal`, and
  `point_laplace`.

## Improvements

- Updated `eb_smoother()` and spatial scoring outputs to use reference-family
  terminology instead of the older nonspatial terminology where multiple
  reference families are supported.
- Added `loglik_reference` and `reference_fits` outputs while retaining
  `loglik_nonspatial` and `nonspatial_fits` as compatibility aliases.
- Added softplus backend selection support for observed and Fisher Laplace
  curvature where applicable.
- Updated package documentation and internal workflow vignettes for the new
  backend and reference-family semantics.

## Fixes

- Fixed softplus posterior variance reporting by replacing first-order
  delta-method approximations with deterministic Gaussian-transform moments.
- Fixed softplus Matern Laplace observation terms so the R implementation uses
  the softplus mean, gradient, and curvature instead of log-link formulas.
- Fixed the TMB softplus implementation to use an AD-safe stable expression.
- Added tests comparing softplus Gaussian moments against high-accuracy
  numerical integration and posterior sampler estimates.
