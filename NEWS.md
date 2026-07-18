# EBSmoothr 0.2.6

## New Features

- Added `BinaryMarkov()` and `ebnm_binary_markov()` for an exact ordered
  intron/exon model with symmetric state transitions. The solver supports a
  fixed transition probability or empirical-Bayes numerical search using exact
  marginal-likelihood evaluations, and returns forward-backward state
  marginals, pairwise transition posteriors, a Viterbi path, and exact posterior
  samples. Because the likelihood need not be concave, the numerical search is
  not guaranteed to find the global maximizer. A flip probability of `0.5`
  reduces exactly to the iid Bernoulli(0.5) binary baseline, so the same solver
  covers both the symmetric Markov model and its equal-probability iid limit.
- Added opt-in PC prior support for L-GP smoothers through
  `pc.penalty = list(scale = c(anchor, alpha))` or the explicit alias
  `pc.penalty = list(latent_scale = c(anchor, alpha))`. The prior is placed on
  the latent process standard deviation `sigma_u = exp(-theta / 2)` with
  `P(sigma_u > anchor) = alpha`; omitting `alpha` defaults to `0.5`.
- The L-GP PC prior is available through both `ebnm_LGP_generator()` and
  `eb_smoother(family = "lgp")`, including fixed-noise and learned-noise L-GP
  fits. The prior contributes to the optimized L-GP objective and returned fits
  include `pc_penalty`, `log_likelihood_pc_prior_theta`, and PC-specific
  `prior_family` labels.

## Validation

- Added L-GP tests confirming that the PC prior affects fitted latent scale and
  posterior smoothness in the expected direction across anchor and alpha
  choices, and that `ebnm_LGP_generator()` and `eb_smoother()` agree for the
  same prior specification.

# EBSmoothr 0.2.5

## Performance

- Sped up the sparse SPD factorization used across Matern and L-GP Laplace
  objectives. `.factorize_spd()` now lets CHOLMOD choose supernodal
  factorization (several times faster than the previous simplicial-only path
  on 2D meshes), reuses cached symbolic analyses through numeric-only
  refactorization when the sparsity pattern repeats across hyperparameter
  evaluations, and computes log-determinants directly from the factor instead
  of materializing the triangular factor with `Matrix::expand()`. The cache
  size is controlled by `options(EBSmoothr.spd_factor_cache_max = 30L)`.
- Stationary Matern SPDE precisions are now assembled directly from the
  template's `M0`/`M1`/`M2` matrices with a precomputed aligned sparsity
  pattern instead of calling `INLA::inla.spde2.precision()` on every
  hyperparameter evaluation (~75x faster per precision build). Non-stationary
  or non-standard templates automatically fall back to INLA.
- Exact Gaussian Matern objectives now precompute the data-side sufficient
  statistics (`A'WA`, `A'Wx`, `A'W1`, and scalar moments) once per
  hyperparameter optimization instead of once per objective evaluation.
  Learned-noise fits cache unit-noise statistics and rescale them per
  evaluation.
- The Fisher-PQL loop now stops early when the linear predictor change drops
  below `pql_tol` instead of always running `pql_inner_iter` passes.
  `fisher_pql_diagnostics$inner_iterations` now reports the number of passes
  actually executed and the new `max_inner_iter` field reports the cap.

## Validation

- Added regression tests asserting the direct SPDE precision assembly matches
  `INLA::inla.spde2.precision()` and that the factorization cache returns
  correct log-determinants and solves across repeated and distinct sparsity
  patterns.
- Updated stale test expectations that still assumed non-identity
  `backend = "auto"` resolves to `laplace_fisher`; since 0.2.4 the auto policy
  routes log and softplus links to `fisher_pql` (these tests were failing
  before the 0.2.5 changes).
- Benchmarked the fisher_pql path against TMB `laplace_fisher` on 2D log and
  softplus problems; fisher_pql remains several times faster at all tested
  sizes, so the auto backend policy is unchanged.

# EBSmoothr 0.2.4

## Improvements

- Changed non-identity Matern `backend = "auto"` to use `backend =
  "fisher_pql"` for log and softplus links.
- Exposed the public `pql_inner_iter` argument for Matern Fisher-PQL fits in
  both `ebnm_Matern_generator()` and `eb_smoother()`. The default remains
  three pseudo-Gaussian exact Step A updates.
- Fixed non-identity initialization for Matern and L-GP fits so log and
  softplus links initialize beta by matching the response-scale mean, keep
  latent prior scales on the linear-predictor scale, and initialize learned
  observation noise on the raw response scale.
- Made the experimental softplus Matern `backend = "inlabru"` use explicit
  prior semantics: known-noise fits no longer synthesize a PC prior when
  `pc.penalty = NULL`, while learned-noise fits require an explicit
  `pc.penalty` with `range`, `sigma`, and `noise` entries.
- Changed inlabru primary `log_likelihood` reporting to the package's manual
  Laplace objective evaluated at the inlabru fitted parameters, retaining raw
  inlabru/INLA marginal likelihood as a diagnostic field.
- Added internal exact-Gaussian learned-noise Matérn support for fixed
  per-observation noise scales, enabling future Fisher-PQL pseudo-response
  fits with `s_i = noise_sd * noise_scale_i`.
- Added non-identity Matern `backend = "fisher_pql"`, which uses three
  Fisher/PQL pseudo-Gaussian exact Matern Step A updates by default and reports
  a Fisher/Laplace score evaluated at the final PQL mode. This score is an
  approximate PQL-mode diagnostic, not a true re-optimized original-model
  Laplace marginal likelihood.
- Added latent linear-predictor posterior summaries as `posterior_eta` for
  Matern Laplace-style fits, enabling direct diagnostics on the latent
  `eta` scale.

## Validation

- Added targeted tests for no-PC known-noise inlabru fits, explicit PC-prior
  handling, learned-noise prior validation, and comparable likelihood
  reporting.
- Added regression tests for the exact learned-noise noise-scale path,
  including unit-scale equivalence and heteroskedastic objective checks.
- Added Fisher-PQL backend tests covering pseudo-response construction,
  known-noise and learned-noise dispatch, Step A mode checks against an R
  sparse reference solve, PQL-mode likelihood reporting, and small-case
  agreement with Laplace backends.

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
