# EBSmoothr

`EBSmoothr` provides empirical Bayes smoothers for Gaussian observations with
L-GP, Matern, and binary Markov priors. It supports known observation
standard errors through `ebnm`-compatible interfaces and can also estimate one
common observation noise standard deviation.

## Installation

Copy and paste the following block into an R console:

```r
install.packages("pak")

pak::repo_add(
  INLA = "https://inla.r-inla-download.org/R/stable"
)

pak::pak("AgueroZZ/EBSmoothr")
```

To install a specific release, append its tag:

```r
pak::pak("AgueroZZ/EBSmoothr@v0.2.6")
```

## Binary Markov normal means

For ordered binary states, `ebnm_binary_markov()` fits

$$
x_i \mid \theta_i, s_i \sim N(\theta_i, s_i^2),
\qquad \theta_i \in \{0, 1\},
$$

under a stationary binary Markov prior with two parameters: the marginal
probability `pi1 = P(theta_i = 1)` and the switching probability
`flip_prob = P(theta_i != theta_{i-1})`. By default both are estimated by
empirical Bayes from exact marginal-likelihood evaluations; either can be fixed.

```r
library(EBSmoothr)

x <- c(0.05, 0.20, 0.55, 0.90, 1.10)
s <- rep(0.20, length(x))

fit <- ebnm_binary_markov(x, s)

fit$fitted_g
fit$posterior[, c("prob_zero", "prob_one")]
fit$viterbi_path
```

Set `pi1 = 0.5` for the symmetric prior (equal marginal probabilities and
equal expected run lengths), and add `flip_prob = 0.5` for the
equal-probability iid reference model:

```r
sym_fit <- ebnm_binary_markov(x, s, pi1 = 0.5)
iid_fit <- ebnm_binary_markov(x, s, flip_prob = 0.5, pi1 = 0.5)
```

The function accepts `g_init`, `fix_g` and `output`, so it can be passed
directly to flashier as an `ebnm_fn`; `s = Inf` is handled as an
uninformative observation:

```r
fl <- flashier::flash(Y, ebnm_fn = list(ebnm_binary_markov, ebnm::ebnm_normal))
```

Conditional on fixed prior parameters, forward-backward inference is exact.
The likelihood need not be concave in the parameters, so the default
numerical search is not guaranteed to find the global maximizer.

## Continuous smoothers

The main continuous-model entry points are:

- `ebnm_LGP_generator()` for one-dimensional local Gaussian-process smoothing;
- `ebnm_Matern_generator()` for one- or two-dimensional Matern smoothing;
- `eb_smoother()` for a higher-level interface that can also learn one common
  observation noise standard deviation.

See the function documentation and `NEWS.md` for model details and release
history.

## Development

After cloning the repository, load the development version with:

```r
devtools::load_all()
```

Run the package tests with:

```r
devtools::test()
```
