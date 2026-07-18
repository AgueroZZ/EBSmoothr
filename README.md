# EBSmoothr

`EBSmoothr` provides empirical Bayes smoothers for Gaussian observations with
L-GP, Matern, and symmetric binary Markov priors. It supports known observation
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

under a symmetric binary Markov prior. The switching probability can be fixed
or estimated by a one-dimensional numerical search using exact marginal-
likelihood evaluations.

```r
library(EBSmoothr)

x <- c(0.05, 0.20, 0.55, 0.90, 1.10)
s <- rep(0.20, length(x))

fit <- ebnm_binary_markov(x, s)

fit$fitted_g
fit$posterior[, c("prob_zero", "prob_one")]
fit$viterbi_path
```

Set `flip_prob = 0.5` for the equal-probability iid reference model:

```r
iid_fit <- ebnm_binary_markov(x, s, flip_prob = 0.5)
```

Conditional on a fixed switching probability, forward-backward inference is
exact. The likelihood need not be concave in the switching probability, so the
default numerical search is not guaranteed to find the global maximizer.

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
