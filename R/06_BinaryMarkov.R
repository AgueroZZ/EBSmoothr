# Symmetric-prior helpers from EBSmoothr <= 0.3.5, kept for backward compatibility.
.normalize_binary_s <- function(s, n) {
  if (!is.numeric(s) || !(length(s) %in% c(1L, n))) {
    stop("`s` must be numeric and either scalar or have length equal to `x`.")
  }
  s <- rep(as.numeric(s), length.out = n)
  if (any(!is.finite(s)) || any(s <= 0)) {
    stop("`s` values must be finite and strictly positive.")
  }
  s
}

.binary_logspace_add <- function(a, b) {
  max_log <- pmax(a, b)
  out <- max_log + log1p(exp(-abs(a - b)))
  both_negative_infinite <- is.infinite(a) & a < 0 & is.infinite(b) & b < 0
  out[both_negative_infinite] <- -Inf
  out
}

.validate_binary_markov_flip_prob <- function(flip_prob, name = "flip_prob") {
  if (!is.numeric(flip_prob) || length(flip_prob) != 1L) {
    stop("`", name, "` must be a numeric scalar.")
  }
  if (!is.finite(flip_prob)) {
    stop("`", name, "` must be a finite numeric scalar.")
  }
  if (flip_prob < 0 || flip_prob > 0.5) {
    stop("`", name, "` must be between 0 and 0.5.")
  }
  as.numeric(flip_prob)
}

.validate_binary_markov_bounds <- function(bounds) {
  if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds))) {
    stop("`flip_prob_bounds` must contain two finite numeric values.")
  }
  bounds <- as.numeric(bounds)
  if (bounds[1] < 0 || bounds[2] > 0.5) {
    stop("`flip_prob_bounds` must stay inside [0, 0.5].")
  }
  if (bounds[1] >= bounds[2]) {
    stop("`flip_prob_bounds` must be in increasing order.")
  }
  bounds
}

.binary_markov_log_transition <- function(flip_prob) {
  log_stay <- log1p(-flip_prob)
  log_flip <- if (flip_prob == 0) -Inf else log(flip_prob)
  matrix(
    c(log_stay, log_flip, log_flip, log_stay),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(from = c("0", "1"), to = c("0", "1"))
  )
}

.binary_markov_forward <- function(log_emission, flip_prob) {
  n <- nrow(log_emission)
  log_transition <- .binary_markov_log_transition(flip_prob)
  log_alpha <- matrix(-Inf, nrow = n, ncol = 2L)
  log_alpha[1L, ] <- log(0.5) + log_emission[1L, ]

  if (n > 1L) {
    for (i in 2:n) {
      log_alpha[i, 1L] <- log_emission[i, 1L] + .binary_logspace_add(
        log_alpha[i - 1L, 1L] + log_transition[1L, 1L],
        log_alpha[i - 1L, 2L] + log_transition[2L, 1L]
      )
      log_alpha[i, 2L] <- log_emission[i, 2L] + .binary_logspace_add(
        log_alpha[i - 1L, 1L] + log_transition[1L, 2L],
        log_alpha[i - 1L, 2L] + log_transition[2L, 2L]
      )
    }
  }

  list(
    log_alpha = log_alpha,
    log_transition = log_transition,
    log_likelihood = as.numeric(.binary_logspace_add(log_alpha[n, 1L], log_alpha[n, 2L]))
  )
}

.binary_markov_forward_backward <- function(log_emission, flip_prob) {
  forward <- .binary_markov_forward(log_emission, flip_prob)
  log_alpha <- forward$log_alpha
  log_transition <- forward$log_transition
  log_likelihood <- forward$log_likelihood
  n <- nrow(log_emission)

  log_beta <- matrix(0, nrow = n, ncol = 2L)
  if (n > 1L) {
    for (i in (n - 1L):1L) {
      log_beta[i, 1L] <- .binary_logspace_add(
        log_transition[1L, 1L] + log_emission[i + 1L, 1L] + log_beta[i + 1L, 1L],
        log_transition[1L, 2L] + log_emission[i + 1L, 2L] + log_beta[i + 1L, 2L]
      )
      log_beta[i, 2L] <- .binary_logspace_add(
        log_transition[2L, 1L] + log_emission[i + 1L, 1L] + log_beta[i + 1L, 1L],
        log_transition[2L, 2L] + log_emission[i + 1L, 2L] + log_beta[i + 1L, 2L]
      )
    }
  }

  state_probability <- exp(log_alpha + log_beta - log_likelihood)
  state_probability <- state_probability / rowSums(state_probability)
  colnames(state_probability) <- c("prob_zero", "prob_one")

  posterior_transitions <- if (n == 1L) {
    data.frame(
      from = integer(),
      to = integer(),
      prob_00 = numeric(),
      prob_01 = numeric(),
      prob_10 = numeric(),
      prob_11 = numeric(),
      prob_flip = numeric()
    )
  } else {
    transition_rows <- vector("list", n - 1L)
    for (i in seq_len(n - 1L)) {
      log_pair <- c(
        log_alpha[i, 1L] + log_transition[1L, 1L] + log_emission[i + 1L, 1L] + log_beta[i + 1L, 1L],
        log_alpha[i, 1L] + log_transition[1L, 2L] + log_emission[i + 1L, 2L] + log_beta[i + 1L, 2L],
        log_alpha[i, 2L] + log_transition[2L, 1L] + log_emission[i + 1L, 1L] + log_beta[i + 1L, 1L],
        log_alpha[i, 2L] + log_transition[2L, 2L] + log_emission[i + 1L, 2L] + log_beta[i + 1L, 2L]
      )
      pair_probability <- exp(log_pair - max(log_pair))
      pair_probability <- pair_probability / sum(pair_probability)
      transition_rows[[i]] <- data.frame(
        from = i,
        to = i + 1L,
        prob_00 = pair_probability[1L],
        prob_01 = pair_probability[2L],
        prob_10 = pair_probability[3L],
        prob_11 = pair_probability[4L],
        prob_flip = pair_probability[2L] + pair_probability[3L]
      )
    }
    do.call(rbind, transition_rows)
  }

  list(
    log_likelihood = log_likelihood,
    log_alpha = log_alpha,
    log_beta = log_beta,
    log_transition = log_transition,
    state_probability = state_probability,
    posterior_transitions = posterior_transitions
  )
}

.binary_markov_viterbi <- function(log_emission, log_transition) {
  n <- nrow(log_emission)
  log_delta <- matrix(-Inf, nrow = n, ncol = 2L)
  back_pointer <- matrix(1L, nrow = n, ncol = 2L)
  log_delta[1L, ] <- log(0.5) + log_emission[1L, ]

  if (n > 1L) {
    for (i in 2:n) {
      for (state in 1:2) {
        candidate <- log_delta[i - 1L, ] + log_transition[, state]
        best_previous <- which.max(candidate)
        log_delta[i, state] <- candidate[best_previous] + log_emission[i, state]
        back_pointer[i, state] <- best_previous
      }
    }
  }

  path_index <- integer(n)
  path_index[n] <- which.max(log_delta[n, ])
  if (n > 1L) {
    for (i in (n - 1L):1L) {
      path_index[i] <- back_pointer[i + 1L, path_index[i + 1L]]
    }
  }
  as.integer(path_index - 1L)
}

.binary_markov_posterior_sampler <- function(log_alpha, log_transition, state_probability) {
  n <- nrow(log_alpha)
  function(nsamp) {
    if (!is.numeric(nsamp) || length(nsamp) != 1L || !is.finite(nsamp) ||
        nsamp < 1 || nsamp != floor(nsamp)) {
      stop("`nsamp` must be a positive integer.")
    }
    nsamp <- as.integer(nsamp)
    draws <- matrix(0L, nrow = nsamp, ncol = n)

    for (sample_index in seq_len(nsamp)) {
      draws[sample_index, n] <- as.integer(
        stats::runif(1) < state_probability[n, 2L]
      )
      if (n > 1L) {
        for (i in (n - 1L):1L) {
          next_state <- draws[sample_index, i + 1L] + 1L
          log_weight <- log_alpha[i, ] + log_transition[, next_state]
          weight <- exp(log_weight - max(log_weight))
          prob_one <- weight[2L] / sum(weight)
          draws[sample_index, i] <- as.integer(stats::runif(1) < prob_one)
        }
      }
    }
    draws
  }
}

# ---------------------------------------------------------------------------
# Binary Markov prior with an estimated marginal probability (EBSmoothr 0.4.0).
# The symmetric helpers above are kept for backward compatibility.
# ---------------------------------------------------------------------------

.validate_binary_markov_pi1 <- function(pi1, name = "pi1") {
  if (!is.numeric(pi1) || length(pi1) != 1L) {
    stop("`", name, "` must be a numeric scalar.")
  }
  if (!is.finite(pi1)) {
    stop("`", name, "` must be a finite numeric scalar.")
  }
  if (pi1 <= 0 || pi1 >= 1) {
    stop("`", name, "` must be strictly between 0 and 1.")
  }
  as.numeric(pi1)
}

# Largest switching probability with non-negative lag-one autocorrelation.
.binary_markov_flip_upper <- function(pi1) 2 * pi1 * (1 - pi1)

.validate_binary_markov_flip_given_pi1 <- function(flip_prob, pi1, name = "flip_prob") {
  if (!is.numeric(flip_prob) || length(flip_prob) != 1L) {
    stop("`", name, "` must be a numeric scalar.")
  }
  if (!is.finite(flip_prob)) {
    stop("`", name, "` must be a finite numeric scalar.")
  }
  upper <- .binary_markov_flip_upper(pi1)
  if (flip_prob < 0 || flip_prob > upper * (1 + 1e-12)) {
    stop("`", name, "` must be between 0 and ", format(signif(upper, 6)),
         " (= 2 * pi1 * (1 - pi1)).")
  }
  min(as.numeric(flip_prob), upper)
}

.validate_binary_markov_pi1_bounds <- function(bounds) {
  if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds))) {
    stop("`pi1_bounds` must contain two finite numeric values.")
  }
  bounds <- as.numeric(bounds)
  if (bounds[1] <= 0 || bounds[2] >= 1) {
    stop("`pi1_bounds` must stay inside (0, 1).")
  }
  if (bounds[1] >= bounds[2]) {
    stop("`pi1_bounds` must be in increasing order.")
  }
  bounds
}

.binary_markov_transition <- function(flip_prob, pi1) {
  p01 <- flip_prob / (2 * (1 - pi1))
  p10 <- flip_prob / (2 * pi1)
  # unnamed: indexing a matrix with dimnames is much slower inside the recursions
  matrix(c(1 - p01, p01, p10, 1 - p10), nrow = 2L, byrow = TRUE)
}

# Accepts s = Inf (flat emission) and s = 0 (theta fixed at x, which must be 0 or 1).
.binary_markov_emission <- function(x, s) {
  n <- length(x)
  if (!is.numeric(s) || !(length(s) %in% c(1L, n))) {
    stop("`s` must be numeric and either scalar or have length equal to `x`.")
  }
  s <- rep(as.numeric(s), length.out = n)
  if (anyNA(s) || any(s < 0)) {
    stop("`s` values must be non-negative (Inf is allowed).")
  }
  log_emission <- matrix(0, n, 2L)
  ok <- is.finite(s) & s > 0
  log_emission[ok, 1L] <- stats::dnorm(x[ok], mean = 0, sd = s[ok], log = TRUE)
  log_emission[ok, 2L] <- stats::dnorm(x[ok], mean = 1, sd = s[ok], log = TRUE)
  pinned <- s == 0
  if (any(pinned)) {
    state <- round(x[pinned])
    if (any(!state %in% c(0, 1)) || any(abs(x[pinned] - state) > 1e-8)) {
      stop("Entries with `s = 0` must have `x` equal to 0 or 1.")
    }
    log_emission[pinned, ] <- cbind(ifelse(state == 0, 0, -Inf), ifelse(state == 1, 0, -Inf))
  }
  log_max <- pmax(log_emission[, 1L], log_emission[, 2L])
  list(s = s, log_emission = log_emission, log_max = log_max,
       E = exp(log_emission - log_max))
}

# Scaled forward pass. `alpha` holds filtered probabilities P(theta_i | x_1..i).
.binary_markov_scaled_forward <- function(em, P, init) {
  E <- em$E
  n <- nrow(E)
  alpha <- matrix(0, n, 2L)
  scale <- numeric(n)
  a <- init * E[1L, ]
  scale[1L] <- sum(a)
  alpha[1L, ] <- a / scale[1L]
  if (n > 1L) {
    for (i in 2:n) {
      a1 <- alpha[i - 1L, 1L]
      a2 <- alpha[i - 1L, 2L]
      a <- c(a1 * P[1L, 1L] + a2 * P[2L, 1L], a1 * P[1L, 2L] + a2 * P[2L, 2L]) * E[i, ]
      scale[i] <- sum(a)
      alpha[i, ] <- a / scale[i]
    }
  }
  list(alpha = alpha, scale = scale,
       log_likelihood = sum(log(scale)) + sum(em$log_max))
}

.binary_markov_scaled_backward <- function(em, P, scale) {
  E <- em$E
  n <- nrow(E)
  beta <- matrix(1, n, 2L)
  if (n > 1L) {
    for (i in (n - 1L):1L) {
      e <- E[i + 1L, ] * beta[i + 1L, ]
      beta[i, ] <- c(P[1L, 1L] * e[1L] + P[1L, 2L] * e[2L],
                     P[2L, 1L] * e[1L] + P[2L, 2L] * e[2L]) / scale[i + 1L]
    }
  }
  beta
}

.binary_markov_log_likelihood <- function(em, flip_prob, pi1) {
  .binary_markov_scaled_forward(em, .binary_markov_transition(flip_prob, pi1),
                                c(1 - pi1, pi1))$log_likelihood
}

.binary_markov_pair_probabilities <- function(em, P, alpha, beta, scale) {
  n <- nrow(alpha)
  if (n == 1L) {
    return(data.frame(from = integer(), to = integer(), prob_00 = numeric(),
                      prob_01 = numeric(), prob_10 = numeric(), prob_11 = numeric(),
                      prob_flip = numeric()))
  }
  i <- seq_len(n - 1L)
  w0 <- em$E[i + 1L, 1L] * beta[i + 1L, 1L] / scale[i + 1L]
  w1 <- em$E[i + 1L, 2L] * beta[i + 1L, 2L] / scale[i + 1L]
  p00 <- alpha[i, 1L] * P[1L, 1L] * w0
  p01 <- alpha[i, 1L] * P[1L, 2L] * w1
  p10 <- alpha[i, 2L] * P[2L, 1L] * w0
  p11 <- alpha[i, 2L] * P[2L, 2L] * w1
  total <- p00 + p01 + p10 + p11
  data.frame(from = i, to = i + 1L, prob_00 = p00 / total, prob_01 = p01 / total,
             prob_10 = p10 / total, prob_11 = p11 / total,
             prob_flip = (p01 + p10) / total)
}

.binary_markov_viterbi_general <- function(log_emission, P, init) {
  n <- nrow(log_emission)
  log_P <- log(P)
  delta <- matrix(-Inf, n, 2L)
  back <- matrix(1L, n, 2L)
  delta[1L, ] <- log(init) + log_emission[1L, ]
  if (n > 1L) {
    for (i in 2:n) {
      for (state in 1:2) {
        candidate <- delta[i - 1L, ] + log_P[, state]
        best <- which.max(candidate)
        delta[i, state] <- candidate[best] + log_emission[i, state]
        back[i, state] <- best
      }
    }
  }
  path <- integer(n)
  path[n] <- which.max(delta[n, ])
  if (n > 1L) {
    for (i in (n - 1L):1L) path[i] <- back[i + 1L, path[i + 1L]]
  }
  as.integer(path - 1L)
}

# Forward-filtering backward-sampling, vectorized over draws.
.binary_markov_ffbs_sampler <- function(alpha, P) {
  n <- nrow(alpha)
  function(nsamp) {
    if (!is.numeric(nsamp) || length(nsamp) != 1L || !is.finite(nsamp) ||
        nsamp < 1 || nsamp != floor(nsamp)) {
      stop("`nsamp` must be a positive integer.")
    }
    nsamp <- as.integer(nsamp)
    draws <- matrix(0L, nrow = nsamp, ncol = n)
    draws[, n] <- as.integer(stats::runif(nsamp) < alpha[n, 2L])
    if (n > 1L) {
      for (i in (n - 1L):1L) {
        nxt <- draws[, i + 1L] + 1L
        w0 <- alpha[i, 1L] * P[1L, nxt]
        w1 <- alpha[i, 2L] * P[2L, nxt]
        draws[, i] <- as.integer(stats::runif(nsamp) < w1 / (w0 + w1))
      }
    }
    draws
  }
}

#' Define the binary Markov prior
#'
#' @description
#' Creates a prior-state object for an ordered binary sequence with states zero
#' and one. The sequence is a stationary two-state Markov chain whose marginal
#' probability of state one is `pi1` at every position, and whose probability of
#' changing state between adjacent positions is `flip_prob`:
#' \deqn{\pi_1 = P(\theta_i = 1), \qquad q = P(\theta_i \ne \theta_{i-1}).}
#' The transition probabilities are
#' \eqn{p_{01} = P(1 \mid 0) = q / \{2(1-\pi_1)\}} and
#' \eqn{p_{10} = P(0 \mid 1) = q / (2\pi_1)}, and the chain starts at its
#' stationary distribution \eqn{(1-\pi_1, \pi_1)}.
#'
#' `pi1 = 0.5` gives the symmetric prior, with transition matrix
#' \deqn{T(q) = \begin{pmatrix}1-q & q \\ q & 1-q\end{pmatrix}.}
#' Restricting \eqn{q \le 2\pi_1(1-\pi_1)} represents non-negative adjacent
#' association; at the upper bound the positions are iid Bernoulli(`pi1`).
#' The expected run lengths are \eqn{1/p_{01}} for state zero,
#' \eqn{1/p_{10}} for state one and \eqn{1/q} over all runs.
#'
#' @param flip_prob Probability of changing state between adjacent positions.
#'   Must be between `0` and `2 * pi1 * (1 - pi1)`.
#' @param pi1 Marginal probability of state one, strictly between 0 and 1.
#'   The default `0.5` gives the symmetric prior.
#'
#' @return A one-row data frame with class `"BinaryMarkov"` containing
#'   `flip_prob`, `pi1`, the transition probabilities `p01` and `p10`, and the
#'   expected run lengths `expected_run_length` (all runs),
#'   `expected_run_length_zero` and `expected_run_length_one`.
#'
#' @examples
#' BinaryMarkov()
#' BinaryMarkov(flip_prob = 0.05)
#' BinaryMarkov(flip_prob = 0.04, pi1 = 0.2)
#'
#' @export
BinaryMarkov <- function(flip_prob = 0.1, pi1 = 0.5) {
  pi1 <- .validate_binary_markov_pi1(pi1)
  flip_prob <- .validate_binary_markov_flip_given_pi1(flip_prob, pi1)
  P <- .binary_markov_transition(flip_prob, pi1)
  inv <- function(v) if (v == 0) Inf else 1 / v
  out <- data.frame(
    flip_prob = flip_prob,
    pi1 = pi1,
    p01 = P[1L, 2L],
    p10 = P[2L, 1L],
    expected_run_length = inv(flip_prob),
    expected_run_length_zero = inv(P[1L, 2L]),
    expected_run_length_one = inv(P[2L, 1L])
  )
  class(out) <- c("BinaryMarkov", "data.frame")
  out
}

#' Exact binary Markov normal-means solver
#'
#' @description
#' Fits an ordered binary normal-means model with state zero and state one
#' under the binary Markov prior defined by [BinaryMarkov()]:
#' \deqn{x_i\mid\theta_i,s_i\sim N(\theta_i,s_i^2), \qquad \theta_i \in \{0, 1\}.}
#'
#' By default both prior parameters, the marginal probability `pi1` and the
#' switching probability `flip_prob`, are estimated by empirical Bayes. Supply
#' `pi1 = 0.5` for the symmetric prior (only `flip_prob` is then estimated),
#' a numeric `flip_prob` to hold it fixed, or both to fix the prior entirely.
#'
#' Conditional on the prior parameters, inference and each marginal-likelihood
#' evaluation are exact, using a scaled forward-backward algorithm. A single
#' free parameter is searched with [stats::optimize()], and the interval
#' endpoints are evaluated separately. With both parameters free, the search
#' uses [stats::optim()] (L-BFGS-B) from several starts, including the
#' symmetric empirical-Bayes fit, so the result is never worse than the
#' symmetric fit. The likelihood need not be concave, so the numerical search
#' is not guaranteed to find the global maximizer.
#'
#' The function can be used directly as an `ebnm_fn` in flashier: it accepts
#' `g_init`, `fix_g` and `output`, treats `s = Inf` as an uninformative
#' observation (for example a time point with no data) and `s = 0` as an
#' observed state.
#'
#' @param x Numeric vector of ordered observations. At least two are needed to
#'   estimate a parameter.
#' @param s Observation standard errors: one non-negative value or a vector the
#'   same length as `x`. `Inf` marks an uninformative observation; `0` fixes
#'   \eqn{\theta_i = x_i}, which must then be 0 or 1.
#' @param flip_prob Optional fixed switching probability. If `NULL`, it is
#'   estimated.
#' @param pi1 Optional fixed marginal probability of state one. If `NULL` (the
#'   default), it is estimated; `pi1 = 0.5` gives the symmetric prior.
#' @param flip_prob_bounds Two increasing values inside `[0, 0.5]` bounding the
#'   search for `flip_prob`. The search is further limited to
#'   `flip_prob <= 2 * pi1 * (1 - pi1)`.
#' @param pi1_bounds Two increasing values inside `(0, 1)` bounding the search
#'   for `pi1`.
#' @param g_init Optional [BinaryMarkov()] object. Used as the fixed prior when
#'   `fix_g = TRUE` and as a starting value otherwise. An object without `pi1`
#'   (from an earlier version) is taken as symmetric.
#' @param fix_g If `TRUE`, the prior is fixed at `g_init`.
#' @param output `NULL` (default) returns every component listed below. A
#'   character vector, as passed by flashier, returns only the requested
#'   components among `"posterior_mean"`, `"posterior_second_moment"`,
#'   `"posterior_sd"`, `"lfsr"`, `"fitted_g"`, `"log_likelihood"`,
#'   `"posterior_sampler"` and `"data"`; `lfsr` is the posterior probability
#'   of state zero.
#'
#' @return An object with class `c("list", "ebnm")`. With `output = NULL`,
#'   important fields include:
#'   \itemize{
#'     \item `posterior`: exact marginal posterior state probabilities and
#'       moments;
#'     \item `posterior_transitions`: exact posterior probabilities for each of
#'       the four adjacent state pairs and their sum `prob_flip`;
#'     \item `fitted_g`: the fitted or fixed [BinaryMarkov()] state;
#'     \item `log_likelihood`: exact Markov marginal log likelihood, with `df`
#'       equal to the number of estimated parameters;
#'     \item `viterbi_path`: most probable joint binary state sequence;
#'     \item `posterior_sampler`: exact forward-filtering backward-sampling
#'       function returning an `nsamp`-by-`length(x)` matrix.
#'   }
#'
#' @examples
#' set.seed(1)
#' theta <- rep(c(0, 1, 0), times = c(20, 6, 14))
#' x <- stats::rnorm(length(theta), mean = theta, sd = 0.3)
#' fit <- ebnm_binary_markov(x, s = 0.3)
#' fit$fitted_g
#' # symmetric prior
#' ebnm_binary_markov(x, s = 0.3, pi1 = 0.5)$fitted_g
#'
#' @export
ebnm_binary_markov <- function(x,
                               s,
                               flip_prob = NULL,
                               pi1 = NULL,
                               flip_prob_bounds = c(0, 0.5),
                               pi1_bounds = c(1e-4, 1 - 1e-4),
                               g_init = NULL,
                               fix_g = FALSE,
                               output = NULL) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop("`x` must be a non-empty numeric vector.")
  }
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop("`x` values must be finite.")
  }
  n <- length(x)
  em <- .binary_markov_emission(x, s)
  loglik <- function(q, p) .binary_markov_log_likelihood(em, q, p)

  profile <- list(estimated = c(flip_prob = FALSE, pi1 = FALSE))
  if (isTRUE(fix_g)) {
    if (is.null(g_init) || is.null(g_init$flip_prob)) {
      stop("`fix_g = TRUE` requires a `BinaryMarkov` `g_init`.")
    }
    fit_pi1 <- if (is.null(g_init$pi1)) 0.5 else .validate_binary_markov_pi1(g_init$pi1)
    fit_flip <- .validate_binary_markov_flip_given_pi1(g_init$flip_prob, fit_pi1)
  } else {
    if (!is.null(pi1)) pi1 <- .validate_binary_markov_pi1(pi1)
    if (!is.null(flip_prob)) {
      flip_prob <- .validate_binary_markov_flip_given_pi1(
        flip_prob, if (is.null(pi1)) 0.5 else pi1)
    }
    if ((is.null(flip_prob) || is.null(pi1)) && n < 2L) {
      stop("Estimating prior parameters requires at least two observations.")
    }
    fbounds <- .validate_binary_markov_bounds(flip_prob_bounds)
    pbounds <- .validate_binary_markov_pi1_bounds(pi1_bounds)

    # 1-D search for flip_prob at fixed pi1 (exact profile; endpoints evaluated)
    search_flip <- function(p) {
      interval <- c(fbounds[1L], min(fbounds[2L], .binary_markov_flip_upper(p)))
      opt <- stats::optimize(function(q) -loglik(q, p), interval = interval)
      cand <- unique(c(interval[1L], opt$minimum, interval[2L]))
      cand_ll <- vapply(cand, function(q) loglik(q, p), numeric(1))
      list(flip_prob = cand[which.max(cand_ll)], log_likelihood = max(cand_ll),
           interval = interval, optimizer_minimum = opt$minimum,
           optimizer_objective = opt$objective, candidate_flip_prob = cand,
           candidate_log_likelihood = cand_ll)
    }

    if (!is.null(flip_prob) && !is.null(pi1)) {
      fit_flip <- flip_prob
      fit_pi1 <- pi1
    } else if (is.null(flip_prob) && !is.null(pi1)) {
      fit_pi1 <- pi1
      sf <- search_flip(pi1)
      fit_flip <- sf$flip_prob
      profile <- c(list(estimated = c(flip_prob = TRUE, pi1 = FALSE)),
                   sf[c("interval", "optimizer_minimum", "optimizer_objective",
                        "candidate_flip_prob", "candidate_log_likelihood")])
    } else if (!is.null(flip_prob) && is.null(pi1)) {
      # pi1 must satisfy flip_prob <= 2 * pi1 * (1 - pi1)
      half_width <- sqrt(max(1 - 2 * flip_prob, 0)) / 2
      interval <- c(max(pbounds[1L], 0.5 - half_width), min(pbounds[2L], 0.5 + half_width))
      if (interval[1L] >= interval[2L]) interval <- c(0.5, 0.5)
      fit_flip <- flip_prob
      if (interval[1L] == interval[2L]) {
        fit_pi1 <- interval[1L]
      } else {
        opt <- stats::optimize(function(p) -loglik(flip_prob, p), interval = interval)
        cand <- unique(c(interval[1L], opt$minimum, interval[2L]))
        cand_ll <- vapply(cand, function(p) loglik(flip_prob, p), numeric(1))
        fit_pi1 <- cand[which.max(cand_ll)]
        profile <- list(estimated = c(flip_prob = FALSE, pi1 = TRUE), interval = interval,
                        optimizer_minimum = opt$minimum, optimizer_objective = opt$objective,
                        candidate_pi1 = cand, candidate_log_likelihood = cand_ll)
      }
    } else {
      # Both free: optimize on an unconstrained scale that respects
      # pi1_bounds and flip_prob <= min(upper bound, 2 * pi1 * (1 - pi1)).
      to_par <- function(z) {
        p <- pbounds[1L] + (pbounds[2L] - pbounds[1L]) * stats::plogis(z[1L])
        hi <- min(fbounds[2L], .binary_markov_flip_upper(p))
        q <- if (hi <= fbounds[1L]) hi else fbounds[1L] + (hi - fbounds[1L]) * stats::plogis(z[2L])
        c(flip_prob = q, pi1 = p)
      }
      to_z <- function(q, p) {
        clamp <- function(u) min(max(u, 1e-6), 1 - 1e-6)
        hi <- min(fbounds[2L], .binary_markov_flip_upper(p))
        c(stats::qlogis(clamp((p - pbounds[1L]) / (pbounds[2L] - pbounds[1L]))),
          stats::qlogis(clamp((q - fbounds[1L]) / max(hi - fbounds[1L], 1e-12))))
      }
      nll <- function(z) {
        v <- to_par(z)
        -loglik(v[["flip_prob"]], v[["pi1"]])
      }
      symmetric <- if (pbounds[1L] <= 0.5 && pbounds[2L] >= 0.5) search_flip(0.5) else NULL
      clamp_p <- function(p) min(max(p, pbounds[1L]), pbounds[2L])
      p_data <- clamp_p(min(max(mean(x[is.finite(em$s)] > 0.5), 0.05), 0.95))
      q_start <- if (is.null(symmetric)) 0.05 else max(symmetric$flip_prob, 1e-3)
      starts <- list(c(q_start, clamp_p(0.5)), c(min(q_start, .binary_markov_flip_upper(p_data) / 2), p_data))
      if (!is.null(g_init) && !is.null(g_init$flip_prob)) {
        starts <- c(list(c(g_init$flip_prob, clamp_p(if (is.null(g_init$pi1)) 0.5 else g_init$pi1))),
                    starts)
      }
      runs <- lapply(starts, function(st) tryCatch(
        stats::optim(to_z(st[1L], st[2L]), nll, method = "L-BFGS-B",
                     lower = c(-20, -20), upper = c(20, 20)),
        error = function(e) NULL))
      runs <- Filter(Negate(is.null), runs)
      cand <- lapply(runs, function(r) to_par(r$par))
      cand_ll <- vapply(cand, function(v) loglik(v[["flip_prob"]], v[["pi1"]]), numeric(1))
      if (!is.null(symmetric)) {
        cand <- c(cand, list(c(flip_prob = symmetric$flip_prob, pi1 = 0.5)))
        cand_ll <- c(cand_ll, symmetric$log_likelihood)
      }
      best <- which.max(cand_ll)
      fit_flip <- cand[[best]][["flip_prob"]]
      fit_pi1 <- cand[[best]][["pi1"]]
      profile <- list(estimated = c(flip_prob = TRUE, pi1 = TRUE), method = "L-BFGS-B",
                      candidates = data.frame(flip_prob = vapply(cand, `[[`, 0, "flip_prob"),
                                              pi1 = vapply(cand, `[[`, 0, "pi1"),
                                              log_likelihood = cand_ll),
                      symmetric_log_likelihood = if (is.null(symmetric)) NA_real_ else symmetric$log_likelihood)
    }
  }

  P <- .binary_markov_transition(fit_flip, fit_pi1)
  init <- c(1 - fit_pi1, fit_pi1)
  fwd <- .binary_markov_scaled_forward(em, P, init)
  beta <- .binary_markov_scaled_backward(em, P, fwd$scale)
  gamma <- fwd$alpha * beta
  prob_one <- gamma[, 2L] / rowSums(gamma)
  fitted_g <- BinaryMarkov(fit_flip, fit_pi1)
  n_estimated <- sum(profile$estimated)

  if (!is.null(output)) {
    out <- list()
    if (any(c("posterior_mean", "posterior_second_moment", "posterior_sd", "lfsr") %in% output)) {
      post <- data.frame(mean = prob_one, second_moment = prob_one,
                         sd = sqrt(prob_one * (1 - prob_one)))
      if ("lfsr" %in% output) post$lfsr <- 1 - prob_one
      out$posterior <- post
    }
    if ("fitted_g" %in% output) out$fitted_g <- fitted_g
    if ("log_likelihood" %in% output) out$log_likelihood <- fwd$log_likelihood
    if ("posterior_sampler" %in% output) out$posterior_sampler <- .binary_markov_ffbs_sampler(fwd$alpha, P)
    if ("data" %in% output) out$data <- data.frame(x = x, s = em$s)
    return(structure(out, class = c("list", "ebnm")))
  }

  symmetric_prior <- isTRUE(all.equal(fit_pi1, 0.5)) && !isTRUE(profile$estimated[["pi1"]])
  family <- if (symmetric_prior) "symmetric" else "asymmetric"
  structure(
    list(
      data = data.frame(x = x, s = em$s),
      posterior = data.frame(
        mean = prob_one,
        var = prob_one * (1 - prob_one),
        second_moment = prob_one,
        prob_zero = 1 - prob_one,
        prob_one = prob_one
      ),
      posterior_transitions = .binary_markov_pair_probabilities(em, P, fwd$alpha, beta, fwd$scale),
      fitted_g = fitted_g,
      initial_probability = c(intron = 1 - fit_pi1, exon = fit_pi1),
      transition_matrix = `dimnames<-`(P, list(from = c("0", "1"), to = c("0", "1"))),
      log_likelihood = structure(fwd$log_likelihood, class = "logLik",
                                 df = as.integer(n_estimated), nobs = n),
      log_likelihood_semantics = paste0("exact_binary_markov_", family,
                                        if (n_estimated > 0L) "_eb" else "_fixed"),
      viterbi_path = .binary_markov_viterbi_general(em$log_emission, P, init),
      posterior_sampler = .binary_markov_ffbs_sampler(fwd$alpha, P),
      profile_optimization = profile,
      prior_family = paste0("binary_markov_", family)
    ),
    class = c("list", "ebnm")
  )
}
