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

#' Define the symmetric binary Markov prior
#'
#' @description
#' Creates a prior-state object for an ordered binary sequence with states zero
#' and one, initial probabilities `(0.5, 0.5)`, and symmetric transition matrix
#' \deqn{T(q) = \begin{pmatrix}1-q & q \\ q & 1-q\end{pmatrix}.}
#' The parameter `flip_prob` is \eqn{q=P(\theta_i\ne\theta_{i-1})}. Restricting
#' \eqn{q\le 0.5} represents non-negative adjacent association. For \eqn{q>0},
#' the expected run length for either state is \eqn{1/q}; \eqn{q=0.5} reduces
#' exactly to the iid binary prior with probability one half.
#'
#' @param flip_prob Symmetric probability of changing state between adjacent
#'   positions. Must be between `0` and `0.5`.
#'
#' @return A one-row data frame with class `"BinaryMarkov"` containing
#'   `flip_prob` and `expected_run_length`.
#'
#' @examples
#' BinaryMarkov()
#' BinaryMarkov(flip_prob = 0.05)
#'
#' @export
BinaryMarkov <- function(flip_prob = 0.1) {
  flip_prob <- .validate_binary_markov_flip_prob(flip_prob)
  out <- data.frame(
    flip_prob = flip_prob,
    expected_run_length = if (flip_prob == 0) Inf else 1 / flip_prob
  )
  class(out) <- c("BinaryMarkov", "data.frame")
  out
}

#' Exact symmetric binary Markov normal-means solver
#'
#' @description
#' Fits an ordered binary normal-means model with intron state zero, exon state
#' one, fixed initial probabilities `(0.5, 0.5)`, and the symmetric transition
#' prior defined by [BinaryMarkov()]. Observations follow
#' \deqn{x_i\mid\theta_i,s_i\sim N(\theta_i,s_i^2).}
#'
#' Conditional on a fixed transition probability, inference and each marginal-
#' likelihood evaluation are exact and use a log-space forward-backward
#' algorithm. If `flip_prob = NULL`, the function numerically searches
#' `flip_prob_bounds` with [stats::optimize()] and separately evaluates both
#' interval endpoints. The likelihood need not be concave, so this numerical
#' search is not guaranteed to find the global maximizer. If a numeric
#' `flip_prob` is supplied, the transition parameter is held fixed.
#'
#' @param x Numeric vector of at least two ordered observations when estimating
#'   `flip_prob`. A single observation is allowed for a fixed-parameter fit.
#' @param s Known observation standard errors. Supply one positive value or a
#'   positive vector with the same length as `x`.
#' @param flip_prob Optional fixed symmetric transition probability. If `NULL`,
#'   search numerically using exact marginal-likelihood evaluations.
#' @param flip_prob_bounds Two increasing values defining the empirical-Bayes
#'   numerical-search interval. Both must lie inside `[0, 0.5]`.
#'
#' @return An object with class `c("list", "ebnm")`. Important fields include:
#'   \itemize{
#'     \item `posterior`: exact marginal posterior state probabilities and
#'       moments;
#'     \item `posterior_transitions`: exact posterior probabilities for each of
#'       the four adjacent state pairs and their sum `prob_flip`;
#'     \item `fitted_g`: the fitted or fixed [BinaryMarkov()] state;
#'     \item `log_likelihood`: exact Markov marginal log likelihood;
#'     \item `viterbi_path`: most probable joint binary state sequence;
#'     \item `posterior_sampler`: exact forward-filtering backward-sampling
#'       function returning an `nsamp`-by-`length(x)` matrix.
#'   }
#'
#' @examples
#' set.seed(1)
#' theta <- rep(c(0, 1), each = 6)
#' x <- stats::rnorm(length(theta), mean = theta, sd = 0.2)
#' fit <- ebnm_binary_markov(x, s = 0.2)
#' fit$fitted_g
#' fit$posterior
#'
#' @export
ebnm_binary_markov <- function(x,
                               s,
                               flip_prob = NULL,
                               flip_prob_bounds = c(0, 0.5)) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop("`x` must be a non-empty numeric vector.")
  }
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop("`x` values must be finite.")
  }
  s <- .normalize_binary_s(s, length(x))

  log_emission <- cbind(
    stats::dnorm(x, mean = 0, sd = s, log = TRUE),
    stats::dnorm(x, mean = 1, sd = s, log = TRUE)
  )

  estimated <- is.null(flip_prob)
  if (estimated) {
    if (length(x) < 2L) {
      stop("Estimating `flip_prob` requires at least two observations.")
    }
    bounds <- .validate_binary_markov_bounds(flip_prob_bounds)

    objective <- function(q) {
      -.binary_markov_forward(log_emission, q)$log_likelihood
    }
    optimization <- stats::optimize(objective, interval = bounds)
    candidate_flip_prob <- unique(c(
      bounds[1L],
      optimization$minimum,
      bounds[2L]
    ))
    candidate_log_likelihood <- vapply(
      candidate_flip_prob,
      function(q) .binary_markov_forward(log_emission, q)$log_likelihood,
      numeric(1)
    )
    best <- which.max(candidate_log_likelihood)
    fitted_flip_prob <- candidate_flip_prob[best]
    profile_optimization <- list(
      estimated = TRUE,
      interval = bounds,
      optimizer_minimum = optimization$minimum,
      optimizer_objective = optimization$objective,
      candidate_flip_prob = candidate_flip_prob,
      candidate_log_likelihood = candidate_log_likelihood
    )
  } else {
    fitted_flip_prob <- .validate_binary_markov_flip_prob(flip_prob)
    profile_optimization <- list(
      estimated = FALSE,
      interval = c(fitted_flip_prob, fitted_flip_prob),
      candidate_flip_prob = fitted_flip_prob
    )
  }

  inference <- .binary_markov_forward_backward(log_emission, fitted_flip_prob)
  prob_one <- as.numeric(inference$state_probability[, 2L])
  posterior <- data.frame(
    mean = prob_one,
    var = prob_one * (1 - prob_one),
    second_moment = prob_one,
    prob_zero = 1 - prob_one,
    prob_one = prob_one
  )
  fitted_g <- BinaryMarkov(fitted_flip_prob)
  transition_matrix <- matrix(
    c(1 - fitted_flip_prob, fitted_flip_prob,
      fitted_flip_prob, 1 - fitted_flip_prob),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(from = c("0", "1"), to = c("0", "1"))
  )
  log_likelihood <- structure(
    inference$log_likelihood,
    class = "logLik",
    df = if (estimated) 1L else 0L,
    nobs = length(x)
  )

  structure(
    list(
      data = data.frame(x = x, s = s),
      posterior = posterior,
      posterior_transitions = inference$posterior_transitions,
      fitted_g = fitted_g,
      initial_probability = c(intron = 0.5, exon = 0.5),
      transition_matrix = transition_matrix,
      log_likelihood = log_likelihood,
      log_likelihood_semantics = if (estimated) {
        "exact_binary_markov_symmetric_eb"
      } else {
        "exact_binary_markov_symmetric_fixed"
      },
      viterbi_path = .binary_markov_viterbi(
        log_emission,
        inference$log_transition
      ),
      posterior_sampler = .binary_markov_posterior_sampler(
        inference$log_alpha,
        inference$log_transition,
        inference$state_probability
      ),
      profile_optimization = profile_optimization,
      prior_family = "binary_markov_symmetric"
    ),
    class = c("list", "ebnm")
  )
}
