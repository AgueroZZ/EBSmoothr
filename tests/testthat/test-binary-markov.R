.binary_markov_enumeration <- function(x, s, flip_prob) {
  n <- length(x)
  states <- as.matrix(expand.grid(rep(list(0:1), n)))
  storage.mode(states) <- "integer"

  log_prior <- rep(log(0.5), nrow(states))
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      flipped <- states[, i] != states[, i + 1L]
      log_prior <- log_prior + ifelse(
        flipped,
        log(flip_prob),
        log1p(-flip_prob)
      )
    }
  }

  log_weight <- log_prior
  for (i in seq_len(n)) {
    log_weight <- log_weight + stats::dnorm(
      x[i],
      mean = states[, i],
      sd = s[i],
      log = TRUE
    )
  }

  max_log_weight <- max(log_weight)
  weight <- exp(log_weight - max_log_weight)
  normalized_weight <- weight / sum(weight)
  log_likelihood <- max_log_weight + log(sum(weight))

  posterior_prob_one <- as.numeric(crossprod(normalized_weight, states))
  transitions <- if (n == 1L) {
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
    do.call(rbind, lapply(seq_len(n - 1L), function(i) {
      state_i <- states[, i]
      state_next <- states[, i + 1L]
      prob_00 <- sum(normalized_weight[state_i == 0L & state_next == 0L])
      prob_01 <- sum(normalized_weight[state_i == 0L & state_next == 1L])
      prob_10 <- sum(normalized_weight[state_i == 1L & state_next == 0L])
      prob_11 <- sum(normalized_weight[state_i == 1L & state_next == 1L])
      data.frame(
        from = i,
        to = i + 1L,
        prob_00 = prob_00,
        prob_01 = prob_01,
        prob_10 = prob_10,
        prob_11 = prob_11,
        prob_flip = prob_01 + prob_10
      )
    }))
  }

  list(
    log_likelihood = log_likelihood,
    posterior_prob_one = posterior_prob_one,
    transitions = transitions,
    viterbi_path = as.integer(states[which.max(log_weight), ])
  )
}

test_that("fixed symmetric binary Markov inference matches enumeration", {
  x <- c(-0.1, 0.25, 0.82, 1.15, 0.55, 0.05)
  s <- c(0.25, 0.4, 0.3, 0.2, 0.5, 0.35)
  flip_prob <- 0.17

  fit <- ebnm_binary_markov(x, s, flip_prob = flip_prob)
  oracle <- .binary_markov_enumeration(x, s, flip_prob)

  expect_s3_class(fit, "ebnm")
  expect_s3_class(fit$fitted_g, "BinaryMarkov")
  expect_equal(fit$fitted_g$flip_prob, flip_prob)
  expect_equal(as.numeric(fit$log_likelihood), oracle$log_likelihood, tolerance = 1e-10)
  expect_equal(fit$posterior$prob_one, oracle$posterior_prob_one, tolerance = 1e-10)
  expect_equal(fit$posterior$mean, oracle$posterior_prob_one, tolerance = 1e-10)
  expect_equal(
    fit$posterior$var,
    oracle$posterior_prob_one * (1 - oracle$posterior_prob_one),
    tolerance = 1e-10
  )
  expect_equal(fit$posterior_transitions, oracle$transitions, tolerance = 1e-10)
  expect_equal(fit$viterbi_path, oracle$viterbi_path)
  expect_equal(attr(fit$log_likelihood, "df"), 0L)
  expect_equal(fit$log_likelihood_semantics, "exact_binary_markov_symmetric_fixed")
})

test_that("symmetric binary Markov model reduces to iid at flip probability one half", {
  x <- c(-0.3, 0.1, 0.45, 0.7, 1.2)
  s <- c(0.3, 0.25, 0.5, 0.35, 0.2)

  fit_markov <- ebnm_binary_markov(x, s, flip_prob = 0.5)
  expected_prob_one <- stats::plogis((2 * x - 1) / (2 * s^2))
  expected_log_likelihood <- sum(log(
    0.5 * stats::dnorm(x, mean = 0, sd = s) +
      0.5 * stats::dnorm(x, mean = 1, sd = s)
  ))
  expected_posterior <- data.frame(
    mean = expected_prob_one,
    var = expected_prob_one * (1 - expected_prob_one),
    second_moment = expected_prob_one,
    prob_zero = 1 - expected_prob_one,
    prob_one = expected_prob_one
  )

  expect_equal(
    as.numeric(fit_markov$log_likelihood),
    expected_log_likelihood,
    tolerance = 1e-12
  )
  expect_equal(fit_markov$posterior, expected_posterior, tolerance = 1e-12)
  expect_equal(fit_markov$viterbi_path, as.integer(expected_prob_one >= 0.5))
  expect_false("Binary" %in% getNamespaceExports("EBSmoothr"))
  expect_false("ebnm_binary" %in% getNamespaceExports("EBSmoothr"))
})

test_that("symmetric binary Markov empirical Bayes fit searches the exact likelihood", {
  set.seed(401)

  n <- 180L
  true_flip_prob <- 0.08
  theta <- integer(n)
  theta[1] <- stats::rbinom(1, size = 1, prob = 0.5)
  for (i in 2:n) {
    theta[i] <- if (stats::runif(1) < true_flip_prob) 1L - theta[i - 1L] else theta[i - 1L]
  }
  s <- 0.3
  x <- stats::rnorm(n, mean = theta, sd = s)

  fit <- ebnm_binary_markov(x, s)
  grid <- seq(0.001, 0.5, length.out = 101)
  grid_log_likelihood <- vapply(
    grid,
    function(q) as.numeric(ebnm_binary_markov(x, s, flip_prob = q)$log_likelihood),
    numeric(1)
  )

  expect_true(fit$fitted_g$flip_prob >= 0 && fit$fitted_g$flip_prob <= 0.5)
  expect_gte(as.numeric(fit$log_likelihood), max(grid_log_likelihood) - 1e-6)
  expect_lt(
    abs(fit$fitted_g$flip_prob - grid[which.max(grid_log_likelihood)]),
    0.01
  )
  expect_equal(attr(fit$log_likelihood, "df"), 1L)
  expect_equal(fit$log_likelihood_semantics, "exact_binary_markov_symmetric_eb")
  expect_true(isTRUE(fit$profile_optimization$estimated))
  expect_false("initial" %in% names(fit$profile_optimization))
  expect_equal(
    fit$profile_optimization$candidate_flip_prob,
    unique(c(
      fit$profile_optimization$interval[1L],
      fit$profile_optimization$optimizer_minimum,
      fit$profile_optimization$interval[2L]
    ))
  )
})

test_that("symmetric binary Markov posterior sampler matches exact marginals", {
  x <- c(0.05, 0.25, 0.8, 1.1, 0.65, 0.1)
  fit <- ebnm_binary_markov(x, s = 0.35, flip_prob = 0.12)

  set.seed(402)
  draws <- fit$posterior_sampler(20000)

  expect_equal(dim(draws), c(20000L, length(x)))
  expect_true(all(draws %in% c(0, 1)))
  expect_equal(colMeans(draws), fit$posterior$prob_one, tolerance = 0.015)
  expect_error(fit$posterior_sampler(0), "positive integer")
  expect_error(fit$posterior_sampler(1.5), "positive integer")
})

test_that("symmetric binary Markov solver handles boundaries and validates inputs", {
  expect_equal(BinaryMarkov()$flip_prob, 0.1)
  expect_error(BinaryMarkov(c(0.1, 0.2)), "numeric scalar")
  expect_error(BinaryMarkov(NA_real_), "finite numeric scalar")
  expect_error(BinaryMarkov(-0.1), "between 0 and 0.5")
  expect_error(BinaryMarkov(0.6), "between 0 and 0.5")

  fit_constant <- ebnm_binary_markov(c(0.1, 0.2, 0.9), s = 0.3, flip_prob = 0)
  oracle_constant <- .binary_markov_enumeration(c(0.1, 0.2, 0.9), rep(0.3, 3), 0)
  expect_equal(
    as.numeric(fit_constant$log_likelihood),
    oracle_constant$log_likelihood,
    tolerance = 1e-12
  )
  expect_equal(fit_constant$posterior$prob_one, oracle_constant$posterior_prob_one, tolerance = 1e-12)

  fit_single <- ebnm_binary_markov(0.7, s = 0.2, flip_prob = 0.1)
  expected_prob_one <- stats::plogis((2 * 0.7 - 1) / (2 * 0.2^2))
  expect_equal(fit_single$posterior$prob_one, expected_prob_one, tolerance = 1e-12)
  expect_equal(nrow(fit_single$posterior_transitions), 0L)

  expect_error(ebnm_binary_markov(0.7, s = 0.2), "at least two observations")
  expect_error(ebnm_binary_markov(c(0, 1), s = 0.2, flip_prob = -0.1), "between 0 and 0.5")
  expect_error(ebnm_binary_markov(c(0, 1), s = 0.2, flip_prob = 0.6), "between 0 and 0.5")
  expect_error(
    ebnm_binary_markov(c(0, 1), s = 0.2, flip_prob_bounds = c(0.2, 0.1)),
    "in increasing order"
  )
  expect_error(
    ebnm_binary_markov(c(0, 1), s = 0.2, flip_prob_bounds = c(-0.1, 0.5)),
    "inside"
  )
})
