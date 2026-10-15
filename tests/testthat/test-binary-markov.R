.binary_markov_enumeration <- function(x, s, flip_prob, pi1 = 0.5) {
  n <- length(x)
  states <- as.matrix(expand.grid(rep(list(0:1), n)))
  storage.mode(states) <- "integer"

  p01 <- flip_prob / (2 * (1 - pi1))
  p10 <- flip_prob / (2 * pi1)
  log_P <- log(matrix(c(1 - p01, p01, p10, 1 - p10), 2L, byrow = TRUE))
  log_prior <- ifelse(states[, 1L] == 1L, log(pi1), log(1 - pi1))
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      log_prior <- log_prior + log_P[cbind(states[, i] + 1L, states[, i + 1L] + 1L)]
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

  fit <- ebnm_binary_markov(x, s, flip_prob = flip_prob, pi1 = 0.5)
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

  fit_markov <- ebnm_binary_markov(x, s, flip_prob = 0.5, pi1 = 0.5)
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

  fit <- ebnm_binary_markov(x, s, pi1 = 0.5)
  grid <- seq(0.001, 0.5, length.out = 101)
  grid_log_likelihood <- vapply(
    grid,
    function(q) as.numeric(ebnm_binary_markov(x, s, flip_prob = q, pi1 = 0.5)$log_likelihood),
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
  expect_equal(fit$profile_optimization$estimated, c(flip_prob = TRUE, pi1 = FALSE))
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
  fit <- ebnm_binary_markov(x, s = 0.35, flip_prob = 0.12, pi1 = 0.5)

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

  fit_constant <- ebnm_binary_markov(c(0.1, 0.2, 0.9), s = 0.3, flip_prob = 0, pi1 = 0.5)
  oracle_constant <- .binary_markov_enumeration(c(0.1, 0.2, 0.9), rep(0.3, 3), 0)
  expect_equal(
    as.numeric(fit_constant$log_likelihood),
    oracle_constant$log_likelihood,
    tolerance = 1e-12
  )
  expect_equal(fit_constant$posterior$prob_one, oracle_constant$posterior_prob_one, tolerance = 1e-12)

  fit_single <- ebnm_binary_markov(0.7, s = 0.2, flip_prob = 0.1, pi1 = 0.5)
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


test_that("asymmetric binary Markov inference matches enumeration", {
  x <- c(-0.1, 0.25, 0.82, 1.15, 0.55, 0.05, 0.3)
  s <- c(0.25, 0.4, 0.3, 0.2, 0.5, 0.35, 0.3)
  for (par in list(c(0.17, 0.2), c(0.05, 0.7), c(2 * 0.3 * 0.7, 0.3))) {
    fit <- ebnm_binary_markov(x, s, flip_prob = par[1], pi1 = par[2])
    oracle <- .binary_markov_enumeration(x, s, par[1], par[2])
    expect_equal(as.numeric(fit$log_likelihood), oracle$log_likelihood, tolerance = 1e-10)
    expect_equal(fit$posterior$prob_one, oracle$posterior_prob_one, tolerance = 1e-10)
    expect_equal(fit$posterior_transitions, oracle$transitions, tolerance = 1e-10)
    expect_equal(fit$viterbi_path, oracle$viterbi_path)
    expect_equal(fit$log_likelihood_semantics, "exact_binary_markov_asymmetric_fixed")
    expect_equal(unname(fit$initial_probability), c(1 - par[2], par[2]))
  }
})

test_that("the symmetric fixed prior matches the EBSmoothr 0.3.5 implementation", {
  set.seed(403)
  x <- c(rep(0, 30), rep(1, 12), rep(0, 18)) + stats::rnorm(60, sd = 0.45)
  for (q in c(0, 0.03, 0.2, 0.5)) {
    fit <- ebnm_binary_markov(x, 0.45, flip_prob = q, pi1 = 0.5)
    log_emission <- cbind(stats::dnorm(x, 0, 0.45, log = TRUE), stats::dnorm(x, 1, 0.45, log = TRUE))
    legacy <- .binary_markov_forward_backward(log_emission, q)
    expect_equal(as.numeric(fit$log_likelihood), legacy$log_likelihood, tolerance = 1e-10)
    expect_equal(fit$posterior$prob_one, unname(legacy$state_probability[, 2L]), tolerance = 1e-10)
  }
})

test_that("default fit estimates pi1 and is never worse than the symmetric fit", {
  set.seed(404)
  n <- 400L
  z <- integer(n)
  z[1] <- 0L
  for (i in 2:n) z[i] <- if (z[i - 1L] == 0L) stats::rbinom(1, 1, 0.02) else 1L - stats::rbinom(1, 1, 0.12)
  x <- z + stats::rnorm(n, sd = 0.4)

  fit <- ebnm_binary_markov(x, 0.4)
  sym <- ebnm_binary_markov(x, 0.4, pi1 = 0.5)
  expect_equal(fit$profile_optimization$estimated, c(flip_prob = TRUE, pi1 = TRUE))
  expect_equal(attr(fit$log_likelihood, "df"), 2L)
  expect_equal(fit$log_likelihood_semantics, "exact_binary_markov_asymmetric_eb")
  expect_gte(as.numeric(fit$log_likelihood), as.numeric(sym$log_likelihood) - 1e-8)
  expect_lt(fit$fitted_g$pi1, 0.4)
  expect_lt(abs(fit$fitted_g$pi1 - mean(z)), 0.1)

  # a coarse grid over (flip_prob, pi1) does not beat the fitted likelihood
  grid <- expand.grid(q = seq(0.005, 0.2, length.out = 25), p = seq(0.05, 0.6, length.out = 23))
  grid <- grid[grid$q <= 2 * grid$p * (1 - grid$p), ]
  grid_ll <- mapply(function(q, p) as.numeric(ebnm_binary_markov(x, 0.4, flip_prob = q, pi1 = p)$log_likelihood),
                    grid$q, grid$p)
  expect_gte(as.numeric(fit$log_likelihood), max(grid_ll) - 1e-6)
})

test_that("pi1 can be estimated with flip_prob fixed", {
  set.seed(405)
  x <- c(rep(0, 70), rep(1, 15), rep(0, 40), rep(1, 10)) + stats::rnorm(135, sd = 0.4)
  fit <- ebnm_binary_markov(x, 0.4, flip_prob = 0.02)
  expect_equal(fit$profile_optimization$estimated, c(flip_prob = FALSE, pi1 = TRUE))
  expect_equal(fit$fitted_g$flip_prob, 0.02)
  grid <- seq(0.02, 0.98, by = 0.01)
  grid_ll <- vapply(grid, function(p) as.numeric(ebnm_binary_markov(x, 0.4, flip_prob = 0.02, pi1 = p)$log_likelihood), 0)
  expect_gte(as.numeric(fit$log_likelihood), max(grid_ll) - 1e-6)
})

test_that("asymmetric posterior sampler matches exact marginals", {
  x <- c(0.05, 0.25, 0.8, 1.1, 0.65, 0.1, -0.2)
  fit <- ebnm_binary_markov(x, s = 0.35, flip_prob = 0.1, pi1 = 0.25)
  set.seed(406)
  draws <- fit$posterior_sampler(20000)
  expect_equal(colMeans(draws), fit$posterior$prob_one, tolerance = 0.015)
  pair_flip <- colMeans(draws[, -1L] != draws[, -ncol(draws)])
  expect_equal(pair_flip, fit$posterior_transitions$prob_flip, tolerance = 0.015)
})

test_that("BinaryMarkov parametrization and validation", {
  g <- BinaryMarkov(0.04, pi1 = 0.2)
  expect_equal(g$p01, 0.04 / (2 * 0.8))
  expect_equal(g$p10, 0.04 / (2 * 0.2))
  expect_equal((1 - g$pi1) * g$p01 + g$pi1 * g$p10, g$flip_prob)
  expect_equal(g$expected_run_length_one, 1 / g$p10)
  sym <- BinaryMarkov(0.1)
  expect_equal(c(sym$p01, sym$p10), c(0.1, 0.1))
  expect_equal(sym$expected_run_length, 10)
  expect_error(BinaryMarkov(0.33, pi1 = 0.2), "between 0 and 0.32")
  expect_error(BinaryMarkov(0.1, pi1 = 1), "strictly between 0 and 1")
  expect_error(ebnm_binary_markov(c(0, 1), 0.2, pi1_bounds = c(0, 1)), "inside \\(0, 1\\)")
})

test_that("infinite and zero standard errors and output subsets", {
  set.seed(407)
  x <- c(rep(0, 20), rep(1, 20)) + stats::rnorm(40, sd = 0.3)
  s <- rep(0.3, 40)
  s[18:23] <- Inf
  fit <- ebnm_binary_markov(x, s)
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
  expect_true(all(fit$posterior$prob_one[18:23] > 0 & fit$posterior$prob_one[18:23] < 1))
  expect_equal(fit$data$s[18], Inf)
  # no information at all: zero log likelihood and posterior equal to the prior
  none <- ebnm_binary_markov(x, Inf, flip_prob = 0.05, pi1 = 0.3)
  expect_equal(as.numeric(none$log_likelihood), 0, tolerance = 1e-12)
  expect_equal(none$posterior$prob_one, rep(0.3, 40), tolerance = 1e-12)

  s0 <- rep(0.3, 40); s0[c(3, 30)] <- 0
  x0 <- x; x0[3] <- 0; x0[30] <- 1
  fit0 <- ebnm_binary_markov(x0, s0, flip_prob = 0.05, pi1 = 0.5)
  expect_equal(fit0$posterior$prob_one[c(3, 30)], c(0, 1))
  expect_error(ebnm_binary_markov(c(0.5, 1), c(0, 0.3), flip_prob = 0.1, pi1 = 0.5), "equal to 0 or 1")
  expect_error(ebnm_binary_markov(c(0, 1), c(-1, 0.3)), "non-negative")

  out <- ebnm_binary_markov(x, 0.3, output = c("posterior_mean", "posterior_second_moment",
                                               "fitted_g", "log_likelihood"))
  full <- ebnm_binary_markov(x, 0.3)
  expect_setequal(names(out), c("posterior", "fitted_g", "log_likelihood"))
  expect_equal(out$posterior$mean, full$posterior$mean)
  expect_equal(out$log_likelihood, as.numeric(full$log_likelihood))
  lf <- ebnm_binary_markov(x, 0.3, g_init = full$fitted_g, fix_g = TRUE, output = "lfsr")
  expect_equal(lf$posterior$lfsr, full$posterior$prob_zero)
  expect_s3_class(out, "ebnm")
})

test_that("fix_g accepts g_init from EBSmoothr <= 0.3.5 (no pi1)", {
  x <- c(0.1, 0.9, 1.1, 0.2)
  old_g <- structure(data.frame(flip_prob = 0.1, expected_run_length = 10),
                     class = c("BinaryMarkov", "data.frame"))
  fit <- ebnm_binary_markov(x, 0.3, g_init = old_g, fix_g = TRUE)
  ref <- ebnm_binary_markov(x, 0.3, flip_prob = 0.1, pi1 = 0.5)
  expect_equal(as.numeric(fit$log_likelihood), as.numeric(ref$log_likelihood))
  expect_error(ebnm_binary_markov(x, 0.3, fix_g = TRUE), "requires")
})

test_that("ebnm_binary_markov works as a flashier ebnm_fn", {
  skip_if_not_installed("flashier")
  set.seed(408)
  n <- 120L; p <- 40L
  l <- integer(n); l[1] <- 0L
  for (i in 2:n) l[i] <- if (stats::runif(1) < 0.04) 1L - l[i - 1L] else l[i - 1L]
  if (all(l == l[1])) l[60:80] <- 1L - l[1]
  f <- stats::rnorm(p)
  Y <- outer(l, f) + matrix(stats::rnorm(n * p, sd = 0.5), n, p)
  Y[sample(length(Y), 0.5 * length(Y))] <- NA
  Y[50:52, ] <- NA  # fully missing time points: flashier passes s = Inf
  fl <- flashier::flash(Y, ebnm_fn = list(ebnm_binary_markov, ebnm::ebnm_normal),
                        greedy_Kmax = 1, backfit = FALSE, verbose = 0)
  expect_equal(fl$n_factors, 1L)
  expect_s3_class(fl$L_ghat[[1]], "BinaryMarkov")
  expect_gt(max(abs(stats::cor(fl$L_pm[, 1], l))), 0.95)
  expect_true(all(fl$L_pm >= 0 & fl$L_pm <= 1))
})
