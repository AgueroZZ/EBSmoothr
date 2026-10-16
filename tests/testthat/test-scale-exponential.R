## ebnm_scale_exponential() (R/08_ebnm_scale_exponential.R).

# One exponential component by numerical integration over u = theta / s.
integrate_component <- function(x, s, sigma) {
  f <- function(u) dexp(u * s, 1 / sigma) * s * dnorm(x, u * s, s)
  brk <- c(0, 50 * min(sigma / s, 1), max(x / s, 0) + 40 + 40 * sigma / s)
  I <- function(g) sum(vapply(1:2, function(i) {
    stats::integrate(g, brk[i], brk[i + 1], rel.tol = 1e-12, abs.tol = 0)$value
  }, numeric(1)))
  p0 <- I(f)
  c(ll = log(p0), m1 = I(function(u) u * s * f(u)) / p0,
    m2 = I(function(u) (u * s)^2 * f(u)) / p0)
}

test_that("component likelihood and moments match numerical integration, also for s >> scale", {
  for (ratio in c(1e-3, 1, 4.9, 5.1, 1e3, 1e9)) {
    for (x in c(-3, 0, 2, 10)) {
      sigma <- 1 / ratio
      ref <- integrate_component(x, 1, sigma)
      ll <- .expmix_loglik_matrix(x, 1, 0, sigma)[, 2]
      mom <- .expmix_trunc_moments(1 / sigma - x)
      expect_equal(ll, unname(ref["ll"]), tolerance = 1e-10)
      expect_equal(mom$h, unname(ref["m1"]), tolerance = 1e-10)
      expect_equal(mom$m2, unname(ref["m2"]), tolerance = 1e-10)
    }
  }
})

test_that("with a one-scale prior it reproduces ebnm_point_exponential()", {
  set.seed(1)
  x <- c(rep(0, 150), rexp(100, 1 / 3)) + rnorm(250)
  pe <- ebnm::ebnm_point_exponential(x, 1, output = ebnm::ebnm_output_all())
  fit <- ebnm_scale_exponential(x, 1, g_init = pe$fitted_g, fix_g = TRUE,
                                output = ebnm::ebnm_output_all())
  expect_equal(as.numeric(fit$log_likelihood), as.numeric(pe$log_likelihood), tolerance = 1e-10)
  expect_equal(fit$posterior$mean, pe$posterior$mean, tolerance = 1e-10)
  expect_equal(fit$posterior$sd, pe$posterior$sd, tolerance = 1e-8)
  expect_equal(fit$posterior$lfsr, pe$posterior$lfsr, tolerance = 1e-10)
})

test_that("the estimated prior is non-negative, normalized and beats one scale on a two-scale truth", {
  set.seed(3)
  x <- c(rep(0, 300), rexp(200, 1), rexp(100, 1 / 20)) + rnorm(600)
  fit <- ebnm_scale_exponential(x, 1)
  expect_s3_class(fit, "ebnm")
  expect_s3_class(fit$fitted_g, "gammamix")
  expect_true(all(fit$fitted_g$pi >= 0))
  expect_equal(sum(fit$fitted_g$pi), 1)
  expect_true(all(fit$posterior$mean >= 0))
  expect_gt(as.numeric(fit$log_likelihood),
            as.numeric(ebnm::ebnm_point_exponential(x, 1)$log_likelihood) + 10)
  # A finer grid nests the default one, so the likelihood cannot go down.
  finer <- ebnm_scale_exponential(x, 1, mult = 2^(1 / 4))
  expect_gte(as.numeric(finer$log_likelihood), as.numeric(fit$log_likelihood) - 1e-6)
})

test_that("flashier's KL term stays <= 0 when s >> scale", {
  set.seed(4)
  s <- rep(1e4, 50)
  x <- rnorm(50, 0, 1e4)
  g <- ebnm::gammamix(pi = c(0.5, 0.5), shape = c(1, 1), scale = c(0, 1e-4), shift = c(0, 0))
  fit <- ebnm_scale_exponential(x, s, g_init = g, fix_g = TRUE,
                                output = ebnm::ebnm_output_all())
  nm_ll <- sum(-0.5 * log(2 * pi * s^2) -
                 0.5 * (x^2 - 2 * x * fit$posterior$mean + fit$posterior$second_moment) / s^2)
  kl <- as.numeric(fit$log_likelihood) - nm_ll
  expect_lte(kl, 1e-8)
  expect_gt(kl, -1e-6)
})

test_that("mode shifts the support and the posterior sampler matches the posterior", {
  set.seed(5)
  x <- 2 + c(rep(0, 100), rexp(100, 1 / 2)) + rnorm(200, 0, 0.5)
  fit <- ebnm_scale_exponential(x, 0.5, mode = 2, output = ebnm::ebnm_output_all())
  expect_true(all(fit$posterior$mean >= 2))
  draws <- fit$posterior_sampler(5000)
  expect_true(all(draws >= 2))
  z <- (colMeans(draws) - fit$posterior$mean) / (fit$posterior$sd / sqrt(5000))
  expect_lt(max(abs(z)), 6)
})

test_that("works as a flashier ebnm_fn", {
  skip_if_not_installed("flashier")
  set.seed(6)
  L <- matrix(rexp(60 * 2) * rbinom(60 * 2, 1, 0.5), 60, 2)
  F <- matrix(rexp(80 * 2) * rbinom(80 * 2, 1, 0.5), 80, 2)
  Y <- L %*% t(F) + matrix(rnorm(60 * 80, 0, 0.5), 60, 80)
  fl <- flashier::flash(Y, ebnm_fn = ebnm_scale_exponential, greedy_Kmax = 3,
                        backfit = TRUE, verbose = 0)
  expect_gte(fl$n_factors, 1)
  expect_true(all(fl$L_pm >= 0))
  expect_true(all(fl$F_pm >= 0))
  expect_true(all(unlist(fl$flash_fit$KL) <= 1e-6))
})
