# Exact Gaussian Matern backend: range bound and optional restart (EBSmoothr 0.3.6).

trap_data <- function() {
  set.seed(1)
  t <- seq(0, 10, length.out = 100)
  s <- rep(0.08, length(t))
  # A level far from 0 with a small dip, fitted with the mean fixed at 0: the likelihood
  # is flat along a ridge towards large ranges.
  x <- 5 - 0.1 * exp(-(t - 2)^2) + rnorm(length(t), sd = s)
  list(t = t, s = s, x = x, g_trap = Matern(theta = log(2e4), sigma = 0.6, beta = 0))
}
out <- c("fitted_g", "log_likelihood")

test_that("a warm start far beyond the domain no longer traps the exact optimizer", {
  d <- trap_data()
  m <- ebnm_Matern_generator(setup = Matern_setup(d$t, alpha = 2))
  warm <- m(d$x, d$s, g_init = d$g_trap, beta_fixed = 0, output = out)
  cold <- m(d$x, d$s, beta_fixed = 0, output = out)
  expect_equal(as.numeric(warm$log_likelihood), as.numeric(cold$log_likelihood), tolerance = 1e-6)
  expect_lte(exp(warm$fitted_g$theta), 20 * diff(range(d$t)) * (1 + 1e-8))

  old <- ebnm_Matern_generator(setup = Matern_setup(d$t, alpha = 2), range_max = Inf, restart = FALSE)
  warm_old <- old(d$x, d$s, g_init = d$g_trap, beta_fixed = 0, output = out)
  expect_gt(as.numeric(warm$log_likelihood), as.numeric(warm_old$log_likelihood))
})

test_that("restart keeps the better of the warm and default starts", {
  d <- trap_data()
  m1 <- ebnm_Matern_generator(setup = Matern_setup(d$t, alpha = 2), range_max = Inf, restart = TRUE)
  m0 <- ebnm_Matern_generator(setup = Matern_setup(d$t, alpha = 2), range_max = Inf, restart = FALSE)
  r1 <- m1(d$x, d$s, g_init = d$g_trap, beta_fixed = 0, output = out)
  r0 <- m0(d$x, d$s, g_init = d$g_trap, beta_fixed = 0, output = out)
  expect_gte(as.numeric(r1$log_likelihood), as.numeric(r0$log_likelihood) - 1e-8)
  expect_true(r1$exact_optimization$restart_used || isTRUE(all.equal(r1$fitted_g, r0$fitted_g)))
})

test_that("the range is never fitted above range_max", {
  d <- trap_data()
  for (rmax in c(15, 50)) {
    m <- ebnm_Matern_generator(setup = Matern_setup(d$t, alpha = 2), range_max = rmax)
    r <- m(d$x, d$s, beta_fixed = 0, output = out)
    expect_lte(exp(r$fitted_g$theta), rmax * (1 + 1e-8))
  }
})

test_that("an interior optimum is unchanged by the bound and by restart without g_init", {
  set.seed(2)
  t <- seq(0, 10, length.out = 80)
  s <- rep(0.2, length(t))
  x <- sin(2 * pi * t / 5) + rnorm(length(t), sd = s)
  new <- ebnm_Matern_generator(setup = Matern_setup(t, alpha = 2))
  old <- ebnm_Matern_generator(setup = Matern_setup(t, alpha = 2), range_max = Inf, restart = FALSE)
  a <- new(x, s, output = out)
  b <- old(x, s, output = out)
  expect_lt(exp(a$fitted_g$theta), 20 * diff(range(t)))
  expect_equal(as.numeric(a$log_likelihood), as.numeric(b$log_likelihood), tolerance = 1e-5)
  expect_equal(a$fitted_g$theta, b$fitted_g$theta, tolerance = 1e-3)
  expect_false(a$exact_optimization$restart_used)
})

test_that("range_max and restart are validated", {
  t <- seq(0, 1, length.out = 10)
  expect_error(ebnm_Matern_generator(setup = Matern_setup(t, alpha = 2), range_max = -1), "range_max")
  expect_error(ebnm_Matern_generator(setup = Matern_setup(t, alpha = 2), restart = NA), "restart")
  expect_silent(ebnm_Matern_generator(setup = Matern_setup(t, alpha = 2), range_max = Inf))
})
