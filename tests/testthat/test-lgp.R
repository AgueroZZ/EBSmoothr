test_that("LGP generator follows the new public beta semantics", {
  set.seed(1)

  t <- seq(0, 1, length.out = 40)
  s <- rep(0.1, length(t))
  x <- 0.5 + sin(2 * pi * t) + rnorm(length(t), sd = s)
  setup <- LGP_setup(t = t, betaprec = 0)
  fit_fun <- ebnm_LGP_generator(setup)

  fit_eb <- fit_fun(x, s, beta_prec = NULL)
  expect_equal(fit_eb$beta_mode, "empirical_bayes")
  expect_null(fit_eb$fitted_g$beta_prec)

  fit_flat <- fit_fun(x, s, beta_prec = 0)
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_equal(fit_flat$fitted_g$beta_prec, 0, tolerance = 1e-10)

  fit_proper <- fit_fun(x, s, beta_prec = 2)
  expect_equal(fit_proper$beta_mode, "prior_proper")
  expect_equal(fit_proper$fitted_g$beta_prec, 2, tolerance = 1e-10)
})

test_that("LGP fixed beta can be used even when fix_g is FALSE", {
  set.seed(2)

  t <- seq(0, 1, length.out = 50)
  s <- rep(0.1, length(t))
  x <- 0.4 + sin(2 * pi * t) + rnorm(length(t), sd = s)
  setup <- LGP_setup(t = t, betaprec = 0)
  fit_fun <- ebnm_LGP_generator(setup)

  beta_fixed <- c(0.4, 0)
  fit <- fit_fun(x, s, beta_fixed = beta_fixed)

  expect_equal(fit$beta_mode, "fixed")
  expect_equal(unname(fit$fitted_beta), beta_fixed, tolerance = 1e-8)
  expect_equal(unname(fit$fitted_g$beta), beta_fixed, tolerance = 1e-8)
})

test_that("LGP g_init beta is accepted and stored in the resolved fit state", {
  set.seed(3)

  t <- seq(0, 1, length.out = 45)
  s <- rep(0.1, length(t))
  x <- 0.2 + sin(2 * pi * t) + rnorm(length(t), sd = s)
  setup <- LGP_setup(t = t, betaprec = 0)
  fit_fun <- ebnm_LGP_generator(setup)

  g_init <- LGP(scale = 0, beta = c(0.1, 0), beta_prec = 0)
  fit <- fit_fun(x, s, g_init = g_init, beta_prec = 0)

  expect_s3_class(fit$g_init, "LGP")
  expect_equal(unname(fit$g_init$beta), c(0.1, 0), tolerance = 1e-10)
  expect_equal(unname(fit$g_init$beta_prec), 0, tolerance = 1e-10)
})

test_that("LGP legacy LGP_setup$betaprec remains a fallback when beta_prec is missing", {
  set.seed(4)

  t <- seq(0, 1, length.out = 35)
  s <- rep(0.1, length(t))
  x <- sin(2 * pi * t) + rnorm(length(t), sd = s)

  fit_legacy_eb <- ebnm_LGP_generator(LGP_setup(t = t, betaprec = -1))(x, s)
  expect_equal(fit_legacy_eb$beta_mode, "empirical_bayes")

  fit_legacy_flat <- ebnm_LGP_generator(LGP_setup(t = t, betaprec = 0))(x, s)
  expect_equal(fit_legacy_flat$beta_mode, "prior_flat")

  fit_override <- ebnm_LGP_generator(LGP_setup(t = t, betaprec = -1))(x, s, beta_prec = 2)
  expect_equal(fit_override$beta_mode, "prior_proper")
  expect_equal(fit_override$fitted_g$beta_prec, 2, tolerance = 1e-10)
})

test_that("LGP wrapper supports known and learned noise with the new beta modes", {
  set.seed(5)

  t <- seq(0, 1, length.out = 40)
  x <- 0.3 + sin(2 * pi * t) + rnorm(length(t), sd = 0.1)
  setup <- LGP_setup(t = t, betaprec = 0)

  fit_known <- eb_smoother(x, s = 0.1, family = "lgp", setup = setup, beta_prec = 0)
  expect_equal(fit_known$beta_mode, "prior_flat")
  expect_equal(unname(fit_known$fitted_g$beta_prec), 0, tolerance = 1e-10)

  fit_learned <- eb_smoother(x, s = NULL, family = "lgp", setup = setup, beta_prec = 2)
  expect_equal(fit_learned$beta_mode, "prior_proper")
  expect_equal(unname(fit_learned$fitted_g$beta_prec), 2, tolerance = 1e-10)
  expect_true(is.finite(fit_learned$fitted_noise_sd))
  expect_gt(fit_learned$fitted_noise_sd, 0)
})

test_that("LGP log-link auto backend uses Fisher Laplace", {
  set.seed(6)

  t <- seq(0, 1, length.out = 16)
  s <- rep(0.08, length(t))
  x <- exp(0.1 + 0.15 * sin(2 * pi * t)) + rnorm(length(t), sd = s)
  setup <- LGP_setup(t = t, num_knots = 8, betaprec = 0, link = "log")

  fit_auto <- ebnm_LGP_generator(setup, link = "log")(x, s)
  fit_observed <- ebnm_LGP_generator(setup, link = "log", backend = "laplace")(x, s)
  fit_wrap <- eb_smoother(x, s = s, family = "lgp", setup = setup, link = "log")
  fit_learned <- eb_smoother(x, s = NULL, family = "lgp", setup = setup, link = "log")

  expect_equal(fit_auto$backend, "laplace_fisher")
  expect_equal(fit_auto$laplace_curvature, "fisher")
  expect_equal(fit_auto$log_likelihood_semantics, "laplace_fisher_prior_flat")
  expect_equal(fit_observed$backend, "laplace")
  expect_equal(fit_observed$laplace_curvature, "observed")
  expect_lt(max(abs(fit_auto$posterior$mean - fit_observed$posterior$mean)), 0.01)
  expect_gt(abs(as.numeric(fit_auto$log_likelihood) - as.numeric(fit_observed$log_likelihood)), 1e-5)
  expect_equal(fit_wrap$backend, "laplace_fisher")
  expect_equal(fit_wrap$raw_fit$laplace_curvature, "fisher")
  expect_equal(fit_learned$backend, "laplace_fisher")
  expect_true(is.finite(fit_learned$fitted_noise_sd))
  expect_gt(fit_learned$fitted_noise_sd, 0)
})

test_that("LGP softplus link defaults to observed Laplace and allows Fisher", {
  skip_if_not_installed("TMB")
  set.seed(9)
  n <- 40
  t <- seq(0, 1, length.out = n)
  eta <- seq(-6, 6, length.out = n)
  x <- log1p(exp(eta)) + rnorm(n, sd = 0.08)
  s <- rep(0.08, n)
  setup <- LGP_setup(t = t, num_knots = 8, betaprec = 0, link = "softplus")

  fit_auto <- ebnm_LGP_generator(setup, link = "softplus")(x, s)
  fit_fisher <- ebnm_LGP_generator(setup, link = "softplus", backend = "laplace_fisher")(x, s)

  expect_equal(fit_auto$backend, "laplace")
  expect_equal(fit_auto$laplace_curvature, "observed")
  expect_equal(fit_fisher$backend, "laplace_fisher")
  expect_equal(fit_fisher$laplace_curvature, "fisher")
})

test_that("LGP generator accepts softplus as a public link", {
  t <- seq(0, 1, length.out = 8)
  setup <- LGP_setup(t = t, num_knots = 6, betaprec = 0, link = "softplus")

  fit_fun <- ebnm_LGP_generator(setup, link = "softplus", backend = "laplace")

  expect_type(fit_fun, "closure")
})

test_that("LGP softplus posterior moments match empirical moments from posterior draws", {
  skip_if_not_installed("TMB")
  set.seed(10)
  n <- 36
  t <- seq(0, 1, length.out = n)
  s <- rep(0.07, n)
  eta <- 0.3 + 1.4 * sin(2 * pi * t)
  x <- log1p(exp(eta)) + rnorm(n, sd = s)
  fit <- ebnm_LGP_generator(
    LGP_setup(t = t, num_knots = 10, betaprec = 0, link = "softplus"),
    link = "softplus",
    backend = "laplace"
  )(x, s)

  nsamp <- 12000
  draws <- fit$posterior_sampler(nsamp)
  draw_mean <- colMeans(draws)
  draw_var <- apply(draws, 2, var)

  mean_tol <- 5 * max(apply(draws, 2, stats::sd) / sqrt(nsamp)) + 1e-4
  var_tol <- 8 * max(sqrt(2 / (nsamp - 1)) * pmax(draw_var, fit$posterior$var)) + 1e-4
  expect_lt(max(abs(draw_mean - fit$posterior$mean)), mean_tol)
  expect_lt(max(abs(draw_var - fit$posterior$var)), var_tol)
})

test_that("LGP fix_params supports scale and beta fixing", {
  set.seed(7)

  t <- seq(0, 1, length.out = 18)
  s <- rep(0.08, length(t))
  x <- exp(0.1 + 0.12 * sin(2 * pi * t)) + rnorm(length(t), sd = s)
  setup <- LGP_setup(t = t, num_knots = 8, betaprec = -1, link = "log")
  g_init <- LGP(scale = 0.25, beta = c(0.1, 0), beta_prec = NULL)

  fit_scale <- ebnm_LGP_generator(setup, link = "log")(x, s, g_init = g_init, fix_params = "scale")
  fit_fix_g <- ebnm_LGP_generator(setup, link = "log")(x, s, g_init = g_init, fix_g = TRUE)
  fit_beta_arg <- ebnm_LGP_generator(setup, link = "log")(x, s, g_init = g_init, beta_fixed = c(0.2, 0.1), fix_params = "beta")
  fit_beta_g <- ebnm_LGP_generator(setup, link = "log")(x, s, g_init = g_init, fix_params = "beta")

  expect_equal(fit_scale$fitted_g$scale, 0.25, tolerance = 1e-12)
  expect_equal(fit_fix_g$fitted_g$scale, 0.25, tolerance = 1e-12)
  expect_equal(unname(fit_beta_arg$fitted_beta), c(0.2, 0.1), tolerance = 1e-12)
  expect_equal(unname(fit_beta_g$fitted_beta), c(0.1, 0), tolerance = 1e-12)

  expect_error(
    ebnm_LGP_generator(setup, link = "log")(x, s, fix_params = "scale"),
    "requires `g_init\\$scale`"
  )
  expect_error(
    ebnm_LGP_generator(setup, link = "log")(x, s, fix_params = "sigma"),
    "unsupported value"
  )
  expect_error(
    ebnm_LGP_generator(setup, link = "log")(x, s, g_init = g_init, fix_params = "beta", beta_prec = 0),
    "cannot be combined with `beta_prec`"
  )
  expect_error(
    ebnm_LGP_generator(LGP_setup(t = t, num_knots = 8, link = "identity"), link = "identity", backend = "laplace_fisher"),
    "only available for `link = \"log\"`"
  )
})
