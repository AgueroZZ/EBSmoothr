test_that("Matern generator supports log-link positive smoothing", {
  set.seed(11)

  loc <- seq(0, 1, length.out = 8)
  s <- rep(0.08, length(loc))
  x <- exp(0.2 + 0.25 * sin(2 * pi * loc)) + rnorm(length(loc), sd = s)

  fit <- ebnm_Matern_generator(locations = loc, link = "log")(x, s)

  expect_equal(fit$backend, "fisher_pql")
  expect_equal(fit$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit$laplace_curvature, "fisher")
  expect_equal(fit$log_likelihood_semantics, "fisher_laplace_at_fisher_pql_mode_empirical_bayes")
  expect_equal(fit$fisher_pql_diagnostics$inner_iterations, 3L)
  expect_equal(fit$link, "log")
  expect_true(all(fit$posterior$mean > 0))
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
})

test_that("Matern non-identity initialization keeps latent and raw scales separate", {
  x <- c(0.5, 1, 2, 4)
  s <- c(1, 2, 1, 2)
  w <- 1 / s^2
  x_floor <- pmax(x, min(x[x > 0]) / 2)

  init_log <- EBSmoothr:::.resolve_matern_g_init(
    x = x,
    s = s,
    penalty_range0 = 1,
    link = "log"
  )
  beta_log <- log(stats::weighted.mean(x, w = w))
  eta_log <- log(x_floor)
  w_log <- x_floor^2 / s^2
  sigma_log <- sqrt(sum(w_log * (eta_log - beta_log)^2) / sum(w_log))
  expect_equal(init_log$beta_init, beta_log, tolerance = 1e-12)
  expect_equal(init_log$sigma_data, sigma_log, tolerance = 1e-12)

  init_softplus <- EBSmoothr:::.resolve_matern_g_init(
    x = x,
    s = s,
    penalty_range0 = 1,
    link = "softplus"
  )
  beta_softplus <- .inverse_softplus_stable(stats::weighted.mean(x, w = w))
  eta_softplus <- .inverse_softplus_stable(x_floor)
  deriv_softplus <- .sigmoid_stable(eta_softplus)
  w_softplus <- deriv_softplus^2 / s^2
  sigma_softplus <- sqrt(sum(w_softplus * (eta_softplus - beta_softplus)^2) / sum(w_softplus))
  expect_equal(init_softplus$beta_init, beta_softplus, tolerance = 1e-12)
  expect_equal(init_softplus$sigma_data, sigma_softplus, tolerance = 1e-12)

  x_learn <- c(0.2, 0.9, 2.5, 7.5, 9)
  raw_scale <- stats::sd(x_learn - mean(x_learn))
  init_log_learn <- EBSmoothr:::.resolve_matern_g_init(
    x = x_learn,
    s = NULL,
    penalty_range0 = 1,
    allow_noise = TRUE,
    link = "log"
  )
  init_softplus_learn <- EBSmoothr:::.resolve_matern_g_init(
    x = x_learn,
    s = NULL,
    penalty_range0 = 1,
    allow_noise = TRUE,
    link = "softplus"
  )
  expect_equal(init_log_learn$beta_init, log(mean(x_learn)), tolerance = 1e-12)
  expect_equal(init_log_learn$noise_sd_init, raw_scale, tolerance = 1e-12)
  expect_gt(abs(init_log_learn$noise_sd_init - init_log_learn$sigma_data), 0.1)
  expect_equal(init_softplus_learn$beta_init, .inverse_softplus_stable(mean(x_learn)), tolerance = 1e-12)
  expect_equal(init_softplus_learn$noise_sd_init, raw_scale, tolerance = 1e-12)
  expect_gt(abs(init_softplus_learn$noise_sd_init - init_softplus_learn$sigma_data), 0.1)
})

test_that("Matern auto backend uses Fisher-PQL for log-link fits", {
  set.seed(111)

  loc <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 4),
    y = seq(0, 1, length.out = 4)
  ))
  s <- rep(0.1, nrow(loc))
  eta <- 0.1 + 0.2 * sin(2 * pi * loc[, 1]) + 0.1 * cos(2 * pi * loc[, 2])
  x <- exp(eta) + rnorm(nrow(loc), sd = s)
  setup <- Matern_setup(locations = loc, max.edge = 0.6)

  fit_auto_fixed <- ebnm_Matern_generator(setup = setup, link = "log")(x, s, beta_fixed = 0)
  fit_laplace_fixed <- ebnm_Matern_generator(
    setup = setup,
    link = "log",
    backend = "laplace_tmb"
  )(x, s, beta_fixed = 0)
  fit_auto_eb <- ebnm_Matern_generator(setup = setup, link = "log")(x, s)
  fit_auto_fixed_g <- ebnm_Matern_generator(setup = setup, link = "log")(
    x,
    s,
    g_init = Matern(theta = log(0.4), sigma = 0.25),
    fix_g = TRUE,
    beta_fixed = 0
  )

  expect_equal(fit_auto_fixed$backend, "fisher_pql")
  expect_equal(fit_auto_fixed$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit_auto_fixed$laplace_curvature, "fisher")
  expect_equal(fit_auto_fixed$log_likelihood_semantics, "fisher_laplace_at_fisher_pql_mode_fixed")
  expect_equal(fit_auto_fixed$fisher_pql_diagnostics$inner_iterations, 3L)
  expect_equal(fit_laplace_fixed$backend, "laplace")
  expect_equal(fit_laplace_fixed$laplace_implementation, "tmb")
  expect_equal(fit_laplace_fixed$laplace_curvature, "observed")
  expect_true(is.finite(as.numeric(fit_auto_fixed$log_likelihood)))
  expect_true(is.finite(as.numeric(fit_laplace_fixed$log_likelihood)))
  expect_lt(max(abs(fit_auto_fixed$posterior$mean - fit_laplace_fixed$posterior$mean)), 0.1)

  expect_equal(fit_auto_eb$backend, "fisher_pql")
  expect_equal(fit_auto_eb$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit_auto_eb$laplace_curvature, "fisher")
  expect_equal(fit_auto_fixed_g$backend, "fisher_pql")
  expect_equal(fit_auto_fixed_g$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit_auto_fixed_g$laplace_curvature, "fisher")
  expect_equal(fit_auto_fixed_g$fisher_pql_diagnostics$inner_iterations, 3L)
})

test_that("Matern softplus auto backend uses Fisher-PQL; INLA unsupported", {
  set.seed(211)
  loc <- seq(0, 1, length.out = 12)
  s <- rep(0.1, length(loc))
  eta <- seq(-6, 6, length.out = length(loc))
  x <- log1p(exp(eta)) + rnorm(length(loc), sd = s)

  fit_auto <- ebnm_Matern_generator(locations = loc, link = "softplus")(x, s)
  fit_fisher <- ebnm_Matern_generator(locations = loc, link = "softplus", backend = "laplace_fisher")(x, s)
  expect_equal(fit_auto$backend, "fisher_pql")
  expect_equal(fit_auto$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit_auto$laplace_curvature, "fisher")
  expect_equal(fit_auto$fisher_pql_diagnostics$inner_iterations, 3L)
  expect_equal(fit_fisher$backend, "laplace_fisher")
  expect_equal(fit_fisher$laplace_curvature, "fisher")
  expect_error(
    ebnm_Matern_generator(locations = loc, link = "softplus", backend = "inla")(x, s),
    "not currently supported by the INLA Matern backend"
  )
})

test_that("Matern softplus posterior moments match empirical moments from posterior draws", {
  set.seed(212)
  loc <- seq(0, 1, length.out = 16)
  s <- rep(0.08, length(loc))
  eta <- -0.4 + 1.2 * sin(2 * pi * loc)
  x <- log1p(exp(eta)) + rnorm(length(loc), sd = s)

  fit <- ebnm_Matern_generator(locations = loc, link = "softplus", backend = "laplace")(x, s)
  nsamp <- 12000
  draws <- fit$posterior_sampler(nsamp)
  draw_mean <- colMeans(draws)
  draw_var <- apply(draws, 2, var)

  mean_tol <- 5 * max(apply(draws, 2, stats::sd) / sqrt(nsamp)) + 1e-4
  var_tol <- 8 * max(sqrt(2 / (nsamp - 1)) * pmax(draw_var, fit$posterior$var)) + 1e-4
  expect_lt(max(abs(draw_mean - fit$posterior$mean)), mean_tol)
  expect_lt(max(abs(draw_var - fit$posterior$var)), var_tol)
})

test_that("Matern log-link manual and INLA backends agree", {
  set.seed(12)

  loc <- seq(0, 1, length.out = 8)
  s <- rep(0.1, length(loc))
  x <- exp(0.15 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = s)

  fit_laplace_eb <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "laplace"
  )(x, s)
  fit_inla_eb <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "inla"
  )(x, s)

  expect_equal(fit_inla_eb$beta_mode, "empirical_bayes")
  expect_false(is.null(fit_inla_eb$beta_profile_optimization))
  expect_equal(fit_laplace_eb$laplace_implementation, "tmb")
  expect_equal(as.numeric(fit_laplace_eb$log_likelihood), as.numeric(fit_inla_eb$log_likelihood), tolerance = 2e-3)
  expect_equal(fit_laplace_eb$fitted_g$theta, fit_inla_eb$fitted_g$theta, tolerance = 3e-2)
  expect_equal(fit_laplace_eb$fitted_g$sigma, fit_inla_eb$fitted_g$sigma, tolerance = 3e-2)
  expect_lt(max(abs(fit_laplace_eb$posterior$mean - fit_inla_eb$posterior$mean)), 0.04)

  fit_laplace <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "laplace"
  )(x, s, beta_fixed = 0.15)
  fit_inla <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "inla"
  )(x, s, beta_fixed = 0.15)

  expect_equal(fit_laplace$laplace_implementation, "tmb")
  expect_equal(as.numeric(fit_laplace$log_likelihood), as.numeric(fit_inla$log_likelihood), tolerance = 1e-3)
  expect_equal(fit_laplace$fitted_g$theta, fit_inla$fitted_g$theta, tolerance = 1e-2)
  expect_equal(fit_laplace$fitted_g$sigma, fit_inla$fitted_g$sigma, tolerance = 1e-2)
  expect_lt(max(abs(fit_laplace$posterior$mean - fit_inla$posterior$mean)), 0.03)

  fit_laplace_flat <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "laplace"
  )(x, s, beta_prec = 0)
  fit_inla_flat <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "inla"
  )(x, s, beta_prec = 0)

  expect_equal(fit_laplace_flat$laplace_implementation, "tmb")
  expect_equal(as.numeric(fit_laplace_flat$log_likelihood), as.numeric(fit_inla_flat$log_likelihood), tolerance = 1e-3)
  expect_equal(fit_laplace_flat$fitted_g$theta, fit_inla_flat$fitted_g$theta, tolerance = 1e-2)
  expect_equal(fit_laplace_flat$fitted_g$sigma, fit_inla_flat$fitted_g$sigma, tolerance = 1e-2)
  expect_lt(max(abs(fit_laplace_flat$posterior$mean - fit_inla_flat$posterior$mean)), 0.04)

  pc <- list(range = c(0.25, 0.5), sigma = c(0.25, 0.5))
  fit_laplace_pc <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "laplace",
    pc.penalty = pc
  )(x, s, beta_fixed = 0.15)
  fit_inla_pc <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "inla",
    pc.penalty = pc
  )(x, s, beta_fixed = 0.15)

  expect_equal(fit_laplace_pc$laplace_implementation, "tmb")
  expect_equal(as.numeric(fit_laplace_pc$log_likelihood), as.numeric(fit_inla_pc$log_likelihood), tolerance = 1e-3)
  expect_lt(max(abs(fit_laplace_pc$posterior$mean - fit_inla_pc$posterior$mean)), 0.04)
})

test_that("Matern log-link Laplace backend supports beta modes", {
  set.seed(13)

  loc <- seq(0, 1, length.out = 7)
  s <- rep(0.1, length(loc))
  x <- exp(0.1 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = s)
  fit_fun <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace")

  fit_eb <- fit_fun(x, s)
  fit_fixed <- fit_fun(x, s, beta_fixed = 0.1)
  fit_flat <- fit_fun(x, s, beta_prec = 0)
  fit_proper <- fit_fun(x, s, beta_prec = 2)
  fit_tmb <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_tmb")(x, s)
  fit_r <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_r")(x, s)
  fit_tmb_flat <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_tmb")(x, s, beta_prec = 0)
  fit_r_flat <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_r")(x, s, beta_prec = 0)

  expect_equal(fit_eb$beta_mode, "empirical_bayes")
  expect_equal(fit_fixed$beta_mode, "fixed")
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_equal(fit_proper$beta_mode, "prior_proper")
  expect_equal(fit_eb$laplace_implementation, "tmb")
  expect_equal(fit_fixed$laplace_implementation, "tmb")
  expect_equal(fit_flat$laplace_implementation, "tmb")
  expect_equal(fit_proper$laplace_implementation, "tmb")
  expect_equal(fit_fixed$laplace_diagnostics$tmb_mode_source, "last.par.best")
  expect_equal(as.numeric(fit_tmb$log_likelihood), as.numeric(fit_r$log_likelihood), tolerance = 1e-5)
  expect_equal(as.numeric(fit_tmb_flat$log_likelihood), as.numeric(fit_r_flat$log_likelihood), tolerance = 1e-5)
  expect_true(all(fit_eb$posterior$mean > 0))
  expect_true(all(fit_fixed$posterior$mean > 0))
  expect_true(all(fit_flat$posterior$mean > 0))
  expect_true(all(fit_proper$posterior$mean > 0))
})

test_that("Matern log-link Fisher backend supports beta modes and learned noise", {
  set.seed(131)

  loc <- seq(0, 1, length.out = 7)
  s <- rep(0.1, length(loc))
  x <- exp(0.1 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = s)
  fit_fun <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_fisher")

  fit_eb <- fit_fun(x, s)
  fit_fixed <- fit_fun(x, s, beta_fixed = 0.1)
  fit_flat <- fit_fun(x, s, beta_prec = 0)
  fit_proper <- fit_fun(x, s, beta_prec = 2)
  learned_eb <- eb_smoother(x, s = NULL, family = "matern", locations = loc, link = "log")
  learned_fixed <- eb_smoother(x, s = NULL, family = "matern", locations = loc, link = "log", beta_fixed = 0.1)
  learned_flat <- eb_smoother(x, s = NULL, family = "matern", locations = loc, link = "log", beta_prec = 0)
  learned_proper <- eb_smoother(x, s = NULL, family = "matern", locations = loc, link = "log", beta_prec = 2)

  for (fit in list(fit_eb, fit_fixed, fit_flat, fit_proper)) {
    expect_equal(fit$backend, "laplace_fisher")
    expect_equal(fit$laplace_curvature, "fisher")
    expect_match(fit$log_likelihood_semantics, "^laplace_fisher_")
    expect_true(all(is.finite(fit$posterior$mean)))
    expect_true(all(is.finite(fit$posterior$var)))
  }

  # eb_smoother() without an explicit backend follows the auto policy, which
  # routes non-identity links to fisher_pql since 0.2.4.
  for (fit in list(learned_eb, learned_fixed, learned_flat, learned_proper)) {
    expect_equal(fit$backend, "fisher_pql")
    expect_equal(fit$laplace_curvature, "fisher")
    expect_match(fit$log_likelihood_semantics, "^fisher_laplace_at_fisher_pql_mode_")
    expect_true(all(is.finite(fit$posterior$mean)))
    expect_true(all(is.finite(fit$posterior$var)))
  }

  expect_equal(fit_eb$beta_mode, "empirical_bayes")
  expect_equal(fit_fixed$beta_mode, "fixed")
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_equal(fit_proper$beta_mode, "prior_proper")
  expect_equal(learned_eb$beta_mode, "empirical_bayes")
  expect_equal(learned_fixed$beta_mode, "fixed")
  expect_equal(learned_flat$beta_mode, "prior_flat")
  expect_equal(learned_proper$beta_mode, "prior_proper")
  expect_true(is.finite(learned_eb$fitted_noise_sd))
  expect_gt(learned_eb$fitted_noise_sd, 0)
})

test_that("Matern observed and Fisher Laplace are close but distinct on ordinary log-link data", {
  set.seed(132)

  loc <- seq(0, 1, length.out = 7)
  s <- rep(0.1, length(loc))
  x <- exp(0.1 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = s)

  fit_fisher <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_fisher")(x, s)
  fit_observed <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_r")(x, s)

  expect_equal(fit_fisher$backend, "laplace_fisher")
  expect_equal(fit_fisher$laplace_curvature, "fisher")
  expect_equal(fit_observed$backend, "laplace")
  expect_equal(fit_observed$laplace_curvature, "observed")
  expect_lt(max(abs(fit_fisher$posterior$mean - fit_observed$posterior$mean)), 0.01)
  expect_gt(abs(as.numeric(fit_fisher$log_likelihood) - as.numeric(fit_observed$log_likelihood)), 1e-5)
})

test_that("Matern Fisher curvature stays positive on deterministic pseudo-pathology", {
  eta <- rep(log(1), 4)
  x <- c(0.9, 1.1, 3.5, 4.0)
  s <- rep(0.1, length(x))

  observed <- EBSmoothr:::.matern_observation_terms(
    eta = eta,
    x = x,
    s = s,
    link = "log",
    laplace_curvature = "observed"
  )
  fisher <- EBSmoothr:::.matern_observation_terms(
    eta = eta,
    x = x,
    s = s,
    link = "log",
    laplace_curvature = "fisher"
  )

  expect_true(any(observed$hess_diag < 0))
  expect_true(all(fisher$hess_diag > 0))

  loc <- seq(0, 1, length.out = 9)
  x_fit <- c(1, 1.05, 8, 1.1, 0.95, 7, 1, 1.05, 0.9)
  fit <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_fisher")(x_fit, rep(0.08, length(loc)))
  expect_equal(fit$backend, "laplace_fisher")
  expect_equal(fit$laplace_curvature, "fisher")
  expect_true(all(is.finite(fit$posterior$mean)))
  expect_true(all(is.finite(fit$posterior$var)))
})

test_that("Matern Fisher-PQL pseudo-response matches link algebra", {
  eta <- c(-0.4, 0.2, 0.8)
  x <- c(0.7, 1.2, 2.4)
  s <- c(0.08, 0.09, 0.1)

  log_pseudo <- EBSmoothr:::.matern_fisher_pql_pseudo_response(
    x = x,
    eta = eta,
    s = s,
    link = "log"
  )
  log_mu <- exp(eta)
  log_g <- exp(eta)
  expect_equal(log_pseudo$z, eta + (x - log_mu) / log_g, tolerance = 1e-12)
  expect_equal(log_pseudo$s, s / log_g, tolerance = 1e-12)
  expect_equal(log_pseudo$noise_scale, 1 / log_g, tolerance = 1e-12)

  softplus_pseudo <- EBSmoothr:::.matern_fisher_pql_pseudo_response(
    x = x,
    eta = eta,
    s = s,
    link = "softplus"
  )
  softplus_mu <- log1p(exp(eta))
  softplus_g <- 1 / (1 + exp(-eta))
  expect_equal(softplus_pseudo$z, eta + (x - softplus_mu) / softplus_g, tolerance = 1e-12)
  expect_equal(softplus_pseudo$s, s / softplus_g, tolerance = 1e-12)
  expect_equal(softplus_pseudo$noise_scale, 1 / softplus_g, tolerance = 1e-12)
})

test_that("Matern Fisher-PQL rejects unsupported links", {
  loc <- seq(0, 1, length.out = 8)
  s <- rep(0.1, length(loc))
  x <- 0.2 + sin(2 * pi * loc)

  expect_error(
    ebnm_Matern_generator(locations = loc, link = "identity", backend = "fisher_pql")(x, s),
    "only available for `link = \"log\"` or `link = \"softplus\"`"
  )
  expect_error(
    eb_smoother(
      x,
      s = s,
      family = "matern",
      locations = loc,
      link = "identity",
      backend = "fisher_pql"
    ),
    "only available for `link = \"log\"` or `link = \"softplus\"`"
  )
})

test_that("Matern Fisher-PQL Step A mode matches the sparse one-step reference", {
  set.seed(219)
  loc <- seq(0, 1, length.out = 8)
  setup <- Matern_setup(locations = loc)
  s <- rep(0.1, length(loc))
  g_init <- Matern(theta = log(0.4), sigma = 0.3)

  for (link in c("log", "softplus")) {
    eta <- if (identical(link, "log")) {
      0.08 + 0.18 * sin(2 * pi * loc)
    } else {
      0.55 + 0.22 * sin(2 * pi * loc)
    }
    mean <- if (identical(link, "log")) exp(eta) else log1p(exp(eta))
    x <- mean + rnorm(length(loc), sd = 0.04)

    fit <- EBSmoothr:::.fit_matern_fisher_pql_known_noise(
      x = x,
      s = s,
      A = setup$A,
      spde_template = setup$spde_template,
      alpha = setup$alpha,
      d = setup$d,
      theta_init = log(0.4),
      sigma_init = 0.3,
      beta_init = 0,
      beta_mode = "fixed",
      beta_fixed = 0,
      link = link,
      fix_g = TRUE,
      pql_max_iter = 1L
    )

    ref <- EBSmoothr:::.matern_fisher_pql_reference_mode_at_params(
      x = x,
      s = s,
      A = setup$A,
      spde_template = setup$spde_template,
      alpha = setup$alpha,
      d = setup$d,
      log_range = log(0.4),
      log_sigma = log(0.3),
      beta_mode = "fixed",
      beta_value = 0,
      beta_start = 0,
      link = link,
      pql_w_start = rep(0, ncol(setup$A))
    )

    expect_equal(fit$laplace_implementation, "exact_fisher_pql")
    expect_equal(fit$fisher_pql_diagnostics$inner_iterations, 1L)
    expect_equal(fit$fisher_pql_diagnostics$stepA_engine, "exact_gaussian")
    expect_equal(fit$fisher_pql_mode, ref$mode, tolerance = 1e-8)
  }
})

test_that("Matern Fisher-PQL exposes pseudo-Gaussian update count", {
  set.seed(2191)
  loc <- seq(0, 1, length.out = 10)
  s <- rep(0.1, length(loc))
  eta <- 0.5 + 0.2 * sin(2 * pi * loc)
  x <- log1p(exp(eta)) + rnorm(length(loc), sd = 0.04)

  fit <- ebnm_Matern_generator(
    locations = loc,
    link = "softplus",
    backend = "fisher_pql"
  )(x, s)

  expect_equal(fit$fisher_pql_diagnostics$max_inner_iter, 3L)
  expect_gte(fit$fisher_pql_diagnostics$inner_iterations, 1L)
  expect_lte(fit$fisher_pql_diagnostics$inner_iterations, 3L)
  expect_equal(
    length(fit$fisher_pql_diagnostics$eta_change),
    fit$fisher_pql_diagnostics$inner_iterations
  )
  expect_equal(
    length(fit$fisher_pql_diagnostics$stepA_log_marginals),
    fit$fisher_pql_diagnostics$inner_iterations
  )

  fit_two <- ebnm_Matern_generator(
    locations = loc,
    link = "softplus",
    backend = "fisher_pql",
    pql_inner_iter = 2L
  )(x, s)

  expect_equal(fit_two$fisher_pql_diagnostics$max_inner_iter, 2L)
  expect_lte(fit_two$fisher_pql_diagnostics$inner_iterations, 2L)
  expect_equal(
    length(fit_two$fisher_pql_diagnostics$eta_change),
    fit_two$fisher_pql_diagnostics$inner_iterations
  )
  expect_equal(
    length(fit_two$fisher_pql_diagnostics$stepA_log_marginals),
    fit_two$fisher_pql_diagnostics$inner_iterations
  )

  fit_learned_two <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    link = "softplus",
    backend = "fisher_pql",
    pql_inner_iter = 2L
  )

  expect_equal(fit_learned_two$raw_fit$fisher_pql_diagnostics$max_inner_iter, 2L)
  expect_lte(fit_learned_two$raw_fit$fisher_pql_diagnostics$inner_iterations, 2L)
  expect_equal(
    length(fit_learned_two$raw_fit$fisher_pql_diagnostics$eta_change),
    fit_learned_two$raw_fit$fisher_pql_diagnostics$inner_iterations
  )
  expect_equal(
    length(fit_learned_two$raw_fit$fisher_pql_diagnostics$stepA_log_marginals),
    fit_learned_two$raw_fit$fisher_pql_diagnostics$inner_iterations
  )

  expect_error(
    ebnm_Matern_generator(
      locations = loc,
      link = "softplus",
      backend = "fisher_pql",
      pql_inner_iter = 0L
    ),
    "positive integer"
  )
  expect_error(
    eb_smoother(
      x,
      s = s,
      family = "matern",
      locations = loc,
      link = "softplus",
      backend = "fisher_pql",
      pql_inner_iter = 0L
    ),
    "positive integer"
  )
})

test_that("Matern Fisher-PQL routes Step A through exact Gaussian fits, not TMB", {
  expect_false(exists(".fit_matern_fisher_pql_tmb", envir = asNamespace("EBSmoothr"), inherits = FALSE))
  expect_false(exists(".matern_fisher_pql_tmb_data", envir = asNamespace("EBSmoothr"), inherits = FALSE))
  pql_body <- paste(
    deparse(body(EBSmoothr:::.fit_matern_fisher_pql_exact_stepA)),
    deparse(body(EBSmoothr:::.fit_matern_fisher_pql_known_noise)),
    deparse(body(EBSmoothr:::.fit_matern_fisher_pql_unknown_noise)),
    collapse = "\n"
  )
  expect_true(grepl(".fit_matern_exact_known_noise", pql_body, fixed = TRUE))
  expect_true(grepl(".fit_matern_exact_unknown_noise", pql_body, fixed = TRUE))
  expect_false(grepl(".matern_tmb_make_adfun", pql_body, fixed = TRUE))
  src_path <- testthat::test_path("../../src/EBSmoothr.cpp")
  if (file.exists(src_path)) {
    src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
    expect_false(grepl("model_id == 2", src, fixed = TRUE))
  }
})

test_that("Matern Fisher-PQL known-noise backend matches Laplace references", {
  cases <- list(
    list(link = "log", seed = 220, ref_backend = "laplace_fisher"),
    list(link = "softplus", seed = 221, ref_backend = "laplace")
  )

  for (case in cases) {
    set.seed(case$seed)
    loc <- seq(0, 1, length.out = 12)
    eta <- if (identical(case$link, "log")) {
      0.05 + 0.18 * sin(2 * pi * loc)
    } else {
      0.6 + 0.35 * sin(2 * pi * loc)
    }
    mean <- if (identical(case$link, "log")) exp(eta) else log1p(exp(eta))
    s <- rep(0.08, length(loc))
    x <- mean + rnorm(length(loc), sd = 0.04)

    fit <- ebnm_Matern_generator(
      locations = loc,
      link = case$link,
      backend = "fisher_pql"
    )(x, s)
    ref <- ebnm_Matern_generator(
      locations = loc,
      link = case$link,
      backend = case$ref_backend
    )(x, s)

    expect_equal(fit$backend, "fisher_pql")
    expect_equal(fit$laplace_implementation, "exact_fisher_pql")
    expect_equal(fit$laplace_curvature, "fisher")
    expect_equal(fit$log_likelihood_semantics, "fisher_laplace_at_fisher_pql_mode_empirical_bayes")
    expect_equal(as.numeric(fit$log_likelihood), fit$log_likelihood_original_at_fisher_pql_mode, tolerance = 1e-12)
    expect_false("log_likelihood_laplace_at_fisher_pql_params" %in% names(fit))
    expect_true(is.finite(fit$log_likelihood_fisher_pql_stepA))
    expect_true(isTRUE(fit$fisher_pql_diagnostics$converged))
    expect_equal(fit$fisher_pql_diagnostics$inner_iterations, 3L)
    expect_equal(fit$fisher_pql_diagnostics$stepA_engine, "exact_gaussian")
    expect_equal(fit$fisher_pql_diagnostics$g_floor, 1e-6)
    expect_true(all(is.finite(fit$fisher_pql_eta_mode)))
    expect_true(all(is.finite(fit$fitted_s)))
    expect_gt(min(fit$fitted_s), 0)
    expect_true(all(is.finite(fit$posterior$mean)))
    expect_true(all(is.finite(fit$posterior$var)))
    expect_lt(max(abs(fit$posterior$mean - ref$posterior$mean)), 0.1)
    expect_lt(sqrt(mean((fit$posterior$mean - ref$posterior$mean)^2)), 0.05)
    ref_w <- ref$posterior_spatial_field$mode
    if (is.null(ref_w)) ref_w <- ref$posterior_spatial_field$mean
    eta_ref <- as.numeric(ref$matern_objective_context$A %*% as.numeric(ref_w)) + ref$fitted_beta
    expect_lt(max(abs(fit$fisher_pql_eta_mode - eta_ref)), 0.15)
    expect_lt(sqrt(mean((fit$fisher_pql_eta_mode - eta_ref)^2)), 0.05)
    expect_lt(abs(fit$fitted_g$theta - ref$fitted_g$theta), 0.15)
    expect_lt(abs(log(fit$fitted_g$sigma) - log(ref$fitted_g$sigma)), 0.15)
    expect_true(is.finite(as.numeric(fit$log_likelihood) - as.numeric(ref$log_likelihood)))
  }
})

test_that("Matern Fisher-PQL learned-noise backend matches Laplace references", {
  cases <- list(
    list(link = "log", seed = 302, ref_backend = "laplace_fisher"),
    list(link = "softplus", seed = 306, ref_backend = "laplace")
  )

  for (case in cases) {
    set.seed(case$seed)
    loc <- seq(0, 1, length.out = 20)
    eta <- if (identical(case$link, "log")) {
      0.05 + 0.25 * sin(2 * pi * loc) + 0.1 * cos(4 * pi * loc)
    } else {
      0.8 + 0.45 * sin(2 * pi * loc) + 0.1 * cos(4 * pi * loc)
    }
    mean <- if (identical(case$link, "log")) exp(eta) else log1p(exp(eta))
    x <- mean + rnorm(length(loc), sd = 0.1)

    fit <- eb_smoother(
      x,
      s = NULL,
      family = "matern",
      locations = loc,
      link = case$link,
      backend = "fisher_pql"
    )
    ref <- eb_smoother(
      x,
      s = NULL,
      family = "matern",
      locations = loc,
      link = case$link,
      backend = case$ref_backend
    )

    expect_equal(fit$backend, "fisher_pql")
    expect_equal(fit$raw_fit$laplace_implementation, "exact_fisher_pql")
    expect_equal(fit$raw_fit$laplace_curvature, "fisher")
    expect_equal(fit$log_likelihood_semantics, "fisher_laplace_at_fisher_pql_mode_empirical_bayes")
    expect_equal(as.numeric(fit$log_likelihood), fit$raw_fit$log_likelihood_original_at_fisher_pql_mode, tolerance = 1e-12)
    expect_false("log_likelihood_laplace_at_fisher_pql_params" %in% names(fit$raw_fit))
    expect_true(is.finite(fit$raw_fit$log_likelihood_fisher_pql_stepA))
    expect_true(isTRUE(fit$raw_fit$fisher_pql_diagnostics$converged))
    diagnostics <- summary(fit)$diagnostics
    expect_true("log_likelihood_fisher_pql_stepA" %in% names(diagnostics))
    expect_true("log_likelihood_original_at_fisher_pql_mode" %in% names(diagnostics))
    expect_true("fisher_pql_diagnostics" %in% names(diagnostics))
    expect_equal(fit$raw_fit$fisher_pql_diagnostics$inner_iterations, 3L)
    expect_equal(fit$raw_fit$fisher_pql_diagnostics$stepA_engine, "exact_gaussian")
    expect_gt(fit$fitted_noise_sd, 0)
    expect_true(all(is.finite(fit$fitted_s)))
    expect_gt(min(fit$fitted_s), 0)
    expect_equal(fit$fitted_s, fit$raw_fit$fitted_s, tolerance = 1e-12)
    A <- fit$raw_fit$matern_objective_context$A
    eta_mode <- fit$raw_fit$fisher_pql_eta_mode
    expect_equal(eta_mode, as.numeric(A %*% fit$raw_fit$fisher_pql_mode) + fit$fitted_beta, tolerance = 1e-10)
    final_pseudo <- EBSmoothr:::.matern_fisher_pql_pseudo_response(
      x = x,
      eta = eta_mode,
      s = rep(fit$fitted_noise_sd, length(x)),
      link = case$link
    )
    expect_equal(fit$fitted_s, final_pseudo$s, tolerance = 1e-12)
    expect_true(all(is.finite(fit$posterior$mean)))
    expect_true(all(is.finite(fit$posterior$var)))
    expect_lt(max(abs(fit$posterior$mean - ref$posterior$mean)), 0.1)
    expect_lt(sqrt(mean((fit$posterior$mean - ref$posterior$mean)^2)), 0.05)
    ref_w <- ref$raw_fit$posterior_spatial_field$mode
    if (is.null(ref_w)) ref_w <- ref$raw_fit$posterior_spatial_field$mean
    eta_ref <- as.numeric(ref$raw_fit$matern_objective_context$A %*% as.numeric(ref_w)) + ref$fitted_beta
    expect_lt(max(abs(fit$raw_fit$fisher_pql_eta_mode - eta_ref)), 0.15)
    expect_lt(sqrt(mean((fit$raw_fit$fisher_pql_eta_mode - eta_ref)^2)), 0.05)
    expect_lt(abs(fit$fitted_g$theta - ref$fitted_g$theta), 0.15)
    expect_lt(abs(log(fit$fitted_g$sigma) - log(ref$fitted_g$sigma)), 0.15)
    expect_lt(abs(log(fit$fitted_noise_sd) - log(ref$fitted_noise_sd)), 0.15)
    expect_true(is.finite(as.numeric(fit$log_likelihood) - as.numeric(ref$log_likelihood)))
  }
})

test_that("Matern TMB Laplace validation catches unusable fits", {
  valid_fit <- list(
    posterior = data.frame(
      mean = c(1, 2),
      var = c(0.1, 0.2),
      second_moment = c(1.1, 4.2)
    ),
    log_likelihood = structure(1, class = "logLik"),
    laplace_diagnostics = list(stepB_joint_nll = 2)
  )
  expect_null(EBSmoothr:::.matern_laplace_tmb_invalid_reason(valid_fit))

  invalid_posterior <- valid_fit
  invalid_posterior$posterior$mean[1] <- Inf
  expect_match(
    EBSmoothr:::.matern_laplace_tmb_invalid_reason(invalid_posterior),
    "non-finite posterior"
  )

  invalid_stepB <- valid_fit
  invalid_stepB$laplace_diagnostics$stepB_joint_nll <- Inf
  expect_match(
    EBSmoothr:::.matern_laplace_tmb_invalid_reason(invalid_stepB),
    "Step B"
  )

  skipped_stepB <- valid_fit
  skipped_stepB$laplace_diagnostics$stepB_joint_nll <- NA_real_
  expect_null(EBSmoothr:::.matern_laplace_tmb_invalid_reason(skipped_stepB))
})

test_that("Public Matern Laplace backends do not fall back to R reference", {
  loc <- seq(0, 1, length.out = 6)
  s <- rep(0.1, length(loc))
  x <- exp(0.1 + 0.2 * sin(2 * pi * loc))

  expect_error(
    ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace", alpha = 1)(x, s),
    "require a successful TMB fit"
  )
  expect_error(
    ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_fisher", alpha = 1)(x, s),
    "require a successful TMB fit"
  )

  fit_r <- ebnm_Matern_generator(locations = loc, link = "log", backend = "laplace_r", alpha = 1)(x, s)
  expect_equal(fit_r$laplace_implementation, "r")
})

test_that("Selected-inverse variance helper matches solve reference", {
  old_options <- options(EBSmoothr.qinv_min_n = 1, EBSmoothr.qinv_max_row_nnz = 8)
  on.exit(options(old_options), add = TRUE)

  Q <- Matrix::bandSparse(
    8,
    k = c(-1, 0, 1),
    diagonals = list(rep(-0.1, 7), rep(2, 8), rep(-0.1, 7))
  )
  Q_dense_inv <- solve(as.matrix(Q))

  A_diag <- Matrix::Diagonal(8)
  diag_var <- EBSmoothr:::.compute_diag_A_Qinv_At(A_diag, Q)
  expect_equal(diag_var, diag(Q_dense_inv), tolerance = 1e-10)

  A_covered <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3, 3),
    j = c(1, 2, 2, 3, 5, 6),
    x = 1,
    dims = c(3, 8)
  )
  qinv <- INLA::inla.qinv(Q)
  expect_true(EBSmoothr:::.qinv_pattern_covers(A_covered, qinv))
  covered_var <- EBSmoothr:::.compute_diag_A_Qinv_At(A_covered, Q)
  covered_ref <- diag(as.matrix(A_covered) %*% Q_dense_inv %*% t(as.matrix(A_covered)))
  expect_equal(covered_var, covered_ref, tolerance = 1e-10)

  Q_diag <- Matrix::Diagonal(8, x = 2)
  A_uncovered <- Matrix::sparseMatrix(i = c(1, 1), j = c(1, 3), x = 1, dims = c(1, 8))
  expect_false(EBSmoothr:::.qinv_pattern_covers(A_uncovered, INLA::inla.qinv(Q_diag)))
  uncovered_var <- EBSmoothr:::.compute_diag_A_Qinv_At(A_uncovered, Q_diag)
  uncovered_ref <- diag(as.matrix(A_uncovered) %*% solve(as.matrix(Q_diag)) %*% t(as.matrix(A_uncovered)))
  expect_equal(uncovered_var, uncovered_ref, tolerance = 1e-10)

  options(EBSmoothr.qinv_max_row_nnz = 1)
  expect_false(EBSmoothr:::.qinv_pattern_covers(A_covered, qinv))
})

test_that("Matern wrapper supports known and learned-noise log links", {
  set.seed(14)

  loc <- seq(0, 1, length.out = 7)
  s <- rep(0.1, length(loc))
  x <- exp(0.1 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = s)

  fit <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    link = "log"
  )

  expect_equal(fit$backend, "fisher_pql")
  expect_equal(fit$raw_fit$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit$raw_fit$laplace_curvature, "fisher")
  expect_equal(fit$log_likelihood_semantics, "fisher_laplace_at_fisher_pql_mode_empirical_bayes")
  expect_equal(fit$raw_fit$fisher_pql_diagnostics$inner_iterations, 3L)
  expect_equal(fit$link, "log")
  expect_true(all(fit$posterior$mean > 0))

  fit_learned <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    link = "log"
  )
  expect_equal(fit_learned$backend, "fisher_pql")
  expect_equal(fit_learned$raw_fit$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit_learned$raw_fit$laplace_curvature, "fisher")
  expect_equal(fit_learned$log_likelihood_semantics, "fisher_laplace_at_fisher_pql_mode_empirical_bayes")
  expect_equal(fit_learned$raw_fit$fisher_pql_diagnostics$max_inner_iter, 3L)
  expect_gte(fit_learned$raw_fit$fisher_pql_diagnostics$inner_iterations, 1L)
  expect_lte(fit_learned$raw_fit$fisher_pql_diagnostics$inner_iterations, 3L)
  expect_equal(fit_learned$noise_mode, "estimated")
  expect_gt(fit_learned$fitted_noise_sd, 0)
  expect_true(all(fit_learned$posterior$mean > 0))
})

test_that("Matern log-link learned-noise Laplace implementations agree across beta modes", {
  set.seed(142)

  loc <- seq(0, 1, length.out = 8)
  x <- exp(0.1 + 0.18 * sin(2 * pi * loc)) + rnorm(length(loc), sd = 0.08)
  base_args <- list(x = x, s = NULL, family = "matern", locations = loc, link = "log")
  beta_cases <- list(
    list(beta_mode = "empirical_bayes"),
    list(beta_mode = "fixed", beta_fixed = 0.1),
    list(beta_mode = "prior_flat", beta_prec = 0),
    list(beta_mode = "prior_proper", beta_prec = 2)
  )

  for (case in beta_cases) {
    fit_tmb <- do.call(eb_smoother, c(base_args, list(backend = "laplace_tmb"), case[setdiff(names(case), "beta_mode")]))
    fit_r <- do.call(eb_smoother, c(base_args, list(backend = "laplace_r"), case[setdiff(names(case), "beta_mode")]))

    expect_equal(fit_tmb$beta_mode, case$beta_mode)
    expect_equal(fit_tmb$raw_fit$laplace_implementation, "tmb")
    expect_equal(as.numeric(fit_tmb$log_likelihood), as.numeric(fit_r$log_likelihood), tolerance = 2e-3)
    expect_equal(unname(fit_tmb$fitted_beta), unname(fit_r$fitted_beta), tolerance = 2e-3)
    expect_lt(max(abs(fit_tmb$posterior$mean - fit_r$posterior$mean)), 3e-3)
  }

  pc <- list(range = c(0.3, 0.5), sigma = c(0.2, 0.5), noise = c(0.08, 0.5))
  fit_pc_auto <- do.call(eb_smoother, c(base_args, list(pc.penalty = pc, beta_fixed = 0.1)))
  expect_equal(fit_pc_auto$backend, "fisher_pql")
  expect_equal(fit_pc_auto$raw_fit$laplace_curvature, "fisher")
  expect_equal(fit_pc_auto$raw_fit$fisher_pql_diagnostics$inner_iterations, 3L)
  expect_true(is.finite(as.numeric(fit_pc_auto$log_likelihood)))

  fit_pc_tmb <- do.call(eb_smoother, c(base_args, list(backend = "laplace_tmb", pc.penalty = pc, beta_fixed = 0.1)))
  fit_pc_r <- do.call(eb_smoother, c(base_args, list(backend = "laplace_r", pc.penalty = pc, beta_fixed = 0.1)))
  expect_equal(fit_pc_tmb$backend, "laplace")
  expect_equal(fit_pc_tmb$raw_fit$laplace_implementation, "tmb")
  expect_equal(as.numeric(fit_pc_tmb$log_likelihood), as.numeric(fit_pc_r$log_likelihood), tolerance = 2e-3)
})

test_that("Matern log-link learned-noise INLA validates against package Laplace", {
  set.seed(143)

  loc <- seq(0, 1, length.out = 10)
  x <- exp(0.1 + 0.15 * sin(2 * pi * loc)) + rnorm(length(loc), sd = 0.08)

  fit_or_error <- tryCatch(
    eb_smoother(
      x,
      s = NULL,
      family = "matern",
      locations = loc,
      link = "log",
      backend = "inla"
    ),
    error = function(e) e
  )
  if (inherits(fit_or_error, "error")) {
    expect_match(conditionMessage(fit_or_error), "INLA Matern validation failed|INLA learned-noise Matern noise profiling failed")
    expect_false(grepl("not currently supported", conditionMessage(fit_or_error), fixed = TRUE))
  } else {
    expect_equal(fit_or_error$backend, "inla")
    expect_equal(fit_or_error$raw_fit$inla_validation$status, "passed")
  }
})

test_that("Matern INLA log-link likelihood uses spatial-mode warm starts", {
  set.seed(141)

  loc <- seq(0, 1, length.out = 120)
  s <- rep(0.1, length(loc))
  x <- exp(0.15 + 0.25 * sin(2 * pi * loc) + 0.08 * cos(6 * pi * loc)) +
    rnorm(length(loc), sd = s)
  setup <- Matern_setup(loc)

  fit_tmb <- ebnm_Matern_generator(setup = setup, link = "log", backend = "laplace_tmb")(x, s, beta_fixed = 0)
  fit_inla <- ebnm_Matern_generator(setup = setup, link = "log", backend = "inla")(x, s, beta_fixed = 0)

  expect_equal(as.numeric(fit_inla$log_likelihood), fit_inla$log_likelihood_laplace_at_inla_mode, tolerance = 1e-10)
  expect_equal(as.numeric(fit_tmb$log_likelihood), as.numeric(fit_inla$log_likelihood), tolerance = 1e-2)
  expect_lt(max(abs(fit_tmb$posterior$mean - fit_inla$posterior$mean)), 0.01)
})

test_that("Identity-link Laplace backend matches exact Matern backend", {
  set.seed(15)

  loc <- seq(0, 1, length.out = 8)
  s <- rep(0.1, length(loc))
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_exact <- ebnm_Matern_generator(locations = loc, backend = "exact")(x, s)
  fit_laplace <- ebnm_Matern_generator(locations = loc, backend = "laplace")(x, s)

  expect_equal(as.numeric(fit_exact$log_likelihood), as.numeric(fit_laplace$log_likelihood), tolerance = 1e-5)
  expect_equal(fit_exact$posterior$mean, fit_laplace$posterior$mean, tolerance = 1e-4)
  expect_equal(fit_laplace$laplace_implementation, "tmb")
})

test_that("Quiet INLA defaults provide an explicit thread specification", {
  opts <- EBSmoothr:::.with_quiet_inla_defaults(INLA:::inla.getOption.default())

  expect_match(opts$num.threads, "^[0-9]+:1$")
  expect_false(grepl("^NA", opts$num.threads))
})

test_that("Exact known-noise Matern fit uses empirical-Bayes beta by default", {
  set.seed(1)

  loc <- seq(0, 1, length.out = 25)
  s <- rep(0.1, length(loc))
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit <- ebnm_Matern_generator(locations = loc)(x, s)

  expect_equal(fit$backend, "exact")
  expect_equal(fit$beta_mode, "empirical_bayes")
  expect_s3_class(fit$fitted_g, "Matern")
  expect_equal(unname(fit$fitted_g$beta), unname(fit$fitted_beta), tolerance = 1e-8)
  expect_null(fit$fitted_g$beta_prec)
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
  expect_equal(nrow(fit$posterior), length(loc))
})

test_that("Identity-link INLA backends report the package exact objective", {
  set.seed(1011)

  loc <- seq(0, 1, length.out = 18)
  s <- rep(0.1, length(loc))
  x <- 0.2 + 0.4 * sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_exact <- ebnm_Matern_generator(locations = loc, backend = "exact")(x, s)
  fit_inla <- ebnm_Matern_generator(locations = loc, backend = "inla")(x, s)
  expect_equal(fit_inla$log_likelihood_semantics, "exact_at_inla_mode_empirical_bayes")
  expect_equal(as.numeric(fit_inla$log_likelihood), fit_inla$log_likelihood_exact_at_stepA_mode, tolerance = 1e-10)
  expect_equal(as.numeric(fit_exact$log_likelihood), as.numeric(fit_inla$log_likelihood), tolerance = 1e-3)
  expect_lt(max(abs(fit_exact$posterior$mean - fit_inla$posterior$mean)), 1e-3)

  fit_wrap <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    backend = "inla"
  )
  expect_equal(fit_wrap$log_likelihood_semantics, "exact_at_inla_mode_empirical_bayes")
  expect_equal(as.numeric(fit_wrap$log_likelihood), fit_wrap$log_likelihood_exact_at_stepA_mode, tolerance = 1e-10)
})

test_that("Matern generator with PC prior stays exact under empirical-Bayes beta", {
  set.seed(2)

  loc <- seq(0, 1, length.out = 20)
  s <- rep(0.12, length(loc))
  x <- sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = 0.2, sigma = 0.3)
  )(x, s)

  expect_equal(fit$backend, "exact")
  expect_equal(fit$beta_mode, "empirical_bayes")
  expect_equal(fit$prior_family, "identity_Matern_pc_exact")
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
  expect_true(is.finite(fit$log_likelihood_pc_prior_theta))
})

test_that("Matern wrapper dispatches PC-prior fits by beta mode", {
  set.seed(3)

  loc <- seq(0, 1, length.out = 20)
  x <- sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)
  pc <- list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))

  fit_eb <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    pc.penalty = pc
  )
  expect_equal(fit_eb$backend, "exact")
  expect_equal(fit_eb$beta_mode, "empirical_bayes")

  fit_flat <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    pc.penalty = pc,
    beta_prec = 0
  )
  expect_equal(fit_flat$backend, "exact")
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_equal(unname(fit_flat$fitted_g$beta_prec), 0, tolerance = 1e-10)

  fit_fixed <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    pc.penalty = pc,
    beta_fixed = 0
  )
  expect_equal(fit_fixed$backend, "exact")
  expect_equal(fit_fixed$beta_mode, "fixed")
  expect_equal(unname(fit_fixed$fitted_beta), 0, tolerance = 1e-10)
})

test_that("Matern exact backend supports fixed and prior beta modes with PC prior", {
  set.seed(4)

  loc <- seq(0, 1, length.out = 20)
  x <- 0.3 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)
  pc <- list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))

  fit_fixed <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    backend = "exact",
    pc.penalty = pc,
    beta_fixed = 0.3
  )
  expect_equal(fit_fixed$backend, "exact")
  expect_equal(fit_fixed$beta_mode, "fixed")
  expect_equal(unname(fit_fixed$fitted_beta), 0.3, tolerance = 1e-10)

  fit_flat <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    backend = "exact",
    pc.penalty = pc,
    beta_prec = 0
  )
  expect_equal(fit_flat$backend, "exact")
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_equal(unname(fit_flat$fitted_g$beta_prec), 0, tolerance = 1e-10)

  fit_proper <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    backend = "exact",
    pc.penalty = pc,
    beta_prec = 2
  )
  expect_equal(fit_proper$backend, "exact")
  expect_equal(fit_proper$beta_mode, "prior_proper")
  expect_equal(unname(fit_proper$fitted_g$beta_prec), 2, tolerance = 1e-10)
})

test_that("Matern INLA backends support empirical-Bayes beta", {
  loc <- seq(0, 1, length.out = 8)
  x <- sin(2 * pi * loc)

  fit_gen <- ebnm_Matern_generator(
    locations = loc,
    link = "log",
    backend = "inla"
  )(exp(0.1 + x), 0.1)
  expect_equal(fit_gen$backend, "inla")
  expect_equal(fit_gen$beta_mode, "empirical_bayes")
  expect_true(is.finite(fit_gen$fitted_beta))
  expect_false(is.null(fit_gen$beta_profile_optimization))

  fit_wrap <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    backend = "inla",
    pc.penalty = list(range = 0.2, sigma = 0.3)
  )
  expect_equal(fit_wrap$backend, "inla")
  expect_equal(fit_wrap$beta_mode, "empirical_bayes")
  expect_true(is.finite(fit_wrap$fitted_beta))
  expect_false(is.null(fit_wrap$beta_profile_optimization))
})

test_that("Learned-noise Matern supports exact empirical-Bayes beta with PC prior", {
  set.seed(5)

  loc <- seq(0, 1, length.out = 20)
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.08)

  fit <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    pc.penalty = list(range = 0.2, sigma = 0.3, noise = 0.1)
  )

  expect_equal(fit$backend, "exact")
  expect_equal(fit$beta_mode, "empirical_bayes")
  expect_true(is.finite(fit$fitted_noise_sd))
  expect_gt(fit$fitted_noise_sd, 0)
})

test_that("Exact learned-noise Matern supports fixed observation noise scales", {
  set.seed(316)

  loc <- seq(0, 1, length.out = 14)
  setup <- Matern_setup(locations = loc, max.edge = 0.35)
  x <- 0.15 + 0.4 * sin(2 * pi * loc) + rnorm(length(loc), sd = 0.08)
  base_args <- list(
    x = x,
    A = setup$A,
    spde_template = setup$spde_template,
    alpha = setup$alpha,
    d = setup$d,
    theta_init = log(0.3),
    sigma_init = 0.35,
    noise_sd_init = 0.08,
    suppress_warnings = TRUE
  )

  fit_null <- do.call(
    .fit_matern_exact_unknown_noise,
    c(base_args, list(beta_mode = "empirical_bayes"))
  )
  fit_unit <- do.call(
    .fit_matern_exact_unknown_noise,
    c(base_args, list(beta_mode = "empirical_bayes", noise_scale = rep(1, length(x))))
  )
  expect_equal(as.numeric(fit_unit$log_likelihood), as.numeric(fit_null$log_likelihood), tolerance = 1e-8)
  expect_equal(fit_unit$fitted_noise_sd, fit_null$fitted_noise_sd, tolerance = 1e-8)
  expect_equal(fit_unit$fitted_s, rep(fit_unit$fitted_noise_sd, length(x)), tolerance = 1e-10)

  noise_scale <- seq(0.75, 1.25, length.out = length(x))
  cases <- list(
    list(beta_mode = "empirical_bayes"),
    list(beta_mode = "fixed", beta_fixed = 0.1),
    list(beta_mode = "prior_flat", beta_prec = 0),
    list(beta_mode = "prior_proper", beta_prec = 2)
  )

  for (case_args in cases) {
    fit <- do.call(
      .fit_matern_exact_unknown_noise,
      c(base_args, case_args, list(noise_scale = noise_scale))
    )
    expect_true(is.finite(fit$fitted_noise_sd))
    expect_true(all(is.finite(fit$fitted_s)))
    expect_equal(fit$fitted_s, fit$fitted_noise_sd * noise_scale, tolerance = 1e-10)

    objective_args <- list(
      x = x,
      s = fit$fitted_s,
      A = setup$A,
      spde_template = setup$spde_template,
      alpha = setup$alpha,
      d = setup$d,
      log_range = fit$fitted_g$theta,
      log_sigma = log(fit$fitted_g$sigma),
      beta_mode = case_args$beta_mode,
      suppress_warnings = TRUE
    )
    if (!is.null(case_args$beta_fixed)) objective_args$beta_fixed <- case_args$beta_fixed
    if (!is.null(case_args$beta_prec)) objective_args$beta_prec <- case_args$beta_prec
    expected_objective <- do.call(.exact_matern_known_noise_objective_at_params, objective_args)
    expect_equal(
      as.numeric(fit$log_likelihood),
      as.numeric(expected_objective$log_marginal),
      tolerance = 1e-6
    )
  }

  expect_error(
    do.call(
      .fit_matern_exact_unknown_noise,
      c(base_args, list(beta_mode = "empirical_bayes", noise_scale = noise_scale[-1]))
    ),
    "noise_scale"
  )
  expect_error(
    do.call(
      .fit_matern_exact_unknown_noise,
      c(base_args, list(beta_mode = "empirical_bayes", noise_scale = replace(noise_scale, 1, -1)))
    ),
    "positive finite"
  )
})

test_that("Learned-noise Matern INLA path supports flat and proper beta priors", {
  set.seed(6)

  loc <- seq(0, 1, length.out = 20)
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.08)
  pc <- list(range = 0.2, sigma = 0.3, noise = 0.1)

  fit_flat <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    backend = "inla",
    pc.penalty = pc,
    beta_prec = 0,
    compute_exact_diagnostic = TRUE
  )
  expect_equal(fit_flat$backend, "inla")
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_true(is.finite(fit_flat$log_likelihood_stepA_penalized))
  expect_true(is.finite(fit_flat$log_likelihood_exact_at_stepA_mode))

  fit_proper <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    backend = "inla",
    pc.penalty = pc,
    beta_prec = 2
  )
  expect_equal(fit_proper$backend, "inla")
  expect_equal(fit_proper$beta_mode, "prior_proper")
  expect_equal(unname(fit_proper$fitted_g$beta_prec), 2, tolerance = 1e-10)
})

test_that("matern_objective_breakdown follows the new beta mode semantics", {
  set.seed(7)

  loc <- seq(0, 1, length.out = 18)
  x <- 0.3 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)

  fit_exact <- ebnm_Matern_generator(locations = loc)(x, 0.1)
  br_exact <- matern_objective_breakdown(fit_exact)
  expect_equal(br_exact$matched_beta_mode, "empirical_bayes")
  expect_true(is.finite(br_exact$matched_exact_objective))

  fit_pc_exact <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    backend = "exact",
    pc.penalty = list(range = 0.2, sigma = 0.3),
    beta_prec = 0
  )$raw_fit
  class(fit_pc_exact) <- c("ebnm", "list")
  br_flat <- matern_objective_breakdown(fit_pc_exact)
  expect_equal(br_flat$matched_beta_mode, "prior_flat")
  expect_true(is.finite(br_flat$matched_exact_objective))
})

test_that("Matern softplus inlabru rejects non-softplus links", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("inlabru")

  set.seed(310)
  loc <- seq(0, 1, length.out = 12)
  s <- rep(0.1, length(loc))
  x_log <- exp(0.1 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = s)
  x_id <- 0.2 * sin(2 * pi * loc) + rnorm(length(loc), sd = s)
  expect_error(
    ebnm_Matern_generator(locations = loc, link = "log", backend = "inlabru")(x_log, s),
    "only supported for `link = \"softplus\"`"
  )
  expect_error(
    ebnm_Matern_generator(locations = loc, link = "identity", backend = "inlabru")(x_id, s),
    "only supported for `link = \"softplus\"`"
  )
})

test_that("Matern softplus inlabru backend agrees with Laplace (s known, ebnm path)", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("inlabru")

  set.seed(311)
  loc <- seq(0, 1, length.out = 60)
  s <- rep(0.08, length(loc))
  eta <- -0.4 + 1.2 * sin(2 * pi * loc)
  x <- log1p(exp(eta)) + rnorm(length(loc), sd = s)

  fit_lap <- ebnm_Matern_generator(
    locations = loc, link = "softplus", backend = "laplace"
  )(x, s)
  fit_bru <- tryCatch(
    ebnm_Matern_generator(
      locations = loc, link = "softplus", backend = "inlabru"
    )(x, s),
    error = function(e) e
  )
  if (inherits(fit_bru, "error")) {
    skip(paste("inlabru iterative INLA did not converge:",
               conditionMessage(fit_bru)))
  }

  expect_equal(fit_bru$backend, "inlabru")
  expect_equal(fit_bru$link, "softplus")
  expect_null(fit_bru$pc_penalty)
  expect_equal(fit_bru$log_likelihood_semantics, "laplace_at_inlabru_params_empirical_bayes")
  expect_true(all(is.finite(fit_bru$posterior$mean)))
  expect_true(all(fit_bru$posterior$var > 0))
  ctx <- fit_bru$matern_objective_context
  expected_objective <- .matern_laplace_known_noise_objective_at_params(
    x = x,
    s = s,
    A = ctx$A,
    spde_template = ctx$spde_template,
    alpha = ctx$alpha,
    d = ctx$d,
    log_range = fit_bru$fitted_g$theta,
    log_sigma = log(fit_bru$fitted_g$sigma),
    beta_mode = fit_bru$beta_mode,
    beta_init = fit_bru$fitted_beta,
    beta_prec = fit_bru$beta_prec,
    link = "softplus",
    pc_penalty = fit_bru$pc_penalty,
    initial_mode = .matern_inla_laplace_initial_mode(fit_bru, fit_bru$beta_mode),
    compute_posterior = FALSE,
    optimize_beta = FALSE
  )
  expect_equal(
    as.numeric(fit_bru$log_likelihood),
    as.numeric(expected_objective$log_marginal),
    tolerance = 1e-8
  )
  expect_true(is.finite(fit_bru$log_likelihood_inlabru_mlik_integration))
  expect_gt(
    abs(as.numeric(fit_bru$log_likelihood) -
          as.numeric(fit_bru$log_likelihood_inlabru_mlik_integration)),
    1e-4
  )
  # Cross-method posterior agreement: tolerances reflect that bru uses
  # iterative linearisation while Laplace uses observed-Hessian curvature
  # at the joint mode, so per-observation moments and hyperparameters need
  # not match bitwise.
  expect_lt(max(abs(fit_bru$posterior$mean - fit_lap$posterior$mean)), 0.10)
  expect_lt(max(abs(fit_bru$posterior$var - fit_lap$posterior$var)), 0.10)
  expect_lt(abs(fit_bru$fitted_beta - fit_lap$fitted_beta), 0.30)
  expect_lt(abs(log(fit_bru$fitted_g$sigma) - log(fit_lap$fitted_g$sigma)), 1.0)
  expect_lt(abs(fit_bru$fitted_g$theta - fit_lap$fitted_g$theta), 1.0)
})

test_that("Matern softplus inlabru uses only explicit PC priors", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("inlabru")

  set.seed(313)
  loc <- seq(0, 1, length.out = 36)
  s <- rep(0.08, length(loc))
  eta <- -0.2 + 0.8 * sin(2 * pi * loc)
  x <- log1p(exp(eta)) + rnorm(length(loc), sd = s)
  pc <- list(range = c(0.15, 0.4), sigma = c(0.5, 0.45))

  fit_bru <- tryCatch(
    ebnm_Matern_generator(
      locations = loc, link = "softplus", backend = "inlabru",
      pc.penalty = pc
    )(x, s),
    error = function(e) e
  )
  if (inherits(fit_bru, "error")) {
    skip(paste("inlabru iterative INLA did not converge:",
               conditionMessage(fit_bru)))
  }

  expect_equal(fit_bru$backend, "inlabru")
  expect_equal(unname(fit_bru$pc_penalty$range), unname(pc$range))
  expect_equal(unname(fit_bru$pc_penalty$sigma), unname(pc$sigma))
  expect_equal(fit_bru$log_likelihood_semantics, "laplace_at_inlabru_params_empirical_bayes")
})

test_that("Matern softplus inlabru backend agrees with Laplace (s = NULL, eb_smoother path)", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("inlabru")

  # Use signal in the near-linear regime of softplus (eta well above 0)
  # because iterative INLA can otherwise wander into NaN regions when the
  # inferred noise precision diverges with default flat priors.
  set.seed(312)
  loc <- seq(0, 1, length.out = 50)
  eta <- 1.2 + 0.6 * sin(2 * pi * loc)
  x <- log1p(exp(eta)) + rnorm(length(loc), sd = 0.08)
  pc <- list(range = c(0.1, 0.5), sigma = c(0.35, 0.5), noise = c(0.2, 0.5))

  expect_error(
    eb_smoother(
      x, s = NULL, family = "matern", locations = loc,
      link = "softplus", backend = "inlabru"
    ),
    "requires an explicit `pc.penalty`"
  )
  expect_error(
    eb_smoother(
      x, s = NULL, family = "matern", locations = loc,
      link = "softplus", backend = "inlabru",
      pc.penalty = list(range = c(0.1, 0.5), sigma = c(0.35, 0.5))
    ),
    "requires explicit `pc.penalty` entries"
  )

  fit_lap <- eb_smoother(
    x, s = NULL, family = "matern", locations = loc,
    link = "softplus", backend = "laplace", pc.penalty = pc
  )
  fit_bru <- tryCatch(
    eb_smoother(
      x, s = NULL, family = "matern", locations = loc,
      link = "softplus", backend = "inlabru", pc.penalty = pc
    ),
    error = function(e) e
  )
  if (inherits(fit_bru, "error")) {
    skip(paste("inlabru iterative INLA did not converge:",
               conditionMessage(fit_bru)))
  }

  expect_equal(fit_bru$backend, "inlabru")
  expect_equal(fit_bru$log_likelihood_semantics, "laplace_at_inlabru_params_learned_noise_empirical_bayes")
  expect_equal(unname(fit_bru$pc_penalty$noise), unname(pc$noise))
  expect_true(is.finite(fit_bru$fitted_noise_sd))
  ctx <- fit_bru$matern_objective_context
  expected_objective <- .matern_laplace_unknown_noise_objective_at_params(
    x = x,
    A = ctx$A,
    spde_template = ctx$spde_template,
    alpha = ctx$alpha,
    d = ctx$d,
    log_range = fit_bru$fitted_g$theta,
    log_sigma = log(fit_bru$fitted_g$sigma),
    log_noise_sd = log(fit_bru$fitted_noise_sd),
    beta_mode = fit_bru$beta_mode_internal,
    beta_init = fit_bru$fitted_beta,
    beta_prec = fit_bru$beta_prec,
    link = "softplus",
    pc_penalty = fit_bru$pc_penalty,
    initial_mode = .matern_inla_laplace_initial_mode(fit_bru, fit_bru$beta_mode_internal),
    compute_posterior = FALSE,
    optimize_beta = FALSE
  )
  expect_equal(
    as.numeric(fit_bru$log_likelihood),
    as.numeric(expected_objective$log_marginal),
    tolerance = 1e-8
  )
  expect_true(is.finite(fit_bru$log_likelihood_inlabru_mlik_integration))
  expect_lt(
    abs(fit_bru$fitted_noise_sd - fit_lap$fitted_noise_sd) /
      fit_lap$fitted_noise_sd,
    0.30
  )
  expect_lt(max(abs(fit_bru$posterior$mean - fit_lap$posterior$mean)), 0.20)
})

test_that("direct SPDE precision assembly matches INLA across hyperparameters", {
  skip_if_not_installed("INLA")

  loc <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 5),
    y = seq(0, 1, length.out = 5)
  ))
  setup <- Matern_setup(locations = loc, max.edge = 0.6)
  ctx <- attr(setup$spde_template, "EBSmoothr_direct_ctx", exact = TRUE)
  expect_false(is.null(ctx))

  for (pars in list(c(log(0.3), log(0.5)), c(log(1.2), log(2)))) {
    pr <- EBSmoothr:::.matern_precision_from_log_params(
      spde_template = setup$spde_template,
      alpha = setup$alpha,
      d = setup$d,
      log_range = pars[1],
      log_sigma = pars[2]
    )
    theta_spde <- EBSmoothr:::.matern_spde_theta_from_log_range_log_sigma(
      pars[1], pars[2], setup$alpha, setup$d
    )
    Q_ref <- Matrix::forceSymmetric(Matrix::Matrix(
      INLA::inla.spde2.precision(setup$spde_template, theta = theta_spde),
      sparse = TRUE
    ))
    expect_lt(
      max(abs(pr$Q - Q_ref)),
      1e-8 * max(abs(Q_ref@x))
    )
    expect_equal(
      pr$Q_factor$logdet,
      as.numeric(Matrix::determinant(Q_ref, logarithm = TRUE)$modulus),
      tolerance = 1e-8
    )
  }

  # Templates without the precomputed context fall back to INLA.
  template_plain <- setup$spde_template
  attr(template_plain, "EBSmoothr_direct_ctx") <- NULL
  pr_fallback <- EBSmoothr:::.matern_precision_from_log_params(
    spde_template = template_plain,
    alpha = setup$alpha,
    d = setup$d,
    log_range = log(0.7),
    log_sigma = log(1.1)
  )
  pr_direct <- EBSmoothr:::.matern_precision_from_log_params(
    spde_template = setup$spde_template,
    alpha = setup$alpha,
    d = setup$d,
    log_range = log(0.7),
    log_sigma = log(1.1)
  )
  expect_lt(max(abs(pr_fallback$Q - pr_direct$Q)), 1e-8 * max(abs(pr_direct$Q@x)))
})

test_that("sparse SPD factorization cache reuses symbolic analyses correctly", {
  set.seed(404)
  n <- 60
  B <- Matrix::rsparsematrix(n, n, density = 0.08)
  Q1 <- Matrix::forceSymmetric(Matrix::crossprod(B) + Matrix::Diagonal(n, 1))
  Q2 <- Matrix::forceSymmetric(Matrix::crossprod(B) + Matrix::Diagonal(n, 4))

  f1 <- EBSmoothr:::.factorize_spd(Q1)
  f2 <- EBSmoothr:::.factorize_spd(Q2)

  expect_equal(
    f1$logdet,
    as.numeric(Matrix::determinant(Q1, logarithm = TRUE)$modulus),
    tolerance = 1e-10
  )
  expect_equal(
    f2$logdet,
    as.numeric(Matrix::determinant(Q2, logarithm = TRUE)$modulus),
    tolerance = 1e-10
  )

  rhs <- rnorm(n)
  expect_lt(
    max(abs(
      as.numeric(EBSmoothr:::.solve_spd_factor(f2, rhs)) -
        as.numeric(Matrix::solve(Q2, rhs))
    )),
    1e-8
  )

  # A different pattern with the same dimension must not reuse the analysis.
  Q3 <- Matrix::forceSymmetric(
    Matrix::crossprod(Matrix::rsparsematrix(n, n, density = 0.12)) +
      Matrix::Diagonal(n, 1)
  )
  f3 <- EBSmoothr:::.factorize_spd(Q3)
  expect_equal(
    f3$logdet,
    as.numeric(Matrix::determinant(Q3, logarithm = TRUE)$modulus),
    tolerance = 1e-10
  )

  cache <- EBSmoothr:::.spd_factor_cache
  old_cache_max <- getOption("EBSmoothr.spd_factor_cache_max")
  on.exit({
    options(EBSmoothr.spd_factor_cache_max = old_cache_max)
    rm(list = ls(cache), envir = cache)
  }, add = TRUE)

  rm(list = ls(cache), envir = cache)
  options(EBSmoothr.spd_factor_cache_max = 2L)
  invisible(EBSmoothr:::.factorize_spd(Q1))
  invisible(EBSmoothr:::.factorize_spd(Q3))
  invisible(EBSmoothr:::.factorize_spd(Matrix::Diagonal(n, 2)))
  n_cached <- sum(vapply(ls(cache), function(k) length(cache[[k]]), integer(1)))
  expect_equal(n_cached, 2L)
})
