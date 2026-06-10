test_that("known-noise Matern wrapper matches the ebnm generator under default beta semantics", {
  set.seed(101)

  loc <- seq(0, 1, length.out = 25)
  s <- 0.12
  x <- 0.3 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_gen <- ebnm_Matern_generator(locations = loc)(x, s)
  fit_wrap <- eb_smoother(x, s = s, family = "matern", locations = loc)

  expect_s3_class(fit_wrap, "eb_smoother_fit")
  expect_equal(fit_wrap$backend, "exact")
  expect_equal(fit_wrap$beta_mode, "empirical_bayes")
  expect_equal(fit_wrap$posterior$mean, fit_gen$posterior$mean, tolerance = 1e-8)
  expect_equal(fit_wrap$fitted_g$theta, fit_gen$fitted_g$theta, tolerance = 1e-8)
  expect_equal(fit_wrap$fitted_g$beta, fit_gen$fitted_g$beta, tolerance = 1e-8)
})

test_that("constant wrapper supports empirical-Bayes, fixed, flat, and proper beta modes", {
  x <- c(1, 2, 3, 4)
  s <- c(0.5, 1, 1.5, 2)

  fit_eb <- eb_smoother(x, s = s, family = "constant")
  expect_equal(fit_eb$beta_mode, "empirical_bayes")
  expect_equal(unname(fit_eb$fitted_g$beta), unname(fit_eb$fitted_beta), tolerance = 1e-10)
  expect_null(fit_eb$fitted_g$beta_prec)

  fit_fixed <- eb_smoother(x, s = s, family = "constant", beta_fixed = 0)
  expect_equal(fit_fixed$beta_mode, "fixed")
  expect_equal(unname(fit_fixed$fitted_beta), 0, tolerance = 1e-10)

  fit_flat <- eb_smoother(x, s = s, family = "constant", beta_prec = 0)
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_equal(unname(fit_flat$fitted_g$beta_prec), 0, tolerance = 1e-10)
  expect_true(all(fit_flat$posterior$var > 0))

  fit_proper <- eb_smoother(x, s = s, family = "constant", beta_prec = 2)
  expect_equal(fit_proper$beta_mode, "prior_proper")
  expect_equal(unname(fit_proper$fitted_g$beta_prec), 2, tolerance = 1e-10)
  expect_true(all(fit_proper$posterior$var > 0))
})

test_that("learned-noise constant wrapper supports flat and proper beta priors", {
  set.seed(102)

  x <- 0.3 + rnorm(30, sd = 0.15)

  fit_flat <- eb_smoother(x, s = NULL, family = "constant", beta_prec = 0)
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_true(is.finite(fit_flat$fitted_noise_sd))
  expect_gt(fit_flat$fitted_noise_sd, 0)
  expect_true(all(fit_flat$posterior$var > 0))

  fit_proper <- eb_smoother(x, s = NULL, family = "constant", beta_prec = 1.5)
  expect_equal(fit_proper$beta_mode, "prior_proper")
  expect_equal(unname(fit_proper$fitted_g$beta_prec), 1.5, tolerance = 1e-10)
  expect_true(is.finite(fit_proper$fitted_noise_sd))
  expect_gt(fit_proper$fitted_noise_sd, 0)
})

test_that("constant log-link supports known and learned noise", {
  x <- c(0.95, 1.15, 1.28, 1.05, 1.18)
  s <- c(0.08, 0.1, 0.09, 0.12, 0.1)
  beta0 <- log(1.1)

  fit_fixed <- eb_smoother(
    x,
    s = s,
    family = "constant",
    link = "log",
    beta_fixed = beta0
  )
  direct_ll <- sum(stats::dnorm(x, mean = exp(beta0), sd = s, log = TRUE))
  expect_equal(as.numeric(fit_fixed$log_likelihood), direct_ll, tolerance = 1e-10)
  expect_equal(unique(fit_fixed$posterior$mean), exp(beta0), tolerance = 1e-10)

  fit_eb <- eb_smoother(x, s = s, family = "constant", link = "log")
  beta_mle <- log(sum(x / s^2) / sum(1 / s^2))
  expect_equal(unname(fit_eb$fitted_beta), beta_mle, tolerance = 1e-10)
  expect_equal(fit_eb$beta_mode, "empirical_bayes")

  fit_flat <- eb_smoother(x, s = s, family = "constant", link = "log", beta_prec = 0)
  expect_equal(fit_flat$beta_mode, "prior_flat")
  expect_equal(fit_flat$log_likelihood_semantics, "laplace_constant_log_prior_flat")
  expect_true(is.finite(as.numeric(fit_flat$log_likelihood)))
  expect_true(all(fit_flat$posterior$mean > 0))

  fit_learned <- eb_smoother(x, s = NULL, family = "constant", link = "log")
  expected_noise <- sqrt(mean((x - exp(fit_learned$fitted_beta))^2))
  expect_equal(fit_learned$fitted_noise_sd, expected_noise, tolerance = 1e-10)
  expect_true(all(fit_learned$posterior$mean > 0))

  fit_learned_flat <- eb_smoother(x, s = NULL, family = "constant", link = "log", beta_prec = 0)
  expect_equal(fit_learned_flat$beta_mode, "prior_flat")
  expect_true(is.finite(as.numeric(fit_learned_flat$log_likelihood)))
  expect_gt(fit_learned_flat$fitted_noise_sd, 0)
})

test_that("eb_smoother validates beta arguments and backend dispatch", {
  loc <- seq(0, 1, length.out = 10)
  x <- sin(2 * pi * loc)
  setup_lgp <- LGP_setup(t = loc)

  expect_error(
    eb_smoother(
      x,
      s = 0.1,
      family = "matern",
      locations = loc,
      beta_fixed = 0,
      beta_prec = 1
    ),
    "cannot be supplied together"
  )

  expect_error(
    eb_smoother(
      x,
      s = 0.1,
      family = "matern",
      locations = loc,
      backend = "inla_pc"
    ),
    "requires `pc.penalty`"
  )

  expect_error(
    eb_smoother(x, s = 0.1, family = "constant", backend = "inla_pc"),
    "only supports the exact backend"
  )

  expect_error(
    eb_smoother(x, s = 0.1, family = "nonspatial"),
    "should be one of"
  )

  expect_error(
    eb_smoother(x, s = 0.1, family = "lgp", backend = "inla_pc", setup = setup_lgp),
    "only supported for `family = \"matern\"`"
  )
})

test_that("point-exponential wrapper supports fixed and profiled noise", {
  x_sparse <- c(0, 0, 0.15, 0, 1.2, 0, 2.1)

  fit_fixed <- eb_smoother(x_sparse, s = 0.1, family = "point_exponential")
  fit_direct <- ebnm::ebnm_point_exponential(x_sparse, s = rep(0.1, length(x_sparse)))
  expect_s3_class(fit_fixed, "eb_smoother_fit")
  expect_equal(fit_fixed$family, "point_exponential")
  expect_equal(fit_fixed$noise_mode, "fixed")
  expect_equal(as.numeric(fit_fixed$log_likelihood), as.numeric(fit_direct$log_likelihood), tolerance = 1e-10)
  expect_true(all(is.finite(fit_fixed$posterior$mean)))
  expect_true(all(is.finite(fit_fixed$posterior$var)))
  expect_true(is.finite(fit_fixed$point_exponential_pi0))
  expect_true(is.finite(fit_fixed$point_exponential_scale))

  fit_zero <- eb_smoother(rep(0, 6), s = NULL, family = "point_exponential")
  expect_equal(fit_zero$noise_mode, "estimated")
  expect_true(is.finite(fit_zero$fitted_noise_sd))
  expect_true(is.finite(as.numeric(fit_zero$log_likelihood)))
  expect_true(all(fit_zero$posterior$mean >= 0))
  expect_equal(fit_zero$fitted_noise_sd, fit_zero$profile_optimization$selected_s, tolerance = 1e-14)

  x_smooth <- exp(0.2 + 0.3 * sin(seq(0, 2 * pi, length.out = 12)))
  fit_smooth <- eb_smoother(x_smooth, s = NULL, family = "point_exponential")
  expect_true(is.finite(fit_smooth$fitted_noise_sd))
  expect_true(is.finite(as.numeric(fit_smooth$log_likelihood)))
  expect_true(fit_smooth$point_exponential_pi0 >= 0 && fit_smooth$point_exponential_pi0 <= 1)
})

test_that("point-normal and point-laplace wrappers support fixed and profiled noise", {
  x <- c(-1.1, 0, 0.05, 0.7, 1.4)

  for (family in c("point_normal", "point_laplace")) {
    ebnm_fun <- switch(
      family,
      point_normal = ebnm::ebnm_point_normal,
      point_laplace = ebnm::ebnm_point_laplace
    )

    fit_fixed <- eb_smoother(x, s = 0.2, family = family)
    fit_direct <- ebnm_fun(x, s = rep(0.2, length(x)))
    expect_s3_class(fit_fixed, "eb_smoother_fit")
    expect_equal(fit_fixed$family, family)
    expect_equal(as.numeric(fit_fixed$log_likelihood), as.numeric(fit_direct$log_likelihood), tolerance = 1e-10)
    expect_true(is.finite(fit_fixed$point_mass_probability))
    expect_true(is.finite(fit_fixed$point_nonzero_scale))

    fit_profile <- eb_smoother(x, s = NULL, family = family, profile_s_upper = 2)
    expect_equal(fit_profile$noise_mode, "estimated")
    expect_true(is.finite(fit_profile$fitted_noise_sd))
    expect_true(is.finite(as.numeric(fit_profile$log_likelihood)))
  }
})

test_that("point-exponential wrapper validates unsupported smoother arguments", {
  x <- c(0, 0.2, 1)

  expect_error(
    eb_smoother(x, s = 0.1, family = "point_exponential", link = "log"),
    "identity observation scale"
  )
  expect_error(
    eb_smoother(x, s = 0.1, family = "point_exponential", locations = seq_along(x)),
    "not used"
  )
  expect_error(
    eb_smoother(x, s = 0.1, family = "point_exponential", backend = "inla"),
    "only supports the exact backend"
  )
  expect_error(
    eb_smoother(x, s = 0.1, family = "point_exponential", beta_fixed = 0),
    "not used"
  )
})

test_that("high-level Matern fix_g holds fitted range and sigma fixed", {
  set.seed(107)

  loc <- seq(0, 1, length.out = 12)
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.08)

  fit <- eb_smoother(x, s = 0.08, family = "matern", locations = loc)
  fixed_fit <- eb_smoother(
    x,
    s = 0.08,
    family = "matern",
    locations = loc,
    g_init = fit$fitted_g,
    fix_g = TRUE
  )

  expect_equal(fixed_fit$fitted_g$theta, fit$fitted_g$theta, tolerance = 1e-12)
  expect_equal(fixed_fit$fitted_g$sigma, fit$fitted_g$sigma, tolerance = 1e-12)
})

test_that("Matern fix_params supports partial hyperparameter and beta fixing", {
  set.seed(108)

  loc <- seq(0, 1, length.out = 12)
  s <- rep(0.08, length(loc))
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)
  g_init <- Matern(theta = log(0.3), sigma = 1, beta = 0.15)

  sigma_fixed <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    g_init = g_init,
    fix_params = "sigma"
  )
  expect_equal(sigma_fixed$fitted_g$sigma, 1, tolerance = 1e-12)
  expect_true(is.finite(sigma_fixed$fitted_g$theta))
  expect_true(is.finite(sigma_fixed$fitted_beta))

  both_fixed <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    g_init = g_init,
    fix_params = c("range", "sigma")
  )
  fix_g_fit <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    g_init = g_init,
    fix_g = TRUE
  )
  expect_equal(both_fixed$fitted_g$theta, fix_g_fit$fitted_g$theta, tolerance = 1e-12)
  expect_equal(both_fixed$fitted_g$sigma, fix_g_fit$fitted_g$sigma, tolerance = 1e-12)
  expect_equal(as.numeric(both_fixed$log_likelihood), as.numeric(fix_g_fit$log_likelihood), tolerance = 1e-10)

  beta_from_argument <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    g_init = g_init,
    beta_fixed = 0.2,
    fix_params = "beta"
  )
  beta_from_g <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    g_init = g_init,
    fix_params = "beta"
  )
  expect_equal(unname(beta_from_argument$fitted_beta), 0.2, tolerance = 1e-12)
  expect_equal(unname(beta_from_g$fitted_beta), 0.15, tolerance = 1e-12)

  expect_error(
    eb_smoother(x, s = s, family = "matern", locations = loc, fix_params = "sigma"),
    "requires `g_init\\$sigma`"
  )
  expect_error(
    eb_smoother(x, s = s, family = "matern", locations = loc, fix_params = "scale"),
    "unsupported value"
  )
  expect_error(
    eb_smoother(x, s = s, family = "matern", locations = loc, g_init = g_init, fix_params = "beta", beta_prec = 0),
    "cannot be combined with `beta_prec`"
  )
  expect_error(
    eb_smoother(x, s = s, family = "matern", locations = loc, backend = "inla", g_init = g_init, fix_params = "sigma"),
    "not supported with the INLA Matern backend"
  )
})

test_that("Matern TMB path respects configurable start count", {
  set.seed(105)

  loc <- seq(0, 1, length.out = 8)
  s <- rep(0.04, length(loc))
  x <- exp(0.1 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = 0.04)

  fit_known_one <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    link = "log",
    backend = "laplace"
  )
  expect_equal(fit_known_one$raw_fit$laplace_diagnostics$stepA_n_starts, 1)

  fit_known_five <- eb_smoother(
    x,
    s = s,
    family = "matern",
    locations = loc,
    link = "log",
    backend = "laplace",
    matern_n_starts = 5
  )
  expect_equal(fit_known_five$raw_fit$laplace_diagnostics$stepA_n_starts, 5)

  fit_learned_one <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    link = "log",
    backend = "laplace"
  )
  expect_equal(fit_learned_one$raw_fit$laplace_diagnostics$stepA_n_starts, 1)

  fit_learned_five <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    link = "log",
    backend = "laplace",
    matern_n_starts = 5
  )
  expect_equal(fit_learned_five$raw_fit$laplace_diagnostics$stepA_n_starts, 5)
})

test_that("Matern backend aliases canonicalize to simpler public backends", {
  set.seed(106)

  loc <- seq(0, 1, length.out = 8)
  x <- exp(0.1 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = 0.04)

  fit_laplace <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    link = "log",
    backend = "laplace"
  )
  expect_equal(fit_laplace$backend, "laplace")
  expect_equal(fit_laplace$raw_fit$laplace_implementation, "tmb")

  fit_laplace_r <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    link = "log",
    backend = "laplace_r"
  )
  expect_equal(fit_laplace_r$backend, "laplace")
  expect_equal(fit_laplace_r$raw_fit$laplace_implementation, "r")

  fit_laplace_alias <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    link = "log",
    backend = "laplace_tmb"
  )
  expect_equal(fit_laplace_alias$backend, "laplace")
  expect_equal(fit_laplace_alias$raw_fit$laplace_implementation, "tmb")

  pc <- list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))
  x_identity <- sin(2 * pi * loc)
  fit_inla <- eb_smoother(
    x_identity,
    s = 0.1,
    family = "matern",
    locations = loc,
    backend = "inla",
    pc.penalty = pc
  )
  fit_inla_alias <- eb_smoother(
    x_identity,
    s = 0.1,
    family = "matern",
    locations = loc,
    backend = "inla_pc",
    pc.penalty = pc
  )
  expect_equal(fit_inla$backend, "inla")
  expect_equal(fit_inla_alias$backend, "inla")
})

test_that("default 2D Matern mesh includes observed locations as vertices", {
  loc <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 4),
    y = seq(0, 1, length.out = 4)
  ))
  setup <- Matern_setup(locations = loc)
  mesh_loc <- setup$mesh$loc[, seq_len(2), drop = FALSE]
  min_sq_dist <- apply(loc, 1, function(z) min(rowSums((sweep(mesh_loc, 2, z))^2)))

  expect_true(all(min_sq_dist < 1e-20))
})

test_that("Matern wrapper mirrors generator Fisher auto backend policy for log-link fits", {
  set.seed(1021)

  loc <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 4),
    y = seq(0, 1, length.out = 4)
  ))
  s <- rep(0.1, nrow(loc))
  eta <- 0.05 + 0.18 * sin(2 * pi * loc[, 1]) + 0.08 * cos(2 * pi * loc[, 2])
  x <- exp(eta) + rnorm(nrow(loc), sd = s)
  setup <- Matern_setup(locations = loc, max.edge = 0.6)

  fit_fixed <- eb_smoother(
    x,
    s = s,
    family = "matern",
    setup = setup,
    link = "log",
    beta_fixed = 0
  )
  fit_eb <- eb_smoother(
    x,
    s = s,
    family = "matern",
    setup = setup,
    link = "log"
  )
  fit_learned_fixed <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    setup = setup,
    link = "log",
    beta_fixed = 0
  )

  expect_equal(fit_fixed$backend, "fisher_pql")
  expect_equal(fit_fixed$raw_fit$backend, "fisher_pql")
  expect_equal(fit_fixed$raw_fit$laplace_curvature, "fisher")
  expect_equal(fit_eb$backend, "fisher_pql")
  expect_equal(fit_eb$raw_fit$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit_eb$raw_fit$laplace_curvature, "fisher")
  expect_equal(fit_learned_fixed$backend, "fisher_pql")
  expect_equal(fit_learned_fixed$raw_fit$laplace_implementation, "exact_fisher_pql")
  expect_equal(fit_learned_fixed$raw_fit$laplace_curvature, "fisher")
})

test_that("constant and Matern wrappers preserve resolved g_init information", {
  set.seed(103)

  loc <- seq(0, 1, length.out = 20)
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)

  fit_matern <- eb_smoother(
    x,
    s = 0.1,
    family = "matern",
    locations = loc,
    g_init = Matern(theta = log(0.3), sigma = 0.2, beta = 0.1, beta_prec = 0)
  )
  expect_s3_class(fit_matern$raw_fit$g_init, "Matern")
  expect_equal(unname(fit_matern$raw_fit$g_init$beta), 0.1, tolerance = 1e-10)

  fit_constant <- eb_smoother(
    x,
    s = NULL,
    family = "constant",
    g_init = Constant(beta = 0.25, beta_prec = 0)
  )
  expect_s3_class(fit_constant$raw_fit$g_init, "Constant")
  expect_equal(unname(fit_constant$raw_fit$g_init$beta), 0.25, tolerance = 1e-10)
})

test_that("print and summary methods reflect the new public beta mode names", {
  set.seed(104)

  loc <- seq(0, 1, length.out = 20)
  x <- 0.25 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)
  fit <- eb_smoother(x, s = 0.1, family = "matern", locations = loc)

  expect_output(print(fit), "beta mode: empirical_bayes")
  expect_output(print(summary(fit)), "beta mode: empirical_bayes")
})
