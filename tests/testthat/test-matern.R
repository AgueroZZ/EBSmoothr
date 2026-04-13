test_that("Matern generator only supports the identity link", {
  loc <- seq(0, 1, length.out = 10)
  expect_error(
    ebnm_Matern_generator(locations = loc, link = "log"),
    "identity"
  )
})

test_that("Exact Matern fit returns finite 1D posterior summaries", {
  set.seed(3)

  loc <- seq(0, 1, length.out = 30)
  s <- rep(0.12, length(loc))
  x <- 0.25 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_fun <- ebnm_Matern_generator(locations = loc)
  fit <- fit_fun(x, s)

  expect_equal(nrow(fit$posterior), length(loc))
  expect_true(all(is.finite(fit$posterior$mean)))
  expect_true(all(is.finite(fit$posterior$var)))
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
  expect_s3_class(fit$fitted_g, "Matern")
  expect_length(fit$fitted_beta, 1)
  expect_equal(fit$beta_mode, "profile")
  expect_null(fit$inla_result)
  expect_true(nrow(fit$posterior_spatial_field) > 0)
})

test_that("Matern PC penalty parsing supports defaults and validates inputs", {
  set.seed(21)

  loc <- seq(0, 1, length.out = 20)
  s <- rep(0.12, length(loc))
  x <- sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_fun_partial <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = 0.2)
  )
  fit_partial <- fit_fun_partial(x, s)

  expect_equal(fit_partial$backend, "inla_pc")
  expect_equal(unname(fit_partial$pc_penalty$range["anchor"]), 0.2, tolerance = 1e-10)
  expect_equal(unname(fit_partial$pc_penalty$range["alpha"]), 0.5, tolerance = 1e-10)
  expect_equal(unname(fit_partial$pc_penalty$sigma["alpha"]), 0.5, tolerance = 1e-10)
  expect_true(unname(fit_partial$pc_penalty$sigma["anchor"]) > 0)

  fit_fun_bad_alpha <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = c(0.2, 1))
  )
  expect_error(
    fit_fun_bad_alpha(x, s),
    "alpha must satisfy 0 < alpha < 1"
  )

  fit_fun_bad_names <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(foo = c(0.2, 0.5))
  )
  expect_error(
    fit_fun_bad_names(x, s),
    "only supports `range` and `sigma`"
  )
})

test_that("Exact Matern marginal likelihood matches dense Gaussian evaluation", {
  build_mesh_A <- getFromNamespace(".build_mesh_A", "EBSmoothr")
  exact_matern_state <- getFromNamespace(".exact_matern_state", "EBSmoothr")

  set.seed(11)

  loc <- seq(0, 1, length.out = 12)
  meshA <- build_mesh_A(matrix(loc, ncol = 1), max.edge = 0.1)
  A <- meshA$A
  spde <- INLA::inla.spde2.matern(meshA$mesh, alpha = 2)

  x <- rnorm(length(loc))
  s <- rep(0.15, length(loc))
  par_grid <- rbind(
    c(log(0.2), log(0.4), 0.0),
    c(log(0.3), log(0.7), 0.5),
    c(log(0.12), log(0.25), -0.2)
  )

  for (k in seq_len(nrow(par_grid))) {
    st <- exact_matern_state(
      x = x,
      s = s,
      A = A,
      spde_template = spde,
      alpha = 2,
      d = 1,
      log_range = par_grid[k, 1],
      log_sigma = par_grid[k, 2],
      beta0 = par_grid[k, 3]
    )

    Sigma <- as.matrix(A) %*% solve(as.matrix(st$Q)) %*% t(as.matrix(A)) + diag(s^2)
    r <- x - par_grid[k, 3]
    dense_ll <- -0.5 * (
      length(x) * log(2 * pi) +
        as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) +
        drop(t(r) %*% solve(Sigma, r))
    )

    expect_equal(st$log_marginal, dense_ll, tolerance = 1e-8)
  }
})

test_that("Exact Matern fixed-parameter mode uses beta_fixed and fixed g_init", {
  set.seed(4)

  loc <- seq(0, 1, length.out = 30)
  s <- rep(0.15, length(loc))
  x <- 1 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_fun <- ebnm_Matern_generator(locations = loc)
  g_fix <- Matern(theta = log(0.2), sigma = 0.05)

  fit_default_beta <- fit_fun(x, s, g_init = g_fix, fix_g = TRUE)
  fit_custom_beta <- fit_fun(x, s, g_init = g_fix, fix_g = TRUE, beta_fixed = 1)

  expect_equal(unname(fit_default_beta$fitted_beta), 0, tolerance = 1e-10)
  expect_equal(unname(fit_custom_beta$fitted_beta), 1, tolerance = 1e-10)
  expect_equal(fit_default_beta$fitted_g$theta, log(0.2), tolerance = 1e-10)
  expect_equal(fit_default_beta$fitted_g$sigma, 0.05, tolerance = 1e-10)
  expect_gt(
    mean(fit_custom_beta$posterior$mean - fit_default_beta$posterior$mean),
    0.2
  )
  expect_gt(
    abs(as.numeric(fit_custom_beta$log_likelihood) - as.numeric(fit_default_beta$log_likelihood)),
    1
  )

  expect_error(
    fit_fun(x, s, beta_fixed = 1),
    "`beta_fixed` can only be supplied"
  )
})

test_that("INLA PC-prior Matern mode returns Step A diagnostics and backend metadata", {
  set.seed(22)

  loc <- seq(0, 1, length.out = 25)
  s <- rep(0.15, length(loc))
  x <- 0.3 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_fun_pc <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))
  )
  fit_pc <- fit_fun_pc(x, s)

  expect_equal(fit_pc$backend, "inla_pc")
  expect_equal(fit_pc$prior_family, "identity_Matern_pc")
  expect_equal(fit_pc$beta_mode, "integrated_flat")
  expect_false(is.null(fit_pc$inla_result))
  expect_equal(
    as.numeric(fit_pc$log_likelihood),
    fit_pc$log_likelihood_stepA_penalized,
    tolerance = 1e-8
  )
  expect_true(is.finite(fit_pc$log_likelihood_stepA_mlik_integration))
  expect_true(is.finite(fit_pc$log_likelihood_stepA_mlik_gaussian))
  expect_true(is.finite(fit_pc$log_likelihood_stepA_joint_log_posterior))
  expect_null(fit_pc$log_likelihood_exact_at_stepA_mode)
  expect_true(all(is.finite(fit_pc$posterior$mean)))
  expect_true(all(is.finite(fit_pc$posterior$var)))
  expect_warning(
    draws <- fit_pc$posterior_sampler(2),
    "not implemented"
  )
  expect_equal(dim(draws), c(2, length(loc)))
  expect_true(all(is.na(draws)))
})

test_that("INLA PC-prior Matern can optionally compute the exact Step A diagnostic", {
  set.seed(2201)

  loc <- seq(0, 1, length.out = 20)
  s <- rep(0.15, length(loc))
  x <- 0.3 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_fun_pc <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = c(0.2, 0.5), sigma = c(0.3, 0.5)),
    compute_exact_diagnostic = TRUE
  )
  fit_pc <- fit_fun_pc(x, s)

  expect_true(is.finite(fit_pc$log_likelihood_exact_at_stepA_mode))
})

test_that("INLA PC-prior fixed mode uses penalized pseudo-objective", {
  set.seed(23)

  loc <- seq(0, 1, length.out = 25)
  s <- rep(0.14, length(loc))
  x <- 1 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_fun_pc <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))
  )

  fit_pc_fixed <- fit_fun_pc(
    x,
    s,
    g_init = Matern(theta = log(0.2), sigma = 0.3),
    fix_g = TRUE,
    beta_fixed = 1
  )

  expect_equal(fit_pc_fixed$backend, "inla_pc")
  expect_equal(unname(fit_pc_fixed$fitted_beta), 1, tolerance = 1e-10)
  expect_equal(fit_pc_fixed$beta_mode, "fixed")
  expect_equal(
    as.numeric(fit_pc_fixed$log_likelihood),
    fit_pc_fixed$log_likelihood_stepB + fit_pc_fixed$log_likelihood_pc_prior_theta,
    tolerance = 1e-8
  )
  expect_null(fit_pc_fixed$log_likelihood_stepA_penalized)
  expect_true(is.finite(fit_pc_fixed$log_likelihood_stepB))
  expect_true(is.finite(fit_pc_fixed$log_likelihood_pc_prior_theta))
  expect_warning(
    draws <- fit_pc_fixed$posterior_sampler(3),
    "not implemented"
  )
  expect_equal(dim(draws), c(3, length(loc)))
  expect_true(all(is.na(draws)))
})

test_that("INLA Step B conditional marginal likelihood is prior-invariant when theta is fixed", {
  set.seed(24)

  loc <- seq(0, 1, length.out = 20)
  s <- rep(0.15, length(loc))
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)
  g_fix <- Matern(theta = log(0.2), sigma = 0.3)

  fit_fun_balanced <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))
  )
  fit_fun_smooth <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = c(0.2, 0.1), sigma = c(0.3, 0.1))
  )

  fit_balanced <- fit_fun_balanced(x, s, g_init = g_fix, fix_g = TRUE, beta_fixed = 0.2)
  fit_smooth <- fit_fun_smooth(x, s, g_init = g_fix, fix_g = TRUE, beta_fixed = 0.2)

  expect_equal(
    fit_balanced$log_likelihood_stepB,
    fit_smooth$log_likelihood_stepB,
    tolerance = 1e-8
  )
  expect_gt(
    abs(as.numeric(fit_balanced$log_likelihood) - as.numeric(fit_smooth$log_likelihood)),
    1e-4
  )
})

test_that("Exact and weak-PC Matern fits remain close on a small two-dimensional example", {
  set.seed(25)

  grid <- expand.grid(
    x = seq(0, 1, length.out = 5),
    y = seq(0, 1, length.out = 4)
  )
  loc <- as.matrix(grid)
  s <- rep(0.18, nrow(loc))
  x <- with(grid, sin(2 * pi * x) * cos(2 * pi * y)) + rnorm(nrow(loc), sd = s)

  fit_exact <- ebnm_Matern_generator(locations = loc, max.edge = c(0.3, 0.45))(x, s)
  fit_pc <- ebnm_Matern_generator(
    locations = loc,
    max.edge = c(0.3, 0.45),
    pc.penalty = list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))
  )(x, s)

  expect_gt(stats::cor(fit_exact$posterior$mean, fit_pc$posterior$mean), 0.99)
  expect_lt(sqrt(mean((fit_exact$posterior$mean - fit_pc$posterior$mean)^2)), 0.05)
})

test_that("matern_objective_breakdown matches recorded exact and fixed-PC objectives", {
  set.seed(2501)

  loc <- seq(0, 1, length.out = 25)
  s <- rep(0.15, length(loc))
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_exact <- ebnm_Matern_generator(locations = loc)(x, s)
  bd_exact <- matern_objective_breakdown(fit_exact)

  expect_equal(bd_exact$matched_beta_mode, "profile")
  expect_equal(
    bd_exact$recorded_log_likelihood,
    bd_exact$matched_exact_objective,
    tolerance = 1e-8
  )
  expect_true(is.null(bd_exact$log_pc_prior_theta))

  fit_pc_fixed <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))
  )(
    x,
    s,
    g_init = Matern(theta = log(0.2), sigma = 0.3),
    fix_g = TRUE,
    beta_fixed = 0.2
  )
  bd_pc_fixed <- matern_objective_breakdown(fit_pc_fixed)

  expect_equal(bd_pc_fixed$matched_beta_mode, "fixed")
  expect_equal(
    bd_pc_fixed$objective_fixed_plus_prior,
    bd_pc_fixed$matched_exact_objective,
    tolerance = 1e-8
  )
  expect_true(is.finite(bd_pc_fixed$log_pc_prior_theta))
})

test_that("matern_objective_breakdown exposes exact integrated-flat comparator for Step A PC fits", {
  set.seed(2502)

  loc <- seq(0, 1, length.out = 25)
  s <- rep(0.15, length(loc))
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_pc <- ebnm_Matern_generator(
    locations = loc,
    pc.penalty = list(range = c(0.2, 0.5), sigma = c(0.3, 0.5))
  )(x, s)

  bd_pc <- matern_objective_breakdown(fit_pc)

  expect_equal(bd_pc$matched_beta_mode, "integrated_flat")
  expect_true(is.finite(bd_pc$objective_integrated_flat_plus_prior))
  expect_true(is.finite(bd_pc$recorded_minus_matched_exact))
})

test_that("Exact Matern optimization improves over the supplied initialization", {
  set.seed(5)

  loc <- seq(0, 1, length.out = 35)
  s <- rep(0.12, length(loc))
  x <- 0.5 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

  fit_fun <- ebnm_Matern_generator(locations = loc)
  g_init <- Matern(theta = log(0.05), sigma = 0.2)
  beta0_init <- stats::weighted.mean(x, 1 / (s^2))

  fit_init <- fit_fun(
    x,
    s,
    g_init = g_init,
    fix_g = TRUE,
    beta_fixed = beta0_init
  )
  fit_opt <- fit_fun(x, s, g_init = g_init)

  expect_gte(
    as.numeric(fit_opt$log_likelihood),
    as.numeric(fit_init$log_likelihood) - 1e-6
  )
})

test_that("Exact Matern supports two-dimensional locations", {
  set.seed(6)

  grid <- expand.grid(
    x = seq(0, 1, length.out = 6),
    y = seq(0, 1, length.out = 5)
  )
  loc <- as.matrix(grid)
  s <- rep(0.15, nrow(loc))
  x <- with(grid, sin(2 * pi * x) * cos(2 * pi * y)) + rnorm(nrow(loc), sd = s)

  fit_fun <- ebnm_Matern_generator(locations = loc, max.edge = c(0.25, 0.4))
  fit <- fit_fun(x, s)

  expect_equal(nrow(fit$posterior), nrow(loc))
  expect_true(all(is.finite(fit$posterior$mean)))
  expect_true(all(is.finite(fit$posterior$var)))
  expect_true(nrow(fit$posterior_spatial_field) > 0)
})
