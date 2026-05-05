test_that("Matern generator supports log-link positive smoothing", {
  set.seed(11)

  loc <- seq(0, 1, length.out = 8)
  s <- rep(0.08, length(loc))
  x <- exp(0.2 + 0.25 * sin(2 * pi * loc)) + rnorm(length(loc), sd = s)

  fit <- ebnm_Matern_generator(locations = loc, link = "log")(x, s)

  expect_equal(fit$backend, "laplace_fisher")
  expect_equal(fit$laplace_implementation, "tmb")
  expect_equal(fit$laplace_curvature, "fisher")
  expect_equal(fit$log_likelihood_semantics, "laplace_fisher_empirical_bayes")
  expect_equal(fit$link, "log")
  expect_true(all(fit$posterior$mean > 0))
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
})

test_that("Matern auto backend uses Fisher Laplace for log-link fits", {
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

  expect_equal(fit_auto_fixed$backend, "laplace_fisher")
  expect_equal(fit_auto_fixed$laplace_implementation, "tmb")
  expect_equal(fit_auto_fixed$laplace_curvature, "fisher")
  expect_equal(fit_auto_fixed$log_likelihood_semantics, "laplace_fisher_fixed")
  expect_equal(fit_laplace_fixed$backend, "laplace")
  expect_equal(fit_laplace_fixed$laplace_implementation, "tmb")
  expect_equal(fit_laplace_fixed$laplace_curvature, "observed")
  expect_true(is.finite(as.numeric(fit_auto_fixed$log_likelihood)))
  expect_true(is.finite(as.numeric(fit_laplace_fixed$log_likelihood)))
  expect_lt(max(abs(fit_auto_fixed$posterior$mean - fit_laplace_fixed$posterior$mean)), 0.05)

  expect_equal(fit_auto_eb$backend, "laplace_fisher")
  expect_equal(fit_auto_eb$laplace_implementation, "tmb")
  expect_equal(fit_auto_eb$laplace_curvature, "fisher")
  expect_equal(fit_auto_fixed_g$backend, "laplace_fisher")
  expect_equal(fit_auto_fixed_g$laplace_implementation, "tmb")
  expect_equal(fit_auto_fixed_g$laplace_curvature, "fisher")
})

test_that("Matern softplus auto backend uses observed Laplace; INLA unsupported", {
  set.seed(211)
  loc <- seq(0, 1, length.out = 12)
  s <- rep(0.1, length(loc))
  eta <- seq(-6, 6, length.out = length(loc))
  x <- log1p(exp(eta)) + rnorm(length(loc), sd = s)

  fit_auto <- ebnm_Matern_generator(locations = loc, link = "softplus")(x, s)
  fit_fisher <- ebnm_Matern_generator(locations = loc, link = "softplus", backend = "laplace_fisher")(x, s)
  expect_equal(fit_auto$backend, "laplace")
  expect_equal(fit_auto$laplace_curvature, "observed")
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

  for (fit in list(fit_eb, fit_fixed, fit_flat, fit_proper, learned_eb, learned_fixed, learned_flat, learned_proper)) {
    expect_equal(fit$backend, "laplace_fisher")
    expect_equal(fit$laplace_curvature, "fisher")
    expect_match(fit$log_likelihood_semantics, "^laplace_fisher_")
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

  expect_equal(fit$backend, "laplace_fisher")
  expect_equal(fit$raw_fit$laplace_implementation, "tmb")
  expect_equal(fit$raw_fit$laplace_curvature, "fisher")
  expect_equal(fit$log_likelihood_semantics, "laplace_fisher_empirical_bayes")
  expect_equal(fit$link, "log")
  expect_true(all(fit$posterior$mean > 0))

  fit_learned <- eb_smoother(
    x,
    s = NULL,
    family = "matern",
    locations = loc,
    link = "log"
  )
  expect_equal(fit_learned$backend, "laplace_fisher")
  expect_equal(fit_learned$raw_fit$laplace_implementation, "tmb")
  expect_equal(fit_learned$raw_fit$laplace_curvature, "fisher")
  expect_equal(fit_learned$log_likelihood_semantics, "laplace_fisher_empirical_bayes")
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
  expect_equal(fit_pc_auto$backend, "laplace_fisher")
  expect_equal(fit_pc_auto$raw_fit$laplace_curvature, "fisher")
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
