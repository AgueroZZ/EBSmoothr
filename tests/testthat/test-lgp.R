test_that("LGP fix_g uses fixed beta coefficients", {
  set.seed(1)

  t <- seq(0, 1, length.out = 40)
  s <- rep(0.1, length(t))
  x <- 0.5 + sin(2 * pi * t) + rnorm(length(t), sd = s)

  setup_diffuse <- LGP_setup(t = t, betaprec = 0)
  fit_fun_diffuse <- ebnm_LGP_generator(setup_diffuse)

  fit_zero <- fit_fun_diffuse(
    x,
    s,
    g_init = LGP(0),
    fix_g = TRUE
  )
  expect_equal(unname(fit_zero$fitted_beta), c(0, 0), tolerance = 1e-10)

  fit_fixed <- fit_fun_diffuse(
    x,
    s,
    g_init = LGP(0),
    fix_g = TRUE,
    beta_fixed = c(0.5, -0.25)
  )
  expect_equal(unname(fit_fixed$fitted_beta), c(0.5, -0.25), tolerance = 1e-10)

  setup_eb_beta <- LGP_setup(t = t, betaprec = -1)
  fit_fun_eb_beta <- ebnm_LGP_generator(setup_eb_beta)
  fit_fixed_eb <- fit_fun_eb_beta(
    x,
    s,
    g_init = LGP(0),
    fix_g = TRUE,
    beta_fixed = c(0.5, -0.25)
  )
  expect_equal(unname(fit_fixed_eb$fitted_beta), c(0.5, -0.25), tolerance = 1e-10)
})

test_that("LGP beta_fixed validation and log-likelihood diagnostics behave as expected", {
  set.seed(2)

  t <- seq(0, 1, length.out = 60)
  s <- rep(0.1, length(t))
  x <- sin(2 * pi * t) + rnorm(length(t), sd = s)

  setup <- LGP_setup(t = t, betaprec = 0)
  fit_fun <- ebnm_LGP_generator(setup)

  expect_error(
    fit_fun(x, s, beta_fixed = c(0, 0)),
    "`beta_fixed` can only be supplied"
  )

  expect_error(
    fit_fun(x, s, g_init = LGP(0), fix_g = TRUE, beta_fixed = 1),
    "length ncol\\(LGP_setup\\$X\\)"
  )

  fit <- fit_fun(x, s)
  expect_lt(
    abs(as.numeric(fit$log_likelihood) - as.numeric(fit$log_likelihood_stepB_laplace)),
    1e-3
  )
})
