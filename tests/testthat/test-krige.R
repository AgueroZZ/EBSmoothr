set.seed(101)
t_obs <- sort(runif(25))
x_obs <- 2 + sin(2 * pi * t_obs) + rnorm(25, 0, 0.2)
s_obs <- runif(25, 0.1, 0.3)

expect_reproduces_fit <- function(fit, setup, tol_mean, tol_var) {
  for (m in c("analytical", "fem")) {
    k <- krige_GP(fit, setup, t_obs, method = m)
    expect_lt(max(abs(k$mean - fit$posterior$mean)), tol_mean)
    expect_lt(max(abs(k$sd^2 - fit$posterior$var)) / max(fit$posterior$var), tol_var)
  }
}

test_that("krige_GP reproduces L-GP fits at the observed locations, for every beta mode", {
  for (bp in c(-1, 0, 1)) {
    st <- LGP_setup(t_obs, num_knots = 12, betaprec = bp)
    fit <- ebnm_LGP_generator(st)(x_obs, s_obs)
    expect_reproduces_fit(fit, st, tol_mean = 1e-4, tol_var = 1e-4)
  }
})

test_that("krige_GP reproduces Matern fits at the observed locations, for every beta mode", {
  for (alpha in c(2, 3)) {
    st <- Matern_setup(t_obs, alpha = alpha)
    gen <- ebnm_Matern_generator(setup = st)
    expect_reproduces_fit(gen(x_obs, s_obs), st, tol_mean = 1e-7, tol_var = 1e-6)
    expect_reproduces_fit(gen(x_obs, s_obs, beta_prec = 0), st, tol_mean = 1e-7, tol_var = 1e-6)
    expect_reproduces_fit(gen(x_obs, s_obs, beta_prec = 1), st, tol_mean = 1e-7, tol_var = 1e-6)
  }
})

test_that("L-GP analytical and fem predictions converge as the knots get denser", {
  tn <- seq(min(t_obs), max(t_obs), length.out = 120)
  gap <- sapply(c(10, 200), function(nk) {
    st <- LGP_setup(t_obs, num_knots = nk, betaprec = -1)
    fit <- ebnm_LGP_generator(st)(x_obs, s_obs)
    a <- krige_GP(fit, st, tn, "analytical")
    b <- krige_GP(fit, st, tn, "fem")
    c(mean = max(abs(a$mean - b$mean)), sd = max(abs(a$sd - b$sd)))
  })
  expect_lt(gap["mean", 2], gap["mean", 1] / 20)
  expect_lt(gap["mean", 2], 1e-3)
  expect_lt(gap["sd", 2], 5e-3)
})

test_that("the IWP kernel matches the L-GP basis covariance for dense knots", {
  for (kn in list(seq(0.1, 1, length.out = 601), seq(-1, 0, length.out = 601), seq(-1, 1, length.out = 601))) {
    for (p in 2:3) {
      tt <- c(-0.8, -0.3, 0.3, 0.7)
      tt <- tt[tt > min(kn) & tt < max(kn)]
      info <- list(type = "LGP", knots = kn, p = p)
      B <- local_poly_helper(kn, tt, p)
      Pd <- diag(compute_weights_precision_helper(kn))
      theta <- 0.7
      basis_cov <- B %*% diag(1 / (exp(theta) * Pd)) %*% t(B)
      kern <- .krige_lgp_cov(matrix(tt), matrix(tt), list(scale = theta), info)
      expect_lt(max(abs(basis_cov - kern)) / max(abs(kern)), 1e-4)
      expect_equal(.krige_cov_diag(matrix(tt), list(scale = theta), info), diag(kern), tolerance = 1e-10)
    }
  }
})

test_that("the Matern kernel matches its closed forms", {
  h <- c(0, 0.05, 0.3, 1.2)
  a <- matrix(0, 1, 1)
  b <- matrix(h, ncol = 1)
  g <- list(theta = log(0.4), sigma = 1.7)
  for (alpha in 1:3) {
    nu <- alpha - 1 / 2
    kappa <- sqrt(8 * nu) / 0.4
    z <- kappa * h
    closed <- 1.7^2 * switch(alpha, exp(-z), (1 + z) * exp(-z), (1 + z + z^2 / 3) * exp(-z))
    expect_equal(as.numeric(.krige_matern_cov(a, b, g, list(alpha = alpha, d = 1L))), closed, tolerance = 1e-10)
  }
})

test_that("krige_GP returns a consistent summary and full covariance", {
  st <- Matern_setup(t_obs, alpha = 2)
  fit <- ebnm_Matern_generator(setup = st)(x_obs, s_obs)
  tn <- c(-0.2, 0.33, 0.5, 1.3)
  k <- krige_GP(fit, st, tn, return_cov = TRUE, level = 0.9, boundary = "stationary")
  expect_equal(k$location, tn)
  expect_equal(sqrt(diag(attr(k, "cov"))), k$sd, tolerance = 1e-10)
  expect_equal(k$upper - k$mean, stats::qnorm(0.95) * k$sd, tolerance = 1e-10)
  expect_warning(krige_GP(fit, st, tn, method = "fem"), "outside the Matern mesh")
  expect_error(krige_GP(list(x = x_obs, s = s_obs, g = structure(list(pi = 1, mean = 0, sd = 1), class = "normalmix")), st, tn), "EBSmoothr")
  st_lgp <- LGP_setup(t_obs, num_knots = 12, betaprec = -1)
  attr(st_lgp, "lgp_design") <- NULL
  expect_error(krige_GP(ebnm_LGP_generator(LGP_setup(t_obs, num_knots = 12, betaprec = -1))(x_obs, s_obs),
                        st_lgp, tn), "0.3.2")
})

test_that("krige_flash and flash_ebnm_data agree with flashier and with krige_GP", {
  skip_if_not_installed("flashier")
  set.seed(7)
  tt <- seq(0, 1, length.out = 30)
  L0 <- cbind(sin(2 * pi * tt), exp(-(tt - 0.5)^2 / 0.02))
  Y <- L0 %*% t(matrix(rnorm(160), 80, 2)) + matrix(rnorm(30 * 80, 0, 0.3), 30, 80)
  tn <- seq(0, 1, length.out = 101)
  for (pr in c("Matern", "LGP")) {
    st <- if (pr == "Matern") Matern_setup(tt, alpha = 2) else LGP_setup(tt, num_knots = 30, betaprec = -1)
    gen <- if (pr == "Matern") ebnm_Matern_generator(setup = st) else ebnm_LGP_generator(st)
    fl <- suppressWarnings(flashier::flash(Y, ebnm_fn = list(gen, ebnm::ebnm_point_normal), greedy_Kmax = 3,
                                           var_type = 0, backfit = TRUE, verbose = 0))
    expect_gte(fl$n_factors, 1)
    expect_equal(krige_flash(fl, st, tt), fl$L_pm, tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(krige_flash(fl, st, tt, method = "fem"), fl$L_pm, tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(abs(krige_flash(fl, st, tt, normalize = "ldf")),
                 abs(flashier::ldf(fl, type = "2")$L), tolerance = 1e-8, ignore_attr = TRUE)
    grid_means <- krige_flash(fl, st, tn)
    for (k in seq_len(fl$n_factors)) {
      d <- flash_ebnm_data(fl, k)
      expect_lt(max(abs(krige_GP(d, st, tt)$mean - fl$L_pm[, k])) / max(abs(fl$L_pm[, k])), 1e-3)
      expect_lt(max(abs(krige_GP(d, st, tn)$mean - grid_means[, k])) / max(abs(grid_means[, k])), 1e-3)
    }
  }
})

test_that("Neumann Matern: images and cosine series agree, and the diagonal matches", {
  for (nu in c(0.5, 1.5, 2.5)) {
    kappa <- 1.3; L <- 1
    s1 <- c(0, 0.2, 0.55, 1); t1 <- c(0.1, 0.5, 0.9)
    D1 <- outer(s1, t1, "-"); D2 <- outer(s1, t1, "+")
    im <- .krige_neumann_corr(D1, D2, nu, kappa, L, method = "images")
    co <- .krige_neumann_corr(D1, D2, nu, kappa, L, method = "cosine")
    expect_lt(max(abs(im - co)) / max(abs(im)), if (nu == 0.5) 1e-3 else 1e-6)
  }
  st <- Matern_setup(t_obs, alpha = 2)
  info <- .krige_resolve_boundary(NULL, .krige_setup_info(st, "Matern"))
  expect_identical(info$boundary, "neumann")
  g <- list(theta = log(0.3), sigma = 1.4)
  loc <- matrix(c(0.02, 0.4, 0.97))
  expect_equal(.krige_cov_diag(loc, g, info), diag(.krige_cov(loc, loc, g, info)), tolerance = 1e-10)
  # far from the ends the Neumann covariance is the stationary one
  mid <- matrix(c(0.45, 0.5))
  info_s <- .krige_resolve_boundary("stationary", info)
  g_short <- list(theta = log(0.05), sigma = 1.4)
  expect_equal(.krige_cov(mid, mid, g_short, info), .krige_cov(mid, mid, g_short, info_s), tolerance = 1e-8)
})

test_that("Matern fem on a fine mesh converges to the Neumann (unextended) or stationary (extended) kriging", {
  st0 <- Matern_setup(t_obs, alpha = 2)
  fit <- ebnm_Matern_generator(setup = st0)(x_obs, s_obs)
  prob <- list(x = x_obs, s = s_obs, g = fit$fitted_g)
  fine_setup <- function(h, E) {
    st <- st0
    nodes <- sort(unique(round(c(t_obs, seq(min(t_obs) - E, max(t_obs) + E, by = h)), 10)))
    st$mesh <- fmesher::fm_mesh_1d(loc = nodes)
    st$A <- fmesher::fm_basis(st$mesh, loc = t_obs)
    st$fem <- .matern_fem_build(mesh = st$mesh, alpha = 2, d = 1L)
    st$spde_template <- st$fem
    st
  }
  tn <- seq(min(t_obs), max(t_obs), length.out = 150)
  gap <- function(st, boundary) {
    a <- krige_GP(prob, st, tn, "analytical", boundary = boundary)
    f <- krige_GP(prob, st, tn, "fem")
    c(mean = max(abs(a$mean - f$mean)) / max(abs(a$mean)), sd = max(abs(a$sd - f$sd)) / max(a$sd))
  }
  g_neu <- sapply(c(1 / 60, 1 / 480), function(h) gap(fine_setup(h, 0), "neumann"))
  g_sta <- sapply(c(1 / 60, 1 / 480), function(h) gap(fine_setup(h, 1), "stationary"))
  expect_lt(g_neu["mean", 2], g_neu["mean", 1] / 10)
  expect_lt(g_neu["mean", 2], 1e-3)
  expect_lt(g_neu["sd", 2], 5e-3)
  expect_lt(g_sta["mean", 2], g_sta["mean", 1] / 10)
  expect_lt(g_sta["mean", 2], 1e-3)
  expect_lt(g_sta["sd", 2], 5e-3)
  # at the observed points the boundary choice does not matter
  for (b in c("neumann", "stationary")) {
    k <- krige_GP(fit, st0, t_obs, boundary = b)
    expect_lt(max(abs(k$mean - fit$posterior$mean)), 1e-7)
  }
  expect_identical(attr(krige_GP(fit, st0, 0.5), "boundary"), "neumann")
  expect_warning(krige_GP(fit, st0, 1.5), "outside the Matern mesh")
  set.seed(3)
  loc2 <- cbind(runif(20), runif(20))
  st2 <- Matern_setup(loc2, alpha = 2)
  fit2 <- ebnm_Matern_generator(setup = st2)(rnorm(20), rep(0.5, 20))
  expect_identical(attr(krige_GP(fit2, st2, cbind(0.5, 0.5)), "boundary"), "stationary")
  expect_error(krige_GP(fit2, st2, cbind(0.5, 0.5), boundary = "neumann"), "only available for 1-d")
})
