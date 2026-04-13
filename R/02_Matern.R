# ---- helper: normalize locations and build mesh/A for d = 1 or 2 ----
.normalize_locations <- function(locations) {
  if (is.data.frame(locations)) locations <- as.matrix(locations)

  if (is.null(dim(locations))) {
    locations <- matrix(as.numeric(locations), ncol = 1)
  } else {
    locations <- as.matrix(locations)
    storage.mode(locations) <- "double"
  }

  if (!ncol(locations) %in% c(1L, 2L)) {
    stop("For Matern, `locations` must be a vector, an n x 1 matrix, or an n x 2 matrix.")
  }

  list(loc = locations, d = ncol(locations))
}

.default_penalty_range <- function(loc_mat) {
  if (ncol(loc_mat) == 1L) {
    rr <- range(loc_mat[, 1], finite = TRUE)
    out <- diff(rr) / 10
  } else {
    ranges <- apply(loc_mat, 2, function(z) diff(range(z, finite = TRUE)))
    out <- min(ranges) / 10
  }
  if (!is.finite(out) || out <= 0) out <- 1
  out
}

.build_mesh_A <- function(loc_mat, max.edge = NULL) {
  d <- ncol(loc_mat)

  if (d == 1L) {
    loc1 <- as.numeric(loc_mat[, 1])
    if (is.null(max.edge)) {
      rr <- range(loc1, finite = TRUE)
      max.edge <- diff(rr) / 10
      if (!is.finite(max.edge) || max.edge <= 0) max.edge <- 1
    }
    mesh <- INLA::inla.mesh.1d(loc = loc1, max.edge = max.edge)
    A <- INLA::inla.spde.make.A(mesh = mesh, loc = loc1)
    return(list(mesh = mesh, A = Matrix::Matrix(A, sparse = TRUE)))
  }

  if (is.null(max.edge)) {
    ranges <- apply(loc_mat, 2, function(z) diff(range(z, finite = TRUE)))
    max.edge <- min(ranges) / 10
    if (!is.finite(max.edge) || max.edge <= 0) max.edge <- 1
  }
  if (length(max.edge) == 1L) max.edge <- c(max.edge, max.edge)
  if (length(max.edge) != 2L) stop("For d = 2, max.edge must be NULL, length-1, or length-2.")

  mesh <- INLA::inla.mesh.2d(loc = loc_mat, max.edge = max.edge)
  A <- INLA::inla.spde.make.A(mesh = mesh, loc = loc_mat)
  list(mesh = mesh, A = Matrix::Matrix(A, sparse = TRUE))
}

.check_single_numeric <- function(z, nm) {
  if (!is.numeric(z) || length(z) != 1L || is.na(z)) {
    stop(nm, " must be a single non-NA numeric.")
  }
  as.numeric(z)
}

.matern_tau_from_range_sigma <- function(range, sigma, alpha, d) {
  nu <- alpha - d / 2
  if (!is.finite(nu) || nu <= 0) {
    stop("`alpha` must satisfy alpha > d / 2 for the Matern SPDE parameterization.")
  }

  kappa <- sqrt(8 * nu) / range
  tau <- sqrt(
    gamma(nu) /
      (gamma(alpha) * (4 * pi)^(d / 2) * kappa^(2 * nu) * sigma^2)
  )

  list(kappa = kappa, tau = tau, nu = nu)
}

.exact_matern_sufficient_stats <- function(x, s, A, spde_template, alpha, d, log_range, log_sigma) {
  if (any(!is.finite(c(log_range, log_sigma)))) {
    stop("All Matern hyperparameters must be finite.")
  }

  range <- exp(log_range)
  sigma <- exp(log_sigma)

  if (!is.finite(range) || range <= 0) stop("Matern range must be positive.")
  if (!is.finite(sigma) || sigma <= 0) stop("Matern sigma must be positive.")

  ts <- .matern_tau_from_range_sigma(range = range, sigma = sigma, alpha = alpha, d = d)
  theta_spde <- c(log(ts$tau), log(ts$kappa))

  Q <- INLA::inla.spde2.precision(spde_template, theta = theta_spde)
  Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))

  w_prec_diag <- 1 / (s^2)
  W <- Matrix::Diagonal(x = w_prec_diag)
  AtW <- Matrix::t(A) %*% W

  Q_post <- Q + AtW %*% A
  Q_post <- Matrix::forceSymmetric(Matrix::Matrix(Q_post, sparse = TRUE))
  chol_post <- Matrix::Cholesky(Q_post, LDL = FALSE, perm = FALSE)

  solve_post <- function(rhs) as.numeric(Matrix::solve(chol_post, rhs, system = "A"))

  x <- as.numeric(x)
  u <- rep(1, length(x))
  c_x <- as.numeric(AtW %*% x)
  c_u <- as.numeric(AtW %*% u)

  quad_x <- sum(w_prec_diag * x^2) - sum(c_x * solve_post(c_x))
  quad_u <- sum(w_prec_diag * u^2) - sum(c_u * solve_post(c_u))
  cross_ux <- sum(w_prec_diag * x * u) - sum(c_u * solve_post(c_x))

  if (!is.finite(quad_u) || quad_u <= 0) {
    stop("The generalized intercept precision is not positive.")
  }

  beta_profile_hat <- cross_ux / quad_u
  quad_profile <- quad_x - cross_ux^2 / quad_u

  logdet_D <- sum(log(s^2))
  logdet_Q <- .compute_logdet_spd(Q)
  logdet_Q_post <- .compute_logdet_spd(Q_post)
  logdet_Sigma <- logdet_D - logdet_Q + logdet_Q_post

  list(
    range = range,
    sigma = sigma,
    kappa = ts$kappa,
    tau = ts$tau,
    Q = Q,
    Q_post = Q_post,
    logdet_Sigma = as.numeric(logdet_Sigma),
    quad_x = as.numeric(quad_x),
    quad_u = as.numeric(quad_u),
    cross_ux = as.numeric(cross_ux),
    beta_profile_hat = as.numeric(beta_profile_hat),
    quad_profile = as.numeric(quad_profile)
  )
}

.exact_matern_loglik_fixed_beta <- function(stats, beta0, n_obs) {
  beta0 <- .check_single_numeric(beta0, "beta0")
  quad_fixed <- stats$quad_x - 2 * beta0 * stats$cross_ux + beta0^2 * stats$quad_u
  as.numeric(-0.5 * (n_obs * log(2 * pi) + stats$logdet_Sigma + quad_fixed))
}

.exact_matern_loglik_profile_beta <- function(stats, n_obs) {
  as.numeric(-0.5 * (n_obs * log(2 * pi) + stats$logdet_Sigma + stats$quad_profile))
}

.exact_matern_loglik_integrated_flat_beta <- function(stats, n_obs) {
  as.numeric(-0.5 * (
    (n_obs - 1) * log(2 * pi) +
      stats$logdet_Sigma +
      log(stats$quad_u) +
      stats$quad_profile
  ))
}

.exact_matern_state <- function(x, s, A, spde_template, alpha, d, log_range, log_sigma, beta0) {
  if (any(!is.finite(c(log_range, log_sigma, beta0)))) {
    stop("All Matern optimization parameters must be finite.")
  }
  stats <- .exact_matern_sufficient_stats(
    x = x,
    s = s,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = log_range,
    log_sigma = log_sigma
  )

  w_prec_diag <- 1 / (s^2)
  W <- Matrix::Diagonal(x = w_prec_diag)
  AtW <- Matrix::t(A) %*% W
  residual <- as.numeric(x - beta0)
  b <- as.numeric(AtW %*% residual)
  m_w <- as.numeric(Matrix::solve(Matrix::Cholesky(stats$Q_post, LDL = FALSE, perm = FALSE), b, system = "A"))
  log_marginal <- .exact_matern_loglik_fixed_beta(stats, beta0 = beta0, n_obs = length(x))

  obs_mean <- as.numeric(beta0 + A %*% m_w)
  obs_var <- .compute_diag_A_Qinv_At(A, stats$Q_post)

  latent_var <- .compute_diag_A_Qinv_At(Matrix::Diagonal(n = nrow(stats$Q_post)), stats$Q_post)
  latent_sd <- sqrt(pmax(latent_var, 0))
  obs_sd <- sqrt(pmax(obs_var, 0))

  q025 <- stats::qnorm(0.025)
  q975 <- stats::qnorm(0.975)
  posterior_spatial_field <- data.frame(
    ID = seq_along(m_w),
    mean = m_w,
    sd = latent_sd,
    `0.025quant` = m_w + q025 * latent_sd,
    `0.5quant` = m_w,
    `0.975quant` = m_w + q975 * latent_sd,
    mode = m_w,
    kld = rep(0, length(m_w)),
    var = latent_var,
    check.names = FALSE
  )

  list(
    log_marginal = as.numeric(log_marginal),
    fitted_g = Matern(theta = log_range, sigma = stats$sigma),
    fitted_beta = as.numeric(beta0),
    Q = stats$Q,
    Q_post = stats$Q_post,
    post_mean_latent = m_w,
    post_var_latent = latent_var,
    posterior = data.frame(
      mean = obs_mean,
      var = obs_var,
      second_moment = obs_mean^2 + obs_var
    ),
    posterior_spatial_field = posterior_spatial_field,
    range = stats$range,
    sigma = stats$sigma,
    kappa = stats$kappa,
    tau = stats$tau
  )
}

.default_pc_alpha <- function() 0.5

.resolve_pc_penalty_component <- function(component, default_anchor, nm) {
  if (is.null(component)) {
    anchor <- default_anchor
    alpha <- .default_pc_alpha()
  } else {
    if (!is.numeric(component) || anyNA(component) || !length(component) %in% c(1L, 2L)) {
      stop(nm, " must be a numeric vector of length 1 or 2.")
    }
    anchor <- as.numeric(component[1])
    alpha <- if (length(component) == 1L) .default_pc_alpha() else as.numeric(component[2])
  }

  if (!is.finite(anchor) || anchor <= 0) {
    stop(nm, " anchor must be a single positive finite number.")
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop(nm, " alpha must satisfy 0 < alpha < 1.")
  }

  c(anchor = anchor, alpha = alpha)
}

.resolve_matern_pc_penalty <- function(pc.penalty, penalty_range0, sigma_anchor0) {
  if (is.null(pc.penalty)) return(NULL)

  if (!is.list(pc.penalty)) {
    stop("`pc.penalty` must be NULL or a list with optional `range` and `sigma` entries.")
  }

  nms <- names(pc.penalty)
  if (is.null(nms)) nms <- rep("", length(pc.penalty))
  if (any(nms == "")) {
    stop("`pc.penalty` entries must be named `range` and/or `sigma`.")
  }
  if (!all(nms %in% c("range", "sigma"))) {
    stop("`pc.penalty` only supports `range` and `sigma` entries.")
  }

  list(
    range = .resolve_pc_penalty_component(pc.penalty$range, penalty_range0, "pc.penalty$range"),
    sigma = .resolve_pc_penalty_component(pc.penalty$sigma, sigma_anchor0, "pc.penalty$sigma")
  )
}

.log_pc_prior_matern_internal <- function(log_range, log_sigma, range_spec, sigma_spec, d) {
  if (is.null(range_spec) || is.null(sigma_spec)) {
    stop("Both range_spec and sigma_spec must be provided.")
  }

  range <- exp(log_range)
  sigma <- exp(log_sigma)

  alpha_r <- as.numeric(range_spec["alpha"])
  alpha_s <- as.numeric(sigma_spec["alpha"])
  rho0 <- as.numeric(range_spec["anchor"])
  sigma0 <- as.numeric(sigma_spec["anchor"])

  R <- -log(alpha_r) * rho0^(d / 2)
  S <- -log(alpha_s) / sigma0

  log_prior_range <- log(d * R / 2) - (d / 2) * log_range - R * range^(-d / 2)
  log_prior_sigma <- log(S) + log_sigma - S * sigma

  as.numeric(log_prior_range + log_prior_sigma)
}

.build_matern_pc_spde <- function(mesh, alpha, pc_penalty, suppress_warnings = TRUE) {
  build <- function() {
    INLA::inla.spde2.pcmatern(
      mesh = mesh,
      alpha = alpha,
      prior.range = unname(pc_penalty$range[c("anchor", "alpha")]),
      prior.sigma = unname(pc_penalty$sigma[c("anchor", "alpha")])
    )
  }

  if (suppress_warnings) suppressWarnings(build()) else build()
}

.matern_posterior_from_inla_result <- function(res, stack, random_name = "spatial.field") {
  ii <- INLA::inla.stack.index(stack, tag = "est")$data

  posterior <- data.frame(
    mean = res$summary.fitted.values$mean[ii],
    var = res$summary.fitted.values$sd[ii]^2
  )
  posterior$second_moment <- posterior$mean^2 + posterior$var

  list(
    posterior = posterior,
    posterior_spatial_field = res$summary.random[[random_name]]
  )
}

.posterior_sampler_unavailable <- function(nsamp, n_obs) {
  nsamp <- .check_single_numeric(nsamp, "nsamp")
  if (nsamp < 1 || nsamp != floor(nsamp)) {
    stop("`nsamp` must be a positive integer.")
  }
  warning("posterior_sampler is not implemented for the INLA-backed Matern PC-prior mode.")
  matrix(NA_real_, nrow = nsamp, ncol = n_obs)
}

.check_matern_fit_for_objective_breakdown <- function(fit) {
  if (!inherits(fit, "ebnm")) {
    stop("`fit` must inherit from class \"ebnm\".")
  }
  if (is.null(fit$matern_objective_context)) {
    stop("This fit object does not contain the stored Matern objective context.")
  }
  if (is.null(fit$fitted_g) || is.null(fit$data)) {
    stop("This fit object is missing `fitted_g` or `data`.")
  }
  invisible(TRUE)
}

#' Compute an exact objective breakdown for a fitted Matern model
#'
#' @description
#' Re-evaluates the fitted Matern model under the package's manual Gaussian
#' formulas and returns a breakdown of exact objective components at the fitted
#' hyperparameters. This function is useful for post-fit validation and for
#' checking whether the recorded \code{log_likelihood} agrees with the exact
#' objective implied by a chosen beta-handling convention.
#'
#' @param fit An \code{ebnm} object returned by \code{\link{ebnm_Matern_generator}}.
#' @param beta_mode Character string indicating which exact beta-handling
#'   convention should be treated as the comparator. Supported values are
#'   \code{"auto"}, \code{"fixed"}, \code{"profile"}, and
#'   \code{"integrated_flat"}.
#'
#' @return A named list containing the fitted hyperparameters, exact
#'   log-likelihood values under fixed/profile/integrated-flat beta handling,
#'   the PC-prior term when relevant, the corresponding exact penalized
#'   objectives, and the gap between the recorded \code{fit$log_likelihood} and
#'   the selected exact comparator.
#'
#' @export
matern_objective_breakdown <- function(fit,
                                       beta_mode = c("auto", "fixed", "profile", "integrated_flat")) {
  .check_matern_fit_for_objective_breakdown(fit)
  beta_mode <- match.arg(beta_mode)

  ctx <- fit$matern_objective_context
  x <- as.numeric(fit$data$x)
  s <- as.numeric(fit$data$s)
  log_range <- .check_single_numeric(fit$fitted_g$theta, "fit$fitted_g$theta")
  log_sigma <- log(.check_single_numeric(fit$fitted_g$sigma, "fit$fitted_g$sigma"))
  beta_fixed <- .check_single_numeric(fit$fitted_beta, "fit$fitted_beta")

  stats <- .exact_matern_sufficient_stats(
    x = x,
    s = s,
    A = ctx$A,
    spde_template = ctx$spde_template,
    alpha = ctx$alpha,
    d = ctx$d,
    log_range = log_range,
    log_sigma = log_sigma
  )

  loglik_fixed_beta <- .exact_matern_loglik_fixed_beta(stats, beta0 = beta_fixed, n_obs = length(x))
  loglik_profile_beta <- .exact_matern_loglik_profile_beta(stats, n_obs = length(x))
  loglik_integrated_flat_beta <- .exact_matern_loglik_integrated_flat_beta(stats, n_obs = length(x))

  log_pc_prior_theta <- if (is.null(fit$pc_penalty)) {
    NULL
  } else {
    .log_pc_prior_matern_internal(
      log_range = log_range,
      log_sigma = log_sigma,
      range_spec = fit$pc_penalty$range,
      sigma_spec = fit$pc_penalty$sigma,
      d = ctx$d
    )
  }

  objective_fixed_plus_prior <- if (is.null(log_pc_prior_theta)) loglik_fixed_beta else loglik_fixed_beta + log_pc_prior_theta
  objective_profile_plus_prior <- if (is.null(log_pc_prior_theta)) loglik_profile_beta else loglik_profile_beta + log_pc_prior_theta
  objective_integrated_flat_plus_prior <- if (is.null(log_pc_prior_theta)) {
    loglik_integrated_flat_beta
  } else {
    loglik_integrated_flat_beta + log_pc_prior_theta
  }

  fit_beta_mode <- fit$beta_mode
  if (is.null(fit_beta_mode)) {
    fit_beta_mode <- if (identical(fit$backend, "exact")) {
      "profile"
    } else if (!is.null(fit$log_likelihood_stepB)) {
      "fixed"
    } else {
      "integrated_flat"
    }
  }

  matched_beta_mode <- if (beta_mode == "auto") fit_beta_mode else beta_mode
  matched_exact_objective <- switch(
    matched_beta_mode,
    fixed = objective_fixed_plus_prior,
    profile = objective_profile_plus_prior,
    integrated_flat = objective_integrated_flat_plus_prior,
    stop("Unsupported `beta_mode`: ", matched_beta_mode)
  )

  list(
    backend = fit$backend,
    beta_mode = fit_beta_mode,
    requested_beta_mode = beta_mode,
    matched_beta_mode = matched_beta_mode,
    log_range = log_range,
    log_sigma = log_sigma,
    range = stats$range,
    sigma = stats$sigma,
    beta_fixed = beta_fixed,
    beta_profile_hat = stats$beta_profile_hat,
    log_pc_prior_theta = log_pc_prior_theta,
    loglik_fixed_beta = loglik_fixed_beta,
    loglik_profile_beta = loglik_profile_beta,
    loglik_integrated_flat_beta = loglik_integrated_flat_beta,
    objective_fixed_plus_prior = objective_fixed_plus_prior,
    objective_profile_plus_prior = objective_profile_plus_prior,
    objective_integrated_flat_plus_prior = objective_integrated_flat_plus_prior,
    recorded_log_likelihood = as.numeric(fit$log_likelihood),
    matched_exact_objective = matched_exact_objective,
    recorded_minus_matched_exact = as.numeric(fit$log_likelihood) - matched_exact_objective
  )
}

.fit_matern_inla_stepA <- function(x,
                                   s,
                                   A,
                                   mesh,
                                   alpha,
                                   theta_init,
                                   pc_penalty,
                                   suppress_warnings = TRUE) {
  n <- length(x)
  spde <- .build_matern_pc_spde(mesh, alpha, pc_penalty, suppress_warnings = suppress_warnings)
  idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

  stackA <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(A, matrix(1, nrow = n, ncol = 1)),
    effects = list(
      spatial.field = idx$spatial.field,
      beta0 = 1
    ),
    tag = "est"
  )

  formulaA <- Y ~ 0 + beta0 + f(spatial.field, model = spde)

  runA <- function() {
    INLA::inla(
      formulaA,
      scale = (1 / s^2),
      control.inla = list(int.strategy = "eb", strategy = "gaussian"),
      control.family = list(
        control.link = list(model = "identity"),
        hyper = list(prec = list(fixed = TRUE, initial = 0))
      ),
      control.fixed = INLA::control.fixed(prec = 0),
      control.mode = INLA::control.mode(theta = theta_init, restart = TRUE),
      data = INLA::inla.stack.data(stackA),
      control.predictor = list(A = INLA::inla.stack.A(stackA), link = 1, compute = TRUE),
      control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
      silent = TRUE
    )
  }

  resA <- if (suppress_warnings) suppressWarnings(runA()) else runA()
  post <- .matern_posterior_from_inla_result(resA, stackA)
  theta_hat <- as.numeric(resA$mode$theta)

  if (length(theta_hat) != 2L) {
    stop("The INLA Step A fit did not return the expected two Matern hyperparameters.")
  }

  list(
    result = resA,
    posterior = post$posterior,
    posterior_spatial_field = post$posterior_spatial_field,
    fitted_g = Matern(theta = theta_hat[1], sigma = exp(theta_hat[2])),
    fitted_beta = as.numeric(resA$summary.fixed$mean[1]),
    log_likelihood_stepA_penalized = as.numeric(resA$misc$log.posterior.mode),
    log_likelihood_stepA_mlik_integration = as.numeric(resA$mlik["log marginal-likelihood (integration)", 1]),
    log_likelihood_stepA_mlik_gaussian = as.numeric(resA$mlik["log marginal-likelihood (Gaussian)", 1]),
    log_likelihood_stepA_joint_log_posterior = as.numeric(resA$joint.hyper[, "Log posterior density"])
  )
}

.fit_matern_inla_stepB <- function(x,
                                   s,
                                   A,
                                   mesh,
                                   alpha,
                                   theta_fixed,
                                   beta_fixed,
                                   pc_penalty,
                                   suppress_warnings = TRUE) {
  n <- length(x)
  spde <- .build_matern_pc_spde(mesh, alpha, pc_penalty, suppress_warnings = suppress_warnings)
  idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

  stackB <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(A, matrix(beta_fixed, nrow = n, ncol = 1)),
    effects = list(
      spatial.field = idx$spatial.field,
      beta0 = 1
    ),
    tag = "est"
  )

  formulaB <- Y ~ 0 + offset(beta0) + f(spatial.field, model = spde)

  runB <- function() {
    INLA::inla(
      formulaB,
      scale = (1 / s^2),
      control.inla = list(int.strategy = "eb", strategy = "gaussian"),
      control.family = list(
        control.link = list(model = "identity"),
        hyper = list(prec = list(fixed = TRUE, initial = 0))
      ),
      control.fixed = INLA::control.fixed(prec = 0),
      data = INLA::inla.stack.data(stackB),
      control.predictor = list(A = INLA::inla.stack.A(stackB), link = 1, compute = TRUE),
      control.mode = INLA::control.mode(theta = theta_fixed, fixed = TRUE),
      control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
      silent = TRUE
    )
  }

  resB <- if (suppress_warnings) suppressWarnings(runB()) else runB()
  post <- .matern_posterior_from_inla_result(resB, stackB)

  list(
    result = resB,
    posterior = post$posterior,
    posterior_spatial_field = post$posterior_spatial_field,
    fitted_g = Matern(theta = theta_fixed[1], sigma = exp(theta_fixed[2])),
    fitted_beta = as.numeric(beta_fixed),
    log_likelihood_stepB = as.numeric(resB$mlik["log marginal-likelihood (integration)", 1])
  )
}



#' Define the Matern GP Family
#'
#' @description Creates an object representing a stationary Matern prior,
#' parameterized by \code{theta = log(range)} and the latent-field marginal
#' standard deviation \code{sigma}.
#'
#' @param theta Numeric scalar. Log-transformed range parameter. Defaults to
#'   \code{NULL}.
#' @param sigma Numeric scalar. Latent-field marginal standard deviation on the
#'   natural scale. Defaults to \code{1}.
#'
#' @return An object of class \code{"Matern"}.
#'
#' @export
Matern <- function(theta = NULL, sigma = 1) {
  structure(list(theta = theta, sigma = sigma), class = "Matern")
}



#' Generate an `ebnm` function for Matern smoothing
#'
#' @description
#' Returns an \code{ebnm}-compatible function that fits a stationary Matern
#' smoother in one or two spatial dimensions.
#'
#' - If \code{locations} is a numeric vector, it is treated as one-dimensional
#'   locations.
#' - If \code{locations} is a matrix or data frame, it must have one or two
#'   columns.
#'
#' The current implementation supports the Gaussian normal-means case with the
#' identity link only. The public \code{link} argument is retained for backward
#' compatibility, but any value other than \code{"identity"} results in an
#' error.
#'
#' @details
#' Let \eqn{w} denote the latent Matern field on the mesh and \eqn{A} the
#' projector matrix from mesh nodes to observation locations. The observation
#' model is
#'
#' \deqn{x_i \mid w, \beta_0, s_i \sim N(\beta_0 + (Aw)_i, s_i^2).}
#'
#' The latent field prior is Gaussian with precision matrix \eqn{Q(range, sigma)}
#' induced by the SPDE representation of a stationary Matern field with fixed
#' \code{alpha}. For fixed \code{(range, sigma, beta0)}, the exact marginal
#' likelihood is
#'
#' \deqn{
#' \ell = -\frac{1}{2} \left[
#' n \log(2\pi) + \log|D| - \log|Q| + \log|Q_{post}| +
#' r^T D^{-1} r - b^T Q_{post}^{-1} b
#' \right],
#' }
#'
#' where \eqn{D = \mathrm{diag}(s_1^2, \ldots, s_n^2)},
#' \eqn{r = x - \beta_0}, \eqn{b = A^T D^{-1} r}, and
#' \eqn{Q_{post} = Q + A^T D^{-1} A}.
#'
#' If \code{pc.penalty = NULL}, the fitter uses the exact Gaussian
#' empirical-Bayes implementation and maximizes the exact marginal likelihood
#' over \code{(log(range), log(sigma), beta0)} using
#' \code{\link[stats:optim]{optim}} with method \code{"BFGS"}.
#'
#' If \code{pc.penalty} is supplied, the fitter switches to an INLA backend
#' using \code{\link[INLA:inla.spde2.pcmatern]{inla.spde2.pcmatern()}}.
#' With \code{fix_g = FALSE}, the primary \code{log_likelihood} is then the
#' Step A penalized objective \code{resA$misc$log.posterior.mode}. With
#' \code{fix_g = TRUE}, the primary \code{log_likelihood} is a penalized
#' pseudo-objective equal to the Step B conditional marginal likelihood plus
#' the PC prior log-density on the internal \code{(log(range), log(sigma))}
#' scale.
#'
#' By default, the INLA PC-prior path does \emph{not} automatically recompute
#' the exact Gaussian objective at the fitted Step A mode. That cross-check can
#' be enabled with \code{compute_exact_diagnostic = TRUE}, or performed later by
#' calling \code{\link{matern_objective_breakdown}} on the fitted object.
#'
#' @param locations Spatial locations.
#' @param max.edge Mesh maximum edge length.
#' @param alpha Smoothness parameter used in the SPDE representation. It is held
#'   fixed and is not optimized.
#' @param suppress_warnings If \code{TRUE}, suppress warnings from INLA mesh and
#'   SPDE construction.
#' @param penalty_range Initial range anchor used when \code{g_init$theta} is
#'   missing. If \code{NULL}, defaults to roughly one-tenth of the spatial
#'   extent.
#' @param pc.penalty Optional list specifying PC priors for the INLA backend.
#'   Supported entries are \code{range} and \code{sigma}. Each entry may be a
#'   numeric vector of length 1, interpreted as the anchor with default tail
#'   probability \code{0.5}, or length 2, interpreted as \code{c(anchor, alpha)}.
#'   Missing entries default to a data-driven anchor and tail probability
#'   \code{0.5}. If \code{NULL}, the exact Gaussian EB backend is used.
#' @param compute_exact_diagnostic Logical scalar. If \code{TRUE} and
#'   \code{pc.penalty} is supplied with \code{fix_g = FALSE}, also evaluate the
#'   manual exact Gaussian log-likelihood at the fitted INLA Step A mode and
#'   store it in \code{log_likelihood_exact_at_stepA_mode}. Defaults to
#'   \code{FALSE}.
#' @param link Backward-compatible link argument. Only \code{"identity"} is
#'   supported by the current Matern implementation.
#'
#' @return A function with signature
#'   \code{function(x, s, g_init = NULL, fix_g = FALSE, beta_fixed = NULL, output = NULL)}.
#'   The returned object includes posterior summaries, fitted hyperparameters,
#'   backend-specific likelihood diagnostics, the mesh, and a summary of the
#'   latent spatial field.
#'
#' @export
ebnm_Matern_generator <- function(locations,
                                  max.edge = NULL,
                                  alpha = 2,
                                  suppress_warnings = TRUE,
                                  penalty_range = NULL,
                                  pc.penalty = NULL,
                                  compute_exact_diagnostic = FALSE,
                                  link = c("identity", "log", "logit", "probit")) {

  link <- match.arg(link)
  if (link != "identity") {
    stop("The exact Matern implementation currently supports `link = \"identity\"` only.")
  }

  loc_info <- .normalize_locations(locations)
  loc_mat <- loc_info$loc
  d <- loc_info$d

  if (alpha <= d / 2) {
    stop("`alpha` must satisfy alpha > d / 2.")
  }
  if (!is.logical(compute_exact_diagnostic) || length(compute_exact_diagnostic) != 1L || is.na(compute_exact_diagnostic)) {
    stop("`compute_exact_diagnostic` must be TRUE or FALSE.")
  }

  meshA <- if (suppress_warnings) suppressWarnings(.build_mesh_A(loc_mat, max.edge = max.edge)) else .build_mesh_A(loc_mat, max.edge = max.edge)
  mesh <- meshA$mesh
  A <- meshA$A

  penalty_range0 <- if (is.null(penalty_range)) .default_penalty_range(loc_mat) else penalty_range
  if (!is.numeric(penalty_range0) || length(penalty_range0) != 1L || is.na(penalty_range0) || penalty_range0 <= 0) {
    stop("`penalty_range` must be a single positive number (or NULL).")
  }

  spde_template <- if (suppress_warnings) {
    suppressWarnings(INLA::inla.spde2.matern(mesh = mesh, alpha = alpha))
  } else {
    INLA::inla.spde2.matern(mesh = mesh, alpha = alpha)
  }

  ebnm_Matern <- function(x, s, g_init = NULL, fix_g = FALSE, beta_fixed = NULL, output = NULL) {
    n <- nrow(loc_mat)

    if (length(x) != n) {
      warning(
        paste0(
          "The length of x must equal the number of locations.\n",
          "length(x) = ", length(x), ", nrow(locations) = ", n, ".\n"
        )
      )
      if (length(s) == 3 && length(x) == 3) {
        warning("Assume this is just an initialization check. Returning ebnm_flat(x).")
        return(ebnm::ebnm_flat(x))
      }
      stop("The length of x must equal the number of locations.")
    }

    if (!(length(s) == 1L || length(s) == length(x))) stop("`s` must have length 1 or length(x).")
    if (length(s) == 1L) s <- rep(s, length(x))
    if (anyNA(x) || anyNA(s)) stop("`x` and `s` must not contain NA.")
    if (any(s <= 0)) stop("All entries of `s` must be > 0.")

    sigma_data <- stats::sd(x)
    if (!is.finite(sigma_data) || sigma_data <= 0) sigma_data <- 1
    pc_penalty0 <- .resolve_matern_pc_penalty(pc.penalty, penalty_range0, sigma_data)

    g_init_missing <- is.null(g_init)
    if (is.null(g_init)) g_init <- Matern(theta = log(penalty_range0), sigma = 1)

    theta_init <- if (is.null(g_init$theta)) log(penalty_range0) else .check_single_numeric(g_init$theta, "g_init$theta")
    sigma_init_from_g <- if (is.null(g_init$sigma)) 1 else .check_single_numeric(g_init$sigma, "g_init$sigma")

    if (!fix_g && !is.null(beta_fixed)) {
      stop("`beta_fixed` can only be supplied when `fix_g = TRUE`.")
    }
    if (fix_g) {
      if (is.null(beta_fixed)) {
        beta0_fixed <- 0
      } else {
        beta0_fixed <- .check_single_numeric(beta_fixed, "beta_fixed")
      }
    } else {
      beta0_fixed <- NULL
    }

    raw_eval_state <- function(par) {
      .exact_matern_state(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[1],
        log_sigma = par[2],
        beta0 = par[3]
      )
    }

    eval_state <- function(par) {
      if (suppress_warnings) {
        suppressWarnings(raw_eval_state(par))
      } else {
        raw_eval_state(par)
      }
    }

    safe_objective <- function(par) {
      state <- tryCatch(eval_state(par), error = function(e) e)
      if (inherits(state, "error")) return(Inf)
      -state$log_marginal
    }

    if (is.null(pc_penalty0)) {
      if (!fix_g) {
        par0 <- c(
          theta_init,
          if (g_init_missing || is.null(g_init$sigma)) log(sigma_data) else log(sigma_init_from_g),
          stats::weighted.mean(x, 1 / (s^2))
        )

        opt <- optim(
          par = par0,
          fn = safe_objective,
          method = "BFGS"
        )
        if (!is.finite(opt$value)) {
          stop("Exact Matern EB optimization failed to produce a finite objective value.")
        }

        state <- tryCatch(eval_state(opt$par), error = function(e) e)
        if (inherits(state, "error")) {
          stop("Exact Matern EB optimization failed at the reported optimum: ", conditionMessage(state))
        }

        fitted_g <- state$fitted_g
        fitted_beta <- state$fitted_beta
        log_likelihood <- state$log_marginal
      } else {
        par_fixed <- c(theta_init, log(sigma_init_from_g), beta0_fixed)
        state <- tryCatch(eval_state(par_fixed), error = function(e) e)
        if (inherits(state, "error")) {
          stop("Exact Matern fixed-parameter evaluation failed: ", conditionMessage(state))
        }

        fitted_g <- state$fitted_g
        fitted_beta <- state$fitted_beta
        log_likelihood <- state$log_marginal
      }

      class(log_likelihood) <- "logLik"

      posterior_sampler <- function(nsamp) {
        nsamp <- .check_single_numeric(nsamp, "nsamp")
        if (nsamp < 1 || nsamp != floor(nsamp)) {
          stop("`nsamp` must be a positive integer.")
        }

        samps_w <- LaplacesDemon::rmvnp(
          n = nsamp,
          mu = state$post_mean_latent,
          Omega = as.matrix(state$Q_post)
        )
        if (is.null(dim(samps_w))) {
          samps_w <- matrix(samps_w, nrow = 1)
        }

        obs_samps <- t(as.matrix(A %*% t(samps_w)))
        sweep(obs_samps, 2, fitted_beta, `+`)
      }

      out <- list(
        posterior = state$posterior,
        fitted_g = fitted_g,
        fitted_beta = fitted_beta,
        beta_mode = if (fix_g) "fixed" else "profile",
        g_init = g_init,
        log_likelihood = log_likelihood,
        log_likelihood_semantics = if (fix_g) "exact_fixed" else "exact_profile",
        posterior_sampler = posterior_sampler,
        data = data.frame(x = x, s = s),
        prior_family = "identity_Matern",
        posterior_spatial_field = state$posterior_spatial_field,
        mesh = mesh,
        inla_result = NULL,
        backend = "exact",
        pc_penalty = NULL,
        matern_objective_context = list(
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d
        )
      )
    } else if (!fix_g) {
      theta_init_pc <- c(
        theta_init,
        if (g_init_missing || is.null(g_init$sigma)) log(sigma_data) else log(sigma_init_from_g)
      )

      pc_fit <- .fit_matern_inla_stepA(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        mesh = mesh,
        alpha = alpha,
        theta_init = theta_init_pc,
        pc_penalty = pc_penalty0,
        suppress_warnings = suppress_warnings
      )

      if (isTRUE(compute_exact_diagnostic)) {
        state_exact <- tryCatch(
          eval_state(c(pc_fit$fitted_g$theta, log(pc_fit$fitted_g$sigma), pc_fit$fitted_beta)),
          error = function(e) e
        )
        if (inherits(state_exact, "error")) {
          stop("Failed to evaluate the exact Gaussian objective at the INLA Step A mode: ", conditionMessage(state_exact))
        }
        log_likelihood_exact_at_stepA_mode <- state_exact$log_marginal
      } else {
        log_likelihood_exact_at_stepA_mode <- NULL
      }

      log_likelihood <- pc_fit$log_likelihood_stepA_penalized
      class(log_likelihood) <- "logLik"

      out <- list(
        posterior = pc_fit$posterior,
        fitted_g = pc_fit$fitted_g,
        fitted_beta = pc_fit$fitted_beta,
        beta_mode = "integrated_flat",
        g_init = g_init,
        log_likelihood = log_likelihood,
        log_likelihood_semantics = "stepA_penalized",
        posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
        data = data.frame(x = x, s = s),
        prior_family = "identity_Matern_pc",
        posterior_spatial_field = pc_fit$posterior_spatial_field,
        mesh = mesh,
        inla_result = pc_fit$result,
        backend = "inla_pc",
        pc_penalty = pc_penalty0,
        log_likelihood_stepA_penalized = pc_fit$log_likelihood_stepA_penalized,
        log_likelihood_stepA_mlik_integration = pc_fit$log_likelihood_stepA_mlik_integration,
        log_likelihood_stepA_mlik_gaussian = pc_fit$log_likelihood_stepA_mlik_gaussian,
        log_likelihood_stepA_joint_log_posterior = pc_fit$log_likelihood_stepA_joint_log_posterior,
        log_likelihood_exact_at_stepA_mode = log_likelihood_exact_at_stepA_mode,
        log_likelihood_stepB = NULL,
        log_likelihood_pc_prior_theta = NULL,
        matern_objective_context = list(
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d
        )
      )
    } else {
      theta_fixed <- c(theta_init, log(sigma_init_from_g))
      pc_fit <- .fit_matern_inla_stepB(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        mesh = mesh,
        alpha = alpha,
        theta_fixed = theta_fixed,
        beta_fixed = beta0_fixed,
        pc_penalty = pc_penalty0,
        suppress_warnings = suppress_warnings
      )

      log_likelihood_pc_prior_theta <- .log_pc_prior_matern_internal(
        log_range = theta_fixed[1],
        log_sigma = theta_fixed[2],
        range_spec = pc_penalty0$range,
        sigma_spec = pc_penalty0$sigma,
        d = d
      )
      log_likelihood <- pc_fit$log_likelihood_stepB + log_likelihood_pc_prior_theta
      class(log_likelihood) <- "logLik"

      out <- list(
        posterior = pc_fit$posterior,
        fitted_g = pc_fit$fitted_g,
        fitted_beta = pc_fit$fitted_beta,
        beta_mode = "fixed",
        g_init = g_init,
        log_likelihood = log_likelihood,
        log_likelihood_semantics = "stepB_plus_pc_prior",
        posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
        data = data.frame(x = x, s = s),
        prior_family = "identity_Matern_pc",
        posterior_spatial_field = pc_fit$posterior_spatial_field,
        mesh = mesh,
        inla_result = pc_fit$result,
        backend = "inla_pc",
        pc_penalty = pc_penalty0,
        log_likelihood_stepA_penalized = NULL,
        log_likelihood_stepA_mlik_integration = NULL,
        log_likelihood_stepA_mlik_gaussian = NULL,
        log_likelihood_stepA_joint_log_posterior = NULL,
        log_likelihood_exact_at_stepA_mode = NULL,
        log_likelihood_stepB = pc_fit$log_likelihood_stepB,
        log_likelihood_pc_prior_theta = log_likelihood_pc_prior_theta,
        matern_objective_context = list(
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d
        )
      )
    }

    structure(out, class = c("list", "ebnm"))
  }

  if (suppress_warnings) suppressWarnings(ebnm_Matern) else ebnm_Matern
}
