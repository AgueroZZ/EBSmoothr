# ---- helper: normalize locations and build mesh/A for d = 1 or 2 ----
utils::globalVariables(".scale_w")

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

.check_matern_n_starts <- function(matern_n_starts) {
  if (!is.numeric(matern_n_starts) || length(matern_n_starts) != 1L ||
      is.na(matern_n_starts) || !is.finite(matern_n_starts)) {
    stop("`matern_n_starts` must be a single finite integer between 1 and 5.")
  }
  if (abs(matern_n_starts - round(matern_n_starts)) > sqrt(.Machine$double.eps)) {
    stop("`matern_n_starts` must be a single finite integer between 1 and 5.")
  }
  matern_n_starts <- as.integer(round(matern_n_starts))
  if (matern_n_starts < 1L || matern_n_starts > 5L) {
    stop("`matern_n_starts` must be between 1 and 5.")
  }
  matern_n_starts
}

.check_matern_pql_inner_iter <- function(pql_inner_iter) {
  pql_inner_iter <- .check_single_numeric(pql_inner_iter, "pql_inner_iter")
  if (!is.finite(pql_inner_iter) || pql_inner_iter < 1 ||
      pql_inner_iter != floor(pql_inner_iter)) {
    stop("`pql_inner_iter` must be a positive integer for Fisher-PQL.")
  }
  as.integer(pql_inner_iter)
}

.build_mesh_A <- function(loc_mat, max.edge = NULL) {
  .with_quiet_inla_defaults({
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
      inner_edge <- min(ranges) / 10
      if (!is.finite(inner_edge) || inner_edge <= 0) inner_edge <- 1
      max.edge <- c(inner_edge, 3 * inner_edge)
    }
    if (length(max.edge) == 1L) max.edge <- c(max.edge, max.edge)
    if (length(max.edge) != 2L) stop("For d = 2, max.edge must be NULL, length-1, or length-2.")

    mesh <- INLA::inla.mesh.2d(loc = loc_mat, max.edge = max.edge)
    A <- INLA::inla.spde.make.A(mesh = mesh, loc = loc_mat)
    list(mesh = mesh, A = Matrix::Matrix(A, sparse = TRUE))
  })
}

.check_single_numeric <- function(z, nm) {
  if (!is.numeric(z) || length(z) != 1L || is.na(z)) {
    stop(nm, " must be a single non-NA numeric.")
  }
  as.numeric(z)
}

.matern_auto_backend <- function(link,
                                 d,
                                 beta_mode,
                                 pc_penalty = NULL,
                                 fix_g = FALSE,
                                 learn_noise = FALSE) {
  if (identical(link, "identity")) {
    return("exact")
  }

  if (identical(link, "log") || identical(link, "softplus")) {
    return("fisher_pql")
  }

  "laplace"
}

.matern_resolve_backend <- function(backend,
                                    link,
                                    d,
                                    beta_mode,
                                    pc_penalty = NULL,
                                    fix_g = FALSE,
                                    learn_noise = FALSE) {
  if (identical(backend, "auto")) {
    return(.matern_auto_backend(
      link = link,
      d = d,
      beta_mode = beta_mode,
      pc_penalty = pc_penalty,
      fix_g = fix_g,
      learn_noise = learn_noise
    ))
  }
  backend <- .matern_canonical_backend(backend, pc_penalty = pc_penalty)
  if (identical(backend, "laplace_fisher") && !(identical(link, "log") || identical(link, "softplus"))) {
    stop("`backend = \"laplace_fisher\"` is only available for `link = \"log\"` or `link = \"softplus\"`.")
  }
  if (identical(backend, "fisher_pql") && !(identical(link, "log") || identical(link, "softplus"))) {
    stop("`backend = \"fisher_pql\"` is only available for `link = \"log\"` or `link = \"softplus\"`.")
  }
  if (identical(backend, "inlabru") && !identical(link, "softplus")) {
    stop("`backend = \"inlabru\"` is currently only supported for `link = \"softplus\"`.")
  }
  if (identical(backend, "exact") && (identical(link, "log") || identical(link, "softplus"))) {
    return("laplace")
  }
  backend
}

.match_matern_backend_arg <- function(backend) {
  choices <- c("auto", "exact", "laplace", "laplace_fisher", "fisher_pql", "laplace_r", "inla", "laplace_tmb", "inla_pc", "inlabru")
  if (length(backend) > 1L) backend <- backend[1L]
  match.arg(backend, choices = choices)
}

.matern_canonical_backend <- function(backend, pc_penalty = NULL) {
  if (identical(backend, "laplace_tmb")) {
    return("laplace")
  }
  if (identical(backend, "inla_pc")) {
    if (is.null(pc_penalty)) {
      stop("`backend = \"inla_pc\"` requires `pc.penalty`; use `backend = \"inla\"` for unpenalized INLA fits.")
    }
    return("inla")
  }
  backend
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

.matern_range_sigma_from_tau_kappa <- function(tau, kappa, alpha, d) {
  nu <- alpha - d / 2
  if (!is.finite(nu) || nu <= 0) {
    stop("`alpha` must satisfy alpha > d / 2 for the Matern SPDE parameterization.")
  }
  if (!is.finite(tau) || tau <= 0) stop("Matern tau must be positive.")
  if (!is.finite(kappa) || kappa <= 0) stop("Matern kappa must be positive.")

  range <- sqrt(8 * nu) / kappa
  sigma <- sqrt(
    gamma(nu) /
      (gamma(alpha) * (4 * pi)^(d / 2) * kappa^(2 * nu) * tau^2)
  )

  list(range = range, sigma = sigma, nu = nu)
}

.matern_spde_theta_from_log_range_log_sigma <- function(log_range, log_sigma, alpha, d) {
  ts <- .matern_tau_from_range_sigma(
    range = exp(log_range),
    sigma = exp(log_sigma),
    alpha = alpha,
    d = d
  )
  c(log(ts$tau), log(ts$kappa))
}

.matern_precision_from_log_params <- function(spde_template, alpha, d, log_range, log_sigma) {
  ts <- .matern_tau_from_range_sigma(
    range = exp(log_range),
    sigma = exp(log_sigma),
    alpha = alpha,
    d = d
  )
  theta_spde <- c(log(ts$tau), log(ts$kappa))
  Q <- INLA::inla.spde2.precision(spde_template, theta = theta_spde)
  Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))

  list(
    Q = Q,
    Q_factor = .factorize_spd(Q),
    range = exp(log_range),
    sigma = exp(log_sigma),
    kappa = ts$kappa,
    tau = ts$tau
  )
}

.exact_matern_sufficient_stats <- function(x, s, A, spde_template, alpha, d, log_range, log_sigma) {
  .with_quiet_inla_defaults({
    if (any(!is.finite(c(log_range, log_sigma)))) {
      stop("All Matern hyperparameters must be finite.")
    }

    range <- exp(log_range)
    sigma <- exp(log_sigma)

    if (!is.finite(range) || range <= 0) stop("Matern range must be positive.")
    if (!is.finite(sigma) || sigma <= 0) stop("Matern sigma must be positive.")

    precision <- .matern_precision_from_log_params(
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      log_range = log_range,
      log_sigma = log_sigma
    )
    Q <- precision$Q
    Q_factor <- precision$Q_factor

    w_prec_diag <- 1 / (s^2)
    W <- Matrix::Diagonal(x = w_prec_diag)
    AtW <- Matrix::t(A) %*% W

    Q_post <- Q + AtW %*% A
    Q_post <- Matrix::forceSymmetric(Matrix::Matrix(Q_post, sparse = TRUE))
    Q_post_factor <- .factorize_spd(Q_post)

    solve_post <- function(rhs) as.numeric(.solve_spd_factor(Q_post_factor, rhs))

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
    logdet_Q <- Q_factor$logdet
    logdet_Q_post <- Q_post_factor$logdet
    logdet_Sigma <- logdet_D - logdet_Q + logdet_Q_post

    list(
      range = range,
      sigma = sigma,
      kappa = precision$kappa,
      tau = precision$tau,
      Q = Q,
      Q_factor = Q_factor,
      Q_post = Q_post,
      Q_post_factor = Q_post_factor,
      logdet_Sigma = as.numeric(logdet_Sigma),
      quad_x = as.numeric(quad_x),
      quad_u = as.numeric(quad_u),
      cross_ux = as.numeric(cross_ux),
      beta_profile_hat = as.numeric(beta_profile_hat),
      quad_profile = as.numeric(quad_profile)
    )
  })
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

.exact_matern_loglik_prior_beta <- function(stats, beta_prec, n_obs) {
  beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
  if (is.null(beta_prec) || beta_prec <= 0) {
    stop("`beta_prec` must be a single positive number for the proper-beta objective.")
  }

  quad_prior <- stats$quad_x - (stats$cross_ux^2) / (stats$quad_u + beta_prec)
  as.numeric(-0.5 * (
    n_obs * log(2 * pi) +
      stats$logdet_Sigma +
      log(stats$quad_u + beta_prec) -
      log(beta_prec) +
      quad_prior
  ))
}

.exact_matern_profile_objective <- function(x, s, A, spde_template, alpha, d, log_range, log_sigma) {
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

  list(
    log_marginal = .exact_matern_loglik_profile_beta(stats, n_obs = length(x)),
    fitted_beta = stats$beta_profile_hat,
    stats = stats
  )
}

.exact_matern_fixed_objective <- function(x,
                                          s,
                                          A,
                                          spde_template,
                                          alpha,
                                          d,
                                          log_range,
                                          log_sigma,
                                          beta0) {
  beta0 <- .check_single_numeric(beta0, "beta0")
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

  list(
    log_marginal = .exact_matern_loglik_fixed_beta(stats, beta0 = beta0, n_obs = length(x)),
    fitted_beta = beta0,
    stats = stats
  )
}

.exact_matern_prior_objective <- function(x,
                                          s,
                                          A,
                                          spde_template,
                                          alpha,
                                          d,
                                          log_range,
                                          log_sigma,
                                          beta_prec) {
  beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
  if (is.null(beta_prec) || beta_prec <= 0) {
    stop("`beta_prec` must be positive for the proper-beta objective.")
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

  list(
    log_marginal = .exact_matern_loglik_prior_beta(stats, beta_prec = beta_prec, n_obs = length(x)),
    fitted_beta = as.numeric(stats$cross_ux / (stats$quad_u + beta_prec)),
    stats = stats,
    beta_prec = beta_prec
  )
}

.exact_matern_unknown_noise_s <- function(noise_sd, n, noise_scale = NULL) {
  if (is.null(noise_scale)) {
    return(rep(noise_sd, n))
  }

  noise_scale <- as.numeric(noise_scale)
  if (length(noise_scale) != n ||
      anyNA(noise_scale) ||
      any(!is.finite(noise_scale)) ||
      any(noise_scale <= 0)) {
    stop("`noise_scale` must be NULL or a positive finite numeric vector with length equal to `length(x)`.")
  }
  as.numeric(noise_sd * noise_scale)
}

.exact_matern_profile_objective_unknown_noise <- function(x,
                                                          A,
                                                          spde_template,
                                                          alpha,
                                                          d,
                                                          log_range,
                                                          log_sigma,
                                                          log_noise_sd,
                                                          noise_scale = NULL) {
  if (!is.finite(log_noise_sd)) {
    stop("`log_noise_sd` must be finite.")
  }

  noise_sd <- exp(log_noise_sd)
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    stop("The learned observation noise SD must be positive.")
  }

  s_eff <- .exact_matern_unknown_noise_s(noise_sd, length(x), noise_scale)
  objective <- .exact_matern_profile_objective(
    x = x,
    s = s_eff,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = log_range,
    log_sigma = log_sigma
  )
  objective$fitted_noise_sd <- as.numeric(noise_sd)
  objective$s <- s_eff
  if (!is.null(noise_scale)) objective$noise_scale <- as.numeric(noise_scale)
  objective
}

.exact_matern_fixed_objective_unknown_noise <- function(x,
                                                        A,
                                                        spde_template,
                                                        alpha,
                                                        d,
                                                        log_range,
                                                        log_sigma,
                                                        beta0,
                                                        log_noise_sd,
                                                        noise_scale = NULL) {
  beta0 <- .check_single_numeric(beta0, "beta0")
  if (!is.finite(log_noise_sd)) {
    stop("`log_noise_sd` must be finite.")
  }

  noise_sd <- exp(log_noise_sd)
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    stop("The learned observation noise SD must be positive.")
  }

  s_eff <- .exact_matern_unknown_noise_s(noise_sd, length(x), noise_scale)
  objective <- .exact_matern_fixed_objective(
    x = x,
    s = s_eff,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = log_range,
    log_sigma = log_sigma,
    beta0 = beta0
  )
  objective$fitted_noise_sd <- as.numeric(noise_sd)
  objective$s <- s_eff
  if (!is.null(noise_scale)) objective$noise_scale <- as.numeric(noise_scale)
  objective
}

.exact_matern_prior_objective_unknown_noise <- function(x,
                                                        A,
                                                        spde_template,
                                                        alpha,
                                                        d,
                                                        log_range,
                                                        log_sigma,
                                                        beta_prec,
                                                        log_noise_sd,
                                                        noise_scale = NULL) {
  beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
  if (is.null(beta_prec) || beta_prec <= 0) {
    stop("`beta_prec` must be positive for the proper-beta objective.")
  }
  if (!is.finite(log_noise_sd)) {
    stop("`log_noise_sd` must be finite.")
  }

  noise_sd <- exp(log_noise_sd)
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    stop("The learned observation noise SD must be positive.")
  }

  s_eff <- .exact_matern_unknown_noise_s(noise_sd, length(x), noise_scale)
  objective <- .exact_matern_prior_objective(
    x = x,
    s = s_eff,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = log_range,
    log_sigma = log_sigma,
    beta_prec = beta_prec
  )
  objective$fitted_noise_sd <- as.numeric(noise_sd)
  objective$s <- s_eff
  if (!is.null(noise_scale)) objective$noise_scale <- as.numeric(noise_scale)
  objective
}

.optimize_exact_matern_objective <- function(par0,
                                             safe_objective,
                                             eval_objective,
                                             failure_label) {
  methods <- c("BFGS", "Nelder-Mead")
  last_message <- NULL

  for (method in methods) {
    opt <- tryCatch(
      optim(par = par0, fn = safe_objective, method = method),
      error = function(e) e
    )
    if (inherits(opt, "error")) {
      last_message <- conditionMessage(opt)
      next
    }
    if (!is.finite(opt$value) || opt$value >= 1e99) {
      last_message <- sprintf("%s optimization reached the finite penalty under method %s.", failure_label, method)
      next
    }

    objective_opt <- tryCatch(eval_objective(opt$par), error = function(e) e)
    if (inherits(objective_opt, "error")) {
      last_message <- conditionMessage(objective_opt)
      next
    }

    return(list(opt = opt, objective = objective_opt, method = method))
  }

  stop(failure_label, ": ", last_message)
}

.matern_fixed_log_param_names <- function(fix_params) {
  out <- character()
  if ("range" %in% fix_params) out <- c(out, "log_range")
  if ("sigma" %in% fix_params) out <- c(out, "log_sigma")
  out
}

.validate_matern_fixed_param_values <- function(fix_params, g_init) {
  if ("range" %in% fix_params && (is.null(g_init) || is.null(g_init$theta))) {
    stop("`fix_params = \"range\"` requires `g_init$theta`.")
  }
  if ("sigma" %in% fix_params && (is.null(g_init) || is.null(g_init$sigma))) {
    stop("`fix_params = \"sigma\"` requires `g_init$sigma`.")
  }
  invisible(TRUE)
}

.matern_optimize_or_eval_objective <- function(par0,
                                               fixed_names = character(),
                                               safe_objective,
                                               eval_objective,
                                               failure_label) {
  if (is.null(names(par0)) || any(names(par0) == "")) {
    stop("Internal error: Matern optimizer parameters must be named.")
  }
  fixed_names <- intersect(unique(fixed_names), names(par0))
  free_names <- setdiff(names(par0), fixed_names)
  expand_par <- function(par_free) {
    out <- par0
    out[free_names] <- as.numeric(par_free)
    out
  }
  if (!length(free_names)) {
    return(list(
      opt = NULL,
      objective = eval_objective(par0),
      method = "fixed"
    ))
  }
  opt_res <- .optimize_exact_matern_objective(
    par0 = par0[free_names],
    safe_objective = function(par_free) safe_objective(expand_par(par_free)),
    eval_objective = function(par_free) eval_objective(expand_par(par_free)),
    failure_label = failure_label
  )
  opt_res$free_names <- free_names
  opt_res$fixed_names <- fixed_names
  opt_res
}

.matern_exact_optimization_diagnostics <- function(opt_res) {
  opt <- opt_res$opt
  if (is.null(opt)) {
    return(list(
      method = opt_res$method,
      convergence = 0L,
      message = "all exact Gaussian Matern hyperparameters fixed",
      counts = NULL,
      value = NA_real_,
      fixed_names = if (is.null(opt_res$fixed_names)) character() else opt_res$fixed_names,
      free_names = if (is.null(opt_res$free_names)) character() else opt_res$free_names
    ))
  }

  list(
    method = opt_res$method,
    convergence = if (is.null(opt$convergence)) NA_integer_ else as.integer(opt$convergence),
    message = if (is.null(opt$message)) "" else as.character(opt$message),
    counts = opt$counts,
    value = as.numeric(opt$value),
    fixed_names = if (is.null(opt_res$fixed_names)) character() else opt_res$fixed_names,
    free_names = if (is.null(opt_res$free_names)) character() else opt_res$free_names
  )
}

.exact_matern_posterior_from_stats <- function(x, s, A, stats, beta0, log_marginal, beta_prec = NULL) {
  beta0 <- .check_single_numeric(beta0, "beta0")
  beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
  w_prec_diag <- 1 / (s^2)
  W <- Matrix::Diagonal(x = w_prec_diag)
  AtW <- Matrix::t(A) %*% W
  residual <- as.numeric(x - beta0)
  b <- as.numeric(AtW %*% residual)
  Q_post_factor <- if (is.null(stats$Q_post_factor)) .factorize_spd(stats$Q_post) else stats$Q_post_factor
  m_w <- as.numeric(.solve_spd_factor(Q_post_factor, b))

  obs_mean <- as.numeric(beta0 + A %*% m_w)
  obs_var <- .compute_diag_A_Qinv_At(A, Q_post_factor)
  latent_var <- .compute_diag_A_Qinv_At(Matrix::Diagonal(n = nrow(stats$Q_post)), Q_post_factor)
  latent_sd <- sqrt(pmax(latent_var, 0))

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
    fitted_g = Matern(theta = log(stats$range), sigma = stats$sigma, beta = beta0, beta_prec = beta_prec),
    fitted_beta = as.numeric(beta0),
    beta_prec = beta_prec,
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

.exact_matern_joint_posterior_from_stats <- function(x,
                                                     s,
                                                     A,
                                                     stats,
                                                     beta_prec,
                                                     log_marginal) {
  beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
  if (is.null(beta_prec) || beta_prec < 0) {
    stop("`beta_prec` must be non-negative for the joint-beta posterior.")
  }

  x <- as.numeric(x)
  s <- as.numeric(s)
  w_prec_diag <- 1 / (s^2)
  W <- Matrix::Diagonal(x = w_prec_diag)
  AtW <- Matrix::t(A) %*% W
  c_x <- as.numeric(AtW %*% x)
  c_u <- as.numeric(AtW %*% rep(1, length(x)))
  cross_x <- sum(w_prec_diag * x)
  quad_u_raw <- sum(w_prec_diag)
  A_aug <- cbind(A, Matrix::Matrix(rep(1, length(x)), ncol = 1, sparse = TRUE))
  rhs <- c(c_x, cross_x)
  Q_joint <- rbind(
    cbind(stats$Q_post, Matrix::Matrix(c_u, ncol = 1, sparse = TRUE)),
    cbind(
      Matrix::Matrix(t(c_u), sparse = TRUE),
      Matrix::Matrix(quad_u_raw + beta_prec, nrow = 1, ncol = 1, sparse = TRUE)
    )
  )
  Q_joint <- Matrix::forceSymmetric(Matrix::Matrix(Q_joint, sparse = TRUE))
  Q_joint_factor <- .factorize_spd(Q_joint)
  post_mean_joint <- as.numeric(.solve_spd_factor(Q_joint_factor, rhs))

  m_w <- post_mean_joint[seq_len(length(c_x))]
  beta_hat <- post_mean_joint[length(post_mean_joint)]
  obs_mean <- as.numeric(A %*% m_w + beta_hat)
  obs_var <- .compute_diag_A_Qinv_At(A_aug, Q_joint_factor)
  A_latent <- cbind(
    Matrix::Diagonal(n = nrow(stats$Q_post)),
    Matrix::Matrix(0, nrow = nrow(stats$Q_post), ncol = 1, sparse = TRUE)
  )
  latent_var <- .compute_diag_A_Qinv_At(A_latent, Q_joint_factor)
  latent_sd <- sqrt(pmax(latent_var, 0))
  beta_var <- as.numeric(.solve_spd_factor(Q_joint_factor, c(rep(0, length(c_x)), 1))[length(post_mean_joint)])

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
    fitted_g = Matern(theta = log(stats$range), sigma = stats$sigma, beta = beta_hat, beta_prec = beta_prec),
    fitted_beta = as.numeric(beta_hat),
    beta_prec = beta_prec,
    beta_var = beta_var,
    Q = stats$Q,
    Q_post = stats$Q_post,
    Q_joint = Q_joint,
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
  log_marginal <- .exact_matern_loglik_fixed_beta(stats, beta0 = beta0, n_obs = length(x))

  .exact_matern_posterior_from_stats(
    x = x,
    s = s,
    A = A,
    stats = stats,
    beta0 = beta0,
    log_marginal = log_marginal
  )
}

.exact_matern_state_unknown_noise <- function(x,
                                              A,
                                              spde_template,
                                              alpha,
                                              d,
                                              log_range,
                                              log_sigma,
                                              beta0,
                                              log_noise_sd,
                                              noise_scale = NULL) {
  if (!is.finite(log_noise_sd)) {
    stop("`log_noise_sd` must be finite.")
  }

  noise_sd <- exp(log_noise_sd)
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    stop("The learned observation noise SD must be positive.")
  }

  s_eff <- .exact_matern_unknown_noise_s(noise_sd, length(x), noise_scale)
  state <- .exact_matern_state(
    x = x,
    s = s_eff,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = log_range,
    log_sigma = log_sigma,
    beta0 = beta0
  )
  state$fitted_noise_sd <- as.numeric(noise_sd)
  state$fitted_s <- s_eff
  if (!is.null(noise_scale)) state$noise_scale <- as.numeric(noise_scale)
  state
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

.resolve_matern_pc_penalty <- function(pc.penalty,
                                       penalty_range0,
                                       sigma_anchor0,
                                       noise_anchor0 = NULL,
                                       allow_noise = FALSE) {
  if (is.null(pc.penalty)) return(NULL)

  if (!is.list(pc.penalty)) {
    stop("`pc.penalty` must be NULL or a list with optional named entries.")
  }

  nms <- names(pc.penalty)
  if (is.null(nms)) nms <- rep("", length(pc.penalty))
  if (any(nms == "")) {
    stop("`pc.penalty` entries must be named.")
  }

  supported_names <- c("range", "sigma")
  if (isTRUE(allow_noise)) supported_names <- c(supported_names, "noise")
  if (!all(nms %in% supported_names)) {
    stop(
      "`pc.penalty` only supports ",
      paste(sprintf("`%s`", supported_names), collapse = " and "),
      " entries."
    )
  }

  out <- list(
    range = .resolve_pc_penalty_component(pc.penalty$range, penalty_range0, "pc.penalty$range"),
    sigma = .resolve_pc_penalty_component(pc.penalty$sigma, sigma_anchor0, "pc.penalty$sigma")
  )

  if (isTRUE(allow_noise)) {
    if (is.null(noise_anchor0)) {
      stop("`noise_anchor0` must be supplied when `allow_noise = TRUE`.")
    }
    out$noise <- .resolve_pc_penalty_component(pc.penalty$noise, noise_anchor0, "pc.penalty$noise")
  }

  out
}

.matern_centered_scale <- function(x, center, s = NULL, weights = NULL) {
  x <- as.numeric(x)
  center <- as.numeric(center)[1]
  resid <- x - center

  out <- if (!is.null(weights)) {
    w <- as.numeric(weights)
    if (length(w) != length(x) || anyNA(w) || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
      NA_real_
    } else {
      sqrt(sum(w * resid^2) / sum(w))
    }
  } else if (is.null(s)) {
    stats::sd(resid)
  } else {
    w <- 1 / (as.numeric(s)^2)
    sqrt(sum(w * resid^2) / sum(w))
  }

  if (!is.finite(out) || out <= 0) {
    out <- stats::sd(x)
  }
  if (!is.finite(out) || out <= 0) {
    out <- 1
  }
  as.numeric(out)
}

.matern_log_scale_observations <- function(x, s = NULL) {
  x <- as.numeric(x)
  floor0 <- .positive_response_floor(x, s = s)
  log(pmax(x, floor0))
}

.matern_link_scale_observations <- function(x, s = NULL, link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  x <- as.numeric(x)
  if (identical(link, "identity")) {
    return(list(eta = x, deriv = rep(1, length(x))))
  }
  x_floored <- pmax(x, .positive_response_floor(x, s = s))
  if (identical(link, "log")) {
    eta <- log(x_floored)
    deriv <- x_floored
  } else {
    eta <- .inverse_softplus_stable(x_floored)
    deriv <- .sigmoid_stable(eta)
  }
  list(eta = as.numeric(eta), deriv = as.numeric(deriv), response = as.numeric(x_floored))
}

.resolve_matern_beta_init <- function(x,
                                      s = NULL,
                                      g_init = NULL,
                                      beta_fixed = NULL,
                                      link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  if (!is.null(beta_fixed)) {
    return(.check_single_numeric(beta_fixed, "beta_fixed"))
  }
  if (!is.null(g_init) && !is.null(g_init$beta)) {
    return(.check_optional_beta_vector(
      g_init$beta,
      "g_init$beta",
      expected_length = 1L,
      allow_null = FALSE
    ))
  }
  .response_beta_init(x, s = s, link = link)
}

.resolve_matern_g_init <- function(x,
                                   s = NULL,
                                   g_init = NULL,
                                   beta_fixed = NULL,
                                   beta_prec = NULL,
                                   penalty_range0,
                                   pc.penalty = NULL,
                                   allow_noise = FALSE,
                                   link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  beta_init <- .resolve_matern_beta_init(
    x = x,
    s = s,
    g_init = g_init,
    beta_fixed = beta_fixed,
    link = link
  )
  link_scale <- .matern_link_scale_observations(x, s = s, link = link)
  sigma_weights <- if (is.null(s)) NULL else (link_scale$deriv / as.numeric(s))^2
  sigma_data <- .matern_centered_scale(
    x = link_scale$eta,
    center = beta_init,
    weights = sigma_weights
  )
  noise_data <- .response_raw_residual_scale(x, eta = beta_init, link = link)
  pc_penalty0 <- .resolve_matern_pc_penalty(
    pc.penalty = pc.penalty,
    penalty_range0 = penalty_range0,
    sigma_anchor0 = sigma_data,
    noise_anchor0 = noise_data,
    allow_noise = allow_noise
  )

  theta_anchor <- if (is.null(pc_penalty0)) {
    penalty_range0
  } else {
    as.numeric(pc_penalty0$range["anchor"])
  }
  sigma_anchor <- if (is.null(pc_penalty0)) {
    sigma_data
  } else {
    as.numeric(pc_penalty0$sigma["anchor"])
  }

  theta_init <- if (!is.null(g_init) && !is.null(g_init$theta)) {
    .check_single_numeric(g_init$theta, "g_init$theta")
  } else {
    log(theta_anchor)
  }
  sigma_init <- if (!is.null(g_init) && !is.null(g_init$sigma)) {
    .check_single_numeric(g_init$sigma, "g_init$sigma")
  } else {
    sigma_anchor
  }
  if (!is.finite(sigma_init) || sigma_init <= 0) {
    sigma_init <- sigma_data
  }

  noise_sd_init <- noise_data
  if (!is.null(pc_penalty0) && isTRUE(allow_noise) && !is.null(pc_penalty0$noise)) {
    noise_sd_init <- as.numeric(pc_penalty0$noise["anchor"])
  }
  if (!is.finite(noise_sd_init) || noise_sd_init <= 0) {
    noise_sd_init <- noise_data
  }

  list(
    theta_init = theta_init,
    sigma_init = sigma_init,
    beta_init = as.numeric(beta_init),
    noise_sd_init = as.numeric(noise_sd_init),
    sigma_data = as.numeric(sigma_data),
    pc_penalty = pc_penalty0,
    g_init = Matern(
      theta = theta_init,
      sigma = sigma_init,
      beta = beta_init,
      beta_prec = beta_prec
    )
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

.log_pc_prior_noise_internal <- function(log_noise_sd, noise_spec) {
  if (is.null(noise_spec)) {
    return(0)
  }
  sigma <- exp(log_noise_sd)
  alpha_s <- as.numeric(noise_spec["alpha"])
  sigma0 <- as.numeric(noise_spec["anchor"])
  S <- -log(alpha_s) / sigma0
  as.numeric(log(S) + log_noise_sd - S * sigma)
}

.quiet_inla_num_threads <- function() {
  cores <- parallel::detectCores(all.tests = FALSE, logical = FALSE)
  if (!is.finite(cores) || is.na(cores) || cores < 1) {
    cores <- 1L
  }
  paste0(max(1L, min(16L, as.integer(cores))), ":1")
}

.quiet_inla_default_options <- function() {
  list(
    inla.arg = NULL,
    fmesher.arg = "",
    num.threads = .quiet_inla_num_threads(),
    smtp = "default",
    safe = TRUE,
    keep = FALSE,
    verbose = FALSE,
    save.memory = FALSE,
    internal.opt = TRUE,
    working.directory = NULL,
    silent = TRUE,
    debug = FALSE,
    show.warning.graph.file = TRUE,
    scale.model.default = FALSE,
    short.summary = FALSE,
    inla.timeout = 0,
    fmesher.timeout = 0,
    inla.mode = "compact",
    malloc.lib = "mi",
    fmesher.evolution = 2L,
    fmesher.evolution.warn = TRUE,
    fmesher.evolution.verbosity = "default",
    INLAjoint.features = FALSE,
    numa = FALSE
  )
}

.with_quiet_inla_defaults <- function(expr) {
  ns <- asNamespace("INLA")
  default_name <- "inla.getOption.default"
  old_default <- get(default_name, envir = ns)

  unlockBinding(default_name, ns)
  assign(default_name, .quiet_inla_default_options, envir = ns)
  lockBinding(default_name, ns)

  on.exit({
    unlockBinding(default_name, ns)
    assign(default_name, old_default, envir = ns)
    lockBinding(default_name, ns)
  }, add = TRUE)

  force(expr)
}

.build_matern_pc_spde <- function(mesh, alpha, pc_penalty, suppress_warnings = TRUE) {
  .with_quiet_inla_defaults({
    build <- function() {
      INLA::inla.spde2.pcmatern(
        mesh = mesh,
        alpha = alpha,
        prior.range = unname(pc_penalty$range[c("anchor", "alpha")]),
        prior.sigma = unname(pc_penalty$sigma[c("anchor", "alpha")])
      )
    }

    if (suppress_warnings) suppressWarnings(build()) else build()
  })
}

.build_matern_inla_spde <- function(mesh, alpha, pc_penalty = NULL, suppress_warnings = TRUE) {
  if (!is.null(pc_penalty)) {
    return(.build_matern_pc_spde(
      mesh = mesh,
      alpha = alpha,
      pc_penalty = pc_penalty,
      suppress_warnings = suppress_warnings
    ))
  }

  .with_quiet_inla_defaults({
    build <- function() {
      INLA::inla.spde2.matern(
        mesh = mesh,
        alpha = alpha,
        theta.prior.prec = 0
      )
    }

    if (suppress_warnings) suppressWarnings(build()) else build()
  })
}

.matern_inla_theta_init <- function(log_range, log_sigma, alpha, d, pc_penalty = NULL) {
  if (!is.null(pc_penalty)) {
    return(c(log_range, log_sigma))
  }
  .matern_spde_theta_from_log_range_log_sigma(
    log_range = log_range,
    log_sigma = log_sigma,
    alpha = alpha,
    d = d
  )
}

.matern_inla_theta_to_matern <- function(theta, alpha, d, pc_penalty = NULL) {
  theta <- as.numeric(theta)
  if (length(theta) != 2L) {
    stop("The INLA Matern fit did not return the expected two hyperparameters.")
  }
  if (!is.null(pc_penalty)) {
    return(list(log_range = theta[1], sigma = exp(theta[2])))
  }
  rs <- .matern_range_sigma_from_tau_kappa(
    tau = exp(theta[1]),
    kappa = exp(theta[2]),
    alpha = alpha,
    d = d
  )
  list(log_range = log(rs$range), sigma = rs$sigma)
}

.matern_pc_noise_hyper <- function(noise_spec, noise_sd_init) {
  noise_sd_init <- .check_single_numeric(noise_sd_init, "noise_sd_init")
  if (noise_sd_init <= 0) {
    stop("`noise_sd_init` must be positive.")
  }

  if (is.null(noise_spec)) {
    return(list(
      prior = "logflat",
      param = numeric(0),
      initial = log(1 / noise_sd_init^2),
      fixed = FALSE
    ))
  }

  list(
    prior = "pc.prec",
    param = unname(noise_spec[c("anchor", "alpha")]),
    initial = log(1 / noise_sd_init^2),
    fixed = FALSE
  )
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

.matern_inla_laplace_initial_mode <- function(fit, beta_mode) {
  field <- fit$posterior_spatial_field
  if (is.null(field)) return(NULL)

  mode <- if (!is.null(field$mode)) {
    field$mode
  } else if (!is.null(field$mean)) {
    field$mean
  } else {
    NULL
  }
  if (is.null(mode)) return(NULL)

  mode <- as.numeric(mode)
  if (!length(mode) || any(!is.finite(mode))) return(NULL)

  if (beta_mode %in% c("prior_flat", "prior_proper")) {
    beta <- fit$fitted_beta
    if (is.null(beta) || length(beta) != 1L || !is.finite(beta)) return(NULL)
    mode <- c(mode, as.numeric(beta))
  }

  mode
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
  if (!is.null(fit$link) && !identical(as.character(fit$link), "identity")) {
    stop("`matern_objective_breakdown()` is only available for identity-link Matern fits.")
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
#'   \code{"auto"}, \code{"fixed"}, \code{"empirical_bayes"},
#'   \code{"prior_flat"}, and \code{"prior_proper"}.
#'
#' @return A named list containing the fitted hyperparameters, exact
#'   log-likelihood values under fixed / empirical-Bayes / flat-prior /
#'   proper-prior beta handling, the PC-prior term when relevant, the
#'   corresponding exact penalized objectives, and the gap between the recorded
#'   \code{fit$log_likelihood} and the selected exact comparator.
#'
#' @export
matern_objective_breakdown <- function(fit,
                                       beta_mode = c("auto", "fixed", "empirical_bayes", "prior_flat", "prior_proper")) {
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
  beta_prec <- if (is.null(fit$fitted_g$beta_prec)) NULL else as.numeric(fit$fitted_g$beta_prec)
  loglik_prior_beta <- if (is.null(beta_prec) || beta_prec <= 0) {
    NULL
  } else {
    .exact_matern_loglik_prior_beta(stats, beta_prec = beta_prec, n_obs = length(x))
  }

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
  objective_prior_proper_plus_prior <- if (is.null(loglik_prior_beta)) {
    NULL
  } else if (is.null(log_pc_prior_theta)) {
    loglik_prior_beta
  } else {
    loglik_prior_beta + log_pc_prior_theta
  }

  fit_beta_mode <- fit$beta_mode
  if (is.null(fit_beta_mode)) {
    fit_beta_mode <- if (identical(fit$backend, "exact")) {
      "empirical_bayes"
    } else if (!is.null(fit$log_likelihood_stepB)) {
      "fixed"
    } else {
      "prior_flat"
    }
  }
  fit_beta_mode <- .eb_smoother_public_beta_mode("matern", fit_beta_mode)

  matched_beta_mode <- if (beta_mode == "auto") fit_beta_mode else beta_mode
  matched_exact_objective <- switch(
    matched_beta_mode,
    fixed = objective_fixed_plus_prior,
    empirical_bayes = objective_profile_plus_prior,
    prior_flat = objective_integrated_flat_plus_prior,
    prior_proper = objective_prior_proper_plus_prior,
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
    beta_prec = beta_prec,
    log_pc_prior_theta = log_pc_prior_theta,
    loglik_fixed_beta = loglik_fixed_beta,
    loglik_empirical_bayes_beta = loglik_profile_beta,
    loglik_prior_flat_beta = loglik_integrated_flat_beta,
    loglik_prior_proper_beta = loglik_prior_beta,
    objective_fixed_plus_prior = objective_fixed_plus_prior,
    objective_empirical_bayes_plus_prior = objective_profile_plus_prior,
    objective_prior_flat_plus_prior = objective_integrated_flat_plus_prior,
    objective_prior_proper_plus_prior = objective_prior_proper_plus_prior,
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
                                   d,
                                   theta_init,
                                   beta_prec = 0,
                                   pc_penalty,
                                   link = c("identity", "log", "softplus"),
                                   suppress_warnings = TRUE) {
  link <- match.arg(link)
  .with_quiet_inla_defaults({
    n <- length(x)
    spde <- .build_matern_inla_spde(mesh, alpha, pc_penalty, suppress_warnings = suppress_warnings)
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
    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    if (is.null(beta_prec)) {
      beta_prec <- 0
    }

    runA <- function() {
      INLA::inla(
        formulaA,
        scale = (1 / s^2),
        control.inla = list(int.strategy = "eb", strategy = "gaussian"),
        control.family = list(
          control.link = list(model = link),
          hyper = list(prec = list(fixed = TRUE, initial = 0))
        ),
        control.fixed = INLA::control.fixed(prec = beta_prec),
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
    matern_theta <- .matern_inla_theta_to_matern(
      theta = theta_hat,
      alpha = alpha,
      d = d,
      pc_penalty = pc_penalty
    )

    list(
      result = resA,
      posterior = post$posterior,
      posterior_spatial_field = post$posterior_spatial_field,
      fitted_g = Matern(theta = matern_theta$log_range, sigma = matern_theta$sigma, beta = as.numeric(resA$summary.fixed$mean[1]), beta_prec = beta_prec),
      fitted_beta = as.numeric(resA$summary.fixed$mean[1]),
      log_likelihood_stepA_penalized = as.numeric(resA$misc$log.posterior.mode),
      log_likelihood_stepA_mlik_integration = as.numeric(resA$mlik["log marginal-likelihood (integration)", 1]),
      log_likelihood_stepA_mlik_gaussian = as.numeric(resA$mlik["log marginal-likelihood (Gaussian)", 1]),
      log_likelihood_stepA_joint_log_posterior = as.numeric(resA$joint.hyper[, "Log posterior density"])
    )
  })
}

.fit_matern_inla_stepA_fixed_beta <- function(x,
                                              s,
                                              A,
                                              mesh,
                                              alpha,
                                              d,
                                              theta_init,
                                              beta_fixed,
                                              pc_penalty,
                                              link = c("identity", "log", "softplus"),
                                              suppress_warnings = TRUE) {
  link <- match.arg(link)
  .with_quiet_inla_defaults({
    n <- length(x)
    beta_fixed <- .check_single_numeric(beta_fixed, "beta_fixed")
    spde <- .build_matern_inla_spde(mesh, alpha, pc_penalty, suppress_warnings = suppress_warnings)
    idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

    stackA <- INLA::inla.stack(
      data = list(Y = as.numeric(x)),
      A = list(A, matrix(beta_fixed, nrow = n, ncol = 1)),
      effects = list(
        spatial.field = idx$spatial.field,
        beta0 = 1
      ),
      tag = "est"
    )

    formulaA <- Y ~ 0 + offset(beta0) + f(spatial.field, model = spde)

    runA <- function() {
      INLA::inla(
        formulaA,
        scale = (1 / s^2),
        control.inla = list(int.strategy = "eb", strategy = "gaussian"),
        control.family = list(
          control.link = list(model = link),
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
      stop("The INLA fixed-beta Step A fit did not return the expected two Matern hyperparameters.")
    }
    matern_theta <- .matern_inla_theta_to_matern(
      theta = theta_hat,
      alpha = alpha,
      d = d,
      pc_penalty = pc_penalty
    )

    list(
      result = resA,
      posterior = post$posterior,
      posterior_spatial_field = post$posterior_spatial_field,
      fitted_g = Matern(theta = matern_theta$log_range, sigma = matern_theta$sigma, beta = beta_fixed, beta_prec = NULL),
      fitted_beta = beta_fixed,
      log_likelihood_stepA_penalized = as.numeric(resA$misc$log.posterior.mode),
      log_likelihood_stepA_mlik_integration = as.numeric(resA$mlik["log marginal-likelihood (integration)", 1]),
      log_likelihood_stepA_mlik_gaussian = as.numeric(resA$mlik["log marginal-likelihood (Gaussian)", 1]),
      log_likelihood_stepA_joint_log_posterior = as.numeric(resA$joint.hyper[, "Log posterior density"])
    )
  })
}

.fit_matern_inla_stepA_unknown_noise <- function(x,
                                                 A,
                                                 mesh,
                                                 alpha,
                                                 d,
                                                 theta_init,
                                                 noise_sd_init,
                                                 beta_prec = 0,
                                                 pc_penalty,
                                                 link = c("identity", "log", "softplus"),
                                                 suppress_warnings = TRUE) {
  link <- match.arg(link)
  .with_quiet_inla_defaults({
    n <- length(x)
    spde <- .build_matern_inla_spde(mesh, alpha, pc_penalty, suppress_warnings = suppress_warnings)
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
    noise_hyper <- .matern_pc_noise_hyper(
      if (is.null(pc_penalty)) NULL else pc_penalty$noise,
      noise_sd_init = noise_sd_init
    )
    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    if (is.null(beta_prec)) {
      beta_prec <- 0
    }

    runA <- function() {
      INLA::inla(
        formulaA,
        scale = rep(1, n),
        control.inla = list(int.strategy = "eb", strategy = "gaussian"),
        control.family = list(
          control.link = list(model = link),
          hyper = list(prec = noise_hyper)
        ),
        control.fixed = INLA::control.fixed(prec = beta_prec),
        control.mode = INLA::control.mode(
          theta = c(noise_hyper$initial, theta_init),
          restart = TRUE
        ),
        data = INLA::inla.stack.data(stackA),
        control.predictor = list(A = INLA::inla.stack.A(stackA), link = 1, compute = TRUE),
        control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
        silent = TRUE
      )
    }

    resA <- if (suppress_warnings) suppressWarnings(runA()) else runA()
    post <- .matern_posterior_from_inla_result(resA, stackA)
    theta_hat <- as.numeric(resA$mode$theta)

    if (length(theta_hat) != 3L) {
      stop("The INLA learned-noise Matern fit did not return the expected three hyperparameters.")
    }
    matern_theta <- .matern_inla_theta_to_matern(
      theta = theta_hat[2:3],
      alpha = alpha,
      d = d,
      pc_penalty = pc_penalty
    )

    list(
      result = resA,
      posterior = post$posterior,
      posterior_spatial_field = post$posterior_spatial_field,
      fitted_g = Matern(theta = matern_theta$log_range, sigma = matern_theta$sigma, beta = as.numeric(resA$summary.fixed$mean[1]), beta_prec = beta_prec),
      fitted_beta = as.numeric(resA$summary.fixed$mean[1]),
      fitted_noise_sd = as.numeric(exp(-0.5 * theta_hat[1])),
      fitted_log_precision = as.numeric(theta_hat[1]),
      log_likelihood_stepA_penalized = as.numeric(resA$misc$log.posterior.mode),
      log_likelihood_stepA_mlik_integration = as.numeric(resA$mlik["log marginal-likelihood (integration)", 1]),
      log_likelihood_stepA_mlik_gaussian = as.numeric(resA$mlik["log marginal-likelihood (Gaussian)", 1]),
      log_likelihood_stepA_joint_log_posterior = as.numeric(resA$joint.hyper[, "Log posterior density"])
    )
  })
}

.fit_matern_inla_stepA_unknown_noise_fixed_beta <- function(x,
                                                            A,
                                                            mesh,
                                                            alpha,
                                                            d,
                                                            theta_init,
                                                            noise_sd_init,
                                                            beta_fixed,
                                                            pc_penalty,
                                                            link = c("identity", "log", "softplus"),
                                                            suppress_warnings = TRUE) {
  link <- match.arg(link)
  .with_quiet_inla_defaults({
    n <- length(x)
    beta_fixed <- .check_single_numeric(beta_fixed, "beta_fixed")
    spde <- .build_matern_inla_spde(mesh, alpha, pc_penalty, suppress_warnings = suppress_warnings)
    idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

    stackA <- INLA::inla.stack(
      data = list(Y = as.numeric(x)),
      A = list(A, matrix(beta_fixed, nrow = n, ncol = 1)),
      effects = list(
        spatial.field = idx$spatial.field,
        beta0 = 1
      ),
      tag = "est"
    )

    formulaA <- Y ~ 0 + offset(beta0) + f(spatial.field, model = spde)
    noise_hyper <- .matern_pc_noise_hyper(
      if (is.null(pc_penalty)) NULL else pc_penalty$noise,
      noise_sd_init = noise_sd_init
    )

    runA <- function() {
      INLA::inla(
        formulaA,
        scale = rep(1, n),
        control.inla = list(int.strategy = "eb", strategy = "gaussian"),
        control.family = list(
          control.link = list(model = link),
          hyper = list(prec = noise_hyper)
        ),
        control.fixed = INLA::control.fixed(prec = 0),
        control.mode = INLA::control.mode(
          theta = c(noise_hyper$initial, theta_init),
          restart = TRUE
        ),
        data = INLA::inla.stack.data(stackA),
        control.predictor = list(A = INLA::inla.stack.A(stackA), link = 1, compute = TRUE),
        control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
        silent = TRUE
      )
    }

    resA <- if (suppress_warnings) suppressWarnings(runA()) else runA()
    post <- .matern_posterior_from_inla_result(resA, stackA)
    theta_hat <- as.numeric(resA$mode$theta)

    if (length(theta_hat) != 3L) {
      stop("The INLA learned-noise fixed-beta Step A fit did not return the expected three hyperparameters.")
    }
    matern_theta <- .matern_inla_theta_to_matern(
      theta = theta_hat[2:3],
      alpha = alpha,
      d = d,
      pc_penalty = pc_penalty
    )

    list(
      result = resA,
      posterior = post$posterior,
      posterior_spatial_field = post$posterior_spatial_field,
      fitted_g = Matern(theta = matern_theta$log_range, sigma = matern_theta$sigma, beta = beta_fixed, beta_prec = NULL),
      fitted_beta = beta_fixed,
      fitted_noise_sd = as.numeric(exp(-0.5 * theta_hat[1])),
      fitted_log_precision = as.numeric(theta_hat[1]),
      log_likelihood_stepA_penalized = as.numeric(resA$misc$log.posterior.mode),
      log_likelihood_stepA_mlik_integration = as.numeric(resA$mlik["log marginal-likelihood (integration)", 1]),
      log_likelihood_stepA_mlik_gaussian = as.numeric(resA$mlik["log marginal-likelihood (Gaussian)", 1]),
      log_likelihood_stepA_joint_log_posterior = as.numeric(resA$joint.hyper[, "Log posterior density"])
    )
  })
}

.fit_matern_inla_profile_beta <- function(beta_init,
                                          fit_at_beta,
                                          objective_at_fit,
                                          failure_label,
                                          beta_search_radius = 2,
                                          beta_step = 0.1) {
  beta_init <- .check_single_numeric(beta_init, "beta_init")
  beta_search_radius <- .check_single_numeric(beta_search_radius, "beta_search_radius")
  if (!is.finite(beta_search_radius) || beta_search_radius <= 0) {
    stop("`beta_search_radius` must be positive.")
  }
  beta_step <- .check_single_numeric(beta_step, "beta_step")
  if (!is.finite(beta_step) || beta_step <= 0) {
    stop("`beta_step` must be positive.")
  }
  beta_step <- min(beta_step, beta_search_radius / 2)
  cache <- new.env(parent = emptyenv())

  eval_profile <- function(beta_value) {
    beta_value <- .check_single_numeric(beta_value, "beta_value")
    key <- format(signif(beta_value, 15), scientific = TRUE)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    fit <- fit_at_beta(beta_value)
    objective_value <- objective_at_fit(fit, beta_value)
    if (!is.finite(objective_value)) {
      stop(failure_label, " produced a non-finite objective.")
    }

    fit$fitted_beta <- beta_value
    fit$fitted_g$beta <- beta_value
    fit$fitted_g$beta_prec <- NULL
    fit$beta_profile_objective <- as.numeric(objective_value)
    out <- list(value = as.numeric(objective_value), fit = fit)
    assign(key, out, envir = cache)
    out
  }

  candidate_betas <- c(
    beta_init,
    beta_init - beta_step,
    beta_init + beta_step
  )
  candidates <- lapply(candidate_betas, function(beta_value) {
    tryCatch(eval_profile(beta_value), error = function(e) NULL)
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) {
    stop(failure_label, ": no finite beta-profile evaluations.")
  }

  values <- vapply(candidates, `[[`, numeric(1), "value")
  if (length(candidates) == 3L) {
    denom <- values[2] - 2 * values[1] + values[3]
    if (is.finite(denom) && denom < 0) {
      vertex <- beta_init + beta_step * (values[2] - values[3]) / (2 * denom)
      if (is.finite(vertex) && abs(vertex - beta_init) <= beta_search_radius) {
        vertex_eval <- tryCatch(eval_profile(vertex), error = function(e) NULL)
        if (!is.null(vertex_eval)) {
          candidates[[length(candidates) + 1L]] <- vertex_eval
          values <- c(values, vertex_eval$value)
        }
      }
    }
  }

  best <- candidates[[which.max(values)]]
  fit <- best$fit
  fit$beta_profile_optimization <- list(
    beta = as.numeric(fit$fitted_beta),
    objective = as.numeric(best$value),
    method = "local_quadratic_profile"
  )
  fit
}

.fit_matern_inla_stepA_profile_beta <- function(x,
                                                s,
                                                A,
                                                mesh,
                                                spde_template,
                                                alpha,
                                                d,
                                                theta_init,
                                                beta_init,
                                                pc_penalty,
                                                link = c("identity", "log", "softplus"),
                                                suppress_warnings = TRUE) {
  link <- match.arg(link)

  fit_at_beta <- function(beta_value) {
    .fit_matern_inla_stepA_fixed_beta(
      x = x,
      s = s,
      A = A,
      mesh = mesh,
      alpha = alpha,
      d = d,
      theta_init = theta_init,
      beta_fixed = beta_value,
      pc_penalty = pc_penalty,
      link = link,
      suppress_warnings = suppress_warnings
    )
  }

  objective_at_fit <- function(fit, beta_value) {
    if (identical(link, "log")) {
      laplace_objective <- .matern_laplace_known_noise_objective_at_params(
        x = x,
        s = s,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = fit$fitted_g$theta,
        log_sigma = log(fit$fitted_g$sigma),
        beta_mode = "empirical_bayes",
        beta_init = beta_value,
        link = link,
        pc_penalty = pc_penalty,
        initial_mode = .matern_inla_laplace_initial_mode(fit, "empirical_bayes"),
        compute_posterior = FALSE,
        optimize_beta = FALSE,
        suppress_warnings = suppress_warnings
      )
      return(as.numeric(laplace_objective$log_marginal))
    }
    exact_objective <- .exact_matern_known_noise_objective_at_params(
      x = x,
      s = s,
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      log_range = fit$fitted_g$theta,
      log_sigma = log(fit$fitted_g$sigma),
      beta_mode = "fixed",
      beta_fixed = beta_value,
      pc_penalty = pc_penalty,
      suppress_warnings = suppress_warnings
    )
    as.numeric(exact_objective$log_marginal)
  }

  .fit_matern_inla_profile_beta(
    beta_init = beta_init,
    fit_at_beta = fit_at_beta,
    objective_at_fit = objective_at_fit,
    failure_label = "INLA Matern beta profiling failed"
  )
}

.fit_matern_inla_stepA_unknown_noise_profile_beta <- function(x,
                                                              A,
                                                              mesh,
                                                              spde_template,
                                                              alpha,
                                                              d,
                                                              theta_init,
                                                              noise_sd_init,
                                                              beta_init,
                                                              pc_penalty,
                                                              link = c("identity", "log", "softplus"),
                                                              suppress_warnings = TRUE) {
  link <- match.arg(link)
  fit_at_beta <- function(beta_value) {
    .fit_matern_inla_stepA_unknown_noise_fixed_beta(
      x = x,
      A = A,
      mesh = mesh,
      alpha = alpha,
      d = d,
      theta_init = theta_init,
      noise_sd_init = noise_sd_init,
      beta_fixed = beta_value,
      pc_penalty = pc_penalty,
      link = link,
      suppress_warnings = suppress_warnings
    )
  }

  objective_at_fit <- function(fit, beta_value) {
    if (identical(link, "log")) {
      laplace_objective <- .matern_laplace_unknown_noise_objective_at_params(
        x = x,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = fit$fitted_g$theta,
        log_sigma = log(fit$fitted_g$sigma),
        log_noise_sd = log(fit$fitted_noise_sd),
        beta_mode = "empirical_bayes",
        beta_init = beta_value,
        link = link,
        pc_penalty = pc_penalty,
        initial_mode = .matern_inla_laplace_initial_mode(fit, "empirical_bayes"),
        compute_posterior = FALSE,
        optimize_beta = FALSE,
        suppress_warnings = suppress_warnings
      )
      return(as.numeric(laplace_objective$log_marginal))
    }
    exact_objective <- .exact_matern_unknown_noise_objective_at_params(
      x = x,
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      log_range = fit$fitted_g$theta,
      log_sigma = log(fit$fitted_g$sigma),
      log_noise_sd = log(fit$fitted_noise_sd),
      beta_mode = "fixed",
      beta_fixed = beta_value,
      pc_penalty = pc_penalty,
      suppress_warnings = suppress_warnings
    )
    as.numeric(exact_objective$log_marginal)
  }

  .fit_matern_inla_profile_beta(
    beta_init = beta_init,
    fit_at_beta = fit_at_beta,
    objective_at_fit = objective_at_fit,
    failure_label = "INLA learned-noise Matern beta profiling failed"
  )
}

.fit_matern_inla_stepA_unknown_noise_log_profile_noise <- function(x,
                                                                   A,
                                                                   mesh,
                                                                   spde_template,
                                                                   alpha,
                                                                   d,
                                                                   theta_init,
                                                                   noise_sd_init,
                                                                   beta_mode,
                                                                   beta_fixed = NULL,
                                                                   beta_prec = NULL,
                                                                   beta_init,
                                                                   pc_penalty,
                                                                   suppress_warnings = TRUE) {
  noise_sd_init <- .check_single_numeric(noise_sd_init, "noise_sd_init")
  if (!is.finite(noise_sd_init) || noise_sd_init <= 0) {
    stop("`noise_sd_init` must be positive.")
  }

  cache <- new.env(parent = emptyenv())
  eval_at_log_noise <- function(log_noise_sd) {
    log_noise_sd <- .check_single_numeric(log_noise_sd, "log_noise_sd")
    key <- format(signif(log_noise_sd, 15), scientific = TRUE)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    noise_sd <- exp(log_noise_sd)
    s <- rep(noise_sd, length(x))
    fit <- if (identical(beta_mode, "fixed")) {
      .fit_matern_inla_stepA_fixed_beta(
        x = x,
        s = s,
        A = A,
        mesh = mesh,
        alpha = alpha,
        d = d,
        theta_init = theta_init,
        beta_fixed = beta_fixed,
        pc_penalty = pc_penalty,
        link = "log",
        suppress_warnings = suppress_warnings
      )
    } else if (identical(beta_mode, "empirical_bayes")) {
      .fit_matern_inla_stepA_profile_beta(
        x = x,
        s = s,
        A = A,
        mesh = mesh,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_init,
        beta_init = beta_init,
        pc_penalty = pc_penalty,
        link = "log",
        suppress_warnings = suppress_warnings
      )
    } else {
      .fit_matern_inla_stepA(
        x = x,
        s = s,
        A = A,
        mesh = mesh,
        alpha = alpha,
        d = d,
        theta_init = theta_init,
        beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec,
        pc_penalty = pc_penalty,
        link = "log",
        suppress_warnings = suppress_warnings
      )
    }

    laplace_objective <- .matern_laplace_unknown_noise_objective_at_params(
      x = x,
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      log_range = fit$fitted_g$theta,
      log_sigma = log(fit$fitted_g$sigma),
      log_noise_sd = log_noise_sd,
      beta_mode = beta_mode,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      beta_init = fit$fitted_beta,
      link = "log",
      pc_penalty = pc_penalty,
      initial_mode = .matern_inla_laplace_initial_mode(fit, beta_mode),
      compute_posterior = FALSE,
      optimize_beta = FALSE,
      suppress_warnings = suppress_warnings
    )
    fit$fitted_noise_sd <- as.numeric(noise_sd)
    fit$log_likelihood_laplace_at_inla_mode <- as.numeric(laplace_objective$log_marginal)
    out <- list(value = as.numeric(laplace_objective$log_marginal), fit = fit)
    assign(key, out, envir = cache)
    out
  }

  log_noise_init <- log(noise_sd_init)
  log_noise_step <- 0.5
  candidate_logs <- c(log_noise_init, log_noise_init - log_noise_step, log_noise_init + log_noise_step)
  candidates <- lapply(candidate_logs, function(log_noise_sd) {
    tryCatch(eval_at_log_noise(log_noise_sd), error = function(e) NULL)
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) {
    stop("INLA learned-noise Matern noise profiling failed: no finite noise-profile evaluations.")
  }

  values <- vapply(candidates, `[[`, numeric(1), "value")
  if (length(candidates) == 3L) {
    denom <- values[2] - 2 * values[1] + values[3]
    if (is.finite(denom) && denom < 0) {
      vertex <- log_noise_init + log_noise_step * (values[2] - values[3]) / (2 * denom)
      if (is.finite(vertex) && abs(vertex - log_noise_init) <= 2 * log_noise_step) {
        vertex_eval <- tryCatch(eval_at_log_noise(vertex), error = function(e) NULL)
        if (!is.null(vertex_eval)) {
          candidates[[length(candidates) + 1L]] <- vertex_eval
          values <- c(values, vertex_eval$value)
        }
      }
    }
  }

  best <- candidates[[which.max(values)]]
  fit <- best$fit
  fit$noise_profile_optimization <- list(
    log_noise_sd = log(fit$fitted_noise_sd),
    noise_sd = as.numeric(fit$fitted_noise_sd),
    objective = as.numeric(best$value),
    method = "local_quadratic_profile"
  )
  fit
}

.matern_inlabru_assert_installed <- function() {
  if (!requireNamespace("inlabru", quietly = TRUE)) {
    stop(
      "`backend = \"inlabru\"` requires the `inlabru` package; install it with ",
      "`install.packages(\"inlabru\")`."
    )
  }
}

.matern_inlabru_assert_softplus <- function(link) {
  if (!identical(link, "softplus")) {
    stop("`backend = \"inlabru\"` is currently only supported for `link = \"softplus\"`.")
  }
}

.matern_inlabru_pc_penalty_policy <- function(pc_penalty_arg,
                                              resolved_pc_penalty,
                                              learn_noise = FALSE) {
  if (is.null(pc_penalty_arg)) {
    if (isTRUE(learn_noise)) {
      stop(
        "`backend = \"inlabru\"` with `s = NULL` requires an explicit ",
        "`pc.penalty` list with `range`, `sigma`, and `noise` entries."
      )
    }
    return(NULL)
  }

  if (!is.list(pc_penalty_arg) || is.null(names(pc_penalty_arg))) {
    stop("`pc.penalty` must be a named list for `backend = \"inlabru\"`.")
  }
  required <- c("range", "sigma", if (isTRUE(learn_noise)) "noise")
  missing <- setdiff(required, names(pc_penalty_arg))
  if (length(missing)) {
    stop(
      "`backend = \"inlabru\"` requires explicit `pc.penalty` entries: ",
      paste(sprintf("`%s`", required), collapse = ", "),
      ". Missing: ",
      paste(sprintf("`%s`", missing), collapse = ", "),
      "."
    )
  }
  null_required <- required[vapply(pc_penalty_arg[required], is.null, logical(1))]
  if (length(null_required)) {
    stop(
      "`backend = \"inlabru\"` requires non-NULL `pc.penalty` entries: ",
      paste(sprintf("`%s`", null_required), collapse = ", "),
      "."
    )
  }

  resolved_pc_penalty
}

.matern_inlabru_known_laplace_objective <- function(fit,
                                                    x,
                                                    s,
                                                    A,
                                                    spde_template,
                                                    alpha,
                                                    d,
                                                    beta_mode,
                                                    beta_fixed = NULL,
                                                    beta_prec = NULL,
                                                    link = "softplus",
                                                    pc_penalty = NULL,
                                                    suppress_warnings = TRUE) {
  .matern_laplace_known_noise_objective_at_params(
    x = x,
    s = s,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = fit$fitted_g$theta,
    log_sigma = log(fit$fitted_g$sigma),
    beta_mode = beta_mode,
    beta_fixed = if (identical(beta_mode, "fixed")) fit$fitted_beta else beta_fixed,
    beta_prec = beta_prec,
    beta_init = fit$fitted_beta,
    link = link,
    pc_penalty = pc_penalty,
    initial_mode = .matern_inla_laplace_initial_mode(fit, beta_mode),
    compute_posterior = FALSE,
    optimize_beta = FALSE,
    suppress_warnings = suppress_warnings
  )
}

.matern_inlabru_unknown_laplace_objective <- function(fit,
                                                      x,
                                                      A,
                                                      spde_template,
                                                      alpha,
                                                      d,
                                                      beta_mode,
                                                      beta_fixed = NULL,
                                                      beta_prec = NULL,
                                                      link = "softplus",
                                                      pc_penalty = NULL,
                                                      suppress_warnings = TRUE) {
  .matern_laplace_unknown_noise_objective_at_params(
    x = x,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = fit$fitted_g$theta,
    log_sigma = log(fit$fitted_g$sigma),
    log_noise_sd = log(fit$fitted_noise_sd),
    beta_mode = beta_mode,
    beta_fixed = if (identical(beta_mode, "fixed")) fit$fitted_beta else beta_fixed,
    beta_prec = beta_prec,
    beta_init = fit$fitted_beta,
    link = link,
    pc_penalty = pc_penalty,
    initial_mode = .matern_inla_laplace_initial_mode(fit, beta_mode),
    compute_posterior = FALSE,
    optimize_beta = FALSE,
    suppress_warnings = suppress_warnings
  )
}

.matern_inlabru_data_frame <- function(x, locations, d) {
  df <- data.frame(.Y = as.numeric(x))
  if (d == 1L) {
    df$.loc <- as.numeric(locations)
  } else {
    df$.loc <- as.matrix(locations)
  }
  df
}

.matern_inlabru_run_iinla <- function(components,
                                      obs,
                                      theta_init,
                                      bru_initial = NULL,
                                      bru_max_iter = 25L,
                                      suppress_warnings = TRUE) {
  options_list <- list(
    control.inla = list(int.strategy = "eb", strategy = "gaussian"),
    control.mode = INLA::control.mode(theta = theta_init, restart = TRUE),
    control.compute = list(mlik = TRUE, hyperpar = TRUE,
                           return.marginals = FALSE, config = TRUE),
    bru_max_iter = bru_max_iter,
    verbose = FALSE
  )
  if (!is.null(bru_initial)) options_list$bru_initial <- bru_initial

  runner <- function() inlabru::bru(components, obs, options = options_list)
  if (suppress_warnings) suppressWarnings(runner()) else runner()
}

.matern_inlabru_decode_theta <- function(fit, alpha, d, pc_penalty,
                                         expect_length, noise_offset) {
  theta_hat <- as.numeric(fit$mode$theta)
  if (length(theta_hat) == 0L && !is.null(fit$bru_info$model$inla$mode$theta)) {
    theta_hat <- as.numeric(fit$bru_info$model$inla$mode$theta)
  }
  if (length(theta_hat) != expect_length) {
    stop("The inlabru Matern fit did not return the expected ", expect_length,
         " hyperparameters.")
  }

  spde_theta <- if (noise_offset > 0L) theta_hat[(noise_offset + 1L):expect_length]
                 else theta_hat
  matern_theta <- .matern_inla_theta_to_matern(
    theta = spde_theta,
    alpha = alpha,
    d = d,
    pc_penalty = pc_penalty
  )
  list(theta_hat = theta_hat, matern_theta = matern_theta)
}

.matern_inlabru_posterior <- function(fit, locations, d,
                                      use_intercept,
                                      beta_offset_value = NULL,
                                      n_samples = 1000L,
                                      seed = 1L) {
  nd <- if (d == 1L) {
    data.frame(.loc = as.numeric(locations))
  } else {
    out <- data.frame(.row = seq_len(nrow(locations)))
    out$.loc <- as.matrix(locations)
    out$.row <- NULL
    out
  }
  if (!use_intercept && !is.null(beta_offset_value)) {
    nd$.beta_offset <- rep(as.numeric(beta_offset_value), nrow(nd))
  }
  predictor <- if (use_intercept) {
    ~ log1p(exp(field + Intercept))
  } else {
    ~ log1p(exp(field + .beta_offset))
  }

  pp <- stats::predict(fit, newdata = nd, formula = predictor,
                       n.samples = n_samples, seed = seed)
  posterior <- data.frame(
    mean = as.numeric(pp$mean),
    var = as.numeric(pp$sd)^2
  )
  posterior$second_moment <- posterior$mean^2 + posterior$var

  list(
    posterior = posterior,
    posterior_spatial_field = fit$summary.random[["field"]]
  )
}

.matern_inlabru_log_marginal <- function(fit) {
  mlik <- fit$mlik
  if (is.null(mlik)) return(NA_real_)
  rownames_mlik <- rownames(mlik)
  if ("log marginal-likelihood (integration)" %in% rownames_mlik) {
    return(as.numeric(mlik["log marginal-likelihood (integration)", 1]))
  }
  as.numeric(mlik[1, 1])
}

.fit_matern_inlabru_stepA <- function(x,
                                      s,
                                      A,
                                      mesh,
                                      alpha,
                                      d,
                                      locations,
                                      theta_init,
                                      beta_init = 0,
                                      beta_prec = 0,
                                      pc_penalty,
                                      link = "softplus",
                                      suppress_warnings = TRUE) {
  .matern_inlabru_assert_softplus(link)
  .matern_inlabru_assert_installed()
  .with_quiet_inla_defaults({
    n <- length(x)
    spde <- .build_matern_inla_spde(mesh, alpha, pc_penalty,
                                     suppress_warnings = suppress_warnings)
    df <- .matern_inlabru_data_frame(x, locations, d)
    df$.scale_w <- as.numeric(1 / s^2)

    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    if (is.null(beta_prec)) beta_prec <- 0
    beta_init <- .check_single_numeric(beta_init, "beta_init")

    components <- ~ Intercept(1, prec.linear = if (beta_prec > 0) beta_prec else 1e-6) +
      field(.loc, model = spde)
    formula <- .Y ~ log1p(exp(field + Intercept))

    obs <- inlabru::bru_obs(
      formula = formula,
      family = "gaussian",
      data = df,
      scale = .scale_w,
      control.family = list(
        hyper = list(prec = list(initial = 0, fixed = TRUE))
      )
    )

    fit <- .matern_inlabru_run_iinla(
      components = components,
      obs = obs,
      theta_init = theta_init,
      bru_initial = list(Intercept = beta_init),
      suppress_warnings = suppress_warnings
    )

    decoded <- .matern_inlabru_decode_theta(
      fit, alpha = alpha, d = d, pc_penalty = pc_penalty,
      expect_length = 2L, noise_offset = 0L
    )
    matern_theta <- decoded$matern_theta
    fitted_beta <- as.numeric(fit$summary.fixed["Intercept", "mean"])
    post <- .matern_inlabru_posterior(
      fit, locations = locations, d = d,
      use_intercept = TRUE
    )

    list(
      result = fit,
      posterior = post$posterior,
      posterior_spatial_field = post$posterior_spatial_field,
      fitted_g = Matern(theta = matern_theta$log_range,
                        sigma = matern_theta$sigma,
                        beta = fitted_beta,
                        beta_prec = beta_prec),
      fitted_beta = fitted_beta,
      log_likelihood_inlabru_mlik_integration = .matern_inlabru_log_marginal(fit),
      log_likelihood_stepA_mlik_integration = .matern_inlabru_log_marginal(fit)
    )
  })
}

.fit_matern_inlabru_stepA_fixed_beta <- function(x,
                                                 s,
                                                 A,
                                                 mesh,
                                                 alpha,
                                                 d,
                                                 locations,
                                                 theta_init,
                                                 beta_fixed,
                                                 pc_penalty,
                                                 link = "softplus",
                                                 suppress_warnings = TRUE) {
  .matern_inlabru_assert_softplus(link)
  .matern_inlabru_assert_installed()
  .with_quiet_inla_defaults({
    n <- length(x)
    beta_fixed <- .check_single_numeric(beta_fixed, "beta_fixed")
    spde <- .build_matern_inla_spde(mesh, alpha, pc_penalty,
                                     suppress_warnings = suppress_warnings)
    df <- .matern_inlabru_data_frame(x, locations, d)
    df$.scale_w <- as.numeric(1 / s^2)
    df$.beta_offset <- rep(beta_fixed, n)

    components <- ~ field(.loc, model = spde)
    formula <- .Y ~ log1p(exp(field + .beta_offset))

    obs <- inlabru::bru_obs(
      formula = formula,
      family = "gaussian",
      data = df,
      scale = .scale_w,
      control.family = list(
        hyper = list(prec = list(initial = 0, fixed = TRUE))
      )
    )

    fit <- .matern_inlabru_run_iinla(
      components = components,
      obs = obs,
      theta_init = theta_init,
      suppress_warnings = suppress_warnings
    )

    decoded <- .matern_inlabru_decode_theta(
      fit, alpha = alpha, d = d, pc_penalty = pc_penalty,
      expect_length = 2L, noise_offset = 0L
    )
    matern_theta <- decoded$matern_theta
    post <- .matern_inlabru_posterior(
      fit, locations = locations, d = d,
      use_intercept = FALSE,
      beta_offset_value = beta_fixed
    )

    list(
      result = fit,
      posterior = post$posterior,
      posterior_spatial_field = post$posterior_spatial_field,
      fitted_g = Matern(theta = matern_theta$log_range,
                        sigma = matern_theta$sigma,
                        beta = beta_fixed,
                        beta_prec = NULL),
      fitted_beta = beta_fixed,
      log_likelihood_inlabru_mlik_integration = .matern_inlabru_log_marginal(fit),
      log_likelihood_stepA_mlik_integration = .matern_inlabru_log_marginal(fit)
    )
  })
}

.fit_matern_inlabru_stepA_unknown_noise <- function(x,
                                                    A,
                                                    mesh,
                                                    alpha,
                                                    d,
                                                    locations,
                                                    theta_init,
                                                    noise_sd_init,
                                                    beta_init = 0,
                                                    beta_prec = 0,
                                                    pc_penalty,
                                                    link = "softplus",
                                                    suppress_warnings = TRUE) {
  .matern_inlabru_assert_softplus(link)
  .matern_inlabru_assert_installed()
  .with_quiet_inla_defaults({
    n <- length(x)
    spde <- .build_matern_inla_spde(mesh, alpha, pc_penalty,
                                     suppress_warnings = suppress_warnings)
    df <- .matern_inlabru_data_frame(x, locations, d)

    noise_hyper <- .matern_pc_noise_hyper(
      if (is.null(pc_penalty)) NULL else pc_penalty$noise,
      noise_sd_init = noise_sd_init
    )
    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    if (is.null(beta_prec)) beta_prec <- 0
    beta_init <- .check_single_numeric(beta_init, "beta_init")

    components <- ~ Intercept(1, prec.linear = if (beta_prec > 0) beta_prec else 1e-6) +
      field(.loc, model = spde)
    formula <- .Y ~ log1p(exp(field + Intercept))

    obs <- inlabru::bru_obs(
      formula = formula,
      family = "gaussian",
      data = df,
      control.family = list(hyper = list(prec = noise_hyper))
    )

    fit <- .matern_inlabru_run_iinla(
      components = components,
      obs = obs,
      theta_init = c(noise_hyper$initial, theta_init),
      bru_initial = list(Intercept = beta_init),
      suppress_warnings = suppress_warnings
    )

    decoded <- .matern_inlabru_decode_theta(
      fit, alpha = alpha, d = d, pc_penalty = pc_penalty,
      expect_length = 3L, noise_offset = 1L
    )
    theta_hat <- decoded$theta_hat
    matern_theta <- decoded$matern_theta
    fitted_beta <- as.numeric(fit$summary.fixed["Intercept", "mean"])
    post <- .matern_inlabru_posterior(
      fit, locations = locations, d = d,
      use_intercept = TRUE
    )

    list(
      result = fit,
      posterior = post$posterior,
      posterior_spatial_field = post$posterior_spatial_field,
      fitted_g = Matern(theta = matern_theta$log_range,
                        sigma = matern_theta$sigma,
                        beta = fitted_beta,
                        beta_prec = beta_prec),
      fitted_beta = fitted_beta,
      fitted_noise_sd = as.numeric(exp(-0.5 * theta_hat[1])),
      fitted_log_precision = as.numeric(theta_hat[1]),
      log_likelihood_inlabru_mlik_integration = .matern_inlabru_log_marginal(fit),
      log_likelihood_stepA_mlik_integration = .matern_inlabru_log_marginal(fit)
    )
  })
}

.fit_matern_inlabru_stepA_unknown_noise_fixed_beta <- function(x,
                                                               A,
                                                               mesh,
                                                               alpha,
                                                               d,
                                                               locations,
                                                               theta_init,
                                                               noise_sd_init,
                                                               beta_fixed,
                                                               pc_penalty,
                                                               link = "softplus",
                                                               suppress_warnings = TRUE) {
  .matern_inlabru_assert_softplus(link)
  .matern_inlabru_assert_installed()
  .with_quiet_inla_defaults({
    n <- length(x)
    beta_fixed <- .check_single_numeric(beta_fixed, "beta_fixed")
    spde <- .build_matern_inla_spde(mesh, alpha, pc_penalty,
                                     suppress_warnings = suppress_warnings)
    df <- .matern_inlabru_data_frame(x, locations, d)
    df$.beta_offset <- rep(beta_fixed, n)

    noise_hyper <- .matern_pc_noise_hyper(
      if (is.null(pc_penalty)) NULL else pc_penalty$noise,
      noise_sd_init = noise_sd_init
    )

    components <- ~ field(.loc, model = spde)
    formula <- .Y ~ log1p(exp(field + .beta_offset))

    obs <- inlabru::bru_obs(
      formula = formula,
      family = "gaussian",
      data = df,
      control.family = list(hyper = list(prec = noise_hyper))
    )

    fit <- .matern_inlabru_run_iinla(
      components = components,
      obs = obs,
      theta_init = c(noise_hyper$initial, theta_init),
      suppress_warnings = suppress_warnings
    )

    decoded <- .matern_inlabru_decode_theta(
      fit, alpha = alpha, d = d, pc_penalty = pc_penalty,
      expect_length = 3L, noise_offset = 1L
    )
    theta_hat <- decoded$theta_hat
    matern_theta <- decoded$matern_theta
    post <- .matern_inlabru_posterior(
      fit, locations = locations, d = d,
      use_intercept = FALSE,
      beta_offset_value = beta_fixed
    )

    list(
      result = fit,
      posterior = post$posterior,
      posterior_spatial_field = post$posterior_spatial_field,
      fitted_g = Matern(theta = matern_theta$log_range,
                        sigma = matern_theta$sigma,
                        beta = beta_fixed,
                        beta_prec = NULL),
      fitted_beta = beta_fixed,
      fitted_noise_sd = as.numeric(exp(-0.5 * theta_hat[1])),
      fitted_log_precision = as.numeric(theta_hat[1]),
      log_likelihood_inlabru_mlik_integration = .matern_inlabru_log_marginal(fit),
      log_likelihood_stepA_mlik_integration = .matern_inlabru_log_marginal(fit)
    )
  })
}

.fit_matern_inlabru_stepA_profile_beta <- function(x,
                                                   s,
                                                   A,
                                                   mesh,
                                                   spde_template,
                                                   alpha,
                                                   d,
                                                   locations,
                                                   theta_init,
                                                   beta_init,
                                                   pc_penalty,
                                                   link = "softplus",
                                                   suppress_warnings = TRUE) {
  .matern_inlabru_assert_softplus(link)

  fit_at_beta <- function(beta_value) {
    .fit_matern_inlabru_stepA_fixed_beta(
      x = x, s = s, A = A, mesh = mesh, alpha = alpha, d = d,
      locations = locations,
      theta_init = theta_init,
      beta_fixed = beta_value,
      pc_penalty = pc_penalty,
      link = link,
      suppress_warnings = suppress_warnings
    )
  }

  objective_at_fit <- function(fit, beta_value) {
    as.numeric(fit$log_likelihood_stepA_mlik_integration)
  }

  .fit_matern_inla_profile_beta(
    beta_init = beta_init,
    fit_at_beta = fit_at_beta,
    objective_at_fit = objective_at_fit,
    failure_label = "inlabru Matern beta profiling failed"
  )
}

.fit_matern_inlabru_stepA_unknown_noise_profile_beta <- function(x,
                                                                 A,
                                                                 mesh,
                                                                 spde_template,
                                                                 alpha,
                                                                 d,
                                                                 locations,
                                                                 theta_init,
                                                                 noise_sd_init,
                                                                 beta_init,
                                                                 pc_penalty,
                                                                 link = "softplus",
                                                                 suppress_warnings = TRUE) {
  .matern_inlabru_assert_softplus(link)

  fit_at_beta <- function(beta_value) {
    .fit_matern_inlabru_stepA_unknown_noise_fixed_beta(
      x = x, A = A, mesh = mesh, alpha = alpha, d = d,
      locations = locations,
      theta_init = theta_init,
      noise_sd_init = noise_sd_init,
      beta_fixed = beta_value,
      pc_penalty = pc_penalty,
      link = link,
      suppress_warnings = suppress_warnings
    )
  }

  objective_at_fit <- function(fit, beta_value) {
    as.numeric(fit$log_likelihood_stepA_mlik_integration)
  }

  .fit_matern_inla_profile_beta(
    beta_init = beta_init,
    fit_at_beta = fit_at_beta,
    objective_at_fit = objective_at_fit,
    failure_label = "inlabru learned-noise Matern beta profiling failed"
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
  .with_quiet_inla_defaults({
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
      fitted_g = Matern(theta = theta_fixed[1], sigma = exp(theta_fixed[2]), beta = as.numeric(beta_fixed), beta_prec = NULL),
      fitted_beta = as.numeric(beta_fixed),
      log_likelihood_stepB = as.numeric(resB$mlik["log marginal-likelihood (integration)", 1])
    )
  })
}

.exact_matern_add_pc_prior <- function(objective, pc_penalty, d) {
  if (is.null(pc_penalty)) {
    objective$log_pc_prior_theta <- NULL
    return(objective)
  }

  log_range <- if (!is.null(objective$stats$theta_log_range)) {
    objective$stats$theta_log_range
  } else {
    log(objective$stats$range)
  }
  log_sigma <- if (!is.null(objective$stats$theta_log_sigma)) {
    objective$stats$theta_log_sigma
  } else {
    log(objective$stats$sigma)
  }

  log_pc_prior_theta <- .log_pc_prior_matern_internal(
    log_range = log_range,
    log_sigma = log_sigma,
    range_spec = pc_penalty$range,
    sigma_spec = pc_penalty$sigma,
    d = d
  )
  log_pc_prior_noise <- 0
  if (!is.null(objective$fitted_noise_sd) && !is.null(pc_penalty$noise)) {
    log_pc_prior_noise <- .log_pc_prior_noise_internal(
      log_noise_sd = log(objective$fitted_noise_sd),
      noise_spec = pc_penalty$noise
    )
  }
  log_pc_prior_theta <- log_pc_prior_theta + log_pc_prior_noise
  objective$log_pc_prior_theta <- log_pc_prior_theta
  objective$log_pc_prior_noise <- log_pc_prior_noise
  objective$log_marginal <- as.numeric(objective$log_marginal + log_pc_prior_theta)
  objective
}

.matern_exact_state_with_beta_mode <- function(x,
                                               s,
                                               A,
                                               objective,
                                               beta_mode,
                                               beta_prec = NULL) {
  if (beta_mode == "fixed" || beta_mode == "empirical_bayes") {
    return(.exact_matern_posterior_from_stats(
      x = x,
      s = s,
      A = A,
      stats = objective$stats,
      beta0 = objective$fitted_beta,
      log_marginal = objective$log_marginal,
      beta_prec = beta_prec
    ))
  }

  .exact_matern_joint_posterior_from_stats(
    x = x,
    s = s,
    A = A,
    stats = objective$stats,
    beta_prec = if (beta_mode == "prior_flat") 0 else beta_prec,
    log_marginal = objective$log_marginal
  )
}

.matern_observation_terms <- function(eta,
                                      x,
                                      s,
                                      link = c("identity", "log", "softplus"),
                                      laplace_curvature = c("observed", "fisher")) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  eta <- as.numeric(eta)
  x <- as.numeric(x)
  s <- as.numeric(s)
  if (identical(link, "identity")) {
    mu <- eta
    grad <- (mu - x) / (s^2)
    hess_diag <- 1 / (s^2)
  } else if (identical(link, "log")) {
    if (any(eta > 40)) {
      stop("The log-link linear predictor overflowed during optimization.")
    }
    mu <- exp(eta)
    grad <- mu * (mu - x) / (s^2)
    hess_diag <- if (identical(laplace_curvature, "fisher")) {
      mu^2 / (s^2)
    } else {
      mu * (2 * mu - x) / (s^2)
    }
  } else {
    mu <- .softplus_stable(eta)
    sigmoid <- .sigmoid_stable(eta)
    grad <- ((mu - x) / (s^2)) * sigmoid
    hess_diag <- if (identical(laplace_curvature, "fisher")) {
      (sigmoid^2) / (s^2)
    } else {
      ((sigmoid^2) + (mu - x) * sigmoid * (1 - sigmoid)) / (s^2)
    }
  }
  resid <- x - mu
  nll <- 0.5 * sum((resid / s)^2 + log(2 * pi * s^2))
  if (!is.finite(nll) || any(!is.finite(grad)) || any(!is.finite(hess_diag))) {
    stop("Non-finite observation terms in Matern Laplace objective.")
  }
  list(nll = as.numeric(nll), grad = as.numeric(grad), hess_diag = as.numeric(hess_diag), mu = mu)
}

.matern_response_moments_from_eta <- function(eta_mean, eta_var, link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  eta_mean <- as.numeric(eta_mean)
  eta_var <- pmax(as.numeric(eta_var), 0)

  if (identical(link, "identity")) {
    return(list(mean = eta_mean, var = eta_var))
  }

  if (identical(link, "log")) {
    mean <- exp(eta_mean + 0.5 * eta_var)
    var <- exp(2 * eta_mean + eta_var) * (exp(eta_var) - 1)
  } else {
    moments <- .softplus_gaussian_moments(eta_mean, eta_var)
    mean <- moments$mean
    var <- moments$var
  }
  list(mean = as.numeric(mean), var = as.numeric(var))
}

.matern_laplace_inner_objective <- function(x,
                                            s,
                                            A,
                                            Q,
                                            Q_factor,
                                            beta_mode,
                                            beta_value = NULL,
                                            beta_prec = NULL,
                                            beta_init = 0,
                                            initial_mode = NULL,
                                            compute_posterior = TRUE,
                                            link = c("identity", "log", "softplus"),
                                            laplace_curvature = c("observed", "fisher")) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  if (!is.logical(compute_posterior) || length(compute_posterior) != 1L || is.na(compute_posterior)) {
    stop("`compute_posterior` must be TRUE or FALSE.")
  }
  x <- as.numeric(x)
  s <- as.numeric(s)
  A <- Matrix::Matrix(A, sparse = TRUE)
  Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))
  n_obs <- length(x)
  n_spde <- ncol(A)
  beta_mode <- as.character(beta_mode)
  integrate_beta <- beta_mode %in% c("prior_flat", "prior_proper")
  beta_prec0 <- if (identical(beta_mode, "prior_proper")) {
    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    if (is.null(beta_prec) || beta_prec <= 0) {
      stop("`beta_prec` must be positive for the proper-beta Laplace objective.")
    }
    beta_prec
  } else {
    0
  }

  if (integrate_beta) {
    A_eta <- cbind(A, Matrix::Matrix(rep(1, n_obs), ncol = 1, sparse = TRUE))
    prior_precision <- rbind(
      cbind(Q, Matrix::Matrix(0, nrow = n_spde, ncol = 1, sparse = TRUE)),
      cbind(
        Matrix::Matrix(0, nrow = 1, ncol = n_spde, sparse = TRUE),
        Matrix::Matrix(beta_prec0, nrow = 1, ncol = 1, sparse = TRUE)
      )
    )
    prior_precision <- Matrix::forceSymmetric(Matrix::Matrix(prior_precision, sparse = TRUE))
    default_z0 <- c(rep(0, n_spde), as.numeric(beta_init)[1])
  } else {
    beta_value <- .check_single_numeric(beta_value, "beta_value")
    A_eta <- A
    prior_precision <- Q
    default_z0 <- rep(0, n_spde)
  }

  z0 <- if (!is.null(initial_mode) &&
            length(initial_mode) == length(default_z0) &&
            all(is.finite(initial_mode))) {
    as.numeric(initial_mode)
  } else {
    default_z0
  }

  objective <- function(z) {
    z <- as.numeric(z)
    eta <- as.numeric(A_eta %*% z)
    if (!integrate_beta) eta <- eta + beta_value
    obs <- .matern_observation_terms(eta = eta, x = x, s = s, link = link, laplace_curvature = laplace_curvature)
    prior_quad <- as.numeric(0.5 * sum(z * as.numeric(prior_precision %*% z)))
    obs$nll + prior_quad
  }

  gradient <- function(z) {
    z <- as.numeric(z)
    eta <- as.numeric(A_eta %*% z)
    if (!integrate_beta) eta <- eta + beta_value
    obs <- .matern_observation_terms(eta = eta, x = x, s = s, link = link, laplace_curvature = laplace_curvature)
    as.numeric(Matrix::t(A_eta) %*% obs$grad + prior_precision %*% z)
  }

  opt <- stats::nlminb(
    start = z0,
    objective = objective,
    gradient = gradient,
    control = list(eval.max = 1000, iter.max = 1000)
  )
  if (!is.finite(opt$objective)) {
    stop("Matern Laplace inner optimization failed.")
  }

  z_mode <- as.numeric(opt$par)
  eta_mode <- as.numeric(A_eta %*% z_mode)
  if (!integrate_beta) eta_mode <- eta_mode + beta_value
  obs <- .matern_observation_terms(eta = eta_mode, x = x, s = s, link = link, laplace_curvature = laplace_curvature)
  W_eta <- Matrix::Diagonal(x = obs$hess_diag)
  H <- prior_precision + Matrix::t(A_eta) %*% W_eta %*% A_eta
  H <- Matrix::forceSymmetric(Matrix::Matrix(H, sparse = TRUE))
  H_factor <- .factorize_spd(H)

  log_prior <- 0.5 * Q_factor$logdet - 0.5 * n_spde * log(2 * pi)
  if (integrate_beta && beta_prec0 > 0) {
    beta_mode_value <- z_mode[length(z_mode)]
    log_prior <- log_prior + 0.5 * log(beta_prec0) - 0.5 * log(2 * pi) -
      0.5 * beta_prec0 * beta_mode_value^2
    w_mode <- z_mode[seq_len(n_spde)]
    log_prior <- log_prior - 0.5 * sum(w_mode * as.numeric(Q %*% w_mode))
  } else if (integrate_beta) {
    w_mode <- z_mode[seq_len(n_spde)]
    log_prior <- log_prior - 0.5 * sum(w_mode * as.numeric(Q %*% w_mode))
  } else {
    log_prior <- log_prior - 0.5 * sum(z_mode * as.numeric(Q %*% z_mode))
  }

  log_joint <- -obs$nll + log_prior
  log_marginal <- log_joint + 0.5 * length(z_mode) * log(2 * pi) - 0.5 * H_factor$logdet

  if (integrate_beta) {
    w_mode <- z_mode[seq_len(n_spde)]
    fitted_beta <- z_mode[length(z_mode)]
    eta_design <- A_eta
  } else {
    w_mode <- z_mode
    fitted_beta <- beta_value
    eta_design <- A
  }

  if (!isTRUE(compute_posterior)) {
    return(list(
      log_marginal = as.numeric(log_marginal),
      mode = z_mode,
      precision = H,
      precision_factor = H_factor,
      fitted_beta = as.numeric(fitted_beta),
      post_mean_latent = as.numeric(w_mode),
      eta_mean = eta_mode,
      integrated_beta = integrate_beta,
      eta_design = eta_design,
      inner_convergence = opt$convergence,
      inner_message = opt$message
    ))
  }

  eta_var <- .compute_diag_A_Qinv_At(eta_design, H_factor)
  response <- .matern_response_moments_from_eta(eta_mode, eta_var, link = link)

  latent_design <- if (integrate_beta) {
    cbind(Matrix::Diagonal(n = n_spde), Matrix::Matrix(0, nrow = n_spde, ncol = 1, sparse = TRUE))
  } else {
    Matrix::Diagonal(n = n_spde)
  }
  latent_var <- .compute_diag_A_Qinv_At(latent_design, H_factor)
  latent_sd <- sqrt(pmax(latent_var, 0))
  q025 <- stats::qnorm(0.025)
  q975 <- stats::qnorm(0.975)

  list(
    log_marginal = as.numeric(log_marginal),
    mode = z_mode,
    precision = H,
    precision_factor = H_factor,
    fitted_beta = as.numeric(fitted_beta),
    post_mean_latent = as.numeric(w_mode),
    post_var_latent = latent_var,
    eta_mean = eta_mode,
    eta_var = eta_var,
    posterior = data.frame(
      mean = response$mean,
      var = response$var,
      second_moment = response$mean^2 + response$var
    ),
    posterior_spatial_field = data.frame(
      ID = seq_along(w_mode),
      mean = w_mode,
      sd = latent_sd,
      `0.025quant` = w_mode + q025 * latent_sd,
      `0.5quant` = w_mode,
      `0.975quant` = w_mode + q975 * latent_sd,
      mode = w_mode,
      kld = rep(0, length(w_mode)),
      var = latent_var,
      check.names = FALSE
    ),
    integrated_beta = integrate_beta,
    eta_design = eta_design,
    inner_convergence = opt$convergence,
    inner_message = opt$message
  )
}

.matern_laplace_known_noise_objective_at_params <- function(x,
                                                            s,
                                                            A,
                                                            spde_template,
                                                            alpha,
                                                            d,
                                                            log_range,
                                                            log_sigma,
                                                            beta_mode,
                                                            beta_fixed = NULL,
                                                            beta_prec = NULL,
                                                            beta_init = 0,
                                                            link = c("identity", "log", "softplus"),
                                                            pc_penalty = NULL,
                                                            initial_mode = NULL,
                                                            compute_posterior = TRUE,
                                                            optimize_beta = TRUE,
                                                            suppress_warnings = TRUE,
                                                            laplace_curvature = c("observed", "fisher")) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  raw_eval <- function(beta_value) {
    precision <- .matern_precision_from_log_params(
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      log_range = log_range,
      log_sigma = log_sigma
    )
    inner <- .matern_laplace_inner_objective(
      x = x,
      s = s,
      A = A,
      Q = precision$Q,
      Q_factor = precision$Q_factor,
      beta_mode = beta_mode,
      beta_value = beta_value,
      beta_prec = beta_prec,
      beta_init = beta_init,
      initial_mode = initial_mode,
      compute_posterior = compute_posterior,
      link = link,
      laplace_curvature = laplace_curvature
    )
    inner$stats <- precision
    inner$stats$theta_log_range <- log_range
    inner$stats$theta_log_sigma <- log_sigma
    inner
  }

  objective <- if (identical(beta_mode, "empirical_bayes") && isTRUE(optimize_beta)) {
    beta_objective <- function(beta_value) {
      out <- tryCatch(raw_eval(beta_value), error = function(e) e)
      if (inherits(out, "error")) return(1e100)
      -out$log_marginal
    }
    opt_beta <- stats::optim(par = as.numeric(beta_init)[1], fn = beta_objective, method = "BFGS")
    if (!is.finite(opt_beta$value) || opt_beta$value >= 1e99) {
      stop("Matern Laplace beta optimization failed.")
    }
    raw_eval(as.numeric(opt_beta$par[1]))
  } else if (identical(beta_mode, "empirical_bayes")) {
    raw_eval(as.numeric(beta_init)[1])
  } else if (identical(beta_mode, "fixed")) {
    raw_eval(.check_single_numeric(beta_fixed, "beta_fixed"))
  } else {
    raw_eval(NULL)
  }

  objective <- if (suppress_warnings) suppressWarnings(objective) else objective
  .exact_matern_add_pc_prior(objective, pc_penalty = pc_penalty, d = d)
}

.matern_laplace_unknown_noise_objective_at_params <- function(x,
                                                              A,
                                                              spde_template,
                                                              alpha,
                                                              d,
                                                              log_range,
                                                              log_sigma,
                                                              log_noise_sd,
                                                              beta_mode,
                                                              beta_fixed = NULL,
                                                              beta_prec = NULL,
                                                              beta_init = 0,
                                                              link = c("identity", "log", "softplus"),
                                                              pc_penalty = NULL,
                                                              initial_mode = NULL,
                                                              compute_posterior = TRUE,
                                                              optimize_beta = TRUE,
                                                              suppress_warnings = TRUE,
                                                              laplace_curvature = c("observed", "fisher")) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  if (!is.finite(log_noise_sd)) {
    stop("`log_noise_sd` must be finite.")
  }
  noise_sd <- exp(log_noise_sd)
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    stop("The learned observation noise SD must be positive.")
  }

  objective <- .matern_laplace_known_noise_objective_at_params(
    x = x,
    s = rep(noise_sd, length(x)),
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = log_range,
    log_sigma = log_sigma,
    beta_mode = beta_mode,
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    beta_init = beta_init,
    link = link,
    pc_penalty = pc_penalty,
    initial_mode = initial_mode,
    compute_posterior = compute_posterior,
    optimize_beta = optimize_beta,
    suppress_warnings = suppress_warnings,
    laplace_curvature = laplace_curvature
  )
  objective$fitted_noise_sd <- as.numeric(noise_sd)
  objective$s <- rep(noise_sd, length(x))
  if (!is.null(pc_penalty) && !is.null(pc_penalty$noise)) {
    log_pc_prior_noise <- .log_pc_prior_noise_internal(
      log_noise_sd = log_noise_sd,
      noise_spec = pc_penalty$noise
    )
    objective$log_pc_prior_noise <- log_pc_prior_noise
    objective$log_pc_prior_theta <- if (is.null(objective$log_pc_prior_theta)) {
      log_pc_prior_noise
    } else {
      objective$log_pc_prior_theta + log_pc_prior_noise
    }
    objective$log_marginal <- as.numeric(objective$log_marginal + log_pc_prior_noise)
  }
  objective
}

.matern_laplace_fit_from_objective <- function(x,
                                               s,
                                               A,
                                               spde_template,
                                               alpha,
                                               d,
                                               objective,
                                               beta_mode,
                                               beta_prec = NULL,
                                               link = c("identity", "log", "softplus"),
                                               pc_penalty = NULL,
                                               fitted_noise_sd = NULL,
                                               laplace_implementation = "r",
                                               laplace_curvature = c("observed", "fisher"),
                                               extra_diagnostics = list()) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  log_likelihood <- objective$log_marginal
  class(log_likelihood) <- "logLik"
  eta_var <- pmax(as.numeric(objective$eta_var), 0)

  posterior_sampler <- function(nsamp) {
    nsamp <- .check_single_numeric(nsamp, "nsamp")
    if (nsamp < 1 || nsamp != floor(nsamp)) {
      stop("`nsamp` must be a positive integer.")
    }
    samps <- LaplacesDemon::rmvnp(
      n = nsamp,
      mu = as.numeric(objective$mode),
      Omega = as.matrix(objective$precision)
    )
    if (is.null(dim(samps))) samps <- matrix(samps, nrow = 1)
    eta_s <- as.matrix(objective$eta_design) %*% t(samps)
    if (!isTRUE(objective$integrated_beta)) {
      eta_s <- sweep(eta_s, 1, objective$fitted_beta, `+`)
    }
    if (identical(link, "identity")) t(eta_s) else if (identical(link, "log")) t(exp(eta_s)) else t(.softplus_stable(eta_s))
  }

  out <- list(
    posterior = objective$posterior,
    posterior_eta = data.frame(
      mean = as.numeric(objective$eta_mean),
      sd = sqrt(eta_var),
      var = eta_var
    ),
    fitted_g = Matern(
      theta = objective$stats$theta_log_range,
      sigma = exp(objective$stats$theta_log_sigma),
      beta = objective$fitted_beta,
      beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec
    ),
    fitted_beta = objective$fitted_beta,
    beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec,
    beta_mode = beta_mode,
    log_likelihood = log_likelihood,
    log_likelihood_semantics = paste0(
      if (identical(laplace_curvature, "fisher")) "laplace_fisher_" else "laplace_",
      beta_mode
    ),
    posterior_sampler = posterior_sampler,
    posterior_spatial_field = objective$posterior_spatial_field,
    backend = if (identical(laplace_curvature, "fisher")) "laplace_fisher" else "laplace",
    laplace_implementation = laplace_implementation,
    laplace_curvature = laplace_curvature,
    link = link,
    prior_family = paste0(
      link,
      "_Matern",
      if (is.null(pc_penalty)) "" else if (identical(laplace_curvature, "fisher")) "_pc_laplace_fisher" else "_pc_laplace",
      if (is.null(fitted_noise_sd)) "" else "_learned_noise"
    ),
    pc_penalty = pc_penalty,
    log_likelihood_pc_prior_theta = objective$log_pc_prior_theta,
    laplace_diagnostics = c(list(
      laplace_implementation = laplace_implementation,
      laplace_curvature = laplace_curvature,
      inner_convergence = objective$inner_convergence,
      inner_message = objective$inner_message
    ), extra_diagnostics),
    matern_objective_context = list(
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d
    )
  )
  if (!is.null(fitted_noise_sd)) {
    out$fitted_noise_sd <- as.numeric(fitted_noise_sd)
  }
  out
}

.matern_laplace_fit_has_finite_posterior <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$posterior)) {
    return(FALSE)
  }

  required_columns <- c("mean", "var", "second_moment")
  if (!all(required_columns %in% names(fit$posterior))) {
    return(FALSE)
  }

  all(vapply(
    fit$posterior[required_columns],
    function(column) all(is.finite(column)),
    logical(1)
  ))
}

.matern_laplace_tmb_invalid_reason <- function(fit) {
  if (inherits(fit, "error")) {
    return(conditionMessage(fit))
  }

  log_likelihood <- if (is.null(fit$log_likelihood)) NA_real_ else as.numeric(fit$log_likelihood)
  if (!length(log_likelihood) || !is.finite(log_likelihood[1])) {
    return("TMB Matern Laplace fit returned a non-finite log-likelihood.")
  }

  stepB_joint_nll <- fit$laplace_diagnostics$stepB_joint_nll
  if (!is.null(stepB_joint_nll) &&
      length(stepB_joint_nll) &&
      any(!is.na(stepB_joint_nll)) &&
      any(!is.finite(stepB_joint_nll))) {
    return("TMB Matern Laplace fit returned a non-finite Step B joint objective.")
  }

  if (!.matern_laplace_fit_has_finite_posterior(fit)) {
    return("TMB Matern Laplace fit returned non-finite posterior moments.")
  }

  NULL
}

.matern_fisher_pql_link_values <- function(eta, link = c("log", "softplus")) {
  link <- match.arg(link)
  eta <- as.numeric(eta)
  if (identical(link, "log")) {
    if (any(eta > 40)) {
      stop("The log-link linear predictor overflowed during Fisher-PQL construction.")
    }
    mu <- exp(eta)
    deriv <- mu
  } else {
    mu <- .softplus_stable(eta)
    deriv <- .sigmoid_stable(eta)
  }
  if (any(!is.finite(mu)) || any(!is.finite(deriv))) {
    stop("Non-finite Fisher-PQL link values.")
  }
  list(mu = as.numeric(mu), deriv = as.numeric(deriv))
}

.matern_fisher_pql_pseudo_response <- function(x,
                                               eta,
                                               s = NULL,
                                               link = c("log", "softplus"),
                                               g_floor = 1e-6,
                                               max_floor_fraction = 0.5) {
  link <- match.arg(link)
  x <- as.numeric(x)
  eta <- as.numeric(eta)
  if (length(eta) != length(x)) {
    stop("Internal error: Fisher-PQL eta and x lengths differ.")
  }
  g_floor <- .check_single_numeric(g_floor, "g_floor")
  if (!is.finite(g_floor) || g_floor <= 0) {
    stop("`g_floor` must be positive for Fisher-PQL.")
  }

  vals <- .matern_fisher_pql_link_values(eta, link = link)
  floor_used <- vals$deriv < g_floor
  floor_fraction <- mean(floor_used)
  if (is.finite(floor_fraction) && floor_fraction > max_floor_fraction) {
    stop("Fisher-PQL derivative floor dominated the pseudo-response construction.")
  }

  deriv_eff <- pmax(vals$deriv, g_floor)
  z <- eta + (x - vals$mu) / deriv_eff
  noise_scale <- 1 / deriv_eff
  if (any(!is.finite(z)) || any(!is.finite(noise_scale))) {
    stop("Non-finite Fisher-PQL pseudo-response.")
  }

  out <- list(
    z = as.numeric(z),
    mu = vals$mu,
    deriv = vals$deriv,
    deriv_eff = as.numeric(deriv_eff),
    noise_scale = as.numeric(noise_scale),
    floor_fraction = as.numeric(floor_fraction)
  )
  if (!is.null(s)) {
    s_eff <- as.numeric(s) * out$noise_scale
    if (length(s_eff) != length(x) || anyNA(s_eff) || any(!is.finite(s_eff)) || any(s_eff <= 0)) {
      stop("Non-finite Fisher-PQL effective observation SD.")
    }
    out$s <- as.numeric(s_eff)
  }
  out
}

.matern_fisher_pql_reference_mode_at_params <- function(x,
                                                       s,
                                                       A,
                                                       spde_template,
                                                       alpha,
                                                       d,
                                                       log_range,
                                                       log_sigma,
                                                       beta_mode,
                                                       beta_value = NULL,
                                                       beta_prec = NULL,
                                                       beta_start = 0,
                                                       link = c("log", "softplus"),
                                                       log_noise_sd = NULL,
                                                       pql_w_start = NULL,
                                                       pql_g_floor = 1e-6) {
  link <- match.arg(link)
  x <- as.numeric(x)
  A <- Matrix::Matrix(A, sparse = TRUE)
  n_obs <- length(x)
  n_spde <- ncol(A)
  if (is.null(pql_w_start)) pql_w_start <- rep(0, n_spde)
  beta_random <- beta_mode %in% c("prior_flat", "prior_proper")
  beta_for_start <- .check_single_numeric(beta_start, "beta_start")
  eta_start <- as.numeric(A %*% as.numeric(pql_w_start)) + beta_for_start
  base_s <- if (is.null(log_noise_sd)) as.numeric(s) else rep(exp(log_noise_sd), n_obs)
  pseudo <- .matern_fisher_pql_pseudo_response(
    x = x,
    eta = eta_start,
    s = base_s,
    link = link,
    g_floor = pql_g_floor
  )
  precision <- .matern_precision_from_log_params(
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = log_range,
    log_sigma = log_sigma
  )
  Q <- precision$Q
  w_prec_diag <- 1 / (pseudo$s^2)
  W <- Matrix::Diagonal(x = w_prec_diag)

  if (isTRUE(beta_random)) {
    beta_prec0 <- if (identical(beta_mode, "prior_proper")) {
      beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
      if (is.null(beta_prec) || beta_prec <= 0) {
        stop("`beta_prec` must be positive for the proper-beta Fisher-PQL objective.")
      }
      beta_prec
    } else {
      0
    }
    A_eta <- cbind(A, Matrix::Matrix(rep(1, n_obs), ncol = 1, sparse = TRUE))
    prior_precision <- rbind(
      cbind(Q, Matrix::Matrix(0, nrow = n_spde, ncol = 1, sparse = TRUE)),
      cbind(
        Matrix::Matrix(0, nrow = 1, ncol = n_spde, sparse = TRUE),
        Matrix::Matrix(beta_prec0, nrow = 1, ncol = 1, sparse = TRUE)
      )
    )
    prior_precision <- Matrix::forceSymmetric(Matrix::Matrix(prior_precision, sparse = TRUE))
    H <- prior_precision + Matrix::t(A_eta) %*% W %*% A_eta
    rhs <- as.numeric(Matrix::t(A_eta) %*% (w_prec_diag * pseudo$z))
    z_mode <- as.numeric(.solve_spd_factor(.factorize_spd(H), rhs))
    w_mode <- z_mode[seq_len(n_spde)]
    fitted_beta <- z_mode[length(z_mode)]
  } else {
    beta_value <- .check_single_numeric(beta_value, "beta_value")
    H <- Q + Matrix::t(A) %*% W %*% A
    rhs <- as.numeric(Matrix::t(A) %*% (w_prec_diag * (pseudo$z - beta_value)))
    w_mode <- as.numeric(.solve_spd_factor(.factorize_spd(H), rhs))
    z_mode <- w_mode
    fitted_beta <- beta_value
  }

  list(
    mode = z_mode,
    w = w_mode,
    beta = as.numeric(fitted_beta),
    pseudo = pseudo
  )
}

.matern_fisher_pql_mode_par_from_exact_fit <- function(fit, beta_mode) {
  w_mode <- as.numeric(fit$posterior_spatial_field$mean)
  if (beta_mode %in% c("prior_flat", "prior_proper")) {
    return(stats::setNames(c(w_mode, as.numeric(fit$fitted_beta)), c(rep("w", length(w_mode)), "beta")))
  }
  stats::setNames(w_mode, rep("w", length(w_mode)))
}

.matern_fisher_pql_exact_stepA_optA <- function(exact_diag) {
  list(
    convergence = if (is.null(exact_diag$convergence)) NA_integer_ else exact_diag$convergence,
    message = if (is.null(exact_diag$message)) "" else exact_diag$message
  )
}

.fit_matern_fisher_pql_exact_stepA <- function(x,
                                              s,
                                              A,
                                              spde_template,
                                              alpha,
                                              d,
                                              theta_init,
                                              sigma_init,
                                              beta_init,
                                              beta_mode,
                                              beta_fixed = NULL,
                                              beta_prec = NULL,
                                              pc_penalty = NULL,
                                              link = c("log", "softplus"),
                                              learn_noise = FALSE,
                                              noise_sd_init = NULL,
                                              fix_g = FALSE,
                                              fix_params = character(),
                                              pql_inner_iter = 3L,
                                              pql_tol = 1e-4,
                                              pql_g_floor = 1e-6,
                                              suppress_warnings = TRUE) {
  link <- match.arg(link)
  pql_inner_iter <- .check_matern_pql_inner_iter(pql_inner_iter)
  pql_tol <- .check_single_numeric(pql_tol, "pql_tol")
  if (!is.finite(pql_tol) || pql_tol < 0) {
    stop("`pql_tol` must be non-negative for Fisher-PQL.")
  }
  if (isTRUE(fix_g)) fix_params <- unique(c(fix_params, "range", "sigma"))
  pql_g_floor <- .check_single_numeric(pql_g_floor, "pql_g_floor")
  if (!is.finite(pql_g_floor) || pql_g_floor <= 0) {
    stop("`pql_g_floor` must be positive for Fisher-PQL.")
  }
  if (isTRUE(learn_noise)) {
    noise_sd_init <- .check_single_numeric(noise_sd_init, "noise_sd_init")
    if (!is.finite(noise_sd_init) || noise_sd_init <= 0) {
      stop("`noise_sd_init` must be positive for learned-noise Fisher-PQL fits.")
    }
  }

  beta_start <- if (identical(beta_mode, "fixed")) {
    .check_single_numeric(beta_fixed, "beta_fixed")
  } else {
    as.numeric(beta_init)[1]
  }
  eta_current <- as.numeric(A %*% rep(0, ncol(A))) + beta_start
  theta_current <- theta_init
  sigma_current <- sigma_init
  beta_current <- beta_start
  noise_current <- noise_sd_init
  pql_floor_fraction <- rep(NA_real_, pql_inner_iter)
  pql_eta_change <- rep(NA_real_, pql_inner_iter)
  pql_stepA_log_marginals <- rep(NA_real_, pql_inner_iter)
  pseudo <- NULL
  stepA_fit <- NULL

  for (pql_iter in seq_len(pql_inner_iter)) {
    pseudo <- .matern_fisher_pql_pseudo_response(
      x = x,
      eta = eta_current,
      s = if (isTRUE(learn_noise)) NULL else as.numeric(s),
      link = link,
      g_floor = pql_g_floor
    )
    pql_floor_fraction[pql_iter] <- pseudo$floor_fraction

    stepA_fit <- if (isTRUE(learn_noise)) {
      .fit_matern_exact_unknown_noise(
        x = pseudo$z,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_current,
        sigma_init = sigma_current,
        noise_sd_init = noise_current,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed,
        beta_prec = beta_prec,
        pc_penalty = pc_penalty,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        fix_params = fix_params,
        noise_scale = pseudo$noise_scale
      )
    } else {
      .fit_matern_exact_known_noise(
        x = pseudo$z,
        s = pseudo$s,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_current,
        sigma_init = sigma_current,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed,
        beta_prec = beta_prec,
        pc_penalty = pc_penalty,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        fix_params = fix_params
      )
    }

    pql_stepA_log_marginals[pql_iter] <- as.numeric(stepA_fit$log_likelihood)
    w_current <- as.numeric(stepA_fit$posterior_spatial_field$mean)
    beta_current <- as.numeric(stepA_fit$fitted_beta)
    eta_next <- as.numeric(A %*% w_current) + beta_current
    pql_eta_change[pql_iter] <- max(abs(eta_next - eta_current))
    eta_current <- eta_next
    theta_current <- as.numeric(stepA_fit$fitted_g$theta)
    sigma_current <- as.numeric(stepA_fit$fitted_g$sigma)
    if (isTRUE(learn_noise)) {
      noise_current <- as.numeric(stepA_fit$fitted_noise_sd)
    }
  }

  exact_diag <- stepA_fit$exact_optimization
  if (is.null(exact_diag)) {
    exact_diag <- list(
      method = "unknown",
      convergence = NA_integer_,
      message = "exact Gaussian Step A diagnostics unavailable",
      counts = NULL,
      value = NA_real_
    )
  }
  exact_convergence <- if (is.null(exact_diag$convergence)) NA_integer_ else as.integer(exact_diag$convergence)
  exact_message <- if (is.null(exact_diag$message)) "" else exact_diag$message
  optA <- .matern_fisher_pql_exact_stepA_optA(exact_diag)
  pql_stepA_log_marginal <- as.numeric(stepA_fit$log_likelihood)
  fitted_log_range <- as.numeric(stepA_fit$fitted_g$theta)
  fitted_log_sigma <- log(as.numeric(stepA_fit$fitted_g$sigma))
  fitted_log_noise <- if (isTRUE(learn_noise)) log(as.numeric(stepA_fit$fitted_noise_sd)) else NULL
  fitted_s_original <- if (isTRUE(learn_noise)) rep(as.numeric(stepA_fit$fitted_noise_sd), length(x)) else as.numeric(s)
  fitted_beta_stepA <- as.numeric(stepA_fit$fitted_beta)
  mode_par <- .matern_fisher_pql_mode_par_from_exact_fit(stepA_fit, beta_mode = beta_mode)

  objective <- .matern_laplace_tmb_objective_from_mode(
    x = as.numeric(x),
    s = as.numeric(fitted_s_original),
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = fitted_log_range,
    log_sigma = fitted_log_sigma,
    log_marginal = NA_real_,
    beta_mode = beta_mode,
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    beta_init = fitted_beta_stepA,
    link = link,
    pc_penalty = pc_penalty,
    log_noise_sd = if (isTRUE(learn_noise)) fitted_log_noise else NULL,
    laplace_curvature = "fisher",
    optA = optA,
    mode_par = mode_par,
    mode_source = "exact_fisher_pql_stepA",
    mode_convergence = exact_convergence,
    mode_message = exact_message
  )
  final_pseudo <- .matern_fisher_pql_pseudo_response(
    x = x,
    eta = objective$eta_mean,
    s = fitted_s_original,
    link = link,
    g_floor = pql_g_floor
  )
  exact_counts <- exact_diag$counts
  exact_outer_iterations <- if (!is.null(exact_counts) && "function" %in% names(exact_counts)) {
    as.integer(exact_counts[["function"]])
  } else {
    NA_integer_
  }
  pql_diagnostics <- list(
    converged = identical(exact_convergence, 0L),
    inner_iterations = as.integer(pql_inner_iter),
    stepA_engine = "exact_gaussian",
    g_floor = as.numeric(pql_g_floor),
    tol = as.numeric(pql_tol),
    initial_floor_fraction = pql_floor_fraction[1L],
    stepA_floor_fraction = pql_floor_fraction[pql_inner_iter],
    final_floor_fraction = final_pseudo$floor_fraction,
    eta_change = pql_eta_change,
    max_eta_change = pql_eta_change[pql_inner_iter],
    eta_converged = pql_eta_change[pql_inner_iter] <= pql_tol,
    final_effective_s_range = range(final_pseudo$s),
    outer_convergence = exact_convergence,
    outer_message = exact_message,
    outer_iterations = exact_outer_iterations,
    outer_evaluations = exact_counts,
    stepA_n_starts = 1L,
    stepA_best_start = 1L,
    stepA_log_marginals = pql_stepA_log_marginals,
    mode_source = "exact_fisher_pql_stepA",
    exact_optimization = exact_diag
  )

  fit <- .matern_laplace_fit_from_objective(
    x = as.numeric(x),
    s = as.numeric(fitted_s_original),
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    objective = objective,
    beta_mode = beta_mode,
    beta_prec = beta_prec,
    link = link,
    pc_penalty = pc_penalty,
    fitted_noise_sd = if (isTRUE(learn_noise)) as.numeric(stepA_fit$fitted_noise_sd) else NULL,
    laplace_implementation = "exact_fisher_pql",
    laplace_curvature = "fisher",
    extra_diagnostics = list(
      fisher_pql = pql_diagnostics,
      stepA_log_marginal = pql_stepA_log_marginal,
      original_at_fisher_pql_mode = as.numeric(objective$log_marginal),
      fitted_log_noise = if (isTRUE(learn_noise)) fitted_log_noise else NULL
    )
  )
  fit$backend <- "fisher_pql"
  fit$log_likelihood_semantics <- paste0("fisher_laplace_at_fisher_pql_mode_", beta_mode)
  fit$log_likelihood_fisher_pql_stepA <- as.numeric(pql_stepA_log_marginal)
  fit$log_likelihood_original_at_fisher_pql_mode <- as.numeric(fit$log_likelihood)
  fit$fisher_pql_diagnostics <- pql_diagnostics
  fit$fisher_pql_mode <- as.numeric(objective$mode)
  fit$fisher_pql_eta_mode <- as.numeric(objective$eta_mean)
  fit$fitted_s <- as.numeric(final_pseudo$s)
  fit$prior_family <- paste0(
    link,
    "_Matern_fisher_pql",
    if (is.null(pc_penalty)) "" else "_pc",
    if (isTRUE(learn_noise)) "_learned_noise" else ""
  )
  fit
}

.fit_matern_fisher_pql_known_noise <- function(x,
                                               s,
                                               A,
                                               spde_template,
                                               alpha,
                                               d,
                                               theta_init,
                                               sigma_init,
                                               beta_init,
                                               beta_mode,
                                               beta_fixed = NULL,
                                               beta_prec = NULL,
                                               pc_penalty = NULL,
                                               link = c("log", "softplus"),
                                               suppress_warnings = TRUE,
                                               fix_g = FALSE,
                                               fix_params = character(),
                                               pql_max_iter = 3L,
                                               pql_tol = 1e-4,
                                               pql_g_floor = 1e-6) {
  link <- match.arg(link)
  .fit_matern_fisher_pql_exact_stepA(
    x = x,
    s = s,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    theta_init = theta_init,
    sigma_init = sigma_init,
    beta_init = beta_init,
    beta_mode = beta_mode,
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    pc_penalty = pc_penalty,
    link = link,
    learn_noise = FALSE,
    fix_g = fix_g,
    fix_params = fix_params,
    pql_inner_iter = pql_max_iter,
    pql_tol = pql_tol,
    pql_g_floor = pql_g_floor,
    suppress_warnings = suppress_warnings
  )
}

.fit_matern_fisher_pql_unknown_noise <- function(x,
                                                 A,
                                                 spde_template,
                                                 alpha,
                                                 d,
                                                 theta_init,
                                                 sigma_init,
                                                 noise_sd_init,
                                                 beta_init,
                                                 beta_mode,
                                                 beta_fixed = NULL,
                                                 beta_prec = NULL,
                                                 pc_penalty = NULL,
                                                 link = c("log", "softplus"),
                                                 suppress_warnings = TRUE,
                                                 fix_g = FALSE,
                                                 fix_params = character(),
                                                 pql_max_iter = 3L,
                                                 pql_tol = 1e-4,
                                                 pql_g_floor = 1e-6) {
  link <- match.arg(link)
  .fit_matern_fisher_pql_exact_stepA(
    x = x,
    s = rep(1, length(x)),
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    theta_init = theta_init,
    sigma_init = sigma_init,
    beta_init = beta_init,
    beta_mode = beta_mode,
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    pc_penalty = pc_penalty,
    link = link,
    learn_noise = TRUE,
    noise_sd_init = noise_sd_init,
    fix_g = fix_g,
    fix_params = fix_params,
    pql_inner_iter = pql_max_iter,
    pql_tol = pql_tol,
    pql_g_floor = pql_g_floor,
    suppress_warnings = suppress_warnings
  )
}

.matern_laplace_tmb_unsupported_reason <- function(alpha, beta_mode) {
  if (!isTRUE(all.equal(as.numeric(alpha), 2))) {
    return("The TMB Matern Laplace implementation currently supports `alpha = 2` only.")
  }
  NULL
}

.matern_laplace_tmb_unavailable_reason <- function(alpha, beta_mode, fixed_log_names = character()) {
  reason <- .matern_laplace_tmb_unsupported_reason(alpha = alpha, beta_mode = beta_mode)
  if (!is.null(reason)) {
    return(reason)
  }
  if (length(fixed_log_names) == 1L) {
    return("The TMB Matern Laplace implementation currently supports fixing both `range` and `sigma`, or neither, but not exactly one of them.")
  }
  NULL
}

.matern_laplace_tmb_supported <- function(alpha, beta_mode) {
  is.null(.matern_laplace_tmb_unsupported_reason(alpha = alpha, beta_mode = beta_mode))
}

.matern_laplace_tmb_failure_message <- function(reason,
                                                link,
                                                backend_use,
                                                learn_noise,
                                                fixed_log_names = character()) {
  noise_label <- if (isTRUE(learn_noise)) "learned" else "known"
  msg <- paste0(
    "Matern TMB Laplace fit failed: ", reason, "\n",
    "  link = \"", link, "\", backend = \"", backend_use, "\", noise = ", noise_label, ".\n",
    "  Public `backend = \"laplace\"` and `backend = \"laplace_fisher\"` require a successful TMB fit; ",
    "they do not fall back to the internal R reference implementation."
  )
  if (identical(link, "log")) {
    msg <- paste0(
      msg,
      "\n  If log-link posterior moments are numerically unstable, consider `link = \"softplus\"`."
    )
  }
  if (length(fixed_log_names) == 1L) {
    msg <- paste0(
      msg,
      "\n  For partial hyperparameter fixing, use explicit `backend = \"laplace_r\"` for internal validation only."
    )
  } else {
    msg <- paste0(
      msg,
      "\n  Explicit `backend = \"laplace_r\"` is available for internal validation only."
    )
  }
  msg
}

.matern_pc_prior_vector <- function(pc_penalty) {
  if (is.null(pc_penalty)) return(rep(1, 4))
  out <- as.numeric(c(
    pc_penalty$range["anchor"],
    pc_penalty$range["alpha"],
    pc_penalty$sigma["anchor"],
    pc_penalty$sigma["alpha"]
  ))
  if (!is.null(pc_penalty$noise)) {
    out <- c(
      out,
      as.numeric(pc_penalty$noise["anchor"]),
      as.numeric(pc_penalty$noise["alpha"])
    )
  }
  out
}

.matern_laplace_tmb_data <- function(x,
                                     s,
                                     A,
                                     spde_template,
                                     alpha,
                                     d,
                                     beta_mode,
                                     beta_prec = NULL,
                                     link = c("identity", "log", "softplus"),
                                     pc_penalty = NULL,
                                     learn_noise = FALSE) {
  link <- match.arg(link)
  reason <- .matern_laplace_tmb_unsupported_reason(alpha = alpha, beta_mode = beta_mode)
  if (!is.null(reason)) stop(reason)

  param <- spde_template$param.inla
  if (is.null(param$M0) || is.null(param$M1) || is.null(param$M2)) {
    stop("The Matern SPDE template does not contain the alpha=2 precision basis matrices.")
  }

  betaprec_internal <- if (identical(beta_mode, "prior_proper")) {
    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    if (is.null(beta_prec) || beta_prec <= 0) {
      stop("`beta_prec` must be positive for proper-prior beta.")
    }
    beta_prec
  } else {
    0
  }

  list(
    model_id = 1L,
    x = as.numeric(x),
    s = as.numeric(s),
    A = as(Matrix::Matrix(A, sparse = TRUE), "TsparseMatrix"),
    M0 = as(Matrix::Matrix(param$M0, sparse = TRUE), "TsparseMatrix"),
    M1 = as(Matrix::Matrix(param$M1, sparse = TRUE), "TsparseMatrix"),
    M2 = as(Matrix::Matrix(param$M2, sparse = TRUE), "TsparseMatrix"),
    betaprec = as.numeric(betaprec_internal),
    matern_alpha = as.numeric(alpha),
    matern_d = as.integer(d),
    link_id = if (identical(link, "identity")) 0L else if (identical(link, "log")) 1L else 2L,
    learn_noise = if (isTRUE(learn_noise)) 1L else 0L,
    use_pc_prior = if (is.null(pc_penalty)) 0L else 1L,
    use_pc_noise_prior = if (!is.null(pc_penalty) && !is.null(pc_penalty$noise)) 1L else 0L,
    pc_prior = .matern_pc_prior_vector(pc_penalty)
  )
}

.matern_tmb_make_adfun <- function(data,
                                   parameters,
                                   random = NULL,
                                   map = list(),
                                   inner.control = NULL,
                                   dll = "EBSmoothr") {
  args <- list(
    data = data,
    parameters = parameters,
    DLL = dll,
    silent = TRUE
  )
  if (!is.null(random)) args$random <- random
  if (length(map) > 0L) args$map <- map
  if (!is.null(inner.control)) args$inner.control <- inner.control
  do.call(TMB::MakeADFun, args)
}

.matern_laplace_tmb_optimize <- function(obj, label, lower = NULL, upper = NULL) {
  if (length(obj$par) == 0L) {
    value <- obj$fn(obj$par)
    return(list(
      par = obj$par,
      value = as.numeric(value),
      convergence = 0L,
      message = "No free fixed-effect parameters."
    ))
  }

  if (is.null(lower)) lower <- rep(-Inf, length(obj$par))
  if (is.null(upper)) upper <- rep(Inf, length(obj$par))

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    lower = lower,
    upper = upper,
    control = list(eval.max = 2000, iter.max = 2000)
  )
  if (!is.finite(opt$objective)) {
    stop(label, " failed with a non-finite objective.")
  }
  list(
    par = opt$par,
    value = as.numeric(opt$objective),
    convergence = opt$convergence,
    message = opt$message,
    iterations = opt$iterations,
    evaluations = opt$evaluations
  )
}

.matern_tmb_named_value <- function(par, name, default) {
  idx <- which(names(par) == name)
  if (length(idx) == 0L) return(default)
  as.numeric(par[idx[1]])
}

.matern_laplace_tmb_extract_random_mode <- function(obj, n_spde, beta_mode) {
  mode <- obj$env$last.par.best
  if (is.null(mode)) return(NULL)

  w_idx <- which(names(mode) == "w")
  if (length(w_idx) != n_spde) return(NULL)
  out <- as.numeric(mode[w_idx])
  names(out) <- rep("w", length(out))

  if (beta_mode %in% c("prior_flat", "prior_proper")) {
    beta_idx <- which(names(mode) == "beta")
    if (length(beta_idx) != 1L) return(NULL)
    out <- c(out, beta = as.numeric(mode[beta_idx]))
  }

  if (!all(is.finite(out))) return(NULL)
  out
}

.matern_laplace_tmb_objective_from_mode <- function(x,
                                                   s,
                                                   A,
                                                   spde_template,
                                                   alpha,
                                                   d,
                                                   log_range,
                                                   log_sigma,
                                                   log_marginal,
                                                   beta_mode,
                                                   beta_fixed = NULL,
                                                   beta_prec = NULL,
                                                   beta_init = 0,
                                                   link = c("identity", "log", "softplus"),
                                                   pc_penalty = NULL,
                                                   log_noise_sd = NULL,
                                                   laplace_curvature = c("observed", "fisher"),
                                                   optA = NULL,
                                                   mode_par,
                                                   mode_source = "unknown",
                                                   mode_convergence = NA_integer_,
                                                   mode_message = NA_character_) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  precision <- .matern_precision_from_log_params(
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = log_range,
    log_sigma = log_sigma
  )
  Q <- precision$Q
  Q_factor <- precision$Q_factor
  n_obs <- length(x)
  n_spde <- ncol(A)
  integrate_beta <- beta_mode %in% c("prior_flat", "prior_proper")
  mode_names <- names(mode_par)
  mode_par <- as.numeric(mode_par)
  names(mode_par) <- mode_names

  if (integrate_beta) {
    beta_prec0 <- if (identical(beta_mode, "prior_proper")) {
      beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
      if (is.null(beta_prec) || beta_prec <= 0) {
        stop("`beta_prec` must be positive for proper-prior beta.")
      }
      beta_prec
    } else {
      0
    }
    w_idx <- which(mode_names == "w")
    beta_idx <- which(mode_names == "beta")
    if (length(w_idx) != n_spde || length(beta_idx) != 1L) {
      stop("The TMB Matern random-effect mode has unexpected parameter names.")
    }
    w_mode <- as.numeric(mode_par[w_idx])
    fitted_beta <- as.numeric(mode_par[beta_idx])
    z_mode <- c(w_mode, fitted_beta)
    A_eta <- cbind(A, Matrix::Matrix(rep(1, n_obs), ncol = 1, sparse = TRUE))
    prior_precision <- rbind(
      cbind(Q, Matrix::Matrix(0, nrow = n_spde, ncol = 1, sparse = TRUE)),
      cbind(
        Matrix::Matrix(0, nrow = 1, ncol = n_spde, sparse = TRUE),
        Matrix::Matrix(beta_prec0, nrow = 1, ncol = 1, sparse = TRUE)
      )
    )
    prior_precision <- Matrix::forceSymmetric(Matrix::Matrix(prior_precision, sparse = TRUE))
    eta_mode <- as.numeric(A_eta %*% z_mode)
  } else {
    w_idx <- which(mode_names == "w")
    if (length(w_idx) != n_spde) {
      stop("The TMB Matern random-effect mode has unexpected latent-field parameter names.")
    }
    w_mode <- as.numeric(mode_par[w_idx])
    fitted_beta <- if (identical(beta_mode, "fixed")) {
      .check_single_numeric(beta_fixed, "beta_fixed")
    } else {
      as.numeric(beta_init)[1]
    }
    z_mode <- w_mode
    A_eta <- A
    prior_precision <- Q
    eta_mode <- as.numeric(A %*% w_mode) + fitted_beta
  }

  obs <- .matern_observation_terms(
    eta = eta_mode,
    x = x,
    s = s,
    link = link,
    laplace_curvature = laplace_curvature
  )
  W_eta <- Matrix::Diagonal(x = obs$hess_diag)
  H <- prior_precision + Matrix::t(A_eta) %*% W_eta %*% A_eta
  H <- Matrix::forceSymmetric(Matrix::Matrix(H, sparse = TRUE))
  H_factor <- .factorize_spd(H)

  log_prior <- 0.5 * Q_factor$logdet - 0.5 * n_spde * log(2 * pi)
  log_prior <- log_prior - 0.5 * sum(w_mode * as.numeric(Q %*% w_mode))
  if (integrate_beta && beta_prec0 > 0) {
    log_prior <- log_prior + 0.5 * log(beta_prec0) - 0.5 * log(2 * pi) -
      0.5 * beta_prec0 * fitted_beta^2
  }
  log_joint <- -obs$nll + log_prior

  eta_var <- .compute_diag_A_Qinv_At(A_eta, H_factor)
  response <- .matern_response_moments_from_eta(eta_mode, eta_var, link = link)

  latent_design <- if (integrate_beta) {
    cbind(Matrix::Diagonal(n = n_spde), Matrix::Matrix(0, nrow = n_spde, ncol = 1, sparse = TRUE))
  } else {
    Matrix::Diagonal(n = n_spde)
  }
  latent_var <- .compute_diag_A_Qinv_At(latent_design, H_factor)
  latent_sd <- sqrt(pmax(latent_var, 0))
  q025 <- stats::qnorm(0.025)
  q975 <- stats::qnorm(0.975)

  log_pc_prior_spde <- if (is.null(pc_penalty)) {
    NULL
  } else {
    .log_pc_prior_matern_internal(
      log_range = log_range,
      log_sigma = log_sigma,
      range_spec = pc_penalty$range,
      sigma_spec = pc_penalty$sigma,
      d = d
    )
  }
  log_pc_prior_noise <- if (!is.null(log_noise_sd) && !is.null(pc_penalty) && !is.null(pc_penalty$noise)) {
    .log_pc_prior_noise_internal(
      log_noise_sd = log_noise_sd,
      noise_spec = pc_penalty$noise
    )
  } else {
    0
  }
  log_pc_prior_theta <- if (is.null(log_pc_prior_spde)) {
    if (log_pc_prior_noise == 0) NULL else log_pc_prior_noise
  } else {
    log_pc_prior_spde + log_pc_prior_noise
  }
  log_pc_prior_total <- if (is.null(log_pc_prior_theta)) 0 else as.numeric(log_pc_prior_theta)
  log_marginal_use <- if (identical(laplace_curvature, "fisher")) {
    log_joint + 0.5 * length(z_mode) * log(2 * pi) - 0.5 * H_factor$logdet +
      log_pc_prior_total
  } else {
    as.numeric(log_marginal)
  }

  precision$theta_log_range <- log_range
  precision$theta_log_sigma <- log_sigma

  list(
    log_marginal = as.numeric(log_marginal_use),
    mode = z_mode,
    precision = H,
    precision_factor = H_factor,
    fitted_beta = as.numeric(fitted_beta),
    post_mean_latent = as.numeric(w_mode),
    post_var_latent = latent_var,
    eta_mean = eta_mode,
    eta_var = eta_var,
    posterior = data.frame(
      mean = response$mean,
      var = response$var,
      second_moment = response$mean^2 + response$var
    ),
    posterior_spatial_field = data.frame(
      ID = seq_along(w_mode),
      mean = w_mode,
      sd = latent_sd,
      `0.025quant` = w_mode + q025 * latent_sd,
      `0.5quant` = w_mode,
      `0.975quant` = w_mode + q975 * latent_sd,
      mode = w_mode,
      kld = rep(0, length(w_mode)),
      var = latent_var,
      check.names = FALSE
    ),
    integrated_beta = integrate_beta,
    eta_design = A_eta,
    stats = precision,
    log_pc_prior_theta = log_pc_prior_theta,
    log_pc_prior_noise = log_pc_prior_noise,
    laplace_curvature = laplace_curvature,
    observed_stepA_log_marginal = as.numeric(log_marginal),
    tmb_mode_source = mode_source,
    inner_convergence = mode_convergence,
    inner_message = mode_message,
    outer_convergence = if (is.null(optA)) NA_integer_ else optA$convergence,
    outer_message = if (is.null(optA)) NA_character_ else optA$message
  )
}

.fit_matern_laplace_tmb_known_noise <- function(x,
                                                s,
                                                A,
                                                spde_template,
                                                alpha,
                                                d,
                                                theta_init,
                                                sigma_init,
                                                beta_init,
                                                beta_mode,
                                                beta_fixed = NULL,
                                                beta_prec = NULL,
                                                pc_penalty = NULL,
                                                link = c("identity", "log", "softplus"),
                                                fix_g = FALSE,
                                                learn_noise = FALSE,
                                                noise_sd_init = NULL,
                                                matern_n_starts = 1L,
                                                laplace_curvature = c("observed", "fisher"),
                                                dll = "EBSmoothr") {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  matern_n_starts <- .check_matern_n_starts(matern_n_starts)
  reason <- .matern_laplace_tmb_unsupported_reason(alpha = alpha, beta_mode = beta_mode)
  if (!is.null(reason)) stop(reason)
  if (isTRUE(learn_noise)) {
    noise_sd_init <- .check_single_numeric(noise_sd_init, "noise_sd_init")
    if (!is.finite(noise_sd_init) || noise_sd_init <= 0) {
      stop("`noise_sd_init` must be positive for learned-noise TMB Matern fits.")
    }
  }

  beta_start <- if (identical(beta_mode, "fixed")) {
    .check_single_numeric(beta_fixed, "beta_fixed")
  } else {
    as.numeric(beta_init)[1]
  }
  tmb_s <- if (isTRUE(learn_noise)) rep(1, length(x)) else s
  tmbdat <- .matern_laplace_tmb_data(
    x = x,
    s = tmb_s,
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    beta_mode = beta_mode,
    beta_prec = beta_prec,
    link = link,
    pc_penalty = pc_penalty,
    learn_noise = learn_noise
  )

  n_spde <- ncol(A)
  par0 <- list(
    log_range = as.numeric(theta_init),
    log_sigma = log(as.numeric(sigma_init)),
    w = rep(0, n_spde),
    beta = as.numeric(beta_start),
    log_noise = if (isTRUE(learn_noise)) log(as.numeric(noise_sd_init)) else 0
  )
  integrate_beta <- beta_mode %in% c("prior_flat", "prior_proper")
  randomA <- if (integrate_beta) c("w", "beta") else "w"
  mapA <- list()
  if (isTRUE(fix_g)) {
    mapA$log_range <- factor(NA)
    mapA$log_sigma <- factor(NA)
  }
  if (identical(beta_mode, "fixed")) {
    mapA$beta <- factor(NA)
  }
  if (!isTRUE(learn_noise)) {
    mapA$log_noise <- factor(NA)
  }

  make_objA <- function(par_start) {
    .matern_tmb_make_adfun(
      data = tmbdat,
      parameters = par_start,
      random = randomA,
      map = mapA,
      dll = dll
    )
  }
  make_bounds <- function(obj) {
    lower <- rep(-Inf, length(obj$par))
    upper <- rep(Inf, length(obj$par))
    range_idx <- which(names(obj$par) == "log_range")
    sigma_idx <- which(names(obj$par) == "log_sigma")
    noise_idx <- which(names(obj$par) == "log_noise")
    if (length(range_idx) > 0L) {
      lower[range_idx] <- as.numeric(theta_init) - log(1000)
      upper[range_idx] <- as.numeric(theta_init) + log(1000)
    }
    if (length(sigma_idx) > 0L) {
      lower[sigma_idx] <- log(as.numeric(sigma_init)) - log(1000)
      upper[sigma_idx] <- log(as.numeric(sigma_init)) + log(1000)
    }
    if (length(noise_idx) > 0L) {
      lower[noise_idx] <- log(as.numeric(noise_sd_init)) - log(1000)
      upper[noise_idx] <- log(as.numeric(noise_sd_init)) + log(1000)
    }
    list(lower = lower, upper = upper)
  }
  shift_start <- function(log_range_offset = 0, log_sigma_offset = 0, log_noise_offset = 0) {
    out <- par0
    if (!isTRUE(fix_g)) {
      out$log_range <- as.numeric(out$log_range) + log_range_offset
      out$log_sigma <- as.numeric(out$log_sigma) + log_sigma_offset
    }
    if (isTRUE(learn_noise)) {
      out$log_noise <- as.numeric(out$log_noise) + log_noise_offset
    }
    out
  }
  stepA_starts <- list(par0)
  if (identical(link, "log") && identical(beta_mode, "empirical_bayes") && !isTRUE(fix_g)) {
    shifted_starts <- if (isTRUE(learn_noise)) {
      list(
        shift_start(log_range_offset = -log(2), log_sigma_offset = log(5), log_noise_offset = -log(50)),
        shift_start(log_range_offset = -log(2), log_sigma_offset = log(10), log_noise_offset = -log(100)),
        shift_start(log_range_offset = 0, log_sigma_offset = log(2), log_noise_offset = -log(10)),
        shift_start(log_range_offset = log(2), log_sigma_offset = -log(2), log_noise_offset = log(2))
      )
    } else {
      list(
        shift_start(log_range_offset = log(5), log_sigma_offset = 0),
        shift_start(log_range_offset = log(2), log_sigma_offset = 0),
        shift_start(log_range_offset = log(10), log_sigma_offset = -log(2)),
        shift_start(log_range_offset = -log(2), log_sigma_offset = log(2))
      )
    }
    n_extra_starts <- matern_n_starts - 1L
    if (n_extra_starts > 0L) {
      stepA_starts <- c(stepA_starts, shifted_starts[seq_len(min(n_extra_starts, length(shifted_starts)))])
    }
  }
  stepA_results <- lapply(seq_along(stepA_starts), function(start_index) {
    obj <- make_objA(stepA_starts[[start_index]])
    bounds <- make_bounds(obj)
    opt <- tryCatch(
      .matern_laplace_tmb_optimize(
        obj,
        label = "TMB Matern Laplace Step A",
        lower = bounds$lower,
        upper = bounds$upper
      ),
      error = function(e) e
    )
    if (inherits(opt, "error")) {
      return(NULL)
    }
    list(obj = obj, opt = opt, start_index = start_index)
  })
  stepA_results <- Filter(Negate(is.null), stepA_results)
  if (!length(stepA_results)) {
    stop("TMB Matern Laplace Step A failed for all starts.")
  }
  stepA_values <- vapply(stepA_results, function(z) as.numeric(z$opt$value), numeric(1))
  best_stepA <- stepA_results[[which.min(stepA_values)]]
  objA <- best_stepA$obj
  optA <- best_stepA$opt
  log_marginal <- -as.numeric(optA$value)

  fitted_log_range <- .matern_tmb_named_value(optA$par, "log_range", as.numeric(theta_init))
  fitted_log_sigma <- .matern_tmb_named_value(optA$par, "log_sigma", log(as.numeric(sigma_init)))
  fitted_log_noise <- .matern_tmb_named_value(optA$par, "log_noise", par0$log_noise)
  fitted_s <- if (isTRUE(learn_noise)) rep(exp(fitted_log_noise), length(x)) else as.numeric(s)
  fitted_beta_stepA <- if (identical(beta_mode, "fixed")) {
    as.numeric(beta_start)
  } else if (identical(beta_mode, "empirical_bayes")) {
    .matern_tmb_named_value(optA$par, "beta", as.numeric(beta_start))
  } else {
    as.numeric(beta_start)
  }

  mode_par <- .matern_laplace_tmb_extract_random_mode(
    obj = objA,
    n_spde = n_spde,
    beta_mode = beta_mode
  )
  mode_source <- "last.par.best"
  mode_convergence <- optA$convergence
  mode_message <- optA$message
  stepB_joint_nll <- NA_real_

  parB <- list(
    log_range = fitted_log_range,
    log_sigma = fitted_log_sigma,
    w = if (is.null(mode_par)) rep(0, n_spde) else as.numeric(mode_par[names(mode_par) == "w"]),
    beta = as.numeric(fitted_beta_stepA),
    log_noise = fitted_log_noise
  )
  if (!is.null(mode_par) && integrate_beta) {
    parB$beta <- as.numeric(mode_par[names(mode_par) == "beta"])
  }
  mapB <- list(
    log_range = factor(NA),
    log_sigma = factor(NA),
    log_noise = factor(NA)
  )
  if (!isTRUE(integrate_beta)) {
    mapB$beta <- factor(NA)
  }

  had_stepA_mode <- !is.null(mode_par)
  refine_mode <- !had_stepA_mode || !identical(as.integer(optA$convergence), 0L)
  if (isTRUE(refine_mode)) {
    objB <- .matern_tmb_make_adfun(
      data = tmbdat,
      parameters = parB,
      map = mapB,
      dll = dll
    )
    optB <- stats::nlminb(
      start = objB$par,
      objective = objB$fn,
      gradient = objB$gr,
      control = list(eval.max = 20000, iter.max = 20000)
    )
    stepB_joint_nll <- as.numeric(optB$objective)
    if (!is.finite(optB$objective)) {
      if (!had_stepA_mode) {
        stop("TMB Matern Laplace Step B failed with a non-finite objective.")
      }
      mode_source <- "last.par.best_stepB_failed"
      mode_convergence <- optB$convergence
      mode_message <- paste("Step B failed with a non-finite objective:", optB$message)
    } else {
      mode_par <- optB$par
      mode_source <- if (had_stepA_mode) "stepB_refine" else "stepB_fallback"
      mode_convergence <- optB$convergence
      mode_message <- optB$message
    }
  }

  objective <- .matern_laplace_tmb_objective_from_mode(
    x = as.numeric(x),
    s = as.numeric(fitted_s),
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    log_range = fitted_log_range,
    log_sigma = fitted_log_sigma,
    log_marginal = log_marginal,
    beta_mode = beta_mode,
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    beta_init = fitted_beta_stepA,
    link = link,
    pc_penalty = pc_penalty,
    log_noise_sd = if (isTRUE(learn_noise)) fitted_log_noise else NULL,
    laplace_curvature = laplace_curvature,
    optA = optA,
    mode_par = mode_par,
    mode_source = mode_source,
    mode_convergence = mode_convergence,
    mode_message = mode_message
  )

  fit_out <- .matern_laplace_fit_from_objective(
    x = as.numeric(x),
    s = as.numeric(fitted_s),
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    objective = objective,
    beta_mode = beta_mode,
    beta_prec = beta_prec,
    link = link,
    pc_penalty = pc_penalty,
    fitted_noise_sd = if (isTRUE(learn_noise)) exp(fitted_log_noise) else NULL,
    laplace_implementation = "tmb",
    laplace_curvature = laplace_curvature,
    extra_diagnostics = list(
      outer_convergence = objective$outer_convergence,
      outer_message = objective$outer_message,
      tmb_mode_source = objective$tmb_mode_source,
      stepA_log_marginal = log_marginal,
      curvature_log_marginal = as.numeric(objective$log_marginal),
      stepA_n_starts = length(stepA_starts),
      stepA_best_start = best_stepA$start_index,
      stepB_joint_nll = stepB_joint_nll,
      fitted_log_noise = if (isTRUE(learn_noise)) fitted_log_noise else NULL
    )
  )

  fit_out
}

.fit_matern_laplace_known_noise <- function(x,
                                            s,
                                            A,
                                            spde_template,
                                            alpha,
                                            d,
                                            theta_init,
                                            sigma_init,
                                            beta_init,
                                            beta_mode,
                                            beta_fixed = NULL,
                                            beta_prec = NULL,
                                            pc_penalty = NULL,
                                            link = c("identity", "log", "softplus"),
                                            suppress_warnings = TRUE,
                                            fix_params = character(),
                                            laplace_curvature = c("observed", "fisher")) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  fixed_names <- .matern_fixed_log_param_names(fix_params)
  last_inner_mode <- NULL
  eval_objective <- function(par, compute_posterior = TRUE) {
    log_range <- par[["log_range"]]
    log_sigma <- par[["log_sigma"]]
    beta_start <- if (identical(beta_mode, "empirical_bayes")) par[["beta"]] else beta_init
    compute_objective <- function() {
      .matern_laplace_known_noise_objective_at_params(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed,
        beta_prec = beta_prec,
        beta_init = beta_start,
        link = link,
        pc_penalty = pc_penalty,
        initial_mode = last_inner_mode,
        compute_posterior = compute_posterior,
        optimize_beta = FALSE,
        suppress_warnings = suppress_warnings,
        laplace_curvature = laplace_curvature
      )
    }
    objective <- if (suppress_warnings) suppressWarnings(compute_objective()) else compute_objective()
    if (!is.null(objective$mode) && all(is.finite(objective$mode))) {
      last_inner_mode <<- as.numeric(objective$mode)
    }
    objective$stats$spde_template <- spde_template
    objective$stats$alpha <- alpha
    objective$stats$d <- d
    objective
  }

  safe_objective <- function(par) {
    objective <- tryCatch(suppressWarnings(eval_objective(par, compute_posterior = FALSE)), error = function(e) e)
    if (inherits(objective, "error")) return(1e100)
    -objective$log_marginal
  }

  par0 <- if (identical(beta_mode, "empirical_bayes")) {
    c(log_range = theta_init, log_sigma = log(sigma_init), beta = beta_init)
  } else {
    c(log_range = theta_init, log_sigma = log(sigma_init))
  }
  opt_res <- .matern_optimize_or_eval_objective(
    par0 = par0,
    fixed_names = fixed_names,
    safe_objective = safe_objective,
    eval_objective = eval_objective,
    failure_label = "Matern Laplace optimization failed"
  )
  objective_opt <- opt_res$objective
  .matern_laplace_fit_from_objective(
    x = as.numeric(x),
    s = as.numeric(s),
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    objective = objective_opt,
    beta_mode = beta_mode,
    beta_prec = beta_prec,
    link = link,
    pc_penalty = pc_penalty,
    laplace_curvature = laplace_curvature
  )
}

.fit_matern_laplace_unknown_noise <- function(x,
                                              A,
                                              spde_template,
                                              alpha,
                                              d,
                                              theta_init,
                                              sigma_init,
                                              noise_sd_init,
                                              beta_init,
                                              beta_mode,
                                              beta_fixed = NULL,
                                              beta_prec = NULL,
                                              pc_penalty = NULL,
                                              link = c("identity", "log", "softplus"),
                                              suppress_warnings = TRUE,
                                              fix_params = character(),
                                              laplace_curvature = c("observed", "fisher")) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  fixed_names <- .matern_fixed_log_param_names(fix_params)
  last_inner_mode <- NULL
  eval_objective <- function(par, compute_posterior = TRUE) {
    log_range <- par[["log_range"]]
    log_sigma <- par[["log_sigma"]]
    log_noise_sd <- par[["log_noise"]]
    beta_start <- if (identical(beta_mode, "empirical_bayes")) par[["beta"]] else beta_init
    compute_objective <- function() {
      .matern_laplace_unknown_noise_objective_at_params(
        x = as.numeric(x),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma,
        log_noise_sd = log_noise_sd,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed,
        beta_prec = beta_prec,
        beta_init = beta_start,
        link = link,
        pc_penalty = pc_penalty,
        initial_mode = last_inner_mode,
        compute_posterior = compute_posterior,
        optimize_beta = FALSE,
        suppress_warnings = suppress_warnings,
        laplace_curvature = laplace_curvature
      )
    }
    objective <- if (suppress_warnings) suppressWarnings(compute_objective()) else compute_objective()
    if (!is.null(objective$mode) && all(is.finite(objective$mode))) {
      last_inner_mode <<- as.numeric(objective$mode)
    }
    objective$stats$spde_template <- spde_template
    objective$stats$alpha <- alpha
    objective$stats$d <- d
    objective
  }

  safe_objective <- function(par) {
    objective <- tryCatch(suppressWarnings(eval_objective(par, compute_posterior = FALSE)), error = function(e) e)
    if (inherits(objective, "error")) return(1e100)
    -objective$log_marginal
  }

  par0 <- if (identical(beta_mode, "empirical_bayes")) {
    c(log_range = theta_init, log_sigma = log(sigma_init), log_noise = log(noise_sd_init), beta = beta_init)
  } else {
    c(log_range = theta_init, log_sigma = log(sigma_init), log_noise = log(noise_sd_init))
  }
  opt_res <- .matern_optimize_or_eval_objective(
    par0 = par0,
    fixed_names = fixed_names,
    safe_objective = safe_objective,
    eval_objective = eval_objective,
    failure_label = "Matern Laplace learned-noise optimization failed"
  )
  objective_opt <- opt_res$objective
  .matern_laplace_fit_from_objective(
    x = as.numeric(x),
    s = as.numeric(objective_opt$s),
    A = A,
    spde_template = spde_template,
    alpha = alpha,
    d = d,
    objective = objective_opt,
    beta_mode = beta_mode,
    beta_prec = beta_prec,
    link = link,
    pc_penalty = pc_penalty,
    fitted_noise_sd = objective_opt$fitted_noise_sd,
    laplace_curvature = laplace_curvature
  )
}

.fit_matern_laplace_dispatch_known_noise <- function(x,
                                                     s,
                                                     A,
                                                     spde_template,
                                                     alpha,
                                                     d,
                                                     theta_init,
                                                     sigma_init,
                                                     beta_init,
                                                     beta_mode,
                                                     beta_fixed = NULL,
                                                     beta_prec = NULL,
                                                     pc_penalty = NULL,
                                                     link = c("identity", "log", "softplus"),
                                                     suppress_warnings = TRUE,
                                                     fix_g = FALSE,
                                                     matern_n_starts = 1L,
                                                     backend_use = c("laplace", "laplace_fisher", "laplace_r", "laplace_tmb"),
                                                     fix_params = character()) {
  link <- match.arg(link)
  backend_use <- match.arg(backend_use)
  if (isTRUE(fix_g)) fix_params <- unique(c(fix_params, "range", "sigma"))
  fixed_log_names <- .matern_fixed_log_param_names(fix_params)
  laplace_curvature <- if (identical(backend_use, "laplace_fisher")) "fisher" else "observed"
  tmb_unavailable_reason <- .matern_laplace_tmb_unavailable_reason(
    alpha = alpha,
    beta_mode = beta_mode,
    fixed_log_names = fixed_log_names
  )
  use_tmb <- switch(
    backend_use,
    laplace_fisher = is.null(tmb_unavailable_reason),
    laplace_tmb = is.null(tmb_unavailable_reason),
    laplace_r = FALSE,
    laplace = is.null(tmb_unavailable_reason)
  )

  fit_r_reference <- function(tmb_fallback_error = NULL) {
    fit <- .fit_matern_laplace_known_noise(
      x = x,
      s = s,
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      theta_init = theta_init,
      sigma_init = sigma_init,
      beta_init = beta_init,
      beta_mode = beta_mode,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      pc_penalty = pc_penalty,
      link = link,
      suppress_warnings = suppress_warnings,
      fix_params = fix_params,
      laplace_curvature = laplace_curvature
    )
    if (!is.null(tmb_fallback_error)) {
      fit$laplace_diagnostics$tmb_fallback_error <- tmb_fallback_error
    }
    fit
  }

  if (isTRUE(use_tmb)) {
    reason <- .matern_laplace_tmb_unsupported_reason(alpha = alpha, beta_mode = beta_mode)
    if (!is.null(reason)) stop(reason)
    fit_tmb <- tryCatch(.fit_matern_laplace_tmb_known_noise(
      x = x,
      s = s,
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      theta_init = theta_init,
      sigma_init = sigma_init,
      beta_init = beta_init,
      beta_mode = beta_mode,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      pc_penalty = pc_penalty,
      link = link,
      fix_g = all(c("range", "sigma") %in% fix_params),
      matern_n_starts = matern_n_starts,
      laplace_curvature = laplace_curvature
    ), error = function(e) e)
    tmb_invalid_reason <- .matern_laplace_tmb_invalid_reason(fit_tmb)
    if (is.null(tmb_invalid_reason)) {
      return(fit_tmb)
    }
    if (identical(backend_use, "laplace_tmb")) {
      stop(tmb_invalid_reason, call. = FALSE)
    }
    stop(.matern_laplace_tmb_failure_message(
      reason = tmb_invalid_reason,
      link = link,
      backend_use = backend_use,
      learn_noise = FALSE,
      fixed_log_names = fixed_log_names
    ), call. = FALSE)
  }

  if (identical(backend_use, "laplace_r")) {
    return(fit_r_reference())
  }
  if (identical(backend_use, "laplace_tmb")) {
    stop(.matern_laplace_tmb_unavailable_reason(
      alpha = alpha,
      beta_mode = beta_mode,
      fixed_log_names = fixed_log_names
    ), call. = FALSE)
  }
  stop(.matern_laplace_tmb_failure_message(
    reason = tmb_unavailable_reason,
    link = link,
    backend_use = backend_use,
    learn_noise = FALSE,
    fixed_log_names = fixed_log_names
  ), call. = FALSE)
}

.fit_matern_laplace_dispatch_unknown_noise <- function(x,
                                                       A,
                                                       spde_template,
                                                       alpha,
                                                       d,
                                                       theta_init,
                                                       sigma_init,
                                                       noise_sd_init,
                                                       beta_init,
                                                       beta_mode,
                                                       beta_fixed = NULL,
                                                       beta_prec = NULL,
                                                       pc_penalty = NULL,
                                                       link = c("identity", "log", "softplus"),
                                                       suppress_warnings = TRUE,
                                                       fix_g = FALSE,
                                                       matern_n_starts = 1L,
                                                       backend_use = c("laplace", "laplace_fisher", "laplace_r", "laplace_tmb"),
                                                       fix_params = character()) {
  link <- match.arg(link)
  backend_use <- match.arg(backend_use)
  if (isTRUE(fix_g)) fix_params <- unique(c(fix_params, "range", "sigma"))
  fixed_log_names <- .matern_fixed_log_param_names(fix_params)
  laplace_curvature <- if (identical(backend_use, "laplace_fisher")) "fisher" else "observed"
  tmb_unavailable_reason <- .matern_laplace_tmb_unavailable_reason(
    alpha = alpha,
    beta_mode = beta_mode,
    fixed_log_names = fixed_log_names
  )

  use_tmb <- switch(
    backend_use,
    laplace_fisher = is.null(tmb_unavailable_reason),
    laplace_tmb = is.null(tmb_unavailable_reason),
    laplace_r = FALSE,
    laplace = is.null(tmb_unavailable_reason)
  )

  fit_r_reference <- function(tmb_fallback_error = NULL) {
    fit <- .fit_matern_laplace_unknown_noise(
      x = x,
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      theta_init = theta_init,
      sigma_init = sigma_init,
      noise_sd_init = noise_sd_init,
      beta_init = beta_init,
      beta_mode = beta_mode,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      pc_penalty = pc_penalty,
      link = link,
      suppress_warnings = suppress_warnings,
      fix_params = fix_params,
      laplace_curvature = laplace_curvature
    )
    if (!is.null(tmb_fallback_error)) {
      fit$laplace_diagnostics$tmb_fallback_error <- tmb_fallback_error
    }
    fit
  }

  if (isTRUE(use_tmb)) {
    reason <- .matern_laplace_tmb_unsupported_reason(alpha = alpha, beta_mode = beta_mode)
    if (!is.null(reason)) stop(reason)
    fit_tmb <- tryCatch(.fit_matern_laplace_tmb_known_noise(
      x = x,
      s = rep(1, length(x)),
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      theta_init = theta_init,
      sigma_init = sigma_init,
      beta_init = beta_init,
      beta_mode = beta_mode,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      pc_penalty = pc_penalty,
      link = link,
      fix_g = all(c("range", "sigma") %in% fix_params),
      learn_noise = TRUE,
      noise_sd_init = noise_sd_init,
      matern_n_starts = matern_n_starts,
      laplace_curvature = laplace_curvature
    ), error = function(e) e)
    tmb_invalid_reason <- .matern_laplace_tmb_invalid_reason(fit_tmb)
    if (is.null(tmb_invalid_reason)) {
      return(fit_tmb)
    }
    if (identical(backend_use, "laplace_tmb")) {
      stop(tmb_invalid_reason, call. = FALSE)
    }
    stop(.matern_laplace_tmb_failure_message(
      reason = tmb_invalid_reason,
      link = link,
      backend_use = backend_use,
      learn_noise = TRUE,
      fixed_log_names = fixed_log_names
    ), call. = FALSE)
  }

  if (identical(backend_use, "laplace_r")) {
    return(fit_r_reference())
  }
  if (identical(backend_use, "laplace_tmb")) {
    stop(.matern_laplace_tmb_unavailable_reason(
      alpha = alpha,
      beta_mode = beta_mode,
      fixed_log_names = fixed_log_names
    ), call. = FALSE)
  }
  stop(.matern_laplace_tmb_failure_message(
    reason = tmb_unavailable_reason,
    link = link,
    backend_use = backend_use,
    learn_noise = TRUE,
    fixed_log_names = fixed_log_names
  ), call. = FALSE)
}

.fit_matern_exact_known_noise <- function(x,
                                          s,
                                          A,
                                          spde_template,
                                          alpha,
                                          d,
                                          theta_init,
                                          sigma_init,
                                          beta_mode,
                                          beta_fixed = NULL,
                                          beta_prec = NULL,
                                          pc_penalty = NULL,
                                          suppress_warnings = TRUE,
                                          fix_g = FALSE,
                                          fix_params = character()) {
  if (isTRUE(fix_g)) fix_params <- unique(c(fix_params, "range", "sigma"))
  fixed_names <- .matern_fixed_log_param_names(fix_params)
  raw_eval_objective <- function(par) {
    if (beta_mode == "empirical_bayes") {
      .exact_matern_profile_objective(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[["log_range"]],
        log_sigma = par[["log_sigma"]]
      )
    } else if (beta_mode == "fixed") {
      .exact_matern_fixed_objective(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[["log_range"]],
        log_sigma = par[["log_sigma"]],
        beta0 = beta_fixed
      )
    } else if (beta_mode == "prior_flat") {
      stats <- .exact_matern_sufficient_stats(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[["log_range"]],
        log_sigma = par[["log_sigma"]]
      )
      list(
        log_marginal = .exact_matern_loglik_integrated_flat_beta(stats, n_obs = length(x)),
        fitted_beta = stats$beta_profile_hat,
        stats = stats,
        beta_prec = 0
      )
    } else {
      .exact_matern_prior_objective(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[["log_range"]],
        log_sigma = par[["log_sigma"]],
        beta_prec = beta_prec
      )
    }
  }

  eval_objective <- function(par) {
    objective <- if (suppress_warnings) {
      suppressWarnings(raw_eval_objective(par))
    } else {
      raw_eval_objective(par)
    }
    objective$stats$theta_log_range <- par[["log_range"]]
    objective$stats$theta_log_sigma <- par[["log_sigma"]]
    .exact_matern_add_pc_prior(objective, pc_penalty = pc_penalty, d = d)
  }

  safe_objective <- function(par) {
    objective <- tryCatch(eval_objective(par), error = function(e) e)
    if (inherits(objective, "error")) return(1e100)
    -objective$log_marginal
  }

  par0 <- c(log_range = theta_init, log_sigma = log(sigma_init))
  opt_res <- .matern_optimize_or_eval_objective(
    par0 = par0,
    fixed_names = fixed_names,
    safe_objective = safe_objective,
    eval_objective = eval_objective,
    failure_label = "Exact Matern optimization failed"
  )
  objective_opt <- opt_res$objective
  state <- .matern_exact_state_with_beta_mode(
    x = as.numeric(x),
    s = as.numeric(s),
    A = A,
    objective = objective_opt,
    beta_mode = beta_mode,
    beta_prec = beta_prec
  )

  log_likelihood <- state$log_marginal
  class(log_likelihood) <- "logLik"

  posterior_sampler <- function(nsamp) {
    nsamp <- .check_single_numeric(nsamp, "nsamp")
    if (nsamp < 1 || nsamp != floor(nsamp)) {
      stop("`nsamp` must be a positive integer.")
    }

    if (beta_mode == "fixed" || beta_mode == "empirical_bayes") {
      samps_w <- LaplacesDemon::rmvnp(
        n = nsamp,
        mu = state$post_mean_latent,
        Omega = as.matrix(state$Q_post)
      )
      if (is.null(dim(samps_w))) {
        samps_w <- matrix(samps_w, nrow = 1)
      }
      obs_samps <- t(as.matrix(A %*% t(samps_w)))
      return(sweep(obs_samps, 2, state$fitted_beta, `+`))
    }

    Q_joint <- state$Q_joint
    A_aug <- cbind(A, Matrix::Matrix(rep(1, nrow(A)), ncol = 1, sparse = TRUE))
    joint_mean <- c(state$post_mean_latent, state$fitted_beta)
    samps_joint <- LaplacesDemon::rmvnp(
      n = nsamp,
      mu = joint_mean,
      Omega = as.matrix(Q_joint)
    )
    if (is.null(dim(samps_joint))) {
      samps_joint <- matrix(samps_joint, nrow = 1)
    }
    t(as.matrix(A_aug %*% t(samps_joint)))
  }

  list(
    posterior = state$posterior,
    fitted_g = state$fitted_g,
    fitted_beta = state$fitted_beta,
    beta_prec = beta_prec,
    beta_mode = beta_mode,
    log_likelihood = log_likelihood,
    log_likelihood_semantics = paste0("exact_", beta_mode),
    posterior_sampler = posterior_sampler,
    posterior_spatial_field = state$posterior_spatial_field,
    backend = "exact",
    link = "identity",
    prior_family = if (is.null(pc_penalty)) "identity_Matern" else "identity_Matern_pc_exact",
    pc_penalty = pc_penalty,
    log_likelihood_pc_prior_theta = objective_opt$log_pc_prior_theta,
    exact_optimization = .matern_exact_optimization_diagnostics(opt_res),
    matern_objective_context = list(
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d
    )
  )
}

.fit_matern_exact_unknown_noise <- function(x,
                                            A,
                                            spde_template,
                                            alpha,
                                            d,
                                            theta_init,
                                            sigma_init,
                                            noise_sd_init,
                                            beta_mode,
                                            beta_fixed = NULL,
                                            beta_prec = NULL,
                                            pc_penalty = NULL,
                                            suppress_warnings = TRUE,
                                            fix_g = FALSE,
                                            fix_params = character(),
                                            noise_scale = NULL) {
  if (isTRUE(fix_g)) fix_params <- unique(c(fix_params, "range", "sigma"))
  fixed_names <- .matern_fixed_log_param_names(fix_params)
  if (!is.null(noise_scale)) {
    noise_scale <- .exact_matern_unknown_noise_s(1, length(x), noise_scale)
  }
  raw_eval_objective <- function(par) {
    if (beta_mode == "empirical_bayes") {
      .exact_matern_profile_objective_unknown_noise(
        x = as.numeric(x),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[["log_range"]],
        log_sigma = par[["log_sigma"]],
        log_noise_sd = par[["log_noise"]],
        noise_scale = noise_scale
      )
    } else if (beta_mode == "fixed") {
      .exact_matern_fixed_objective_unknown_noise(
        x = as.numeric(x),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[["log_range"]],
        log_sigma = par[["log_sigma"]],
        beta0 = beta_fixed,
        log_noise_sd = par[["log_noise"]],
        noise_scale = noise_scale
      )
    } else if (beta_mode == "prior_flat") {
      noise_sd <- exp(par[["log_noise"]])
      s_eff <- .exact_matern_unknown_noise_s(noise_sd, length(x), noise_scale)
      stats <- .exact_matern_sufficient_stats(
        x = as.numeric(x),
        s = s_eff,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[["log_range"]],
        log_sigma = par[["log_sigma"]]
      )
      objective <- list(
        log_marginal = .exact_matern_loglik_integrated_flat_beta(stats, n_obs = length(x)),
        fitted_beta = stats$beta_profile_hat,
        fitted_noise_sd = noise_sd,
        s = s_eff,
        stats = stats,
        beta_prec = 0
      )
      if (!is.null(noise_scale)) objective$noise_scale <- as.numeric(noise_scale)
      objective
    } else {
      .exact_matern_prior_objective_unknown_noise(
        x = as.numeric(x),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = par[["log_range"]],
        log_sigma = par[["log_sigma"]],
        beta_prec = beta_prec,
        log_noise_sd = par[["log_noise"]],
        noise_scale = noise_scale
      )
    }
  }

  eval_objective <- function(par) {
    objective <- if (suppress_warnings) {
      suppressWarnings(raw_eval_objective(par))
    } else {
      raw_eval_objective(par)
    }
    objective$stats$theta_log_range <- par[["log_range"]]
    objective$stats$theta_log_sigma <- par[["log_sigma"]]
    .exact_matern_add_pc_prior(objective, pc_penalty = pc_penalty, d = d)
  }

  safe_objective <- function(par) {
    objective <- tryCatch(eval_objective(par), error = function(e) e)
    if (inherits(objective, "error")) return(1e100)
    -objective$log_marginal
  }

  par0 <- c(log_range = theta_init, log_sigma = log(sigma_init), log_noise = log(noise_sd_init))
  opt_res <- .matern_optimize_or_eval_objective(
    par0 = par0,
    fixed_names = fixed_names,
    safe_objective = safe_objective,
    eval_objective = eval_objective,
    failure_label = "Exact Matern learned-noise optimization failed"
  )
  objective_opt <- opt_res$objective
  state <- .matern_exact_state_with_beta_mode(
    x = as.numeric(x),
    s = as.numeric(objective_opt$s),
    A = A,
    objective = objective_opt,
    beta_mode = beta_mode,
    beta_prec = beta_prec
  )
  state$fitted_noise_sd <- objective_opt$fitted_noise_sd
  state$fitted_s <- as.numeric(objective_opt$s)
  if (!is.null(noise_scale)) state$noise_scale <- as.numeric(noise_scale)

  log_likelihood <- state$log_marginal
  class(log_likelihood) <- "logLik"

  posterior_sampler <- function(nsamp) {
    nsamp <- .check_single_numeric(nsamp, "nsamp")
    if (nsamp < 1 || nsamp != floor(nsamp)) {
      stop("`nsamp` must be a positive integer.")
    }

    if (beta_mode == "fixed" || beta_mode == "empirical_bayes") {
      samps_w <- LaplacesDemon::rmvnp(
        n = nsamp,
        mu = state$post_mean_latent,
        Omega = as.matrix(state$Q_post)
      )
      if (is.null(dim(samps_w))) {
        samps_w <- matrix(samps_w, nrow = 1)
      }
      obs_samps <- t(as.matrix(A %*% t(samps_w)))
      return(sweep(obs_samps, 2, state$fitted_beta, `+`))
    }

    A_aug <- cbind(A, Matrix::Matrix(rep(1, nrow(A)), ncol = 1, sparse = TRUE))
    joint_mean <- c(state$post_mean_latent, state$fitted_beta)
    samps_joint <- LaplacesDemon::rmvnp(
      n = nsamp,
      mu = joint_mean,
      Omega = as.matrix(state$Q_joint)
    )
    if (is.null(dim(samps_joint))) {
      samps_joint <- matrix(samps_joint, nrow = 1)
    }
    t(as.matrix(A_aug %*% t(samps_joint)))
  }

  list(
    posterior = state$posterior,
    fitted_g = state$fitted_g,
    fitted_beta = state$fitted_beta,
    fitted_noise_sd = state$fitted_noise_sd,
    fitted_s = state$fitted_s,
    beta_prec = beta_prec,
    beta_mode = beta_mode,
    log_likelihood = log_likelihood,
    log_likelihood_semantics = paste0("exact_", beta_mode),
    posterior_sampler = posterior_sampler,
    posterior_spatial_field = state$posterior_spatial_field,
    backend = "exact",
    link = "identity",
    prior_family = if (is.null(pc_penalty)) "identity_Matern_learned_noise" else "identity_Matern_pc_exact_learned_noise",
    pc_penalty = pc_penalty,
    log_likelihood_pc_prior_theta = objective_opt$log_pc_prior_theta,
    noise_scale = if (is.null(noise_scale)) NULL else as.numeric(noise_scale),
    exact_optimization = .matern_exact_optimization_diagnostics(opt_res),
    matern_objective_context = list(
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d
    )
  )
}

.exact_matern_known_noise_objective_at_params <- function(x,
                                                          s,
                                                          A,
                                                          spde_template,
                                                          alpha,
                                                          d,
                                                          log_range,
                                                          log_sigma,
                                                          beta_mode,
                                                          beta_fixed = NULL,
                                                          beta_prec = NULL,
                                                          pc_penalty = NULL,
                                                          suppress_warnings = TRUE) {
  raw_eval <- function() {
    if (beta_mode == "empirical_bayes") {
      .exact_matern_profile_objective(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma
      )
    } else if (beta_mode == "fixed") {
      .exact_matern_fixed_objective(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma,
        beta0 = beta_fixed
      )
    } else if (beta_mode == "prior_flat") {
      stats <- .exact_matern_sufficient_stats(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma
      )
      list(
        log_marginal = .exact_matern_loglik_integrated_flat_beta(stats, n_obs = length(x)),
        fitted_beta = stats$beta_profile_hat,
        stats = stats,
        beta_prec = 0
      )
    } else if (beta_mode == "prior_proper") {
      .exact_matern_prior_objective(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma,
        beta_prec = beta_prec
      )
    } else {
      stop("Unsupported beta mode: ", beta_mode)
    }
  }

  objective <- if (suppress_warnings) suppressWarnings(raw_eval()) else raw_eval()
  objective$stats$theta_log_range <- log_range
  objective$stats$theta_log_sigma <- log_sigma
  .exact_matern_add_pc_prior(objective, pc_penalty = pc_penalty, d = d)
}

.exact_matern_unknown_noise_objective_at_params <- function(x,
                                                            A,
                                                            spde_template,
                                                            alpha,
                                                            d,
                                                            log_range,
                                                            log_sigma,
                                                            log_noise_sd,
                                                            beta_mode,
                                                            beta_fixed = NULL,
                                                            beta_prec = NULL,
                                                            pc_penalty = NULL,
                                                            suppress_warnings = TRUE,
                                                            noise_scale = NULL) {
  raw_eval <- function() {
    if (beta_mode == "empirical_bayes") {
      .exact_matern_profile_objective_unknown_noise(
        x = as.numeric(x),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma,
        log_noise_sd = log_noise_sd,
        noise_scale = noise_scale
      )
    } else if (beta_mode == "fixed") {
      .exact_matern_fixed_objective_unknown_noise(
        x = as.numeric(x),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma,
        beta0 = beta_fixed,
        log_noise_sd = log_noise_sd,
        noise_scale = noise_scale
      )
    } else if (beta_mode == "prior_flat") {
      noise_sd <- exp(log_noise_sd)
      s_eff <- .exact_matern_unknown_noise_s(noise_sd, length(x), noise_scale)
      stats <- .exact_matern_sufficient_stats(
        x = as.numeric(x),
        s = s_eff,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma
      )
      objective <- list(
        log_marginal = .exact_matern_loglik_integrated_flat_beta(stats, n_obs = length(x)),
        fitted_beta = stats$beta_profile_hat,
        fitted_noise_sd = noise_sd,
        s = s_eff,
        stats = stats,
        beta_prec = 0
      )
      if (!is.null(noise_scale)) objective$noise_scale <- as.numeric(noise_scale)
      objective
    } else if (beta_mode == "prior_proper") {
      .exact_matern_prior_objective_unknown_noise(
        x = as.numeric(x),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = log_range,
        log_sigma = log_sigma,
        beta_prec = beta_prec,
        log_noise_sd = log_noise_sd,
        noise_scale = noise_scale
      )
    } else {
      stop("Unsupported beta mode: ", beta_mode)
    }
  }

  objective <- if (suppress_warnings) suppressWarnings(raw_eval()) else raw_eval()
  objective$stats$theta_log_range <- log_range
  objective$stats$theta_log_sigma <- log_sigma
  .exact_matern_add_pc_prior(objective, pc_penalty = pc_penalty, d = d)
}

.validate_matern_setup <- function(setup) {
  required_names <- c("locations", "d", "mesh", "A", "spde_template", "alpha", "max.edge", "penalty_range")
  missing_names <- setdiff(required_names, names(setup))
  if (length(missing_names) > 0L) {
    stop("`setup` is missing required Matern components: ", paste(missing_names, collapse = ", "), ".")
  }

  loc_mat <- as.matrix(setup$locations)
  if (!ncol(loc_mat) %in% c(1L, 2L)) {
    stop("`setup$locations` must be an n x 1 or n x 2 numeric matrix.")
  }

  if (!is.numeric(setup$alpha) || length(setup$alpha) != 1L || is.na(setup$alpha)) {
    stop("`setup$alpha` must be a single numeric value.")
  }
  if (setup$alpha <= setup$d / 2) {
    stop("`setup$alpha` must satisfy alpha > d / 2.")
  }

  setup$locations <- loc_mat
  setup$d <- as.integer(setup$d)
  setup$A <- Matrix::Matrix(setup$A, sparse = TRUE)
  if (nrow(setup$A) != nrow(loc_mat)) {
    stop("`setup$A` must have one row per location.")
  }
  setup$penalty_range <- as.numeric(setup$penalty_range)
  class(setup) <- unique(c("Matern_setup", class(setup)))
  setup
}

.resolve_matern_setup <- function(locations = NULL,
                                  setup = NULL,
                                  max.edge = NULL,
                                  alpha = 2,
                                  suppress_warnings = TRUE,
                                  penalty_range = NULL) {
  if (is.null(setup)) {
    if (is.null(locations)) {
      stop("Either `locations` or `setup` must be supplied.")
    }
    out <- Matern_setup(
      locations = locations,
      max.edge = max.edge,
      alpha = alpha,
      suppress_warnings = suppress_warnings
    )
  } else {
    if (!is.null(locations)) {
      stop("Supply either `locations` or `setup`, not both.")
    }
    if (!is.null(max.edge)) {
      stop("`max.edge` cannot be supplied when `setup` is provided.")
    }
    out <- .validate_matern_setup(setup)
    if (!isTRUE(all.equal(as.numeric(alpha), as.numeric(out$alpha)))) {
      stop("When `setup` is provided, `alpha` must match `setup$alpha`.")
    }
  }

  if (!is.null(penalty_range)) {
    if (!is.numeric(penalty_range) || length(penalty_range) != 1L || is.na(penalty_range) || penalty_range <= 0) {
      stop("`penalty_range` must be a single positive number (or NULL).")
    }
    out$penalty_range <- as.numeric(penalty_range)
  }

  out
}

#' Build reusable mesh/projector objects for Matern smoothing
#'
#' @description
#' Precomputes the location-dependent objects shared by the Matern smoothers:
#' the normalized location matrix, projector matrix \code{A}, mesh, SPDE
#' template, and a default penalty-range scale.
#'
#' The returned object can be reused across repeated fits with
#' \code{\link{ebnm_Matern_generator}} or \code{\link{eb_smoother}}, which is
#' particularly useful when smoothing many vectors on the same spatial grid.
#'
#' @param locations Numeric vector, matrix, or data frame of one- or
#'   two-dimensional locations.
#' @param max.edge Optional mesh-resolution control passed to INLA mesh
#'   construction. For one-dimensional locations this may be a scalar; for
#'   two-dimensional locations it may be \code{NULL}, length 1, or length 2.
#'   When \code{NULL}, the two-dimensional default uses a coarser outer mesh
#'   while keeping observed locations as mesh vertices.
#' @param alpha Smoothness order for the SPDE representation. Must satisfy
#'   \code{alpha > d / 2}.
#' @param suppress_warnings Logical scalar. If \code{TRUE}, suppresses warnings
#'   from INLA mesh/SPDE construction.
#'
#' @return A list of class \code{"Matern_setup"} containing the normalized
#'   locations, spatial dimension, mesh, projector matrix, SPDE template,
#'   \code{alpha}, \code{max.edge}, and a default \code{penalty_range}.
#'
#' @export
Matern_setup <- function(locations,
                         max.edge = NULL,
                         alpha = 2,
                         suppress_warnings = TRUE) {
  .with_quiet_inla_defaults({
    loc_info <- .normalize_locations(locations)
    loc_mat <- loc_info$loc
    d <- loc_info$d

    if (alpha <= d / 2) {
      stop("`alpha` must satisfy alpha > d / 2.")
    }

    meshA <- if (suppress_warnings) {
      suppressWarnings(.build_mesh_A(loc_mat, max.edge = max.edge))
    } else {
      .build_mesh_A(loc_mat, max.edge = max.edge)
    }

    spde_template <- if (suppress_warnings) {
      suppressWarnings(INLA::inla.spde2.matern(mesh = meshA$mesh, alpha = alpha))
    } else {
      INLA::inla.spde2.matern(mesh = meshA$mesh, alpha = alpha)
    }

    structure(
      list(
        locations = loc_mat,
        d = d,
        mesh = meshA$mesh,
        A = meshA$A,
        spde_template = spde_template,
        alpha = as.numeric(alpha),
        max.edge = max.edge,
        penalty_range = .default_penalty_range(loc_mat)
      ),
      class = c("Matern_setup", "list")
    )
  })
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
#' @param beta Optional scalar intercept state stored alongside the Matern
#'   hyperparameters.
#' @param beta_prec Optional non-negative scalar beta prior precision stored
#'   alongside the Matern hyperparameters.
#'
#' @return An object of class \code{"Matern"}.
#'
#' @export
Matern <- function(theta = NULL, sigma = 1, beta = NULL, beta_prec = NULL) {
  if (!is.null(theta)) {
    theta <- .check_single_numeric(theta, "theta")
  }
  sigma <- .check_single_numeric(sigma, "sigma")
  if (sigma <= 0) {
    stop("`sigma` must be positive.")
  }

  beta <- .check_optional_beta_vector(beta, "beta", expected_length = 1L)
  beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")

  structure(
    list(
      theta = theta,
      sigma = sigma,
      beta = beta,
      beta_prec = beta_prec
    ),
    class = "Matern"
  )
}



#' Generate an `ebnm` function for Matern smoothing
#'
#' @description
#' Returns an \code{ebnm}-compatible function that fits a stationary Matern
#' smoother in one or two spatial dimensions.
#'
#' This is the recommended public interface when the observation standard
#' errors \code{s} are known and you want a spatial \code{ebnm} workflow. If
#' the standard errors are unknown and you want to learn one common noise SD
#' instead, use \code{\link{eb_smoother}} with \code{family = "matern"} and
#' \code{s = NULL}.
#'
#' - If \code{locations} is a numeric vector, it is treated as one-dimensional
#'   locations.
#' - If \code{locations} is a matrix or data frame, it must have one or two
#'   columns.
#' - Alternatively, callers may supply a reusable \code{setup} object returned
#'   by \code{\link{Matern_setup}}.
#'
#' The identity link uses the exact Gaussian backend by default. The log link
#' supports positive mean functions by fitting
#' \eqn{x_i \sim N(\exp(\beta_0 + (Aw)_i), s_i^2)} with a sparse Fisher
#' Laplace backend by default, and an independent INLA backend is available for
#' cross-checking. INLA empirical-Bayes beta fits are implemented by external
#' profiling over fixed-beta INLA offset fits. The softplus link fits
#' \eqn{x_i \sim N(\log(1 + \exp(\beta_0 + (Aw)_i)), s_i^2)} with the package
#' Laplace backends; INLA does not support this link.
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
#' The public beta semantics are:
#' \itemize{
#'   \item \code{beta_fixed}: treat \eqn{\beta_0} as known.
#'   \item \code{beta_prec = NULL}: estimate \eqn{\beta_0} by empirical Bayes.
#'   \item \code{beta_prec = 0}: use a flat prior on \eqn{\beta_0}.
#'   \item \code{beta_prec > 0}: use a proper zero-mean Gaussian prior on \eqn{\beta_0}.
#' }
#'
#' For identity-link fits, the primary \code{log_likelihood} is the package
#' exact Gaussian objective for every backend. INLA-backed identity fits report
#' that exact objective evaluated at the INLA mode, with INLA's native Step A
#' quantities retained as diagnostics.
#'
#' For log-link and softplus-link fits, \code{backend = "auto"} uses
#' \code{"fisher_pql"}. Fisher-PQL builds repeated pseudo-Gaussian exact Matern
#' Step A fits and uses \code{pql_inner_iter} pseudo-response updates. It
#' reports a Fisher/Laplace score evaluated at the final PQL mode; this is an
#' approximate PQL-mode score, not a true re-optimized original-model Laplace
#' marginal likelihood. Explicit \code{backend = "laplace"} and
#' \code{"laplace_tmb"} retain observed-Hessian Laplace semantics through the
#' package-owned TMB implementation. Explicit \code{backend = "laplace_fisher"}
#' keeps the conditional mode equal to the mode of the original log posterior
#' and replaces only the Laplace posterior precision/log-determinant observation
#' curvature by the Fisher/Gauss-Newton curvature. For log-link normal
#' observations this replaces
#' \eqn{\exp(\eta_i)(2\exp(\eta_i)-x_i)/s_i^2} by
#' \eqn{\exp(2\eta_i)/s_i^2}; softplus uses the analogous squared first
#' derivative curvature. Fisher-Laplace fits report
#' \code{backend = "laplace_fisher"}, \code{laplace_curvature = "fisher"}, and
#' \code{log_likelihood_semantics = "laplace_fisher_<beta_mode>"}.
#' INLA-backed log-link fits report the observed-Hessian Laplace objective
#' evaluated at the INLA mode, with INLA's native marginal likelihood
#' quantities retained as diagnostics.
#' Inlabru-backed softplus fits likewise report the package observed-Hessian
#' Laplace objective evaluated at the inlabru fitted hyperparameters, with
#' inlabru's native marginal likelihood retained as a diagnostic.
#'
#' For non-identity links, posterior moments are reported on the response scale.
#' The log link uses exact log-normal moment formulas under the marginal
#' Gaussian Laplace posterior for \eqn{\eta_i}; the softplus link uses fixed
#' deterministic Gauss-Hermite quadrature for the same marginal Gaussian
#' approximation.
#'
#' @param locations Optional raw locations used to build the Matern mesh and
#'   projector matrix. Supply either \code{locations} or \code{setup}, but not
#'   both.
#' @param setup Optional reusable object returned by
#'   \code{\link{Matern_setup}}. Supply either \code{locations} or
#'   \code{setup}, but not both.
#' @param max.edge Mesh maximum edge length.
#' @param alpha Smoothness parameter used in the SPDE representation. It is held
#'   fixed and is not optimized.
#' @param suppress_warnings If \code{TRUE}, suppress warnings from INLA mesh and
#'   SPDE construction.
#' @param penalty_range Initial range anchor used when \code{g_init$theta} is
#'   missing. If \code{NULL}, defaults to roughly one-tenth of the spatial
#'   extent.
#' @param pc.penalty Optional list specifying PC priors for the Matern
#'   hyperparameters.
#'   Supported entries are \code{range} and \code{sigma}. Each entry may be a
#'   numeric vector of length 1, interpreted as the anchor with default tail
#'   probability \code{0.5}, or length 2, interpreted as \code{c(anchor, alpha)}.
#'   Missing entries default to a data-driven anchor and tail probability
#'   \code{0.5}. If \code{NULL}, the Matern hyperparameters are fit by
#'   unpenalized EB.
#' @param compute_exact_diagnostic Logical scalar retained for compatibility.
#'   Identity-link INLA fits now store the package exact objective at the INLA
#'   mode as their primary \code{log_likelihood}; log-link INLA fits store the
#'   package Laplace objective at the INLA mode.
#' @param backend Backend choice. \code{"auto"} uses \code{"exact"} for the
#'   identity link and \code{"fisher_pql"} for the log and softplus links.
#'   \code{"laplace_fisher"} and \code{"fisher_pql"} are available for
#'   \code{link = "log"} and \code{link = "softplus"}.
#'   \code{"laplace"} uses the observed-Hessian TMB implementation and errors
#'   if that TMB fit is unavailable or invalid. \code{"laplace_fisher"} uses the same TMB mode and
#'   hyperparameter fit when supported, but builds the returned posterior
#'   approximation and Fisher Laplace score with Fisher/Gauss-Newton curvature.
#'   \code{"fisher_pql"} is an approximate pseudo-likelihood backend that uses
#'   \code{pql_inner_iter} Fisher/PQL pseudo-Gaussian exact Matern Step A
#'   updates. It reports a Fisher/Laplace score evaluated at the final PQL mode;
#'   this is an approximate PQL-mode score, not a true re-optimized
#'   original-model Laplace marginal likelihood.
#'   The R reference backend \code{"laplace_r"} is accepted for internal
#'   validation only and is never reached automatically from public backends.
#'   \code{"inla"} is the independent INLA implementation and uses the supplied
#'   \code{pc.penalty}, if any, as the source of PC-prior penalization.
#'   \code{"inlabru"} runs the SPDE model through \code{inlabru}'s
#'   iterative-linearised INLA method and is currently only supported for
#'   \code{link = "softplus"}; it provides an opt-in alternative to the
#'   default Laplace path for cross-validation and benchmarking, not an expected
#'   acceleration path. For known-noise fits, \code{pc.penalty = NULL} keeps the
#'   non-PC SPDE parameterisation so the prior semantics match the manual
#'   Laplace default. For learned-noise fits with \code{s = NULL},
#'   \code{"inlabru"} requires an explicit \code{pc.penalty} list with
#'   \code{range}, \code{sigma}, and \code{noise} entries.
#'   Compatibility aliases \code{"laplace_tmb"} and \code{"inla_pc"} are still
#'   accepted.
#' @param pql_inner_iter Positive integer number of Fisher/PQL pseudo-Gaussian
#'   Step A updates used by \code{backend = "fisher_pql"} and by
#'   \code{backend = "auto"} when it resolves to Fisher-PQL for log or softplus
#'   Matern fits. The default is \code{3}.
#' @param link Link function. \code{"identity"} fits unconstrained Gaussian
#'   means. \code{"log"} and \code{"softplus"} fit positive mean functions and
#'   report posterior summaries on the positive response scale. \code{"logit"}
#'   and \code{"probit"} are reserved and currently unsupported.
#'
#' @return A function with signature
#'   \code{function(x, s, g_init = NULL, fix_g = FALSE, beta_fixed = NULL, beta_prec = NULL, fix_params = character(), output = NULL)}.
#'   The returned object includes posterior summaries, fitted hyperparameters,
#'   fitted/fixed \code{beta}, the resolved \code{g_init}, backend-specific
#'   likelihood diagnostics, the mesh, and a summary of the latent spatial
#'   field.
#'
#' The returned function accepts \code{fix_params} for partial parameter
#' fixing. Allowed values are \code{"range"}, \code{"sigma"}, and
#' \code{"beta"}. Fixed range uses \code{g_init$theta}, fixed sigma uses
#' \code{g_init$sigma}, and fixed beta uses \code{beta_fixed} when supplied or
#' otherwise \code{g_init$beta}. \code{fix_g = TRUE} is retained as a shortcut
#' for \code{fix_params = c("range", "sigma")}. Range/sigma fixing is not
#' supported by the INLA backend.
#'
#' @export
ebnm_Matern_generator <- function(locations = NULL,
                                  setup = NULL,
                                  max.edge = NULL,
                                  alpha = 2,
                                  suppress_warnings = TRUE,
                                  penalty_range = NULL,
                                  pc.penalty = NULL,
                                  compute_exact_diagnostic = FALSE,
                                  backend = c("auto", "exact", "laplace", "laplace_fisher", "fisher_pql", "inla", "inlabru"),
                                  pql_inner_iter = 3L,
                                  link = c("identity", "log", "softplus", "logit", "probit")) {

  backend <- .match_matern_backend_arg(backend)
  link <- match.arg(link)
  if (!link %in% c("identity", "log", "softplus")) {
    stop("The Matern implementation currently supports `link = \"identity\"`, `link = \"log\"`, and `link = \"softplus\"` only.")
  }
  pql_inner_iter <- .check_matern_pql_inner_iter(pql_inner_iter)

  setup0 <- .resolve_matern_setup(
    locations = locations,
    setup = setup,
    max.edge = max.edge,
    alpha = alpha,
    suppress_warnings = suppress_warnings,
    penalty_range = penalty_range
  )
  loc_mat <- setup0$locations
  d <- setup0$d
  alpha <- setup0$alpha
  if (!is.logical(compute_exact_diagnostic) || length(compute_exact_diagnostic) != 1L || is.na(compute_exact_diagnostic)) {
    stop("`compute_exact_diagnostic` must be TRUE or FALSE.")
  }
  mesh <- setup0$mesh
  A <- setup0$A
  spde_template <- setup0$spde_template
  penalty_range0 <- setup0$penalty_range

  ebnm_Matern <- function(x,
                          s,
                          g_init = NULL,
                          fix_g = FALSE,
                          beta_fixed = NULL,
                          beta_prec = NULL,
                          fix_params = character(),
                          output = NULL) {
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
    if (!is.logical(fix_g) || length(fix_g) != 1L || is.na(fix_g)) {
      stop("`fix_g` must be TRUE or FALSE.")
    }
    beta_prec_supplied <- !missing(beta_prec)
    fix_params_use <- .normalize_fix_params(
      fix_params = fix_params,
      allowed = c("range", "sigma", "beta"),
      fix_g = fix_g,
      fix_g_params = c("range", "sigma"),
      nm = "fix_params"
    )
    .validate_matern_fixed_param_values(fix_params_use, g_init)
    beta_fixed <- .resolve_fixed_beta_from_fix_params(
      fix_params = fix_params_use,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      beta_prec_supplied = beta_prec_supplied,
      g_init = g_init,
      expected_length = 1L
    )

    beta_spec <- .eb_smoother_resolve_beta_spec(
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      g_init_beta_prec = if (is.null(g_init)) NULL else g_init$beta_prec
    )
    beta_mode <- beta_spec$mode
    beta0_fixed <- if (is.null(beta_spec$beta_fixed)) {
      NULL
    } else {
      .check_single_numeric(beta_spec$beta_fixed, "beta_fixed")
    }
    beta_prec_use <- beta_spec$beta_prec

    resolved_init <- .resolve_matern_g_init(
      x = x,
      s = s,
      g_init = g_init,
      beta_fixed = beta0_fixed,
      beta_prec = beta_prec_use,
      penalty_range0 = penalty_range0,
      pc.penalty = pc.penalty,
      allow_noise = FALSE,
      link = link
    )
    theta_init <- resolved_init$theta_init
    sigma_init <- resolved_init$sigma_init
    g_init_resolved <- resolved_init$g_init
    pc_penalty0 <- resolved_init$pc_penalty
    backend_use <- .matern_resolve_backend(
      backend = backend,
      link = link,
      d = d,
      beta_mode = beta_mode,
      pc_penalty = pc_penalty0,
      fix_g = all(c("range", "sigma") %in% fix_params_use)
    )
    if (identical(backend_use, "inla") &&
        any(c("range", "sigma") %in% fix_params_use)) {
      stop("`fix_params` values `range` and `sigma` are not supported with the INLA Matern backend.")
    }
    if (identical(backend_use, "inla") && identical(link, "softplus")) {
      stop("`link = \"softplus\"` is not currently supported by the INLA Matern backend; use `backend = \"laplace\"` or `backend = \"laplace_fisher\"`.")
    }

    exact_fit_from_state <- function(state,
                                     log_likelihood_semantics,
                                     prior_family,
                                     log_likelihood_pc_prior_theta = NULL) {
      log_likelihood <- state$log_marginal
      class(log_likelihood) <- "logLik"

      posterior_sampler <- function(nsamp) {
        nsamp <- .check_single_numeric(nsamp, "nsamp")
        if (nsamp < 1 || nsamp != floor(nsamp)) {
          stop("`nsamp` must be a positive integer.")
        }

        if (beta_mode %in% c("fixed", "empirical_bayes")) {
          samps_w <- LaplacesDemon::rmvnp(
            n = nsamp,
            mu = state$post_mean_latent,
            Omega = as.matrix(state$Q_post)
          )
          if (is.null(dim(samps_w))) {
            samps_w <- matrix(samps_w, nrow = 1)
          }
          obs_samps <- t(as.matrix(A %*% t(samps_w)))
          return(sweep(obs_samps, 2, state$fitted_beta, `+`))
        }

        A_aug <- cbind(A, Matrix::Matrix(rep(1, nrow(A)), ncol = 1, sparse = TRUE))
        joint_mean <- c(state$post_mean_latent, state$fitted_beta)
        samps_joint <- LaplacesDemon::rmvnp(
          n = nsamp,
          mu = joint_mean,
          Omega = as.matrix(state$Q_joint)
        )
        if (is.null(dim(samps_joint))) {
          samps_joint <- matrix(samps_joint, nrow = 1)
        }
        t(as.matrix(A_aug %*% t(samps_joint)))
      }

      structure(
        list(
          posterior = state$posterior,
          fitted_g = state$fitted_g,
          fitted_beta = state$fitted_beta,
          beta_prec = beta_prec_use,
          beta_mode = beta_mode,
          g_init = g_init_resolved,
          log_likelihood = log_likelihood,
          log_likelihood_semantics = log_likelihood_semantics,
          posterior_sampler = posterior_sampler,
          data = data.frame(x = x, s = s),
          prior_family = prior_family,
          posterior_spatial_field = state$posterior_spatial_field,
          mesh = mesh,
          inla_result = NULL,
          backend = "exact",
          link = "identity",
          pc_penalty = pc_penalty0,
          log_likelihood_pc_prior_theta = log_likelihood_pc_prior_theta,
          matern_objective_context = list(
            A = A,
            spde_template = spde_template,
            alpha = alpha,
            d = d
          )
        ),
        class = c("list", "ebnm")
      )
    }

    if (identical(backend_use, "fisher_pql")) {
      fit <- .fit_matern_fisher_pql_known_noise(
        x = x,
        s = s,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_init,
        sigma_init = sigma_init,
        beta_init = resolved_init$beta_init,
        beta_mode = beta_mode,
        beta_fixed = beta0_fixed,
        beta_prec = beta_prec_use,
        pc_penalty = pc_penalty0,
        link = link,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        fix_params = fix_params_use,
        pql_max_iter = pql_inner_iter
      )
      fit$data <- data.frame(x = x, s = s)
      fit$mesh <- mesh
      fit$inla_result <- NULL
      fit$g_init <- g_init_resolved
      return(structure(fit, class = c("list", "ebnm")))
    }

    if (backend_use %in% c("laplace", "laplace_fisher", "laplace_r", "laplace_tmb")) {
      fit <- .fit_matern_laplace_dispatch_known_noise(
        x = x,
        s = s,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_init,
        sigma_init = sigma_init,
        beta_init = resolved_init$beta_init,
        beta_mode = beta_mode,
        beta_fixed = beta0_fixed,
        beta_prec = beta_prec_use,
        pc_penalty = pc_penalty0,
        link = link,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        backend_use = backend_use,
        fix_params = fix_params_use
      )
      fit$data <- data.frame(x = x, s = s)
      fit$mesh <- mesh
      fit$inla_result <- NULL
      fit$g_init <- g_init_resolved
      return(structure(fit, class = c("list", "ebnm")))
    }

    if (identical(backend_use, "inlabru")) {
      if (any(c("range", "sigma") %in% fix_params_use)) {
        stop("`fix_params` values `range` and `sigma` are not supported with the inlabru Matern backend.")
      }
      pc_penalty_inlabru <- .matern_inlabru_pc_penalty_policy(
        pc_penalty_arg = pc.penalty,
        resolved_pc_penalty = pc_penalty0,
        learn_noise = FALSE
      )
      theta_init_inlabru <- .matern_inla_theta_init(
        log_range = theta_init,
        log_sigma = log(sigma_init),
        alpha = alpha,
        d = d,
        pc_penalty = pc_penalty_inlabru
      )
      ibru_fit <- if (identical(beta_mode, "fixed")) {
        .fit_matern_inlabru_stepA_fixed_beta(
          x = as.numeric(x),
          s = as.numeric(s),
          A = A,
          mesh = mesh,
          alpha = alpha,
          d = d,
          locations = loc_mat,
          theta_init = theta_init_inlabru,
          beta_fixed = beta0_fixed,
          pc_penalty = pc_penalty_inlabru,
          link = link,
          suppress_warnings = suppress_warnings
        )
      } else {
        .fit_matern_inlabru_stepA(
          x = as.numeric(x),
          s = as.numeric(s),
          A = A,
          mesh = mesh,
          alpha = alpha,
          d = d,
          locations = loc_mat,
          theta_init = theta_init_inlabru,
          beta_init = resolved_init$beta_init,
          beta_prec = if (identical(beta_mode, "prior_proper")) beta_prec_use else 0,
          pc_penalty = pc_penalty_inlabru,
          link = link,
          suppress_warnings = suppress_warnings
        )
      }

      comparable_objective <- .matern_inlabru_known_laplace_objective(
        fit = ibru_fit,
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        beta_mode = beta_mode,
        beta_fixed = beta0_fixed,
        beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
        link = link,
        pc_penalty = pc_penalty_inlabru,
        suppress_warnings = suppress_warnings
      )
      log_likelihood <- as.numeric(comparable_objective$log_marginal)
      class(log_likelihood) <- "logLik"

      out <- list(
        posterior = ibru_fit$posterior,
        fitted_g = ibru_fit$fitted_g,
        fitted_beta = ibru_fit$fitted_beta,
        beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
        beta_mode = beta_mode,
        g_init = g_init_resolved,
        log_likelihood = log_likelihood,
        log_likelihood_semantics = paste0("laplace_at_inlabru_params_", beta_mode),
        posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
        data = data.frame(x = x, s = s),
        prior_family = paste0(link, "_Matern_inlabru"),
        posterior_spatial_field = ibru_fit$posterior_spatial_field,
        mesh = mesh,
        inla_result = ibru_fit$result,
        backend = "inlabru",
        link = link,
        pc_penalty = pc_penalty_inlabru,
        log_likelihood_laplace_at_inlabru_params = as.numeric(comparable_objective$log_marginal),
        log_likelihood_inlabru_mlik_integration = ibru_fit$log_likelihood_inlabru_mlik_integration,
        log_likelihood_stepA_mlik_integration = ibru_fit$log_likelihood_stepA_mlik_integration,
        beta_profile_optimization = ibru_fit$beta_profile_optimization,
        beta_profile_objective = ibru_fit$beta_profile_objective,
        matern_objective_context = list(
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d
        )
      )
      return(structure(out, class = c("list", "ebnm")))
    }

    if (backend_use == "exact") {
      if (!identical(link, "identity")) {
        stop("`backend = \"exact\"` is only available for `link = \"identity\"`; use `backend = \"laplace\"` for log-link Matern fits.")
      }
      fit <- .fit_matern_exact_known_noise(
        x = x,
        s = s,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_init,
        sigma_init = sigma_init,
        beta_mode = beta_mode,
        beta_fixed = beta0_fixed,
        beta_prec = beta_prec_use,
        pc_penalty = pc_penalty0,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        fix_params = fix_params_use
      )
      fit$data <- data.frame(x = x, s = s)
      fit$mesh <- mesh
      fit$inla_result <- NULL
      fit$g_init <- g_init_resolved
      return(structure(fit, class = c("list", "ebnm")))
    }

    theta_init_pc <- .matern_inla_theta_init(
      log_range = theta_init,
      log_sigma = log(sigma_init),
      alpha = alpha,
      d = d,
      pc_penalty = pc_penalty0
    )
    if (beta_mode == "fixed") {
      pc_fit <- .fit_matern_inla_stepA_fixed_beta(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        mesh = mesh,
        alpha = alpha,
        d = d,
        theta_init = theta_init_pc,
        beta_fixed = beta0_fixed,
        pc_penalty = pc_penalty0,
        link = link,
        suppress_warnings = suppress_warnings
      )
    } else if (beta_mode == "empirical_bayes") {
      pc_fit <- .fit_matern_inla_stepA_profile_beta(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        mesh = mesh,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_init_pc,
        beta_init = resolved_init$beta_init,
        pc_penalty = pc_penalty0,
        link = link,
        suppress_warnings = suppress_warnings
      )
    } else {
      pc_fit <- .fit_matern_inla_stepA(
        x = as.numeric(x),
        s = as.numeric(s),
        A = A,
        mesh = mesh,
        alpha = alpha,
        d = d,
        theta_init = theta_init_pc,
        beta_prec = if (beta_mode == "prior_flat") 0 else beta_prec_use,
        pc_penalty = pc_penalty0,
        link = link,
        suppress_warnings = suppress_warnings
      )
    }

    comparable_at_inla_mode <- NULL
    if (identical(link, "log")) {
      comparable_at_inla_mode <- .matern_laplace_known_noise_objective_at_params(
        x = x,
        s = s,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = pc_fit$fitted_g$theta,
        log_sigma = log(pc_fit$fitted_g$sigma),
        beta_mode = beta_mode,
        beta_fixed = beta0_fixed,
        beta_prec = beta_prec_use,
        beta_init = pc_fit$fitted_beta,
        link = link,
        pc_penalty = pc_penalty0,
        initial_mode = .matern_inla_laplace_initial_mode(pc_fit, beta_mode),
        compute_posterior = FALSE,
        optimize_beta = FALSE,
        suppress_warnings = suppress_warnings
      )
    } else {
      comparable_at_inla_mode <- .exact_matern_known_noise_objective_at_params(
        x = x,
        s = s,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = pc_fit$fitted_g$theta,
        log_sigma = log(pc_fit$fitted_g$sigma),
        beta_mode = if (identical(beta_mode, "empirical_bayes")) "fixed" else beta_mode,
        beta_fixed = if (identical(beta_mode, "empirical_bayes")) pc_fit$fitted_beta else beta0_fixed,
        beta_prec = beta_prec_use,
        pc_penalty = pc_penalty0,
        suppress_warnings = suppress_warnings
      )
    }

    if (identical(link, "identity")) {
      log_likelihood_exact_at_stepA_mode <- comparable_at_inla_mode$log_marginal
    } else {
      log_likelihood_exact_at_stepA_mode <- NULL
    }

    log_likelihood <- comparable_at_inla_mode$log_marginal
    class(log_likelihood) <- "logLik"

    out <- list(
      posterior = pc_fit$posterior,
      fitted_g = pc_fit$fitted_g,
      fitted_beta = pc_fit$fitted_beta,
      beta_prec = beta_prec_use,
      beta_mode = beta_mode,
      g_init = g_init_resolved,
      log_likelihood = log_likelihood,
      log_likelihood_semantics = if (identical(link, "log")) {
        paste0("laplace_at_inla_mode_", beta_mode)
      } else {
        paste0("exact_at_inla_mode_", beta_mode)
      },
      posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
      data = data.frame(x = x, s = s),
      prior_family = paste0(link, "_Matern", if (is.null(pc_penalty0)) "_inla" else "_pc_inla"),
      posterior_spatial_field = pc_fit$posterior_spatial_field,
      mesh = mesh,
      inla_result = pc_fit$result,
      backend = "inla",
      link = link,
      pc_penalty = pc_penalty0,
      log_likelihood_stepA_penalized = pc_fit$log_likelihood_stepA_penalized,
      log_likelihood_stepA_mlik_integration = pc_fit$log_likelihood_stepA_mlik_integration,
      log_likelihood_stepA_mlik_gaussian = pc_fit$log_likelihood_stepA_mlik_gaussian,
      log_likelihood_stepA_joint_log_posterior = pc_fit$log_likelihood_stepA_joint_log_posterior,
      log_likelihood_laplace_at_inla_mode = if (identical(link, "log")) comparable_at_inla_mode$log_marginal else NULL,
      log_likelihood_exact_at_stepA_mode = log_likelihood_exact_at_stepA_mode,
      log_likelihood_stepB = NULL,
      log_likelihood_pc_prior_theta = NULL,
      beta_profile_optimization = pc_fit$beta_profile_optimization,
      beta_profile_objective = pc_fit$beta_profile_objective,
      matern_objective_context = list(
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d
      )
    )

    structure(out, class = c("list", "ebnm"))
  }

  if (suppress_warnings) suppressWarnings(ebnm_Matern) else ebnm_Matern
}
