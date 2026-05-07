.normalize_eb_smoother_s <- function(s, n) {
  if (!(length(s) == 1L || length(s) == n)) {
    stop("`s` must have length 1 or length(x).")
  }
  if (length(s) == 1L) s <- rep(s, n)
  if (anyNA(s)) stop("`s` must not contain NA.")
  if (any(s <= 0)) stop("All entries of `s` must be > 0.")
  as.numeric(s)
}

.eb_smoother_data_frame <- function(x, s) {
  data.frame(x = as.numeric(x), s = as.numeric(s))
}

.point_family_profile_s_range <- function(x,
                                          profile_s_lower = 1e-8,
                                          profile_s_upper = NULL) {
  if (!is.numeric(profile_s_lower) || length(profile_s_lower) != 1L ||
      is.na(profile_s_lower) || !is.finite(profile_s_lower)) {
    stop("`profile_s_lower` must be a single finite numeric.")
  }
  profile_s_lower <- max(as.numeric(profile_s_lower), sqrt(.Machine$double.xmin))

  if (is.null(profile_s_upper)) {
    finite_x <- x[is.finite(x)]
    if (!length(finite_x)) stop("`x` must contain at least one finite value.")
    data_scale <- max(stats::sd(finite_x), max(abs(finite_x)), 1, na.rm = TRUE)
    profile_s_upper <- 10 * data_scale
  }
  if (!is.numeric(profile_s_upper) || length(profile_s_upper) != 1L ||
      is.na(profile_s_upper) || !is.finite(profile_s_upper)) {
    stop("`profile_s_upper` must be `NULL` or a single finite numeric.")
  }
  profile_s_upper <- as.numeric(profile_s_upper)
  if (profile_s_upper <= profile_s_lower) {
    stop("`profile_s_upper` must be larger than `profile_s_lower`.")
  }

  c(lower = profile_s_lower, upper = profile_s_upper)
}

.point_family_point_mass <- function(fitted_g) {
  if (is.null(fitted_g$pi)) return(NA_real_)
  spread <- if (!is.null(fitted_g$scale)) fitted_g$scale else fitted_g$sd
  if (is.null(spread)) return(NA_real_)
  zero_index <- which.min(abs(as.numeric(spread)))
  as.numeric(fitted_g$pi[zero_index])
}

.point_family_nonzero_scale <- function(fitted_g) {
  spread <- if (!is.null(fitted_g$scale)) fitted_g$scale else fitted_g$sd
  if (is.null(spread)) return(NA_real_)
  nonzero_scale <- as.numeric(spread)[as.numeric(spread) > 0]
  if (!length(nonzero_scale)) return(NA_real_)
  nonzero_scale[1]
}

.point_family_names <- c("point_exponential", "point_normal", "point_laplace")

.point_family_function <- function(family) {
  switch(
    family,
    point_exponential = ebnm::ebnm_point_exponential,
    point_normal = ebnm::ebnm_point_normal,
    point_laplace = ebnm::ebnm_point_laplace,
    stop("Unsupported point family: ", family, call. = FALSE)
  )
}

.fit_point_family_smoother <- function(x,
                                       s = NULL,
                                       family,
                                       g_init = NULL,
                                       fix_g = FALSE,
                                       mode = 0,
                                       scale = "estimate",
                                       output = ebnm::ebnm_output_default(),
                                       optmethod = NULL,
                                       control = NULL,
                                       profile_s_lower = 1e-8,
                                       profile_s_upper = NULL,
                                       profile_s_tol = .Machine$double.eps^0.25) {
  x <- as.numeric(x)
  if (anyNA(x)) stop("`x` must not contain NA.")
  if (!is.numeric(profile_s_tol) || length(profile_s_tol) != 1L ||
      is.na(profile_s_tol) || !is.finite(profile_s_tol) || profile_s_tol <= 0) {
    stop("`profile_s_tol` must be a single positive finite numeric.")
  }
  if (!is.logical(fix_g) || length(fix_g) != 1L || is.na(fix_g)) {
    stop("`fix_g` must be TRUE or FALSE.")
  }

  ebnm_fun <- .point_family_function(family)

  run_ebnm <- function(s_value) {
    s_use <- .normalize_eb_smoother_s(s_value, length(x))
    ebnm_fun(
      x = x,
      s = s_use,
      mode = mode,
      scale = scale,
      g_init = g_init,
      fix_g = fix_g,
      output = output,
      optmethod = optmethod,
      control = control
    )
  }

  fixed_s <- !is.null(s)
  profile_range <- NULL
  profile_optimization <- NULL

  if (fixed_s) {
    s_vec <- .normalize_eb_smoother_s(s, length(x))
    fit <- run_ebnm(s_vec)
    fitted_noise_sd <- if (length(s) == 1L) as.numeric(s) else NA_real_
    log_likelihood_semantics <- paste0(family, "_fixed_noise")
  } else {
    profile_range <- .point_family_profile_s_range(
      x = x,
      profile_s_lower = profile_s_lower,
      profile_s_upper = profile_s_upper
    )
    log_s_range <- log(profile_range)

    objective <- function(log_s) {
      fit_at_s <- tryCatch(run_ebnm(exp(log_s)), error = function(e) NULL)
      if (is.null(fit_at_s)) return(Inf)
      log_likelihood <- as.numeric(fit_at_s$log_likelihood)
      if (!is.finite(log_likelihood)) return(Inf)
      -log_likelihood
    }

    opt <- stats::optimize(f = objective, interval = log_s_range, tol = profile_s_tol)
    candidate_log_s <- unique(c(log_s_range[1], opt$minimum, log_s_range[2]))
    candidate_fits <- lapply(candidate_log_s, function(log_s) {
      fit_at_s <- tryCatch(run_ebnm(exp(log_s)), error = function(e) NULL)
      if (is.null(fit_at_s)) return(NULL)
      log_likelihood <- as.numeric(fit_at_s$log_likelihood)
      if (!is.finite(log_likelihood)) return(NULL)
      list(fit = fit_at_s, s = exp(log_s), log_likelihood = log_likelihood)
    })
    valid_fits <- candidate_fits[!vapply(candidate_fits, is.null, logical(1))]
    if (!length(valid_fits)) {
      stop("Could not fit the ", family, " model at any profiled `s` value.")
    }

    best <- valid_fits[[which.max(vapply(valid_fits, `[[`, numeric(1), "log_likelihood"))]]
    fit <- best$fit
    fitted_noise_sd <- best$s
    s_vec <- rep(fitted_noise_sd, length(x))
    log_likelihood_semantics <- paste0("profile_", family)
    profile_optimization <- list(
      lower = unname(profile_range[["lower"]]),
      upper = unname(profile_range[["upper"]]),
      optimizer_s = exp(opt$minimum),
      optimizer_objective = opt$objective,
      selected_s = fitted_noise_sd,
      selected_log_likelihood = best$log_likelihood
    )
  }

  posterior <- fit$posterior
  if (is.null(posterior$var)) {
    posterior$var <- as.numeric(posterior$sd)^2
  }
  if (is.null(posterior$second_moment)) {
    posterior$second_moment <- as.numeric(posterior$mean)^2 + as.numeric(posterior$var)
  }

  out <- list(
    posterior = posterior,
    fitted_g = fit$fitted_g,
    fitted_beta = numeric(0),
    beta_mode = NULL,
    fitted_noise_sd = fitted_noise_sd,
    log_likelihood = fit$log_likelihood,
    log_likelihood_semantics = log_likelihood_semantics,
    posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
    data = .eb_smoother_data_frame(x, s_vec),
    prior_family = family,
    backend = "exact",
    profile_s_range = profile_range,
    profile_optimization = profile_optimization,
    point_mass_probability = .point_family_point_mass(fit$fitted_g),
    point_nonzero_scale = .point_family_nonzero_scale(fit$fitted_g)
  )
  if (identical(family, "point_exponential")) {
    out$point_exponential_pi0 <- out$point_mass_probability
    out$point_exponential_scale <- out$point_nonzero_scale
  }
  out
}

.check_optional_beta_prec <- function(beta_prec, nm = "beta_prec") {
  if (is.null(beta_prec)) {
    return(NULL)
  }

  beta_prec <- .check_single_numeric(beta_prec, nm)
  if (beta_prec < 0) {
    stop(nm, " must be NULL or a single non-negative numeric.")
  }
  as.numeric(beta_prec)
}

.check_optional_beta_vector <- function(beta,
                                        nm,
                                        expected_length = NULL,
                                        allow_null = TRUE) {
  if (is.null(beta)) {
    if (allow_null) return(NULL)
    stop(nm, " must not be NULL.")
  }
  if (!is.numeric(beta) || anyNA(beta)) {
    stop(nm, " must be a numeric vector with no NA.")
  }
  beta <- as.numeric(beta)
  if (!is.null(expected_length) && length(beta) != expected_length) {
    stop(nm, " must have length ", expected_length, ".")
  }
  beta
}

.normalize_fix_params <- function(fix_params,
                                  allowed,
                                  fix_g = FALSE,
                                  fix_g_params = character(),
                                  nm = "fix_params") {
  if (is.null(fix_params)) {
    fix_params <- character()
  }
  if (!is.character(fix_params) || anyNA(fix_params)) {
    stop("`", nm, "` must be a character vector.")
  }
  invalid <- setdiff(fix_params, allowed)
  if (length(invalid) > 0L) {
    stop(
      "`", nm, "` contains unsupported value(s): ",
      paste(sprintf('"%s"', invalid), collapse = ", "),
      ". Allowed values are ",
      paste(sprintf('"%s"', allowed), collapse = ", "),
      "."
    )
  }
  if (isTRUE(fix_g)) {
    fix_params <- c(fix_params, fix_g_params)
  }
  unique(fix_params)
}

.resolve_fixed_beta_from_fix_params <- function(fix_params,
                                                beta_fixed = NULL,
                                                beta_prec = NULL,
                                                beta_prec_supplied = FALSE,
                                                g_init = NULL,
                                                expected_length = NULL,
                                                beta_nm = "beta_fixed") {
  if (!("beta" %in% fix_params)) {
    return(beta_fixed)
  }
  if (isTRUE(beta_prec_supplied) && !is.null(beta_prec)) {
    stop("`fix_params = \"beta\"` cannot be combined with `beta_prec`.")
  }
  if (!is.null(beta_fixed)) {
    return(.check_optional_beta_vector(
      beta_fixed,
      beta_nm,
      expected_length = expected_length,
      allow_null = FALSE
    ))
  }
  if (!is.null(g_init) && !is.null(g_init$beta)) {
    return(.check_optional_beta_vector(
      g_init$beta,
      "g_init$beta",
      expected_length = expected_length,
      allow_null = FALSE
    ))
  }
  stop("`fix_params = \"beta\"` requires `beta_fixed` or `g_init$beta`.")
}

.beta_mode_from_prec <- function(beta_prec) {
  if (is.null(beta_prec)) {
    return("empirical_bayes")
  }
  if (beta_prec == 0) {
    return("prior_flat")
  }
  "prior_proper"
}

.eb_smoother_resolve_beta_spec <- function(beta_fixed = NULL,
                                           beta_prec = NULL,
                                           g_init_beta_prec = NULL,
                                           legacy_beta_prec = NULL,
                                           nm = "beta_prec") {
  beta_prec_explicit <- .check_optional_beta_prec(beta_prec, nm = nm)
  beta_prec_from_g <- .check_optional_beta_prec(g_init_beta_prec, nm = "g_init$beta_prec")
  beta_prec_legacy <- if (is.null(legacy_beta_prec)) {
    NULL
  } else {
    .check_single_numeric(legacy_beta_prec, "legacy_beta_prec")
  }

  if (!is.null(beta_fixed) && !is.null(beta_prec_explicit)) {
    stop("`beta_fixed` and `beta_prec` cannot be supplied together.")
  }

  if (!is.null(beta_fixed)) {
    return(list(
      mode = "fixed",
      beta_fixed = beta_fixed,
      beta_prec = NULL,
      beta_prec_source = "fixed"
    ))
  }

  if (!is.null(beta_prec_explicit)) {
    return(list(
      mode = .beta_mode_from_prec(beta_prec_explicit),
      beta_fixed = NULL,
      beta_prec = beta_prec_explicit,
      beta_prec_source = "argument"
    ))
  }

  if (!is.null(beta_prec_from_g)) {
    return(list(
      mode = .beta_mode_from_prec(beta_prec_from_g),
      beta_fixed = NULL,
      beta_prec = beta_prec_from_g,
      beta_prec_source = "g_init"
    ))
  }

  if (!is.null(beta_prec_legacy)) {
    if (beta_prec_legacy < 0) {
      return(list(
        mode = "empirical_bayes",
        beta_fixed = NULL,
        beta_prec = NULL,
        beta_prec_source = "legacy_setup"
      ))
    }
    return(list(
      mode = .beta_mode_from_prec(beta_prec_legacy),
      beta_fixed = NULL,
      beta_prec = beta_prec_legacy,
      beta_prec_source = "legacy_setup"
    ))
  }

  list(
    mode = "empirical_bayes",
    beta_fixed = NULL,
    beta_prec = NULL,
    beta_prec_source = "default"
  )
}

.eb_smoother_public_beta_mode <- function(family, beta_mode) {
  if (is.null(beta_mode)) {
    return(NULL)
  }

  beta_mode <- as.character(beta_mode)
  beta_mode_map <- c(
    profile = "empirical_bayes",
    optimized = "empirical_bayes",
    integrated_flat = "prior_flat",
    integrated = "prior_flat",
    prior_proper = "prior_proper",
    fixed = "fixed",
    empirical_bayes = "empirical_bayes",
    prior_flat = "prior_flat"
  )

  if (!is.na(beta_mode_map[beta_mode])) {
    return(unname(beta_mode_map[beta_mode]))
  }

  beta_mode
}

.as_eb_smoother_fit <- function(fit,
                                family,
                                noise_mode,
                                fitted_noise_sd,
                                backend = NULL) {
  raw_fit <- fit
  out <- fit
  out$fitted_noise_sd <- if (length(fitted_noise_sd) == 0L) NA_real_ else as.numeric(fitted_noise_sd)
  out$family <- family
  out$backend <- if (is.null(backend)) fit$backend else backend
  out$noise_mode <- noise_mode
  out$beta_mode_internal <- fit$beta_mode
  out$beta_mode <- .eb_smoother_public_beta_mode(family, fit$beta_mode)
  out$raw_fit <- raw_fit
  class(out) <- c("eb_smoother_fit", "list")
  out
}

.eb_smoother_num_obs <- function(x) {
  if (!is.null(x$data) && is.data.frame(x$data)) {
    return(nrow(x$data))
  }
  if (!is.null(x$posterior) && is.data.frame(x$posterior)) {
    return(nrow(x$posterior))
  }
  NA_integer_
}

.eb_smoother_beta_mode <- function(x) {
  if (is.null(x$beta_mode)) return(NULL)
  as.character(x$beta_mode)
}

.format_eb_smoother_number <- function(x, digits = 4) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x)) return("NA")
  format(signif(x, digits = digits), trim = TRUE, scientific = FALSE)
}

.format_eb_smoother_vector <- function(x, digits = 4, max_entries = 4L) {
  x <- as.numeric(x)
  if (!length(x)) return("numeric(0)")
  shown <- utils::head(x, max_entries)
  pieces <- vapply(shown, .format_eb_smoother_number, character(1), digits = digits)
  out <- paste(pieces, collapse = ", ")
  if (length(x) > max_entries) out <- paste0(out, ", ...")
  if (length(x) == 1L) out else paste0("c(", out, ")")
}

.eb_smoother_prior_parameters <- function(x) {
  if (is.null(x$fitted_g)) return(NULL)

  if (identical(x$family, "matern")) {
    theta <- x$fitted_g$theta
    sigma <- x$fitted_g$sigma
    out <- list(
      range = if (length(theta) == 1L && is.finite(theta)) exp(theta) else NA_real_,
      sigma = if (length(sigma) == 1L) as.numeric(sigma) else NA_real_
    )
    if (!is.null(x$fitted_g$beta_prec)) {
      out$beta_prec <- as.numeric(x$fitted_g$beta_prec)
    }
    return(out)
  }

  if (identical(x$family, "lgp")) {
    scale <- x$fitted_g$scale
    out <- list(
      scale = if (length(scale) == 1L) as.numeric(scale) else NA_real_
    )
    if (!is.null(x$fitted_g$beta_prec)) {
      out$beta_prec <- as.numeric(x$fitted_g$beta_prec)
    }
    return(out)
  }

  if (identical(x$family, "constant")) {
    out <- list()
    if (!is.null(x$fitted_g$beta_prec)) {
      out$beta_prec <- as.numeric(x$fitted_g$beta_prec)
    }
    return(out)
  }

  if (x$family %in% .point_family_names) {
    return(list(
      point_mass_probability = .point_family_point_mass(x$fitted_g),
      nonzero_scale = .point_family_nonzero_scale(x$fitted_g)
    ))
  }

  NULL
}

.eb_smoother_posterior_summary <- function(x) {
  if (is.null(x$posterior) || !is.data.frame(x$posterior) || nrow(x$posterior) == 0L) {
    return(NULL)
  }

  post_mean <- as.numeric(x$posterior$mean)
  post_sd <- sqrt(pmax(as.numeric(x$posterior$var), 0))
  probs <- c(0, 0.25, 0.5, 0.75, 1)
  nm <- c("min", "q25", "median", "q75", "max")

  list(
    mean_quantiles = stats::setNames(
      as.numeric(stats::quantile(post_mean, probs = probs, names = FALSE)),
      nm
    ),
    sd_quantiles = stats::setNames(
      as.numeric(stats::quantile(post_sd, probs = probs, names = FALSE)),
      nm
    )
  )
}

.eb_smoother_diagnostics <- function(x) {
  diagnostic_fields <- c(
    "log_likelihood_stepA_penalized",
    "log_likelihood_stepA_mlik_integration",
    "log_likelihood_inlabru_mlik_integration",
    "log_likelihood_laplace_at_inlabru_params",
    "log_likelihood_fisher_pql_stepA",
    "log_likelihood_original_at_fisher_pql_mode",
    "log_likelihood_stepA_mlik_gaussian",
    "log_likelihood_stepA_joint_log_posterior",
    "log_likelihood_exact_at_stepA_mode",
    "log_likelihood_stepA",
    "log_likelihood_stepB_joint",
    "log_likelihood_stepB_laplace",
    "fisher_pql_diagnostics",
    "inla_validation"
  )

  out <- list()
  for (nm in diagnostic_fields) {
    value <- x[[nm]]
    if (!is.null(value)) {
      out[[nm]] <- if (inherits(value, "logLik")) as.numeric(value) else value
    }
  }
  out
}

.print_named_numeric_lines <- function(x, indent = "  ", digits = 4) {
  if (!length(x)) return(invisible(NULL))
  for (nm in names(x)) {
    cat(indent, nm, ": ", .format_eb_smoother_number(x[[nm]], digits = digits), "\n", sep = "")
  }
  invisible(NULL)
}

#' Print and summarize `eb_smoother_fit` objects
#'
#' @description
#' `print()` shows a compact one-screen overview of an `eb_smoother_fit`.
#' `summary()` reorganizes the stored fit contents into a richer summary object
#' without recomputing expensive quantities.
#'
#' @param x For `print()`, an `eb_smoother_fit` object or a
#'   `summary.eb_smoother_fit` object.
#' @param object An `eb_smoother_fit` object returned by [eb_smoother()].
#' @param ... Currently unused.
#'
#' @return
#' `print.eb_smoother_fit()` and `print.summary.eb_smoother_fit()` return their
#' input invisibly.
#'
#' `summary.eb_smoother_fit()` returns an object of class
#' `c("summary.eb_smoother_fit", "list")` containing the fitted model metadata,
#' compact posterior summaries, and any stored backend diagnostics.
#'
#' @examples
#' loc <- seq(0, 1, length.out = 30)
#' x <- sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)
#'
#' fit <- eb_smoother(x, s = 0.1, family = "matern", locations = loc)
#' fit
#' summary(fit)
#'
#' @name eb_smoother_fit_methods
NULL

#' @rdname eb_smoother_fit_methods
#' @export
print.eb_smoother_fit <- function(x, ...) {
  n_obs <- .eb_smoother_num_obs(x)
  prior <- .eb_smoother_prior_parameters(x)
  beta_mode <- .eb_smoother_beta_mode(x)

  cat("Empirical-Bayes smoother fit\n")
  cat("  class: eb_smoother_fit\n")
  cat("  family: ", x$family, "\n", sep = "")
  cat("  backend: ", x$backend, "\n", sep = "")
  cat("  noise mode: ", x$noise_mode, "\n", sep = "")
  if (!is.null(beta_mode)) {
    cat("  beta mode: ", beta_mode, "\n", sep = "")
  }
  if (is.finite(n_obs)) {
    cat("  observations: ", n_obs, "\n", sep = "")
  }
  if (!is.null(prior) && length(prior)) {
    prior_text <- paste(
      sprintf(
        "%s = %s",
        names(prior),
        vapply(prior, .format_eb_smoother_number, character(1))
      ),
      collapse = ", "
    )
    cat("  fitted prior: ", prior_text, "\n", sep = "")
  }
  cat("  fitted beta: ", .format_eb_smoother_vector(x$fitted_beta), "\n", sep = "")
  if (is.finite(x$fitted_noise_sd)) {
    cat("  fitted noise SD: ", .format_eb_smoother_number(x$fitted_noise_sd), "\n", sep = "")
  }
  cat("  log-likelihood: ", .format_eb_smoother_number(x$log_likelihood), "\n", sep = "")

  invisible(x)
}

#' @rdname eb_smoother_fit_methods
#' @export
summary.eb_smoother_fit <- function(object, ...) {
  out <- list(
    class = "eb_smoother_fit",
    family = object$family,
    backend = object$backend,
    noise_mode = object$noise_mode,
    beta_mode = object$beta_mode,
    prior_family = object$prior_family,
    observations = .eb_smoother_num_obs(object),
    fitted_prior = .eb_smoother_prior_parameters(object),
    fitted_beta = object$fitted_beta,
    fitted_noise_sd = object$fitted_noise_sd,
    log_likelihood = if (inherits(object$log_likelihood, "logLik")) {
      as.numeric(object$log_likelihood)
    } else {
      object$log_likelihood
    },
    posterior = .eb_smoother_posterior_summary(object),
    diagnostics = .eb_smoother_diagnostics(object)
  )
  class(out) <- c("summary.eb_smoother_fit", "list")
  out
}

#' @rdname eb_smoother_fit_methods
#' @export
print.summary.eb_smoother_fit <- function(x, ...) {
  cat("Empirical-Bayes smoother summary\n")
  cat("  class: ", x$class, "\n", sep = "")
  cat("  family: ", x$family, "\n", sep = "")
  cat("  backend: ", x$backend, "\n", sep = "")
  cat("  noise mode: ", x$noise_mode, "\n", sep = "")
  if (!is.null(x$beta_mode)) {
    cat("  beta mode: ", x$beta_mode, "\n", sep = "")
  }
  if (!is.null(x$prior_family)) {
    cat("  prior family: ", x$prior_family, "\n", sep = "")
  }
  if (is.finite(x$observations)) {
    cat("  observations: ", x$observations, "\n", sep = "")
  }

  if (!is.null(x$fitted_prior) && length(x$fitted_prior)) {
    cat("\nFitted prior\n")
    .print_named_numeric_lines(x$fitted_prior)
  }

  cat("\nFitted coefficients\n")
  cat("  beta: ", .format_eb_smoother_vector(x$fitted_beta), "\n", sep = "")
  if (is.finite(x$fitted_noise_sd)) {
    cat("  noise SD: ", .format_eb_smoother_number(x$fitted_noise_sd), "\n", sep = "")
  }
  cat("  log-likelihood: ", .format_eb_smoother_number(x$log_likelihood), "\n", sep = "")

  if (!is.null(x$posterior) && length(x$posterior)) {
    cat("\nPosterior mean quantiles\n")
    .print_named_numeric_lines(x$posterior$mean_quantiles)
    cat("\nPosterior SD quantiles\n")
    .print_named_numeric_lines(x$posterior$sd_quantiles)
  }

  if (length(x$diagnostics)) {
    cat("\nDiagnostics\n")
    for (nm in names(x$diagnostics)) {
      value <- x$diagnostics[[nm]]
      text <- if (length(value) == 1L) {
        .format_eb_smoother_number(value)
      } else {
        .format_eb_smoother_vector(value)
      }
      cat("  ", nm, ": ", text, "\n", sep = "")
    }
  }

  invisible(x)
}

.build_lgp_A_in_par_order <- function(B, X, par_names) {
  B <- Matrix::Matrix(B, sparse = TRUE)
  X <- Matrix::Matrix(X, sparse = TRUE)

  u_pos <- which(par_names == "U")
  b_pos <- which(par_names == "beta")

  if (length(u_pos) == ncol(B) &&
      length(b_pos) == ncol(X) &&
      (length(b_pos) == 0L || max(u_pos) < min(b_pos))) {
    return(cbind(B, X))
  }

  cols <- vector("list", length(par_names))
  u_counter <- 0L
  b_counter <- 0L
  for (k in seq_along(par_names)) {
    if (par_names[k] == "U") {
      u_counter <- u_counter + 1L
      cols[[k]] <- B[, u_counter, drop = FALSE]
    } else if (par_names[k] == "beta") {
      b_counter <- b_counter + 1L
      cols[[k]] <- X[, b_counter, drop = FALSE]
    } else {
      stop("Unexpected free parameter name: ", par_names[k])
    }
  }
  do.call(cbind, cols)
}

.fit_lgp_unknown_noise <- function(x,
                                   LGP_setup,
                                   g_init = NULL,
                                   beta_fixed = NULL,
                                   beta_prec = NULL,
                                   beta_prec_missing = FALSE,
                                   link = c("identity", "log", "softplus"),
                                   dll = "EBSmoothr") {
  link <- match.arg(link)

  if (is.null(LGP_setup$X) || is.null(LGP_setup$B) || is.null(LGP_setup$P)) {
    stop("`setup` must contain X, B, and P.")
  }

  x <- as.numeric(x)
  if (anyNA(x)) stop("`x` must not contain NA.")

  n <- length(x)
  if (nrow(LGP_setup$X) != n || nrow(LGP_setup$B) != n) {
    stop("length(x) must match nrow(X) and nrow(B) in `setup`.")
  }

  link_id_arg <- if (identical(link, "identity")) 0L else if (identical(link, "log")) 1L else 2L
  tmbdat <- LGP_setup
  tmbdat$x <- x
  tmbdat$s <- rep(1, n)
  tmbdat$link_id <- as.integer(link_id_arg)
  tmbdat$learn_noise <- 1L
  tmbdat$model_id <- 0L

  pB <- ncol(tmbdat$B)
  pX <- ncol(tmbdat$X)

  beta_spec <- .eb_smoother_resolve_beta_spec(
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    g_init_beta_prec = if (is.null(g_init)) NULL else g_init$beta_prec,
    legacy_beta_prec = if (isTRUE(beta_prec_missing)) tmbdat$betaprec else NULL
  )
  beta_mode <- beta_spec$mode
  beta_fixed_use <- if (is.null(beta_spec$beta_fixed)) {
    NULL
  } else {
    .check_optional_beta_vector(
      beta_spec$beta_fixed,
      "beta_fixed",
      expected_length = pX,
      allow_null = FALSE
    )
  }
  beta_prec_use <- beta_spec$beta_prec
  betaprec <- if (beta_mode == "empirical_bayes") -1 else beta_prec_use
  tmbdat$betaprec <- if (is.null(betaprec)) -1 else betaprec

  if (is.null(g_init)) g_init <- LGP(0, beta = rep(0, pX), beta_prec = beta_prec_use)
  theta0 <- .check_single_numeric(g_init$scale, "g_init$scale")
  beta_init <- if (!is.null(beta_fixed_use)) {
    beta_fixed_use
  } else if (!is.null(g_init$beta)) {
    .check_optional_beta_vector(g_init$beta, "g_init$beta", expected_length = pX, allow_null = FALSE)
  } else {
    rep(0, pX)
  }
  noise_sd0 <- stats::sd(x)
  if (!is.finite(noise_sd0) || noise_sd0 <= 0) noise_sd0 <- 1

  par0 <- list(
    theta = theta0,
    U = rep(0, pB),
    beta = as.numeric(beta_init),
    log_noise = log(noise_sd0)
  )

  fitted_theta <- theta0
  fitted_beta <- as.numeric(beta_init)
  fitted_log_noise <- log(noise_sd0)
  ll_stepA <- NA_real_

  if (beta_mode == "fixed") {
    objA <- TMB::MakeADFun(
      data = tmbdat,
      parameters = within(par0, beta <- as.numeric(beta_fixed_use)),
      map = list(beta = factor(rep(NA, pX))),
      DLL = dll,
      random = "U",
      silent = TRUE
    )
    optA <- optim(par = objA$par, fn = objA$fn, gr = objA$gr, method = "BFGS")

    ll_stepA <- -as.numeric(optA$value)
    fitted_theta <- as.numeric(optA$par[["theta"]])
    fitted_log_noise <- as.numeric(optA$par[["log_noise"]])
    fitted_beta <- as.numeric(beta_fixed_use)
  } else if (betaprec < 0) {
    objA <- TMB::MakeADFun(
      data = tmbdat,
      parameters = par0,
      DLL = dll,
      random = "U",
      silent = TRUE
    )
    optA <- optim(par = objA$par, fn = objA$fn, gr = objA$gr, method = "BFGS")

    ll_stepA <- -as.numeric(optA$value)
    fitted_theta <- as.numeric(optA$par[["theta"]])
    fitted_log_noise <- as.numeric(optA$par[["log_noise"]])

    beta_idx <- which(names(optA$par) == "beta")
    if (length(beta_idx) != pX) {
      stop(sprintf("Step A: expected %d beta entries, got %d.", pX, length(beta_idx)))
    }
    fitted_beta <- as.numeric(optA$par[beta_idx])
  } else {
    objA <- TMB::MakeADFun(
      data = tmbdat,
      parameters = par0,
      DLL = dll,
      random = c("U", "beta"),
      silent = TRUE
    )
    optA <- optim(par = objA$par, fn = objA$fn, gr = objA$gr, method = "BFGS")

    ll_stepA <- -as.numeric(optA$value)
    fitted_theta <- as.numeric(optA$par[["theta"]])
    fitted_log_noise <- as.numeric(optA$par[["log_noise"]])
  }

  fitted_noise_sd <- exp(fitted_log_noise)
  fixed_theta_and_beta <- identical(beta_mode, "fixed") || (betaprec < 0)

  if (fixed_theta_and_beta) {
    mapB <- list(
      theta = factor(NA),
      beta = factor(rep(NA, pX)),
      log_noise = factor(NA)
    )
    parB <- list(
      theta = as.numeric(fitted_theta),
      U = rep(0, pB),
      beta = as.numeric(fitted_beta),
      log_noise = as.numeric(fitted_log_noise)
    )

    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = parB,
      map = mapB,
      DLL = dll,
      silent = TRUE
    )

    optB <- nlminb(
      start = ff$par,
      objective = ff$fn,
      gradient = ff$gr,
      control = list(eval.max = 20000, iter.max = 20000)
    )

      H <- numDeriv::hessian(function(w) ff$fn(w), optB$par)
      prec <- Matrix::forceSymmetric(H)

      U_hat <- as.numeric(optB$par)
      beta_hat <- as.numeric(if (beta_mode == "fixed") beta_fixed_use else fitted_beta)
      eta_mean <- as.numeric(tmbdat$B %*% U_hat + tmbdat$X %*% beta_hat)
      eta_var <- .compute_diag_A_Qinv_At(tmbdat$B, prec)
      free_dim <- pB
  } else {
    mapB <- list(
      theta = factor(NA),
      log_noise = factor(NA)
    )
    parB <- list(
      theta = as.numeric(fitted_theta),
      U = rep(0, pB),
      beta = rep(0, pX),
      log_noise = as.numeric(fitted_log_noise)
    )

    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = parB,
      map = mapB,
      DLL = dll,
      silent = TRUE
    )

    optB <- nlminb(
      start = ff$par,
      objective = ff$fn,
      gradient = ff$gr,
      control = list(eval.max = 20000, iter.max = 20000)
    )

    H <- numDeriv::hessian(function(w) ff$fn(w), optB$par)
    prec <- Matrix::forceSymmetric(H)

    par_names <- names(optB$par)
    u_idx <- which(par_names == "U")
    b_idx <- which(par_names == "beta")
    if (length(u_idx) != pB || length(b_idx) != pX) {
      stop("Step B: unexpected names or lengths in optB$par.")
    }

    U_hat <- as.numeric(optB$par[u_idx])
    beta_hat <- as.numeric(optB$par[b_idx])
    eta_mean <- as.numeric(tmbdat$B %*% U_hat + tmbdat$X %*% beta_hat)
    A <- .build_lgp_A_in_par_order(tmbdat$B, tmbdat$X, par_names)
    eta_var <- .compute_diag_A_Qinv_At(A, prec)
    free_dim <- pB + pX
  }

  response_moments <- .lgp_response_moments_from_eta(eta_mean, eta_var, link = link)
  post_mean <- response_moments$mean
  post_var <- response_moments$var

  posterior <- data.frame(mean = as.numeric(post_mean), var = as.numeric(post_var))
  posterior$second_moment <- posterior$mean^2 + posterior$var

  ll_stepB_joint <- -ff$fn(optB$par)
  ll_stepB_laplace <- ll_stepB_joint - 0.5 * .compute_logdet_spd(prec) +
    0.5 * free_dim * log(2 * pi)

  class(ll_stepB_joint) <- "logLik"
  class(ll_stepB_laplace) <- "logLik"

  log_likelihood <- structure(ll_stepA, class = "logLik")

  posterior_sampler <- function(nsamp) {
    nsamp <- .check_single_numeric(nsamp, "nsamp")
    if (nsamp < 1 || nsamp != floor(nsamp)) {
      stop("`nsamp` must be a positive integer.")
    }

    samps <- LaplacesDemon::rmvnp(n = nsamp, mu = as.numeric(optB$par), Omega = as.matrix(prec))
    if (fixed_theta_and_beta) {
      eta_s <- as.matrix(tmbdat$B) %*% t(samps) + as.matrix(tmbdat$X) %*% beta_hat
    } else {
      par_names <- names(optB$par)
      A <- .build_lgp_A_in_par_order(tmbdat$B, tmbdat$X, par_names)
      eta_s <- as.matrix(A) %*% t(samps)
    }

    if (tmbdat$link_id == 0L) t(eta_s) else if (tmbdat$link_id == 1L) t(exp(eta_s)) else t(.softplus_stable(eta_s))
  }

  list(
    posterior = posterior,
    fitted_g = LGP(
      scale = fitted_theta,
      beta = beta_hat,
      beta_prec = if (beta_mode == "prior_flat") 0 else beta_prec_use
    ),
    fitted_beta = beta_hat,
    beta_mode = beta_mode,
    beta_prec = if (beta_mode == "prior_flat") 0 else beta_prec_use,
    fitted_noise_sd = as.numeric(fitted_noise_sd),
    log_likelihood = log_likelihood,
    log_likelihood_semantics = paste0("laplace_", beta_mode),
    log_likelihood_stepA = structure(ll_stepA, class = "logLik"),
    log_likelihood_stepB_joint = ll_stepB_joint,
    log_likelihood_stepB_laplace = ll_stepB_laplace,
    posterior_sampler = posterior_sampler,
    data = .eb_smoother_data_frame(x, rep(fitted_noise_sd, n)),
    prior_family = "LGP_learned_noise",
    backend = "tmb",
    laplace_implementation = "tmb",
    laplace_curvature = "observed",
    link = link,
    g_init = LGP(theta0, beta = beta_init, beta_prec = if (beta_mode == "prior_flat") 0 else beta_prec_use)
  )
}

.constant_posterior_sampler <- function(beta0, n_obs) {
  beta0 <- .check_single_numeric(beta0, "beta0")
  n_obs <- as.integer(.check_single_numeric(n_obs, "n_obs"))

  function(nsamp) {
    nsamp <- .check_single_numeric(nsamp, "nsamp")
    if (nsamp < 1 || nsamp != floor(nsamp)) {
      stop("`nsamp` must be a positive integer.")
    }
    matrix(beta0, nrow = nsamp, ncol = n_obs)
  }
}

.constant_posterior_summary <- function(beta0, n_obs) {
  beta0 <- .check_single_numeric(beta0, "beta0")
  data.frame(
    mean = rep(beta0, n_obs),
    var = rep(0, n_obs),
    second_moment = rep(beta0^2, n_obs)
  )
}

#' Define the constant Gaussian family
#'
#' @description
#' Creates an object representing the constant Gaussian baseline used by
#' [eb_smoother()] when `family = "constant"`.
#'
#' This family state is mainly useful as the constant-reference option inside
#' the high-level `eb_smoother()` and `spatial_scores()` workflows.
#'
#' @param beta Optional scalar intercept state.
#' @param beta_prec Optional non-negative scalar prior precision on `beta`.
#'
#' @return An object of class `"Constant"`.
#'
#' @export
Constant <- function(beta = NULL, beta_prec = NULL) {
  beta <- .check_optional_beta_vector(beta, "beta", expected_length = 1L)
  beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")

  structure(
    list(
      beta = beta,
      beta_prec = beta_prec
    ),
    class = "Constant"
  )
}

.constant_beta_state <- function(beta_mode, beta_hat, beta_prec = NULL) {
  Constant(
    beta = beta_hat,
    beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec
  )
}

.constant_response_moments <- function(beta_mean, beta_var, n, link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  beta_mean <- as.numeric(beta_mean)[1]
  beta_var <- pmax(as.numeric(beta_var)[1], 0)

  if (identical(link, "identity")) {
    mean <- beta_mean
    var <- beta_var
  } else if (identical(link, "log")) {
    mean <- exp(beta_mean + 0.5 * beta_var)
    var <- exp(2 * beta_mean + beta_var) * (exp(beta_var) - 1)
  } else {
    moments <- .softplus_gaussian_moments(beta_mean, beta_var)
    mean <- moments$mean
    var <- moments$var
  }

  data.frame(
    mean = rep(as.numeric(mean), n),
    var = rep(as.numeric(var), n),
    second_moment = rep(as.numeric(mean^2 + var), n)
  )
}

.constant_response_sampler <- function(beta_hat, beta_var, n, link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  beta_hat <- as.numeric(beta_hat)[1]
  beta_var <- pmax(as.numeric(beta_var)[1], 0)
  function(nsamp) {
    nsamp <- .check_single_numeric(nsamp, "nsamp")
    if (nsamp < 1 || nsamp != floor(nsamp)) {
      stop("`nsamp` must be a positive integer.")
    }
    beta_samp <- if (beta_var > 0) {
      stats::rnorm(nsamp, mean = beta_hat, sd = sqrt(beta_var))
    } else {
      rep(beta_hat, nsamp)
    }
    value <- if (identical(link, "identity")) {
      beta_samp
    } else if (identical(link, "log")) {
      exp(beta_samp)
    } else {
      .softplus_stable(beta_samp)
    }
    matrix(value, nrow = nsamp, ncol = n)
  }
}

.constant_log_positive_mean_beta <- function(x, w = NULL, nm = "weighted mean") {
  x <- as.numeric(x)
  if (is.null(w)) {
    mu_hat <- mean(x)
  } else {
    mu_hat <- sum(as.numeric(w) * x) / sum(as.numeric(w))
  }
  if (!is.finite(mu_hat) || mu_hat <= 0) {
    stop("The log-link constant ", nm, " is not positive, so beta has no finite optimum.")
  }
  log(mu_hat)
}

.constant_log_beta_interval <- function(beta_start) {
  beta_start <- as.numeric(beta_start)[1]
  if (!is.finite(beta_start)) beta_start <- 0
  beta_start + c(-30, 30)
}

.constant_log_observation_loglik <- function(beta, x, s) {
  mu <- exp(as.numeric(beta)[1])
  sum(stats::dnorm(as.numeric(x), mean = mu, sd = as.numeric(s), log = TRUE))
}

.constant_log_hessian <- function(beta, x, s, beta_prec = 0) {
  beta <- as.numeric(beta)[1]
  mu <- exp(beta)
  w <- 1 / (as.numeric(s)^2)
  hess <- sum(w * mu * (2 * mu - as.numeric(x))) + as.numeric(beta_prec)
  if (is.finite(hess) && hess > 0) {
    return(as.numeric(hess))
  }

  nll <- function(b) {
    -.constant_log_observation_loglik(b, x = x, s = s) +
      0.5 * as.numeric(beta_prec) * b^2
  }
  eps <- 1e-4
  hess_fd <- (nll(beta + eps) - 2 * nll(beta) + nll(beta - eps)) / eps^2
  if (!is.finite(hess_fd) || hess_fd <= 0) {
    stop("The log-link constant beta Hessian is not positive at the fitted mode.")
  }
  as.numeric(hess_fd)
}

.fit_constant_log_known_noise <- function(x,
                                            s,
                                            beta_mode,
                                            beta_fixed = NULL,
                                            beta_prec = NULL) {
  x <- as.numeric(x)
  s <- as.numeric(s)
  n <- length(x)
  w <- 1 / (s^2)
  beta_hat <- NULL
  beta_var <- 0
  log_likelihood <- NULL
  log_likelihood_semantics <- paste0("exact_constant_log_", beta_mode)

  if (beta_mode == "fixed") {
    beta_hat <- .check_single_numeric(beta_fixed, "beta_fixed")
    log_likelihood <- .constant_log_observation_loglik(beta_hat, x = x, s = s)
  } else if (beta_mode == "empirical_bayes") {
    beta_hat <- .constant_log_positive_mean_beta(x, w = w)
    log_likelihood <- .constant_log_observation_loglik(beta_hat, x = x, s = s)
  } else if (beta_mode %in% c("prior_flat", "prior_proper")) {
    beta_prec0 <- if (identical(beta_mode, "prior_proper")) {
      beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
      if (is.null(beta_prec) || beta_prec <= 0) {
        stop("`beta_prec` must be positive for proper-prior beta.")
      }
      beta_prec
    } else {
      0
    }
    beta_start <- if (beta_prec0 == 0) {
      .constant_log_positive_mean_beta(x, w = w)
    } else {
      positive_x <- x[is.finite(x) & x > 0]
      log(if (length(positive_x)) mean(positive_x) else 1)
    }
    nll <- function(beta) {
      -.constant_log_observation_loglik(beta, x = x, s = s) +
        0.5 * beta_prec0 * beta^2
    }
    opt <- stats::optimize(nll, interval = .constant_log_beta_interval(beta_start))
    beta_hat <- as.numeric(opt$minimum)
    beta_hess <- .constant_log_hessian(beta_hat, x = x, s = s, beta_prec = beta_prec0)
    beta_var <- 1 / beta_hess
    log_prior <- if (beta_prec0 > 0) {
      stats::dnorm(beta_hat, mean = 0, sd = 1 / sqrt(beta_prec0), log = TRUE)
    } else {
      0
    }
    log_likelihood <- .constant_log_observation_loglik(beta_hat, x = x, s = s) +
      log_prior + 0.5 * log(2 * pi) - 0.5 * log(beta_hess)
    log_likelihood_semantics <- paste0("laplace_constant_log_", beta_mode)
  } else {
    stop("Unsupported beta mode: ", beta_mode)
  }

  class(log_likelihood) <- "logLik"
  list(
    posterior = .constant_response_moments(beta_hat, beta_var, n = n, link = "log"),
    fitted_g = .constant_beta_state(beta_mode, beta_hat, beta_prec = beta_prec),
    fitted_beta = beta_hat,
    beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec,
    beta_mode = beta_mode,
    beta_var = beta_var,
    log_likelihood = log_likelihood,
    log_likelihood_semantics = log_likelihood_semantics,
    posterior_sampler = .constant_response_sampler(beta_hat, beta_var, n = n, link = "log"),
    data = .eb_smoother_data_frame(x, s),
    prior_family = "log_constant",
    backend = "exact",
    link = "log"
  )
}

.fit_constant_known_noise <- function(x,
                                        s,
                                        beta_mode,
                                        beta_fixed = NULL,
                                        beta_prec = NULL,
                                        link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  if (identical(link, "log")) {
    return(.fit_constant_log_known_noise(
      x = x,
      s = s,
      beta_mode = beta_mode,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec
    ))
  } else if (identical(link, "softplus")) {
    stop("`family = \"constant\"` with `link = \"softplus\"` is not currently supported.")
  }

  x <- as.numeric(x)
  s <- as.numeric(s)
  n <- length(x)
  w <- 1 / (s^2)
  quad_u <- sum(w)
  cross_ux <- sum(w * x)
  quad_x <- sum(w * x^2)
  logdet_sigma <- sum(log(s^2))

  beta_hat <- NULL
  beta_var <- 0
  log_likelihood <- NULL

  if (beta_mode == "fixed") {
    beta_hat <- .check_single_numeric(beta_fixed, "beta_fixed")
    quad <- quad_x - 2 * beta_hat * cross_ux + beta_hat^2 * quad_u
    log_likelihood <- -0.5 * (n * log(2 * pi) + logdet_sigma + quad)
  } else if (beta_mode == "empirical_bayes") {
    beta_hat <- cross_ux / quad_u
    quad <- quad_x - (cross_ux^2) / quad_u
    log_likelihood <- -0.5 * (n * log(2 * pi) + logdet_sigma + quad)
  } else if (beta_mode == "prior_flat") {
    beta_hat <- cross_ux / quad_u
    beta_var <- 1 / quad_u
    quad <- quad_x - (cross_ux^2) / quad_u
    log_likelihood <- -0.5 * ((n - 1) * log(2 * pi) + logdet_sigma + log(quad_u) + quad)
  } else if (beta_mode == "prior_proper") {
    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    beta_hat <- cross_ux / (quad_u + beta_prec)
    beta_var <- 1 / (quad_u + beta_prec)
    quad <- quad_x - (cross_ux^2) / (quad_u + beta_prec)
    log_likelihood <- -0.5 * (
      n * log(2 * pi) +
        logdet_sigma +
        log(quad_u + beta_prec) -
        log(beta_prec) +
        quad
    )
  } else {
    stop("Unsupported beta mode: ", beta_mode)
  }

  class(log_likelihood) <- "logLik"
  posterior <- .constant_response_moments(beta_hat, beta_var, n = n, link = "identity")

  list(
    posterior = posterior,
    fitted_g = .constant_beta_state(beta_mode, beta_hat, beta_prec = beta_prec),
    fitted_beta = beta_hat,
    beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec,
    beta_mode = beta_mode,
    beta_var = beta_var,
    log_likelihood = log_likelihood,
    log_likelihood_semantics = paste0("exact_constant_", beta_mode),
    posterior_sampler = .constant_response_sampler(beta_hat, beta_var, n = n, link = "identity"),
    data = .eb_smoother_data_frame(x, s),
    prior_family = "identity_constant",
    backend = "exact",
    link = "identity"
  )
}

.fit_constant_log_unknown_noise <- function(x,
                                              beta_mode,
                                              beta_fixed = NULL,
                                              beta_prec = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 1L) {
    stop("`x` must contain at least one observation.")
  }

  beta_hat <- NULL
  beta_var <- 0
  noise_sd <- NULL
  log_likelihood <- NULL
  log_likelihood_semantics <- paste0("exact_constant_log_", beta_mode)

  if (beta_mode == "fixed") {
    beta_hat <- .check_single_numeric(beta_fixed, "beta_fixed")
    rss <- sum((x - exp(beta_hat))^2)
    sigma2_hat <- max(rss / n, .Machine$double.eps)
    noise_sd <- sqrt(sigma2_hat)
    log_likelihood <- -0.5 * n * (log(2 * pi) + 1 + log(sigma2_hat))
  } else if (beta_mode == "empirical_bayes") {
    beta_hat <- .constant_log_positive_mean_beta(x)
    rss <- sum((x - exp(beta_hat))^2)
    sigma2_hat <- max(rss / n, .Machine$double.eps)
    noise_sd <- sqrt(sigma2_hat)
    log_likelihood <- -0.5 * n * (log(2 * pi) + 1 + log(sigma2_hat))
  } else if (beta_mode %in% c("prior_flat", "prior_proper")) {
    beta_prec0 <- if (identical(beta_mode, "prior_proper")) {
      beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
      if (is.null(beta_prec) || beta_prec <= 0) {
        stop("`beta_prec` must be positive for proper-prior beta.")
      }
      beta_prec
    } else {
      0
    }
    scale0 <- stats::sd(x)
    if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
    eval_at_log_noise <- function(log_noise_sd) {
      noise <- exp(log_noise_sd)
      fit_known <- .fit_constant_log_known_noise(
        x = x,
        s = rep(noise, n),
        beta_mode = beta_mode,
        beta_fixed = beta_fixed,
        beta_prec = beta_prec
      )
      list(value = as.numeric(fit_known$log_likelihood), fit = fit_known, noise = noise)
    }
    objective <- function(log_noise_sd) {
      out <- tryCatch(eval_at_log_noise(log_noise_sd), error = function(e) e)
      if (inherits(out, "error")) return(1e100)
      -out$value
    }
    opt <- stats::optimize(objective, interval = log(scale0) + c(-12, 12))
    final <- eval_at_log_noise(opt$minimum)
    beta_hat <- final$fit$fitted_beta
    beta_var <- final$fit$beta_var
    noise_sd <- final$noise
    log_likelihood <- final$value
    log_likelihood_semantics <- paste0("laplace_constant_log_", beta_mode)
  } else {
    stop("Unsupported beta mode: ", beta_mode)
  }

  class(log_likelihood) <- "logLik"
  list(
    posterior = .constant_response_moments(beta_hat, beta_var, n = n, link = "log"),
    fitted_g = .constant_beta_state(beta_mode, beta_hat, beta_prec = beta_prec),
    fitted_beta = beta_hat,
    fitted_noise_sd = noise_sd,
    beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec,
    beta_mode = beta_mode,
    beta_var = beta_var,
    log_likelihood = log_likelihood,
    log_likelihood_semantics = log_likelihood_semantics,
    posterior_sampler = .constant_response_sampler(beta_hat, beta_var, n = n, link = "log"),
    data = .eb_smoother_data_frame(x, rep(noise_sd, n)),
    prior_family = "log_constant_learned_noise",
    backend = "exact",
    link = "log"
  )
}

.fit_constant_unknown_noise <- function(x,
                                          beta_mode,
                                          beta_fixed = NULL,
                                          beta_prec = NULL,
                                          link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  if (identical(link, "log")) {
    return(.fit_constant_log_unknown_noise(
      x = x,
      beta_mode = beta_mode,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec
    ))
  } else if (identical(link, "softplus")) {
    stop("`family = \"constant\"` with `link = \"softplus\"` is not currently supported.")
  }

  x <- as.numeric(x)
  n <- length(x)
  if (n < 1L) {
    stop("`x` must contain at least one observation.")
  }

  beta_hat <- NULL
  beta_var <- 0
  noise_sd <- NULL
  log_likelihood <- NULL

  if (beta_mode == "fixed") {
    beta_hat <- .check_single_numeric(beta_fixed, "beta_fixed")
    rss <- sum((x - beta_hat)^2)
    sigma2_hat <- max(rss / n, .Machine$double.eps)
    noise_sd <- sqrt(sigma2_hat)
    log_likelihood <- -0.5 * n * (log(2 * pi) + 1 + log(sigma2_hat))
  } else if (beta_mode == "empirical_bayes") {
    beta_hat <- mean(x)
    rss <- sum((x - beta_hat)^2)
    sigma2_hat <- max(rss / n, .Machine$double.eps)
    noise_sd <- sqrt(sigma2_hat)
    log_likelihood <- -0.5 * n * (log(2 * pi) + 1 + log(sigma2_hat))
  } else if (beta_mode == "prior_flat") {
    if (n < 2L) {
      stop("`beta_prec = 0` with learned noise requires at least two observations.")
    }
    beta_hat <- mean(x)
    rss <- sum((x - beta_hat)^2)
    sigma2_hat <- max(rss / (n - 1), .Machine$double.eps)
    noise_sd <- sqrt(sigma2_hat)
    beta_var <- sigma2_hat / n
    log_likelihood <- -0.5 * ((n - 1) * (log(2 * pi) + 1 + log(sigma2_hat)) + log(n))
  } else if (beta_mode == "prior_proper") {
    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    scale0 <- stats::sd(x)
    if (!is.finite(scale0) || scale0 <= 0) {
      scale0 <- 1
    }
    objective <- function(log_noise_sd) {
      sigma2 <- exp(2 * log_noise_sd)
      logdet <- n * log(sigma2) + log1p(n / (beta_prec * sigma2))
      quad <- sum(x^2) / sigma2 - (sum(x)^2) / (sigma2 * (beta_prec * sigma2 + n))
      -(-0.5 * (n * log(2 * pi) + logdet + quad))
    }
    opt <- stats::optimize(
      f = objective,
      interval = log(scale0) + c(-12, 12)
    )
    noise_sd <- exp(opt$minimum)
    sigma2_hat <- noise_sd^2
    beta_hat <- sum(x) / (beta_prec * sigma2_hat + n)
    beta_var <- 1 / (beta_prec + n / sigma2_hat)
    log_likelihood <- -opt$objective
  } else {
    stop("Unsupported beta mode: ", beta_mode)
  }

  class(log_likelihood) <- "logLik"
  posterior <- .constant_response_moments(beta_hat, beta_var, n = n, link = "identity")

  list(
    posterior = posterior,
    fitted_g = .constant_beta_state(beta_mode, beta_hat, beta_prec = beta_prec),
    fitted_beta = beta_hat,
    fitted_noise_sd = noise_sd,
    beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec,
    beta_mode = beta_mode,
    beta_var = beta_var,
    log_likelihood = log_likelihood,
    log_likelihood_semantics = paste0("exact_constant_", beta_mode),
    posterior_sampler = .constant_response_sampler(beta_hat, beta_var, n = n, link = "identity"),
    data = .eb_smoother_data_frame(x, rep(noise_sd, n)),
    prior_family = "identity_constant_learned_noise",
    backend = "exact",
    link = "identity"
  )
}

.matern_inla_validation_tolerances <- function(n) {
  n <- as.integer(n)
  list(
    log_likelihood_abs = 0.1 + 1e-4 * max(n, 1L),
    beta_abs = 0.05,
    log_range_abs = 0.15,
    log_sigma_abs = 0.15,
    log_noise_sd_abs = 0.20,
    posterior_mean_max_abs = 0.05
  )
}

.matern_validation_delta <- function(x, y) {
  x <- as.numeric(x)[1]
  y <- as.numeric(y)[1]
  if (!is.finite(x) || !is.finite(y)) return(Inf)
  abs(x - y)
}

.validate_matern_inla_against_laplace <- function(inla_fit,
                                                  reference_fit,
                                                  context,
                                                  tolerances = NULL) {
  if (is.null(tolerances)) {
    tolerances <- .matern_inla_validation_tolerances(nrow(inla_fit$posterior))
  }

  deltas <- list(
    log_likelihood_abs = .matern_validation_delta(inla_fit$log_likelihood, reference_fit$log_likelihood),
    beta_abs = .matern_validation_delta(inla_fit$fitted_beta, reference_fit$fitted_beta),
    log_range_abs = .matern_validation_delta(inla_fit$fitted_g$theta, reference_fit$fitted_g$theta),
    log_sigma_abs = .matern_validation_delta(log(inla_fit$fitted_g$sigma), log(reference_fit$fitted_g$sigma)),
    posterior_mean_max_abs = max(abs(inla_fit$posterior$mean - reference_fit$posterior$mean), na.rm = TRUE),
    posterior_var_max_abs = max(abs(inla_fit$posterior$var - reference_fit$posterior$var), na.rm = TRUE)
  )
  if (!is.null(inla_fit$fitted_noise_sd) || !is.null(reference_fit$fitted_noise_sd)) {
    deltas$log_noise_sd_abs <- .matern_validation_delta(
      log(inla_fit$fitted_noise_sd),
      log(reference_fit$fitted_noise_sd)
    )
  }

  checked <- intersect(names(tolerances), names(deltas))
  pass <- vapply(checked, function(nm) {
    is.finite(deltas[[nm]]) && deltas[[nm]] <= tolerances[[nm]]
  }, logical(1))

  validation <- list(
    status = if (all(pass)) "passed" else "failed",
    context = context,
    reference_backend = reference_fit$backend,
    reference_laplace_implementation = reference_fit$laplace_implementation,
    tolerances = tolerances,
    deltas = deltas,
    checked = checked
  )

  if (!all(pass)) {
    failed <- checked[!pass]
    failed_text <- paste(
      sprintf(
        "%s = %.6g > %.6g",
        failed,
        vapply(failed, function(nm) deltas[[nm]], numeric(1)),
        vapply(failed, function(nm) tolerances[[nm]], numeric(1))
      ),
      collapse = "; "
    )
    stop(
      "INLA Matern validation failed for ", context, ". ",
      "The INLA fit did not match the package Laplace reference: ",
      failed_text,
      ". Use `backend = \"laplace\"` for this case.",
      call. = FALSE
    )
  }

  validation
}

#' Fit an empirical-Bayes smoother with known or learned noise
#'
#' @description
#' Provides a unified front-end for the package's smoothing families.
#' Use `eb_smoother()` as the high-level smoothing interface when you want to
#' fit one of the package's prior families directly.
#'
#' The primary documented use case is learned-noise smoothing: when `s = NULL`,
#' `eb_smoother()` estimates one common observation noise SD together with the
#' prior-family parameters. This is the main package entry point when standard
#' errors are not known in advance.
#'
#' When `s` is supplied, `eb_smoother()` still supports the same smoothing
#' families with known observation standard errors. For known-`s` workflows
#' embedded in `ebnm`-compatible pipelines such as `flashier` and `ebmf`, the
#' lower-level `ebnm_LGP_generator()` and `ebnm_Matern_generator()` remain the
#' recommended interfaces.
#'
#' For the `"matern"` and `"constant"` families, the intercept is optimized
#' by marginal likelihood by default. The exact Gaussian implementations may
#' profile that intercept analytically for efficiency, but this is
#' mathematically equivalent to joint optimization. Supplying `beta_fixed`
#' instead holds the intercept fixed, while `beta_prec` can be used to switch
#' to a flat (`0`) or proper (`>0`) Gaussian prior on the intercept. These beta
#' semantics are on the linear-predictor scale for every link, including
#' `link = "log"` and `link = "softplus"`.
#'
#' @param x Numeric vector of observations.
#' @param s Observation standard errors. Supply a scalar or length-`length(x)`
#'   vector for known noise, or `NULL` to learn one common scalar SD.
#' @param family Smoother family. One of `"matern"`, `"constant"`,
#'   `"point_exponential"`, `"point_normal"`, `"point_laplace"`, or `"lgp"`.
#' @param backend Backend choice. For `family = "matern"`, `"auto"` uses the
#'   exact Gaussian backend for `link = "identity"` and `"fisher_pql"` for
#'   `link = "log"` and `link = "softplus"`. Matern also
#'   supports explicit `"laplace"`, `"laplace_fisher"`, `"fisher_pql"`,
#'   `"inla"`, and `"inlabru"`. The `"fisher_pql"` backend is an approximate
#'   Fisher/PQL backend for non-identity Matern fits. It uses three
#'   pseudo-Gaussian exact Matern Step A updates by default and reports a
#'   Fisher/Laplace score evaluated at the final PQL mode, not a true
#'   re-optimized original-model Laplace marginal likelihood. The `"inlabru"`
#'   backend runs
#'   the SPDE model through
#'   `inlabru`'s iterative-linearised INLA method and is currently only
#'   supported for `link = "softplus"`. It is an experimental validation
#'   backend, not an expected acceleration path. For known-noise fits,
#'   `pc.penalty = NULL` keeps the non-PC SPDE parameterisation; learned-noise
#'   inlabru fits with `s = NULL` require an explicit `pc.penalty` list with
#'   `range`, `sigma`, and `noise` entries.
#'   Compatibility aliases `"laplace_tmb"` and `"inla_pc"` are still accepted.
#'   Explicit `"laplace"`, `"laplace_tmb"`, and `"inla"` retain
#'   observed-Hessian Laplace semantics for log-link and softplus-link fits
#'   where supported. Public Matern Laplace backends use the package-owned TMB
#'   implementation or error; the R reference backend `"laplace_r"` is retained
#'   for internal validation only. For `family = "lgp"`, `"auto"` uses
#'   `"tmb"` for the identity link, `"laplace_fisher"` for the log link, and
#'   `"laplace"` for the softplus link. The `"constant"` and point-mass
#'   reference families support only the exact backend.
#' @param locations Optional raw spatial locations for the Matern family.
#' @param setup Family-specific setup object:
#'   [Matern_setup()] for `family = "matern"` or [LGP_setup()] for
#'   `family = "lgp"`.
#' @param g_init Optional prior-family-specific initialization object.
#'   When supported, `g_init$beta` supplies an initialization value for `beta`
#'   and `g_init$beta_prec` can provide a fallback prior precision.
#' @param fix_g Logical. If `TRUE`, supported prior-family parameters in
#'   `g_init` are held fixed while scoring or fitting. This is mainly useful
#'   for permutation scoring with `family = "matern"` or a point-mass
#'   reference family.
#' @param beta_fixed Optional fixed intercept. For `"matern"` and
#'   `"constant"`, `NULL` means optimize the intercept by marginal
#'   likelihood; a numeric scalar holds it fixed. For `"lgp"`, this should be
#'   a numeric vector of length `ncol(setup$X)`.
#' @param beta_prec Optional non-negative prior precision on `beta`.
#'   `NULL` means empirical-Bayes estimation, `0` requests a flat prior, and a
#'   positive value requests a proper zero-mean Gaussian prior.
#' @param fix_params Character vector naming prior parameters to fix while
#'   estimating the remaining parameters. For `"matern"`, allowed values are
#'   `"range"`, `"sigma"`, and `"beta"`; `"range"` uses `g_init$theta`,
#'   `"sigma"` uses `g_init$sigma`, and `"beta"` uses `beta_fixed` when
#'   supplied or otherwise `g_init$beta`. For `"lgp"`, allowed values are
#'   `"scale"` and `"beta"`; `"scale"` uses `g_init$scale`. `fix_g = TRUE`
#'   is retained as a shortcut for fixing Matern range/sigma or L-GP scale.
#' @param pc.penalty Optional Matern PC prior specification. Supported entries
#'   are `range`, `sigma`, and, when `s = NULL`, `noise`.
#' @param penalty_range Initial range anchor for the Matern family when
#'   `g_init$theta` is missing.
#' @param alpha Matern SPDE smoothness parameter.
#' @param max.edge Optional mesh size control for the Matern family. For
#'   two-dimensional Matern fits with `max.edge = NULL`, the default keeps all
#'   observed locations as mesh vertices and uses a coarser outer mesh for
#'   screening-oriented fits.
#' @param matern_n_starts Number of Step A optimization starts for the
#'   package-owned TMB Matern Laplace learned-noise path. The default is `1`;
#'   use `5` to reproduce the earlier multistart behavior.
#' @param profile_s_lower,profile_s_upper Lower and upper bounds for profiling
#'   the scalar observation noise SD when `family` is one of the point-mass
#'   reference families and `s = NULL`. If `profile_s_upper = NULL`, it is set
#'   from the scale of `x`.
#' @param profile_s_tol Optimization tolerance for the profiled
#'   point-mass-reference noise SD.
#' @param suppress_warnings If `TRUE`, suppress INLA mesh and SPDE warnings.
#' @param compute_exact_diagnostic If `TRUE`, store the exact Gaussian
#'   log-likelihood evaluated at the Matern INLA mode.
#' @param link Link used by the selected family. The Matern family supports
#'   `"identity"`, `"log"`, and `"softplus"` fits with known or learned noise.
#'   The constant family supports `"identity"` and `"log"` constant-baseline
#'   fits. The L-GP family supports `"identity"`, `"log"`, and `"softplus"`.
#'   Softplus posterior moments are deterministic Gauss-Hermite transforms of
#'   each marginal Gaussian Laplace posterior for the linear predictor.
#' @param ... Additional arguments. `dll` is used for the learned-noise L-GP TMB
#'   path. `mode`, `scale`, `output`, `optmethod`, and `control` are forwarded
#'   to the selected `ebnm` point-mass family.
#'
#' @return An object of class `c("eb_smoother_fit", "list")`. Use `print(fit)`
#'   for a compact overview and `summary(fit)` for posterior and diagnostic
#'   summaries of the stored fit.
#'
#' @examples
#' loc <- seq(0, 1, length.out = 40)
#' x <- sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)
#' fit_learned <- eb_smoother(
#'   x,
#'   s = NULL,
#'   family = "matern",
#'   locations = loc,
#'   pc.penalty = list(range = 0.2, sigma = 0.3, noise = 0.1)
#' )
#' fit_learned
#' summary(fit_learned)
#' head(fit_learned$posterior)
#'
#' fit_baseline <- eb_smoother(x, s = NULL, family = "constant", beta_prec = 0)
#' summary(fit_baseline)
#'
#' x_pos <- exp(0.1 + 0.4 * sin(2 * pi * loc)) + rnorm(length(loc), sd = 0.1)
#' fit_positive <- eb_smoother(
#'   x_pos,
#'   s = NULL,
#'   family = "matern",
#'   locations = loc,
#'   link = "log"
#' )
#' fit_positive_baseline <- eb_smoother(x_pos, s = NULL, family = "constant", link = "log")
#'
#' fit_sigma_fixed <- eb_smoother(
#'   x,
#'   s = 0.1,
#'   family = "matern",
#'   locations = loc,
#'   g_init = Matern(sigma = 1),
#'   fix_params = "sigma"
#' )
#'
#' fit_sparse_baseline <- eb_smoother(x_pos, s = NULL, family = "point_exponential")
#' as.numeric(fit_positive$log_likelihood) - as.numeric(fit_sparse_baseline$log_likelihood)
#'
#' @export
eb_smoother <- function(x,
                        s = NULL,
                        family = c("matern", "constant", "point_exponential", "point_normal", "point_laplace", "lgp"),
                        backend = c("auto", "exact", "laplace", "laplace_fisher", "fisher_pql", "inla", "inlabru"),
                        locations = NULL,
                        setup = NULL,
                        g_init = NULL,
                        fix_g = FALSE,
                        beta_fixed = NULL,
                        beta_prec = NULL,
                        fix_params = character(),
                        pc.penalty = NULL,
                        penalty_range = NULL,
                        alpha = 2,
                        max.edge = NULL,
                        matern_n_starts = 1L,
                        profile_s_lower = 1e-8,
                        profile_s_upper = NULL,
                        profile_s_tol = .Machine$double.eps^0.25,
                        suppress_warnings = TRUE,
                        compute_exact_diagnostic = FALSE,
                        link = c("identity", "log", "softplus"),
                        ...) {
  beta_prec_missing <- missing(beta_prec)
  alpha_missing <- missing(alpha)
  family <- match.arg(family)
  backend <- .match_matern_backend_arg(backend)
  link <- match.arg(link)
  if (!is.logical(fix_g) || length(fix_g) != 1L || is.na(fix_g)) {
    stop("`fix_g` must be TRUE or FALSE.")
  }

  x <- as.numeric(x)
  if (anyNA(x)) stop("`x` must not contain NA.")

  extra_args <- list(...)
  point_extra <- c("mode", "scale", "output", "optmethod", "control")
  allowed_extra <- c("dll", point_extra)
  extra_nms <- names(extra_args)
  if (is.null(extra_nms)) extra_nms <- rep("", length(extra_args))
  unknown_extra <- extra_nms[!(extra_nms %in% c("", allowed_extra))]
  if (length(unknown_extra) > 0) {
    stop("Unused arguments: ", paste(sprintf("`%s`", unique(unknown_extra)), collapse = ", "))
  }
  dll <- if (!is.null(extra_args$dll)) extra_args$dll else "EBSmoothr"
  point_extra_supplied <- intersect(extra_nms, point_extra)
  if (!(family %in% .point_family_names) && length(point_extra_supplied) > 0L) {
    stop(
      paste(sprintf("`%s`", unique(point_extra_supplied)), collapse = ", "),
      " can only be supplied for point-mass reference families."
    )
  }
  if (!identical(family, "lgp") && !is.null(extra_args$dll)) {
    stop("`dll` is only used when `family = \"lgp\"`.")
  }

  s_known <- !is.null(s)
  s_vec <- if (s_known) .normalize_eb_smoother_s(s, length(x)) else NULL

  if (family == "matern") {
    if (is.null(locations) && is.null(setup)) {
      stop("Either `locations` or `setup` must be supplied when `family = \"matern\"`.")
    }
  } else if (family %in% c("constant", .point_family_names)) {
    if (!is.null(locations) || !is.null(setup)) {
      stop("`locations` and `setup` are not used when `family = \"", family, "\"`.")
    }
  } else {
    if (is.null(setup)) stop("`setup` must be supplied when `family = \"lgp\"`.")
  }

  if (family == "lgp" && backend %in% c("laplace_r", "laplace_tmb", "inla", "inla_pc")) {
    stop("`backend` values `laplace_r`, `laplace_tmb`, `inla`, and `inla_pc` are only supported for `family = \"matern\"`.")
  }
  if (family != "matern" && !is.null(pc.penalty)) {
    stop("`pc.penalty` is only supported for `family = \"matern\"`.")
  }

  if (family %in% .point_family_names) {
    if (!identical(link, "identity")) {
      stop("`family = \"", family, "\"` is defined on the identity observation scale; use `link = \"identity\"`.")
    }
    backend_use <- if (backend == "auto") "exact" else backend
    if (!identical(backend_use, "exact")) {
      stop("`family = \"", family, "\"` only supports the exact backend.")
    }
    if (isTRUE(fix_g) && is.null(g_init)) {
      stop("`g_init` must be supplied when `fix_g = TRUE` for `family = \"", family, "\"`.")
    }
    if (!is.null(beta_fixed) || (!beta_prec_missing && !is.null(beta_prec))) {
      stop("`beta_fixed` and `beta_prec` are not used when `family = \"", family, "\"`.")
    }
    if (length(fix_params) > 0L) {
      stop("`fix_params` is not used when `family = \"", family, "\"`.")
    }
    if (!is.null(penalty_range)) stop("`penalty_range` is only supported for `family = \"matern\"`.")
    if (!is.null(max.edge)) stop("`max.edge` is only supported for `family = \"matern\"`.")
    if (!alpha_missing) stop("`alpha` is only supported for `family = \"matern\"`.")
    if (isTRUE(compute_exact_diagnostic)) {
      stop("`compute_exact_diagnostic` is only supported for `family = \"matern\"`.")
    }

    fit <- .fit_point_family_smoother(
      x = x,
      s = if (s_known) s_vec else NULL,
      family = family,
      g_init = g_init,
      fix_g = fix_g,
      mode = if (!is.null(extra_args$mode)) extra_args$mode else 0,
      scale = if (!is.null(extra_args$scale)) extra_args$scale else "estimate",
      output = if (!is.null(extra_args$output)) extra_args$output else ebnm::ebnm_output_default(),
      optmethod = extra_args$optmethod,
      control = extra_args$control,
      profile_s_lower = profile_s_lower,
      profile_s_upper = profile_s_upper,
      profile_s_tol = profile_s_tol
    )
    return(
      .as_eb_smoother_fit(
        fit = fit,
        family = family,
        noise_mode = if (s_known) "fixed" else "estimated",
        fitted_noise_sd = if (s_known) {
          if (length(s) == 1L) as.numeric(s) else NA_real_
        } else {
          fit$fitted_noise_sd
        },
        backend = "exact"
      )
    )
  }

  if (family == "matern") {
    matern_n_starts <- .check_matern_n_starts(matern_n_starts)
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
      beta_prec_supplied = !beta_prec_missing,
      g_init = g_init,
      expected_length = 1L
    )
    beta_spec <- .eb_smoother_resolve_beta_spec(
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      g_init_beta_prec = if (is.null(g_init)) NULL else g_init$beta_prec
    )
    beta_mode <- beta_spec$mode
    beta_fixed_use <- beta_spec$beta_fixed
    beta_prec_use <- beta_spec$beta_prec

    if (s_known && !is.null(pc.penalty) && !is.null(names(pc.penalty)) && "noise" %in% names(pc.penalty)) {
      stop("`pc.penalty$noise` is only supported when `s = NULL` on the Matern PC-prior path.")
    }

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
    if (length(x) != nrow(loc_mat)) {
      stop("The length of `x` must equal the number of locations.")
    }

    mesh <- setup0$mesh
    A <- setup0$A
    spde_template <- setup0$spde_template
    penalty_range0 <- setup0$penalty_range

    resolved_init <- .resolve_matern_g_init(
      x = x,
      s = if (s_known) s_vec else NULL,
      g_init = g_init,
      beta_fixed = beta_fixed_use,
      beta_prec = beta_prec_use,
      penalty_range0 = penalty_range0,
      pc.penalty = pc.penalty,
      allow_noise = !s_known,
      link = link
    )

    backend_use <- .matern_resolve_backend(
      backend = backend,
      link = link,
      d = d,
      beta_mode = beta_mode,
      pc_penalty = resolved_init$pc_penalty,
      learn_noise = !s_known,
      fix_g = all(c("range", "sigma") %in% fix_params_use)
    )
    if (identical(backend_use, "inla") &&
        any(c("range", "sigma") %in% fix_params_use)) {
      stop("`fix_params` values `range` and `sigma` are not supported with the INLA Matern backend.")
    }

    if (identical(backend_use, "fisher_pql")) {
      fit <- if (s_known) .fit_matern_fisher_pql_known_noise(
        x = x,
        s = s_vec,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = resolved_init$theta_init,
        sigma_init = resolved_init$sigma_init,
        beta_init = resolved_init$beta_init,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        pc_penalty = resolved_init$pc_penalty,
        link = link,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        fix_params = fix_params_use
      ) else .fit_matern_fisher_pql_unknown_noise(
        x = x,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = resolved_init$theta_init,
        sigma_init = resolved_init$sigma_init,
        noise_sd_init = resolved_init$noise_sd_init,
        beta_init = resolved_init$beta_init,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        pc_penalty = resolved_init$pc_penalty,
        link = link,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        fix_params = fix_params_use
      )
      fit$data <- .eb_smoother_data_frame(x, if (s_known) s_vec else fit$fitted_s)
      fit$g_init <- resolved_init$g_init
      return(
        .as_eb_smoother_fit(
          fit = fit,
          family = "matern",
          noise_mode = if (s_known) "fixed" else "estimated",
          fitted_noise_sd = if (s_known) {
            if (length(s) == 1L) as.numeric(s) else NA_real_
          } else {
            fit$fitted_noise_sd
          },
          backend = "fisher_pql"
        )
      )
    }

    if (backend_use %in% c("laplace", "laplace_fisher", "laplace_r", "laplace_tmb")) {
      fit <- if (s_known) .fit_matern_laplace_dispatch_known_noise(
        x = x,
        s = s_vec,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = resolved_init$theta_init,
        sigma_init = resolved_init$sigma_init,
        beta_init = resolved_init$beta_init,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        pc_penalty = resolved_init$pc_penalty,
        link = link,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        matern_n_starts = matern_n_starts,
        backend_use = backend_use,
        fix_params = fix_params_use
      ) else .fit_matern_laplace_dispatch_unknown_noise(
        x = x,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = resolved_init$theta_init,
        sigma_init = resolved_init$sigma_init,
        noise_sd_init = resolved_init$noise_sd_init,
        beta_init = resolved_init$beta_init,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        pc_penalty = resolved_init$pc_penalty,
        link = link,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        matern_n_starts = matern_n_starts,
        backend_use = backend_use,
        fix_params = fix_params_use
      )
      fit$data <- .eb_smoother_data_frame(x, if (s_known) s_vec else rep(fit$fitted_noise_sd, length(x)))
      fit$g_init <- resolved_init$g_init
      return(
        .as_eb_smoother_fit(
          fit = fit,
          family = "matern",
          noise_mode = if (s_known) "fixed" else "estimated",
          fitted_noise_sd = if (s_known) {
            if (length(s) == 1L) as.numeric(s) else NA_real_
          } else {
            fit$fitted_noise_sd
          },
          backend = fit$backend
        )
      )
    }

    if (identical(backend_use, "inlabru")) {
      if (any(c("range", "sigma") %in% fix_params_use)) {
        stop("`fix_params` values `range` and `sigma` are not supported with the inlabru Matern backend.")
      }
      pc_penalty_inlabru <- .matern_inlabru_pc_penalty_policy(
        pc_penalty_arg = pc.penalty,
        resolved_pc_penalty = resolved_init$pc_penalty,
        learn_noise = !s_known
      )
      theta_init_inlabru <- .matern_inla_theta_init(
        log_range = resolved_init$theta_init,
        log_sigma = log(resolved_init$sigma_init),
        alpha = alpha,
        d = d,
        pc_penalty = pc_penalty_inlabru
      )
      if (s_known) {
        fit <- if (identical(beta_mode, "fixed")) {
          .fit_matern_inlabru_stepA_fixed_beta(
            x = x, s = s_vec, A = A, mesh = mesh, alpha = alpha, d = d,
            locations = loc_mat,
            theta_init = theta_init_inlabru,
            beta_fixed = beta_fixed_use,
            pc_penalty = pc_penalty_inlabru,
            link = link,
            suppress_warnings = suppress_warnings
          )
        } else {
          .fit_matern_inlabru_stepA(
            x = x, s = s_vec, A = A, mesh = mesh, alpha = alpha, d = d,
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
          fit = fit,
          x = x,
          s = s_vec,
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d,
          beta_mode = beta_mode,
          beta_fixed = beta_fixed_use,
          beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
          link = link,
          pc_penalty = pc_penalty_inlabru,
          suppress_warnings = suppress_warnings
        )
        log_likelihood <- as.numeric(comparable_objective$log_marginal)
        class(log_likelihood) <- "logLik"
        out <- list(
          posterior = fit$posterior,
          fitted_g = fit$fitted_g,
          fitted_beta = fit$fitted_beta,
          beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
          beta_mode = beta_mode,
          g_init = resolved_init$g_init,
          log_likelihood = log_likelihood,
          log_likelihood_semantics = paste0("laplace_at_inlabru_params_", beta_mode),
          posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
          data = .eb_smoother_data_frame(x, s_vec),
          prior_family = paste0(link, "_Matern_inlabru"),
          posterior_spatial_field = fit$posterior_spatial_field,
          mesh = mesh,
          inla_result = fit$result,
          backend = "inlabru",
          link = link,
          pc_penalty = pc_penalty_inlabru,
          log_likelihood_laplace_at_inlabru_params = as.numeric(comparable_objective$log_marginal),
          log_likelihood_inlabru_mlik_integration = fit$log_likelihood_inlabru_mlik_integration,
          log_likelihood_stepA_mlik_integration = fit$log_likelihood_stepA_mlik_integration,
          beta_profile_optimization = fit$beta_profile_optimization,
          beta_profile_objective = fit$beta_profile_objective,
          matern_objective_context = list(A = A, spde_template = spde_template,
                                           alpha = alpha, d = d)
        )
        return(
          .as_eb_smoother_fit(
            fit = out,
            family = "matern",
            noise_mode = "fixed",
            fitted_noise_sd = if (length(s) == 1L) as.numeric(s) else NA_real_,
            backend = "inlabru"
          )
        )
      }

      fit <- if (identical(beta_mode, "fixed")) {
        .fit_matern_inlabru_stepA_unknown_noise_fixed_beta(
          x = x, A = A, mesh = mesh, alpha = alpha, d = d,
          locations = loc_mat,
          theta_init = theta_init_inlabru,
          noise_sd_init = resolved_init$noise_sd_init,
          beta_fixed = beta_fixed_use,
          pc_penalty = pc_penalty_inlabru,
          link = link,
          suppress_warnings = suppress_warnings
        )
      } else {
        .fit_matern_inlabru_stepA_unknown_noise(
          x = x, A = A, mesh = mesh, alpha = alpha, d = d,
          locations = loc_mat,
          theta_init = theta_init_inlabru,
          noise_sd_init = resolved_init$noise_sd_init,
          beta_init = resolved_init$beta_init,
          beta_prec = if (identical(beta_mode, "prior_proper")) beta_prec_use else 0,
          pc_penalty = pc_penalty_inlabru,
          link = link,
          suppress_warnings = suppress_warnings
        )
      }

      comparable_objective <- .matern_inlabru_unknown_laplace_objective(
        fit = fit,
        x = x,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
        link = link,
        pc_penalty = pc_penalty_inlabru,
        suppress_warnings = suppress_warnings
      )
      log_likelihood <- as.numeric(comparable_objective$log_marginal)
      class(log_likelihood) <- "logLik"
      out <- list(
        posterior = fit$posterior,
        fitted_g = fit$fitted_g,
        fitted_beta = fit$fitted_beta,
        fitted_noise_sd = fit$fitted_noise_sd,
        beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
        beta_mode = beta_mode,
        g_init = resolved_init$g_init,
        log_likelihood = log_likelihood,
        log_likelihood_semantics = paste0("laplace_at_inlabru_params_learned_noise_", beta_mode),
        posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
        data = .eb_smoother_data_frame(x, rep(fit$fitted_noise_sd, length(x))),
        prior_family = paste0(link, "_Matern_inlabru_learned_noise"),
        posterior_spatial_field = fit$posterior_spatial_field,
        mesh = mesh,
        inla_result = fit$result,
        backend = "inlabru",
        link = link,
        pc_penalty = pc_penalty_inlabru,
        log_likelihood_laplace_at_inlabru_params = as.numeric(comparable_objective$log_marginal),
        log_likelihood_inlabru_mlik_integration = fit$log_likelihood_inlabru_mlik_integration,
        log_likelihood_stepA_mlik_integration = fit$log_likelihood_stepA_mlik_integration,
        beta_profile_optimization = fit$beta_profile_optimization,
        beta_profile_objective = fit$beta_profile_objective,
        matern_objective_context = list(A = A, spde_template = spde_template,
                                         alpha = alpha, d = d)
      )
      return(
        .as_eb_smoother_fit(
          fit = out,
          family = "matern",
          noise_mode = "estimated",
          fitted_noise_sd = fit$fitted_noise_sd,
          backend = "inlabru"
        )
      )
    }

    if (backend_use == "exact") {
      if (!identical(link, "identity")) {
        stop("`backend = \"exact\"` is only available for `link = \"identity\"`; use `backend = \"laplace\"` for log-link Matern fits.")
      }
      if (s_known) {
        fit <- .fit_matern_exact_known_noise(
          x = x,
          s = s_vec,
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d,
          theta_init = resolved_init$theta_init,
          sigma_init = resolved_init$sigma_init,
          beta_mode = beta_mode,
          beta_fixed = beta_fixed_use,
          beta_prec = beta_prec_use,
          pc_penalty = resolved_init$pc_penalty,
          suppress_warnings = suppress_warnings,
          fix_g = fix_g,
          fix_params = fix_params_use
        )
        fit$data <- .eb_smoother_data_frame(x, s_vec)
        fit$g_init <- resolved_init$g_init
        return(
          .as_eb_smoother_fit(
            fit = fit,
            family = "matern",
            noise_mode = "fixed",
            fitted_noise_sd = if (length(s) == 1L) as.numeric(s) else NA_real_,
            backend = "exact"
          )
        )
      }

      fit <- .fit_matern_exact_unknown_noise(
        x = x,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = resolved_init$theta_init,
        sigma_init = resolved_init$sigma_init,
        noise_sd_init = resolved_init$noise_sd_init,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        pc_penalty = resolved_init$pc_penalty,
        suppress_warnings = suppress_warnings,
        fix_g = fix_g,
        fix_params = fix_params_use
      )
      fit$data <- .eb_smoother_data_frame(x, fit$fitted_s)
      fit$g_init <- resolved_init$g_init
      return(
        .as_eb_smoother_fit(
          fit = fit,
          family = "matern",
          noise_mode = "estimated",
          fitted_noise_sd = fit$fitted_noise_sd,
          backend = "exact"
        )
      )
    }

    if (s_known) {
      theta_init_inla <- .matern_inla_theta_init(
        log_range = resolved_init$theta_init,
        log_sigma = log(resolved_init$sigma_init),
        alpha = alpha,
        d = d,
        pc_penalty = resolved_init$pc_penalty
      )
      fit <- if (identical(beta_mode, "fixed")) {
        .fit_matern_inla_stepA_fixed_beta(
          x = x,
          s = s_vec,
          A = A,
          mesh = mesh,
          alpha = alpha,
          d = d,
          theta_init = theta_init_inla,
          beta_fixed = beta_fixed_use,
          pc_penalty = resolved_init$pc_penalty,
          link = link,
          suppress_warnings = suppress_warnings
        )
      } else if (identical(beta_mode, "empirical_bayes")) {
        .fit_matern_inla_stepA_profile_beta(
          x = x,
          s = s_vec,
          A = A,
          mesh = mesh,
          spde_template = spde_template,
          alpha = alpha,
          d = d,
          theta_init = theta_init_inla,
          beta_init = resolved_init$beta_init,
          pc_penalty = resolved_init$pc_penalty,
          link = link,
          suppress_warnings = suppress_warnings
        )
      } else {
        .fit_matern_inla_stepA(
          x = x,
          s = s_vec,
          A = A,
          mesh = mesh,
          alpha = alpha,
          d = d,
          theta_init = theta_init_inla,
          beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
          pc_penalty = resolved_init$pc_penalty,
          link = link,
          suppress_warnings = suppress_warnings
        )
      }

      comparable_at_inla_mode <- NULL
      if (identical(link, "log")) {
        comparable_at_inla_mode <- .matern_laplace_known_noise_objective_at_params(
          x = x,
          s = s_vec,
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d,
          log_range = fit$fitted_g$theta,
          log_sigma = log(fit$fitted_g$sigma),
          beta_mode = beta_mode,
          beta_fixed = beta_fixed_use,
          beta_prec = beta_prec_use,
          beta_init = fit$fitted_beta,
          link = link,
          pc_penalty = resolved_init$pc_penalty,
          initial_mode = .matern_inla_laplace_initial_mode(fit, beta_mode),
          compute_posterior = FALSE,
          optimize_beta = FALSE,
          suppress_warnings = suppress_warnings
        )
      } else {
        comparable_at_inla_mode <- .exact_matern_known_noise_objective_at_params(
          x = x,
          s = s_vec,
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d,
          log_range = fit$fitted_g$theta,
          log_sigma = log(fit$fitted_g$sigma),
          beta_mode = if (identical(beta_mode, "empirical_bayes")) "fixed" else beta_mode,
          beta_fixed = if (identical(beta_mode, "empirical_bayes")) fit$fitted_beta else beta_fixed_use,
          beta_prec = beta_prec_use,
          pc_penalty = resolved_init$pc_penalty,
          suppress_warnings = suppress_warnings
        )
      }

      log_likelihood <- comparable_at_inla_mode$log_marginal
      class(log_likelihood) <- "logLik"

      out <- list(
        posterior = fit$posterior,
        fitted_g = fit$fitted_g,
        fitted_beta = fit$fitted_beta,
        beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
        beta_mode = beta_mode,
        g_init = resolved_init$g_init,
        log_likelihood = log_likelihood,
        log_likelihood_semantics = if (identical(link, "log")) {
          paste0("laplace_at_inla_mode_", beta_mode)
        } else {
          paste0("exact_at_inla_mode_", beta_mode)
        },
        posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
        data = .eb_smoother_data_frame(x, s_vec),
        prior_family = paste0(link, "_Matern", if (is.null(resolved_init$pc_penalty)) "_inla" else "_pc_inla"),
        posterior_spatial_field = fit$posterior_spatial_field,
        mesh = mesh,
        inla_result = fit$result,
        backend = "inla",
        link = link,
        pc_penalty = resolved_init$pc_penalty,
        log_likelihood_stepA_penalized = fit$log_likelihood_stepA_penalized,
        log_likelihood_stepA_mlik_integration = fit$log_likelihood_stepA_mlik_integration,
        log_likelihood_stepA_mlik_gaussian = fit$log_likelihood_stepA_mlik_gaussian,
        log_likelihood_stepA_joint_log_posterior = fit$log_likelihood_stepA_joint_log_posterior,
        log_likelihood_laplace_at_inla_mode = if (identical(link, "log")) comparable_at_inla_mode$log_marginal else NULL,
        beta_profile_optimization = fit$beta_profile_optimization,
        beta_profile_objective = fit$beta_profile_objective,
        matern_objective_context = list(
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d
        )
      )

      if (identical(link, "identity")) {
        out$log_likelihood_exact_at_stepA_mode <- comparable_at_inla_mode$log_marginal
      } else {
        out$log_likelihood_exact_at_stepA_mode <- NULL
      }

      return(
        .as_eb_smoother_fit(
          fit = out,
          family = "matern",
          noise_mode = "fixed",
          fitted_noise_sd = if (length(s) == 1L) as.numeric(s) else NA_real_,
          backend = "inla"
        )
      )
    }

    theta_init_inla <- .matern_inla_theta_init(
      log_range = resolved_init$theta_init,
      log_sigma = log(resolved_init$sigma_init),
      alpha = alpha,
      d = d,
      pc_penalty = resolved_init$pc_penalty
    )
    fit <- if (identical(link, "log")) {
      .fit_matern_inla_stepA_unknown_noise_log_profile_noise(
        x = x,
        A = A,
        mesh = mesh,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_init_inla,
        noise_sd_init = resolved_init$noise_sd_init,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        beta_init = resolved_init$beta_init,
        pc_penalty = resolved_init$pc_penalty,
        suppress_warnings = suppress_warnings
      )
    } else if (identical(beta_mode, "fixed")) {
      .fit_matern_inla_stepA_unknown_noise_fixed_beta(
        x = x,
        A = A,
        mesh = mesh,
        alpha = alpha,
        d = d,
        theta_init = theta_init_inla,
        noise_sd_init = resolved_init$noise_sd_init,
        beta_fixed = beta_fixed_use,
        pc_penalty = resolved_init$pc_penalty,
        link = link,
        suppress_warnings = suppress_warnings
      )
    } else if (identical(beta_mode, "empirical_bayes")) {
      .fit_matern_inla_stepA_unknown_noise_profile_beta(
        x = x,
        A = A,
        mesh = mesh,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = theta_init_inla,
        noise_sd_init = resolved_init$noise_sd_init,
        beta_init = resolved_init$beta_init,
        pc_penalty = resolved_init$pc_penalty,
        link = link,
        suppress_warnings = suppress_warnings
      )
    } else {
      .fit_matern_inla_stepA_unknown_noise(
        x = x,
        A = A,
        mesh = mesh,
        alpha = alpha,
        d = d,
        theta_init = theta_init_inla,
        noise_sd_init = resolved_init$noise_sd_init,
        beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
        pc_penalty = resolved_init$pc_penalty,
        link = link,
        suppress_warnings = suppress_warnings
      )
    }

    comparable_at_inla_mode <- NULL
    if (identical(link, "log")) {
      comparable_at_inla_mode <- .matern_laplace_unknown_noise_objective_at_params(
        x = x,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = fit$fitted_g$theta,
        log_sigma = log(fit$fitted_g$sigma),
        log_noise_sd = log(fit$fitted_noise_sd),
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        beta_init = fit$fitted_beta,
        link = link,
        pc_penalty = resolved_init$pc_penalty,
        initial_mode = .matern_inla_laplace_initial_mode(fit, beta_mode),
        compute_posterior = FALSE,
        optimize_beta = FALSE,
        suppress_warnings = suppress_warnings
      )
    } else {
      comparable_at_inla_mode <- .exact_matern_unknown_noise_objective_at_params(
        x = x,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        log_range = fit$fitted_g$theta,
        log_sigma = log(fit$fitted_g$sigma),
        log_noise_sd = log(fit$fitted_noise_sd),
        beta_mode = if (identical(beta_mode, "empirical_bayes")) "fixed" else beta_mode,
        beta_fixed = if (identical(beta_mode, "empirical_bayes")) fit$fitted_beta else beta_fixed_use,
        beta_prec = beta_prec_use,
        pc_penalty = resolved_init$pc_penalty,
        suppress_warnings = suppress_warnings
      )
    }

    log_likelihood <- comparable_at_inla_mode$log_marginal
    class(log_likelihood) <- "logLik"

    out <- list(
      posterior = fit$posterior,
      fitted_g = fit$fitted_g,
      fitted_beta = fit$fitted_beta,
      fitted_noise_sd = fit$fitted_noise_sd,
      beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
      beta_mode = beta_mode,
      g_init = resolved_init$g_init,
      log_likelihood = log_likelihood,
      log_likelihood_semantics = if (identical(link, "log")) {
        paste0("laplace_at_inla_mode_", beta_mode)
      } else {
        paste0("exact_at_inla_mode_", beta_mode)
      },
      posterior_sampler = function(nsamp) .posterior_sampler_unavailable(nsamp, length(x)),
      data = .eb_smoother_data_frame(x, rep(fit$fitted_noise_sd, length(x))),
      prior_family = paste0(link, "_Matern", if (is.null(resolved_init$pc_penalty)) "_inla_learned_noise" else "_pc_inla_learned_noise"),
      posterior_spatial_field = fit$posterior_spatial_field,
      mesh = mesh,
      inla_result = fit$result,
      backend = "inla",
      link = link,
      pc_penalty = resolved_init$pc_penalty,
      log_likelihood_stepA_penalized = fit$log_likelihood_stepA_penalized,
      log_likelihood_stepA_mlik_integration = fit$log_likelihood_stepA_mlik_integration,
      log_likelihood_stepA_mlik_gaussian = fit$log_likelihood_stepA_mlik_gaussian,
      log_likelihood_stepA_joint_log_posterior = fit$log_likelihood_stepA_joint_log_posterior,
      log_likelihood_laplace_at_inla_mode = if (identical(link, "log")) comparable_at_inla_mode$log_marginal else NULL,
      beta_profile_optimization = fit$beta_profile_optimization,
      beta_profile_objective = fit$beta_profile_objective,
      noise_profile_optimization = fit$noise_profile_optimization,
      matern_objective_context = list(
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d
      )
    )

    if (identical(link, "identity")) {
      out$log_likelihood_exact_at_stepA_mode <- comparable_at_inla_mode$log_marginal
    } else {
      out$log_likelihood_exact_at_stepA_mode <- NULL
    }

    if (identical(link, "log")) {
      reference_fit <- .fit_matern_laplace_dispatch_unknown_noise(
        x = x,
        A = A,
        spde_template = spde_template,
        alpha = alpha,
        d = d,
        theta_init = resolved_init$theta_init,
        sigma_init = resolved_init$sigma_init,
        noise_sd_init = resolved_init$noise_sd_init,
        beta_init = resolved_init$beta_init,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        pc_penalty = resolved_init$pc_penalty,
        link = link,
        suppress_warnings = suppress_warnings,
        matern_n_starts = matern_n_starts,
        backend_use = "laplace"
      )
      out$inla_validation <- .validate_matern_inla_against_laplace(
        inla_fit = out,
        reference_fit = reference_fit,
        context = paste0(
          "log-link learned-noise Matern ",
          if (is.null(resolved_init$pc_penalty)) "without PC prior" else "with PC prior"
        )
      )
    }

    return(
      .as_eb_smoother_fit(
        fit = out,
        family = "matern",
        noise_mode = "estimated",
        fitted_noise_sd = out$fitted_noise_sd,
        backend = "inla"
      )
    )
  }

  if (family == "constant") {
    backend_use <- if (backend == "auto") "exact" else backend
    if (backend_use != "exact") {
      stop("`family = \"constant\"` only supports the exact backend.")
    }
    if (isTRUE(fix_g)) {
      stop("`fix_g` is not used when `family = \"constant\"`.")
    }
    if (length(fix_params) > 0L) {
      stop("`fix_params` is not used when `family = \"constant\"`.")
    }

    beta_spec <- .eb_smoother_resolve_beta_spec(
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      g_init_beta_prec = if (is.null(g_init)) NULL else g_init$beta_prec
    )
    beta_mode <- beta_spec$mode
    beta_fixed_use <- beta_spec$beta_fixed
    beta_prec_use <- beta_spec$beta_prec
    beta_init <- if (!is.null(beta_fixed_use)) {
      beta_fixed_use
    } else if (!is.null(g_init) && !is.null(g_init$beta)) {
      .check_optional_beta_vector(g_init$beta, "g_init$beta", expected_length = 1L, allow_null = FALSE)
    } else if (identical(link, "log") && s_known) {
      .constant_log_positive_mean_beta(x, w = 1 / (s_vec^2), nm = "initial weighted mean")
    } else if (identical(link, "log")) {
      .constant_log_positive_mean_beta(x, nm = "initial mean")
    } else if (s_known) {
      stats::weighted.mean(x, w = 1 / (s_vec^2))
    } else {
      mean(x)
    }
    g_init_resolved <- Constant(
      beta = beta_init,
      beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use
    )

    if (s_known) {
      fit <- .fit_constant_known_noise(
        x = x,
        s = s_vec,
        beta_mode = beta_mode,
        beta_fixed = beta_fixed_use,
        beta_prec = beta_prec_use,
        link = link
      )
      fit$g_init <- g_init_resolved
      return(
        .as_eb_smoother_fit(
          fit = fit,
          family = "constant",
          noise_mode = "fixed",
          fitted_noise_sd = if (length(s) == 1L) as.numeric(s) else NA_real_,
          backend = "exact"
        )
      )
    }

    fit <- .fit_constant_unknown_noise(
      x = x,
      beta_mode = beta_mode,
      beta_fixed = beta_fixed_use,
      beta_prec = beta_prec_use,
      link = link
    )
    fit$g_init <- g_init_resolved
    return(
      .as_eb_smoother_fit(
        fit = fit,
        family = "constant",
        noise_mode = "estimated",
        fitted_noise_sd = fit$fitted_noise_sd,
        backend = "exact"
      )
    )
  }

  fix_params_use <- .normalize_fix_params(
    fix_params = fix_params,
    allowed = c("scale", "beta"),
    fix_g = fix_g,
    fix_g_params = "scale",
    nm = "fix_params"
  )
  if ("scale" %in% fix_params_use && (is.null(g_init) || is.null(g_init$scale))) {
    stop("`fix_params = \"scale\"` requires `g_init$scale`.")
  }
  beta_fixed <- .resolve_fixed_beta_from_fix_params(
    fix_params = fix_params_use,
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    beta_prec_supplied = !beta_prec_missing,
    g_init = g_init,
    expected_length = ncol(setup$X)
  )

  backend_use <- if (backend == "auto" || backend == "exact") {
    if (identical(link, "log")) "laplace_fisher" else if (identical(link, "softplus")) "laplace" else "tmb"
  } else {
    backend
  }
  if (identical(backend_use, "laplace_fisher") && !(identical(link, "log") || identical(link, "softplus"))) {
    stop("`backend = \"laplace_fisher\"` is only available for `link = \"log\"` or `link = \"softplus\"`.")
  }
  if (identical(link, "softplus") && identical(backend_use, "tmb")) {
    stop("`link = \"softplus\"` is not supported with `backend = \"tmb\"`; use `backend = \"laplace\"` or `backend = \"laplace_fisher\"`.")
  }
  if (!(backend_use %in% c("tmb", "laplace", "laplace_fisher"))) {
    stop("Unsupported backend for `family = \"lgp\"`.")
  }

  if (s_known) {
    fit_fun <- ebnm_LGP_generator(
      LGP_setup = setup,
      link = link,
      backend = if (identical(backend_use, "laplace_fisher")) "laplace_fisher" else if (identical(backend_use, "laplace")) "laplace" else "tmb",
      dll = dll
    )
    fit_args <- list(
      x = x,
      s = s_vec,
      g_init = g_init,
      fix_g = fix_g,
      beta_fixed = beta_fixed,
      fix_params = fix_params_use
    )
    if (!beta_prec_missing) {
      fit_args$beta_prec <- beta_prec
    }
    fit <- do.call(fit_fun, fit_args)
    fit$data <- .eb_smoother_data_frame(x, s_vec)

    return(
      .as_eb_smoother_fit(
        fit = fit,
          family = "lgp",
          noise_mode = "fixed",
          fitted_noise_sd = if (length(s) == 1L) as.numeric(s) else NA_real_,
          backend = fit$backend
        )
      )
  }

  fit <- if (identical(backend_use, "laplace_fisher") || identical(backend_use, "laplace") || length(fix_params_use) > 0L) {
    .fit_lgp_laplace_r(
      x = x,
      s = NULL,
      LGP_setup = setup,
      g_init = g_init,
      fix_g = fix_g,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      beta_prec_missing = beta_prec_missing,
      fix_params = fix_params_use,
      link = link,
      learn_noise = TRUE,
      laplace_curvature = if (identical(backend_use, "laplace_fisher")) "fisher" else "observed"
    )
  } else {
    .fit_lgp_unknown_noise(
      x = x,
      LGP_setup = setup,
      g_init = g_init,
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      beta_prec_missing = beta_prec_missing,
      link = link,
      dll = dll
    )
  }

  .as_eb_smoother_fit(
    fit = fit,
    family = "lgp",
    noise_mode = "estimated",
    fitted_noise_sd = fit$fitted_noise_sd,
    backend = fit$backend
  )
}
