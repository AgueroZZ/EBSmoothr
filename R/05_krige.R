# Posterior of a fitted L-GP or Matern smoother at new locations ("kriging").
#
# Model (identity link): x_i ~ N(f(t_i), s_i^2), f(t) = X(t) beta + g(t), where g is the
# zero-mean Gaussian process of the prior and X the prior's fixed-effect design
# (Matern: intercept; L-GP: global polynomial of degree p - 1). The fitted prior's
# finite representation writes g(t) = A(t) u with u ~ N(0, Q^{-1})
# (Matern: SPDE finite elements; L-GP: local polynomial basis).
#
# Two ways to get f at new locations t*:
#   "fem":        f(t*) = A(t*) u + X(t*) beta, using the posterior of (u, beta).
#   "analytical": the continuous GP the representation approximates
#                 (Matern covariance; IWP-p for L-GP), conditioned on the fitted posterior of
#                 g at the observed locations: g(t*) | g(t_obs) ~ N(C g(t_obs), K** - C K_o*),
#                 C = K_*o K_oo^+. At the observed locations both methods return the fit.

## ---- internal helpers ----------------------------------------------------------

.krige_prior_type <- function(g, what = "`g`") {
  if (inherits(g, "Matern")) return("Matern")
  if (inherits(g, "LGP")) return("LGP")
  stop(what, " must be a fitted EBSmoothr `LGP` or `Matern` prior, not class ",
       paste(class(g), collapse = "/"), ".", call. = FALSE)
}

.krige_is_setup <- function(setup) {
  is.list(setup) && (!is.null(attr(setup, "lgp_design")) ||
                       all(c("mesh", "fem", "A", "locations") %in% names(setup)))
}

.krige_setup_info <- function(setup, type) {
  if (identical(type, "Matern")) {
    need <- c("locations", "d", "mesh", "A", "fem", "alpha")
    if (!is.list(setup) || !all(need %in% names(setup))) {
      stop("For a Matern prior, `setup` must be the object returned by `Matern_setup()`.", call. = FALSE)
    }
    return(list(type = "Matern", loc = as.matrix(setup$locations), d = as.integer(setup$d),
                mesh = setup$mesh, A_obs = as.matrix(setup$A), fem = setup$fem,
                alpha = as.integer(setup$alpha)))
  }
  des <- attr(setup, "lgp_design")
  if (is.null(des)) {
    stop("This L-GP setup carries no design information (locations, knots, p). ",
         "Rebuild it with `LGP_setup()` from EBSmoothr >= 0.3.2.", call. = FALSE)
  }
  list(type = "LGP", loc = matrix(as.numeric(des$t), ncol = 1), d = 1L,
       t = as.numeric(des$t), p = as.integer(des$p), knots = as.numeric(des$knots),
       A_obs = as.matrix(setup$B), X_obs = as.matrix(setup$X), P = as.matrix(setup$P))
}

.krige_as_loc <- function(loc, d, nm = "new_locations") {
  if (is.data.frame(loc)) loc <- as.matrix(loc)
  if (is.null(dim(loc))) loc <- matrix(as.numeric(loc), ncol = 1)
  loc <- as.matrix(loc)
  storage.mode(loc) <- "double"
  if (ncol(loc) != d) stop("`", nm, "` must have ", d, " column(s), matching the setup.", call. = FALSE)
  if (anyNA(loc)) stop("`", nm, "` must not contain NA.", call. = FALSE)
  loc
}

## basis A(t) of the finite representation and fixed-effect design X(t)
.krige_design <- function(info, loc, warn_outside = FALSE) {
  if (identical(info$type, "Matern")) {
    A <- as.matrix(fmesher::fm_basis(info$mesh, loc = if (info$d == 1L) loc[, 1] else loc))
    if (warn_outside) {
      outside <- if (info$d == 1L) {
        loc[, 1] < min(info$mesh$loc) | loc[, 1] > max(info$mesh$loc)
      } else {
        rowSums(A) == 0
      }
      if (any(outside)) {
        warning(sum(outside), " location(s) lie outside the Matern mesh; the finite-element ",
                "representation cannot extrapolate there (1-d: held at the boundary node; ",
                "2-d: prior mean). Use method = \"analytical\" instead.", call. = FALSE)
      }
    }
    return(list(A = A, X = matrix(1, nrow(loc), 1)))
  }
  list(A = as.matrix(local_poly_helper(knots = info$knots, refined_x = loc[, 1], p = info$p)),
       X = as.matrix(global_poly_helper(x = loc[, 1], p = info$p)))
}

## prior precision Q of the representation's coefficients u
.krige_prior_precision <- function(g, info) {
  if (identical(info$type, "Matern")) {
    tk <- .matern_tau_from_range_sigma(range = exp(as.numeric(g$theta)), sigma = as.numeric(g$sigma),
                                       alpha = info$alpha, d = info$d)
    return(as.matrix(.matern_fem_precision(info$fem, kappa = tk$kappa, tau = tk$tau)))
  }
  exp(as.numeric(g$scale)) * info$P
}

.krige_distance <- function(a, b) {
  if (ncol(a) == 1L) return(abs(outer(a[, 1], b[, 1], "-")))
  d2 <- 0
  for (j in seq_len(ncol(a))) d2 <- d2 + outer(a[, j], b[, j], "-")^2
  sqrt(d2)
}

## Matern correlation 2^(1-nu) (z)^nu K_nu(z) / Gamma(nu) at z = kappa * distance >= 0
## (any array; keeps its shape)
.krige_matern_corr <- function(z, nu) {
  out <- z
  out[] <- 1
  pos <- z > 0
  out[pos] <- exp((1 - nu) * log(2) - lgamma(nu) + nu * log(z[pos]) +
                    log(besselK(z[pos], nu, expon.scaled = TRUE)) - z[pos])
  out
}

## Fourier transform of the unit-variance 1-d Matern correlation:
## 2 sqrt(pi) Gamma(nu + 1/2) / Gamma(nu) * kappa^(2 nu) / (kappa^2 + omega^2)^(nu + 1/2)
.krige_matern_spec <- function(omega, nu, kappa) {
  exp(log(2) + 0.5 * log(pi) + lgamma(nu + 0.5) - lgamma(nu) + 2 * nu * log(kappa) -
        (nu + 0.5) * log(kappa^2 + omega^2))
}

.krige_matern_par <- function(g, info) {
  nu <- info$alpha - info$d / 2
  list(nu = nu, kappa = sqrt(8 * nu) / exp(as.numeric(g$theta)), sig2 = as.numeric(g$sigma)^2)
}

## continuous stationary Matern covariance, EBSmoothr parameterization: nu = alpha - d/2,
## kappa = sqrt(8 nu) / range, sigma = marginal SD
.krige_matern_cov <- function(a, b, g, info) {
  mp <- .krige_matern_par(g, info)
  mp$sig2 * .krige_matern_corr(mp$kappa * .krige_distance(a, b), mp$nu)
}

## Unit-variance Matern with Neumann boundaries on an interval of length L, as a function of
## D1 = s - t and D2 = s + t - 2 * lower end (arrays of the same shape). EBSmoothr's 1-d
## finite-element Matern has natural (Neumann) boundary conditions at the ends of the mesh, and
## for fine meshes it converges to this covariance. Method of images:
##   c_N = sum_k [ c(D1 + 2kL) + c(D2 + 2kL) ],
## with ~exp(-2 kappa L k) decay; for very long ranges the equivalent cosine series from Poisson
## summation converges faster:
##   c_N = (1/L) [ S(0) + sum_{m >= 1} S(pi m / L) (cos(pi m D1 / L) + cos(pi m D2 / L)) ].
.krige_neumann_corr <- function(D1, D2, nu, kappa, L, method = c("auto", "images", "cosine")) {
  method <- match.arg(method)
  n_img <- ceiling(50 / (2 * kappa * L)) + 1
  if (identical(method, "images") || (identical(method, "auto") && n_img <= 2000)) {
    acc <- 0
    for (k in -n_img:n_img) {
      acc <- acc + .krige_matern_corr(kappa * abs(D1 + 2 * k * L), nu) +
        .krige_matern_corr(kappa * abs(D2 + 2 * k * L), nu)
    }
    return(acc)
  }
  n_cos <- 4000
  acc <- .krige_matern_spec(0, nu, kappa)
  for (m in seq_len(n_cos)) {
    w <- pi * m / L
    acc <- acc + .krige_matern_spec(w, nu, kappa) * (cos(w * D1) + cos(w * D2))
  }
  acc / L
}

.krige_neumann_domain <- function(info) range(as.numeric(info$mesh$loc))

.krige_matern_neumann_cov <- function(a, b, g, info) {
  mp <- .krige_matern_par(g, info)
  dom <- .krige_neumann_domain(info)
  D1 <- outer(a[, 1], b[, 1], "-")
  D2 <- outer(a[, 1], b[, 1], "+") - 2 * dom[1]
  mp$sig2 * .krige_neumann_corr(D1, D2, mp$nu, mp$kappa, diff(dom))
}

## the IWP coordinate system implied by local_poly_helper(): side (+1/-1) and distance
## from the process origin (knots[1] if all knots >= 0; mirrored if all <= 0; two
## independent one-sided processes from 0 otherwise). Distance <= 0 means g = 0 there.
.krige_iwp_coords <- function(t, knots) {
  if (min(knots) >= 0) return(list(side = rep(1, length(t)), u = t - min(knots)))
  if (max(knots) <= 0) return(list(side = rep(-1, length(t)), u = max(knots) - t))
  list(side = ifelse(t < 0, -1, 1), u = abs(t))
}

## covariance of a p-fold integrated Wiener process started at 0 with unit diffusion:
## sum_k choose(p-1, k) (M - m)^k m^(2p-1-k) / (2p-1-k) / ((p-1)!)^2, m = min, M = max
.krige_iwp_unit <- function(u1, u2, p) {
  m <- outer(u1, u2, pmin)
  M <- outer(u1, u2, pmax)
  m[m < 0] <- 0
  out <- 0
  for (k in 0:(p - 1)) out <- out + choose(p - 1, k) * (M - m)^k * m^(2 * p - 1 - k) / (2 * p - 1 - k)
  out / factorial(p - 1)^2
}

## L-GP: theta is the log precision of the local coefficients, so the IWP diffusion
## variance is exp(-theta) on the scale of t
.krige_lgp_cov <- function(a, b, g, info) {
  ca <- .krige_iwp_coords(a[, 1], info$knots)
  cb <- .krige_iwp_coords(b[, 1], info$knots)
  K <- .krige_iwp_unit(ca$u, cb$u, info$p) * outer(ca$side, cb$side, "==")
  K[ca$u <= 0, ] <- 0
  K[, cb$u <= 0] <- 0
  exp(-as.numeric(g$scale)) * K
}

## prior variance of g at each location (diagonal of .krige_cov(loc, loc))
.krige_cov_diag <- function(loc, g, info) {
  if (identical(info$type, "Matern")) {
    if (!identical(info$boundary, "neumann")) return(rep(as.numeric(g$sigma)^2, nrow(loc)))
    mp <- .krige_matern_par(g, info)
    dom <- .krige_neumann_domain(info)
    return(mp$sig2 * .krige_neumann_corr(rep(0, nrow(loc)), 2 * (loc[, 1] - dom[1]),
                                         mp$nu, mp$kappa, diff(dom)))
  }
  u <- pmax(.krige_iwp_coords(loc[, 1], info$knots)$u, 0)
  p <- info$p
  exp(-as.numeric(g$scale)) * u^(2 * p - 1) / (2 * p - 1) / factorial(p - 1)^2
}

.krige_cov <- function(a, b, g, info) {
  if (!identical(info$type, "Matern")) return(.krige_lgp_cov(a, b, g, info))
  if (identical(info$boundary, "neumann")) .krige_matern_neumann_cov(a, b, g, info) else .krige_matern_cov(a, b, g, info)
}

## boundary treatment of the analytical Matern: "neumann" (default in 1-d, matching the natural
## boundary conditions of EBSmoothr's finite elements at the mesh ends) or "stationary"
.krige_resolve_boundary <- function(boundary, info, loc_new = NULL, warn = TRUE) {
  if (!identical(info$type, "Matern")) {
    info$boundary <- NA_character_
    return(info)
  }
  if (is.null(boundary)) boundary <- if (info$d == 1L) "neumann" else "stationary"
  boundary <- match.arg(boundary, c("neumann", "stationary"))
  if (identical(boundary, "neumann") && info$d != 1L) {
    stop("`boundary = \"neumann\"` is only available for 1-d Matern priors; use \"stationary\".", call. = FALSE)
  }
  if (identical(boundary, "neumann") && warn && !is.null(loc_new)) {
    dom <- .krige_neumann_domain(info)
    outside <- loc_new[, 1] < dom[1] | loc_new[, 1] > dom[2]
    if (any(outside)) {
      warning(sum(outside), " location(s) lie outside the Matern mesh [", signif(dom[1], 4), ", ",
              signif(dom[2], 4), "], where the Neumann-boundary Matern is not defined; its values there ",
              "are the reflection about the nearest end. Use boundary = \"stationary\" to extrapolate.",
              call. = FALSE)
    }
  }
  info$boundary <- boundary
  info
}

.krige_fixed_design_obs <- function(info) {
  if (identical(info$type, "Matern")) matrix(1, nrow(info$A_obs), 1) else info$X_obs
}

## Gaussian posterior of the representation (identity link, exact). beta is fixed at
## g$beta when g$beta_prec is NULL (empirical Bayes / fixed), otherwise integrated
## with a N(0, 1 / beta_prec) prior (beta_prec = 0: flat).
.krige_latent_posterior <- function(x, s, g, info) {
  Ao <- info$A_obs
  Xo <- .krige_fixed_design_obs(info)
  if (length(x) != nrow(Ao)) {
    stop("length(x) = ", length(x), " does not match the ", nrow(Ao),
         " observed locations of the setup.", call. = FALSE)
  }
  Q0 <- .krige_prior_precision(g, info)
  q <- ncol(Ao)
  pX <- ncol(Xo)
  w <- 1 / s^2
  integrated <- !is.null(g$beta_prec)
  if (integrated) {
    D <- cbind(Ao, Xo)
    Qprior <- rbind(cbind(Q0, matrix(0, q, pX)),
                    cbind(matrix(0, pX, q), diag(as.numeric(g$beta_prec), pX)))
    y <- as.numeric(x)
    beta_fixed <- NULL
  } else {
    D <- Ao
    Qprior <- Q0
    beta_fixed <- as.numeric(g$beta)
    if (length(beta_fixed) != pX) stop("`g$beta` must have length ", pX, ".", call. = FALSE)
    y <- as.numeric(x - Xo %*% beta_fixed)
  }
  Qpost <- Qprior + crossprod(D, w * D)
  Qpost <- (Qpost + t(Qpost)) / 2
  R <- chol(Qpost)
  Sigma <- chol2inv(R)
  mean <- as.numeric(Sigma %*% crossprod(D, w * y))
  list(mean = mean, Sigma = Sigma, q = q, pX = pX, integrated = integrated, beta_fixed = beta_fixed)
}

## pseudo-inverse of a symmetric PSD matrix (drops directions with zero prior variance,
## e.g. the L-GP origin, where g = 0 exactly)
.krige_pinv_psd <- function(K, tol = 1e-10) {
  e <- eigen((K + t(K)) / 2, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  V <- e$vectors[, keep, drop = FALSE]
  V %*% (t(V) / e$values[keep])
}

## the linear map from g(t_obs) to E[g(t*) | g(t_obs)], and the conditional covariance
.krige_conditional <- function(loc_new, g, info, full_cov) {
  Koo <- .krige_cov(info$loc, info$loc, g, info)
  Kno <- .krige_cov(loc_new, info$loc, g, info)
  Cmat <- Kno %*% .krige_pinv_psd(Koo)
  if (full_cov) {
    Kc <- .krige_cov(loc_new, loc_new, g, info) - Cmat %*% t(Kno)
    Kc <- (Kc + t(Kc)) / 2
  } else {
    Kc <- pmax(.krige_cov_diag(loc_new, g, info) - rowSums(Cmat * Kno), 0)
  }
  list(C = Cmat, Kc = Kc)
}

.krige_fit_parts <- function(fit) {
  if (is.list(fit) && !is.null(fit$fitted_g) && !is.null(fit$data)) {
    link <- if (is.null(fit$link)) "identity" else fit$link
    if (!identical(link, "identity")) {
      stop("krige_GP() currently supports the identity link only (fit has link = \"", link, "\").",
           call. = FALSE)
    }
    if (!is.null(fit$fitted_noise_sd)) {
      stop("krige_GP() needs known standard errors; fits that learned the noise SD are not supported.",
           call. = FALSE)
    }
    return(list(x = as.numeric(fit$data$x), s = as.numeric(fit$data$s), g = fit$fitted_g))
  }
  if (is.list(fit) && all(c("x", "s", "g") %in% names(fit))) {
    return(list(x = as.numeric(fit$x), s = rep_len(as.numeric(fit$s), length(fit$x)), g = fit$g))
  }
  stop("`fit` must be an EBSmoothr ebnm fit (output of the function returned by ",
       "`ebnm_LGP_generator()` or `ebnm_Matern_generator()`), or a list(x, s, g), e.g. from ",
       "`flash_ebnm_data()`.", call. = FALSE)
}

## ---- exported functions ---------------------------------------------------------

#' Posterior of a fitted L-GP or Matern smoother at new locations
#'
#' @description
#' Given an L-GP or Matern smoother fitted at a set of observed locations, returns the
#' posterior mean, standard deviation and a pointwise credible interval of the smooth
#' function at any other locations. The fitted hyperparameters are kept fixed
#' (empirical Bayes plug-in), and the prior's fixed effects (Matern intercept, L-GP
#' global polynomial) are included: fixed at their fitted value under empirical Bayes,
#' integrated over otherwise.
#'
#' Two methods are available:
#' \describe{
#'   \item{`"analytical"`}{Kriging with the continuous Gaussian process that the fitted
#'     prior approximates: the Matern covariance
#'     \eqn{\sigma^2 2^{1-\nu} (\kappa h)^\nu K_\nu(\kappa h) / \Gamma(\nu)} with
#'     \eqn{\nu = \alpha - d/2} and \eqn{\kappa = \sqrt{8\nu} / \mathrm{range}} (in 1-d,
#'     by default with Neumann boundaries at the mesh ends; see `boundary`), or, for
#'     the L-GP, a \eqn{p}-fold integrated Wiener process with diffusion variance
#'     \eqn{\exp(-\theta)}. The fitted posterior of the process at the observed
#'     locations is propagated through
#'     \eqn{g(t^*) \mid g(t) \sim N(C g(t), K_{**} - C K_{o*})},
#'     \eqn{C = K_{*o} K_{oo}^{+}}. Costs \eqn{O(n^3)} in the number of observed
#'     locations.}
#'   \item{`"fem"`}{Evaluates the prior's own finite representation (Matern finite
#'     elements, L-GP local polynomial basis) at the new locations, using the posterior
#'     of its coefficients. Cheap, but it cannot extrapolate outside the Matern mesh, and
#'     between coarse knots it is only as good as the basis.}
#' }
#' At the observed locations both methods reproduce the fitted posterior.
#'
#' Only the identity link with known standard errors is supported.
#'
#' @param fit An EBSmoothr ebnm fit, i.e. the output of the function returned by
#'   [ebnm_LGP_generator()] or [ebnm_Matern_generator()]; or a list with elements `x`,
#'   `s` and `g` (the data and fitted prior of such a fit), as returned by
#'   [flash_ebnm_data()].
#' @param setup The setup the prior was built from: the output of [LGP_setup()] or
#'   [Matern_setup()] (for the L-GP, created with EBSmoothr >= 0.3.2).
#' @param new_locations Locations at which to evaluate the posterior: a numeric vector
#'   (1-d) or a matrix with one column per coordinate.
#' @param method `"analytical"` (default) or `"fem"`; see Description.
#' @param level Coverage of the pointwise credible interval.
#' @param return_cov If `TRUE`, the full posterior covariance at `new_locations` is
#'   attached as attribute `"cov"`.
#' @param boundary Matern only, `method = "analytical"`: which continuous Matern to krige
#'   with. `"neumann"` (the default for 1-d priors) uses the Matern with Neumann boundary
#'   conditions at the ends of the mesh, computed by the method of images. This is the process
#'   EBSmoothr's 1-d finite elements approximate (their natural boundary conditions), so it
#'   matches the fitted prior near the ends. `"stationary"` uses the stationary Matern on the
#'   whole line (the only option, and the default, for 2-d priors). The two differ within
#'   about one range of the ends of the mesh.
#'
#' @return A data frame with the location(s), `mean`, `sd`, `lower` and `upper`.
#'
#' @seealso [krige_flash()] for posterior means of all factors of a flashier fit.
#' @export
krige_GP <- function(fit, setup, new_locations, method = c("analytical", "fem"),
                     level = 0.95, return_cov = FALSE, boundary = NULL) {
  method <- match.arg(method)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  parts <- .krige_fit_parts(fit)
  type <- .krige_prior_type(parts$g, "The fitted prior")
  info <- .krige_setup_info(setup, type)
  loc_new <- .krige_as_loc(new_locations, info$d)
  info <- .krige_resolve_boundary(boundary, info, loc_new, warn = identical(method, "analytical"))
  post <- .krige_latent_posterior(parts$x, parts$s, parts$g, info)
  des <- .krige_design(info, loc_new, warn_outside = identical(method, "fem"))

  if (identical(method, "fem")) {
    Mg <- des$A
    Kc <- if (return_cov) matrix(0, nrow(loc_new), nrow(loc_new)) else rep(0, nrow(loc_new))
  } else {
    cond <- .krige_conditional(loc_new, parts$g, info, full_cov = return_cov)
    Mg <- cond$C %*% info$A_obs       # g(t*) = C A_obs u + independent conditional noise
    Kc <- cond$Kc
  }
  if (post$integrated) {
    M <- cbind(Mg, des$X)
    mean <- as.numeric(M %*% post$mean)
  } else {
    M <- Mg
    mean <- as.numeric(M %*% post$mean + des$X %*% post$beta_fixed)
  }
  if (return_cov) {
    cov <- M %*% post$Sigma %*% t(M) + Kc
    cov <- (cov + t(cov)) / 2
    v <- diag(cov)
  } else {
    v <- rowSums((M %*% post$Sigma) * M) + Kc
  }
  sd <- sqrt(pmax(v, 0))
  zq <- stats::qnorm(1 - (1 - level) / 2)
  loc_df <- if (info$d == 1L) data.frame(location = loc_new[, 1]) else {
    stats::setNames(as.data.frame(loc_new), paste0("location", seq_len(info$d)))
  }
  out <- cbind(loc_df, mean = mean, sd = sd, lower = mean - zq * sd, upper = mean + zq * sd)
  attr(out, "method") <- method
  attr(out, "boundary") <- if (identical(method, "analytical")) info$boundary else NA_character_
  attr(out, "level") <- level
  if (return_cov) attr(out, "cov") <- cov
  out
}

#' Posterior means of flashier loadings (or factors) on a new grid
#'
#' @description
#' A lightweight companion to [krige_GP()] for flashier fits whose loadings (or factors)
#' have an EBSmoothr L-GP or Matern prior. It needs only the flash object, the setup and
#' the new locations, and returns the posterior mean of every factor there. It works
#' from the posterior means and fitted priors that flashier stores, so no data are
#' reconstructed; for standard deviations and intervals use [flash_ebnm_data()] with
#' [krige_GP()].
#'
#' For `method = "analytical"` the mean is
#' \eqn{X(t^*)\beta + C\,(m - X(t_{obs})\beta)} with \eqn{m} the posterior means at the
#' observed locations and \eqn{C} the kriging map of the fitted continuous prior (see
#' [krige_GP()]). For `method = "fem"` the coefficients of the finite representation
#' are recovered from \eqn{m}, which requires that they are determined by it (e.g. a
#' 1-d Matern mesh whose nodes are the observed locations, or L-GP knots at a subset
#' of them).
#'
#' @param fl A flash object (from `flashier`).
#' @param setup The [LGP_setup()] or [Matern_setup()] object used to build the prior of
#'   dimension `n`; or a list of such setups, one per factor.
#' @param new_locations Numeric vector (1-d) or matrix of new locations.
#' @param n Which dimension carries the smoothing prior: `1` for loadings (rows of the
#'   data matrix, `fl$L_pm`), `2` for factors (`fl$F_pm`).
#' @param method `"analytical"` (default) or `"fem"`.
#' @param normalize `"none"` (default) returns the curves on flashier's scale. `"ldf"`
#'   divides each curve by the Euclidean norm of the factor's posterior means at the
#'   observed locations, so that it matches `flashier::ldf(fl, type = "2")` there.
#' @param kset Indices of the factors to evaluate (default: all).
#' @param boundary Matern only: boundary treatment for `method = "analytical"`, as in
#'   [krige_GP()] (default `"neumann"` for 1-d priors).
#'
#' @return A matrix with one row per new location and one column per factor in `kset`.
#'   Factors that are identically zero give a column of zeros.
#' @export
krige_flash <- function(fl, setup, new_locations, n = 1, method = c("analytical", "fem"),
                        normalize = c("none", "ldf"), kset = NULL, boundary = NULL) {
  method <- match.arg(method)
  normalize <- match.arg(normalize)
  if (!n %in% c(1, 2)) stop("`n` must be 1 or 2.", call. = FALSE)
  pm <- if (n == 1) fl$L_pm else fl$F_pm
  ghat <- if (n == 1) fl$L_ghat else fl$F_ghat
  if (is.null(pm) || is.null(ghat)) stop("`fl` does not look like a flash object.", call. = FALSE)
  pm <- as.matrix(pm)
  if (is.null(kset)) kset <- seq_len(ncol(pm))
  setups <- if (.krige_is_setup(setup)) rep(list(setup), ncol(pm)) else setup
  if (length(setups) != ncol(pm)) {
    stop("`setup` must be one setup or a list with one setup per factor.", call. = FALSE)
  }

  out <- NULL
  loc_new <- NULL
  for (j in seq_along(kset)) {
    k <- kset[j]
    m <- pm[, k]
    if (is.null(out)) {
      d0 <- if (!is.null(attr(setups[[k]], "lgp_design"))) 1L else as.integer(setups[[k]]$d)
      loc_new <- .krige_as_loc(new_locations, d0)
      out <- matrix(0, nrow(loc_new), length(kset), dimnames = list(NULL, paste0("factor", kset)))
    }
    if (all(m == 0)) next
    g <- ghat[[k]]
    type <- .krige_prior_type(g, sprintf("The prior of factor %d", k))
    info <- .krige_setup_info(setups[[k]], type)
    info <- .krige_resolve_boundary(boundary, info, loc_new,
                                    warn = identical(method, "analytical") && j == 1L)
    if (length(m) != nrow(info$loc)) {
      stop(sprintf("Factor %d has %d posterior means but the setup has %d locations.",
                   k, length(m), nrow(info$loc)), call. = FALSE)
    }
    beta <- as.numeric(g$beta)
    Xo <- .krige_fixed_design_obs(info)
    des <- .krige_design(info, loc_new, warn_outside = identical(method, "fem") && j == 1L)
    resid <- as.numeric(m - Xo %*% beta)
    if (identical(method, "analytical")) {
      Kno <- .krige_cov(loc_new, info$loc, g, info)
      Koo <- .krige_cov(info$loc, info$loc, g, info)
      curve <- as.numeric(des$X %*% beta + Kno %*% (.krige_pinv_psd(Koo) %*% resid))
    } else {
      u <- qr.solve(info$A_obs, resid)
      if (max(abs(info$A_obs %*% u - resid)) > 1e-8 * max(1, max(abs(m)))) {
        stop(sprintf(paste0("Factor %d: the finite representation is not determined by the ",
                            "posterior means at the observed locations; use method = ",
                            "\"analytical\", or krige_GP() with flash_ebnm_data()."), k), call. = FALSE)
      }
      curve <- as.numeric(des$A %*% u + des$X %*% beta)
    }
    if (identical(normalize, "ldf")) curve <- curve / sqrt(sum(m^2))
    out[, j] <- curve
  }
  attr(out, "locations") <- loc_new
  out
}

#' Reconstruct the ebnm problem of one flashier factor
#'
#' @description
#' Returns the data `x`, standard errors `s` and fitted prior `g` that flashier's last
#' update of factor `k` passes to the prior of dimension `n`, recomputed from the final
#' state of the fit. Pass the result to [krige_GP()] to get posterior standard
#' deviations and intervals on a new grid.
#'
#' The reconstruction uses the final state of the fit, so it matches the data of the
#' factor's last update only once the fit has converged: after `flash_backfit()` (or
#' `flash(..., backfit = TRUE)`). After greedy fitting alone, the early factors were fit
#' before the later ones existed, and their reconstructed data will differ.
#'
#' Only dense data matrices are supported (constant, row-wise, column-wise or full
#' residual precisions).
#'
#' @param fl A flash object (from `flashier`).
#' @param k Index of the factor.
#' @param n Dimension carrying the smoothing prior: `1` for loadings, `2` for factors.
#'
#' @return A list with elements `x`, `s` and `g`.
#' @export
flash_ebnm_data <- function(fl, k, n = 1) {
  ff <- fl$flash_fit
  if (is.null(ff) || is.null(ff$EF) || is.null(ff$Y)) stop("`fl` does not look like a flash object.", call. = FALSE)
  if (!n %in% c(1, 2)) stop("`n` must be 1 or 2.", call. = FALSE)
  Y <- as.matrix(ff$Y)
  EF <- ff$EF
  EF2 <- ff$EF2
  if (k < 1 || k > ncol(EF[[1]])) stop("`k` is out of range.", call. = FALSE)
  n1 <- nrow(Y); n2 <- ncol(Y)
  tau <- ff$tau
  Tm <- if (length(tau) == 1L) {
    matrix(tau, n1, n2)
  } else if (is.matrix(tau)) {
    tau
  } else if (length(tau) == n1 && identical(as.integer(ff$est.tau.dim), 1L)) {
    matrix(tau, n1, n2)
  } else if (length(tau) == n2 && identical(as.integer(ff$est.tau.dim), 2L)) {
    matrix(tau, n1, n2, byrow = TRUE)
  } else {
    stop("Unsupported residual precision structure in the flash object.", call. = FALSE)
  }
  Z <- ff$Z
  W <- if (is.null(Z) || length(Z) == 1L) Tm else Tm * as.matrix(Z)
  W[is.na(Y)] <- 0
  others <- setdiff(seq_len(ncol(EF[[1]])), k)
  R <- Y
  if (length(others)) R <- Y - EF[[1]][, others, drop = FALSE] %*% t(EF[[2]][, others, drop = FALSE])
  R[is.na(R)] <- 0
  if (n == 1) {
    prec <- as.numeric(W %*% EF2[[2]][, k])
    x <- as.numeric((W * R) %*% EF[[2]][, k]) / prec
    g <- fl$L_ghat[[k]]
  } else {
    prec <- as.numeric(crossprod(W, EF2[[1]][, k]))
    x <- as.numeric(crossprod(W * R, EF[[1]][, k])) / prec
    g <- fl$F_ghat[[k]]
  }
  list(x = x, s = 1 / sqrt(prec), g = g)
}
