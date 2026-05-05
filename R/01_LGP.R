# ---- internal: sparse SPD factorization and derived helpers ----
.factorize_spd <- function(Q, perm = TRUE) {
  if (inherits(Q, "ebsmooth_spd_factor")) {
    return(Q)
  }

  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = TRUE)
  Q <- Matrix::forceSymmetric(Q)
  cholQ <- Matrix::Cholesky(Q, LDL = FALSE, perm = perm)
  R <- if (inherits(cholQ, "CHMfactor")) {
    Matrix::expand(cholQ)$L
  } else {
    as(cholQ, "dtrMatrix")
  }

  structure(
    list(
      matrix = Q,
      chol = cholQ,
      logdet = as.numeric(2 * sum(log(abs(Matrix::diag(R))))),
      perm = perm
    ),
    class = "ebsmooth_spd_factor"
  )
}

.solve_spd_factor <- function(factor, rhs) {
  factor <- .factorize_spd(factor)
  Matrix::solve(factor$chol, rhs, system = "A")
}

.qinv_min_n <- function() {
  min_n <- getOption("EBSmoothr.qinv_min_n", 2000L)
  if (!is.numeric(min_n) || length(min_n) != 1L || is.na(min_n) || min_n < 1) {
    return(2000L)
  }
  as.integer(ceiling(min_n))
}

.qinv_max_row_nnz <- function() {
  max_nnz <- getOption("EBSmoothr.qinv_max_row_nnz", 50L)
  if (!is.numeric(max_nnz) || length(max_nnz) != 1L || is.na(max_nnz) || max_nnz < 1) {
    return(50L)
  }
  as.integer(floor(max_nnz))
}

.matrix_pattern_keys <- function(M) {
  M <- Matrix::drop0(Matrix::Matrix(M, sparse = TRUE))
  M <- as(M, "TsparseMatrix")
  i <- M@i + 1
  j <- M@j + 1
  n <- nrow(M)
  unique(pmin(i, j) + (pmax(i, j) - 1) * n)
}

.A_required_qinv_pattern_keys <- function(A, max_row_nnz = .qinv_max_row_nnz()) {
  A <- Matrix::drop0(Matrix::Matrix(A, sparse = TRUE))
  A <- as(A, "TsparseMatrix")
  rows <- A@i + 1
  cols <- A@j + 1
  if (!length(cols)) {
    return(numeric())
  }

  row_nnz <- tabulate(rows, nbins = nrow(A))
  if (any(row_nnz > max_row_nnz)) {
    return(NULL)
  }

  n <- ncol(A)
  by_row <- split(cols, rows)
  keys <- unlist(lapply(by_row, function(idx) {
    idx <- unique(idx)
    if (length(idx) == 1L) {
      return(idx + (idx - 1) * n)
    }
    grid <- expand.grid(idx, idx)
    pmin(grid[[1]], grid[[2]]) + (pmax(grid[[1]], grid[[2]]) - 1) * n
  }), use.names = FALSE)
  unique(keys)
}

.qinv_pattern_covers <- function(A, qinv) {
  required <- .A_required_qinv_pattern_keys(A)
  if (is.null(required)) {
    return(FALSE)
  }
  if (!length(required)) {
    return(TRUE)
  }
  available <- .matrix_pattern_keys(qinv)
  all(required %in% available)
}

.compute_diag_A_Qinv_At <- function(A, Q) {
  factor <- .factorize_spd(Q)
  A <- Matrix::Matrix(A, sparse = TRUE)

  use_selected_inverse <- ncol(A) >= .qinv_min_n() && requireNamespace("INLA", quietly = TRUE)
  if (isTRUE(use_selected_inverse)) {
    qinv <- tryCatch(INLA::inla.qinv(factor$matrix), error = function(e) NULL)
    if (!is.null(qinv)) {
      qinv <- Matrix::Matrix(qinv, sparse = TRUE)
      if (inherits(A, "diagonalMatrix")) {
        diag_A <- as.numeric(Matrix::diag(A))
        return(diag_A^2 * as.numeric(Matrix::diag(qinv)))
      }
      if (.qinv_pattern_covers(A, qinv)) {
        V <- qinv %*% Matrix::t(A)
        return(as.numeric(Matrix::rowSums(A * Matrix::t(V))))
      }
    }
  }

  V <- .solve_spd_factor(factor, Matrix::t(A))
  as.numeric(Matrix::rowSums(A * Matrix::t(V)))
}

.compute_logdet_spd <- function(Q) {
  .factorize_spd(Q)$logdet
}

# internal function to compute the precision matrix for the L-GP, based on the differences between knots
compute_weights_precision_helper <- function (x) {
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}

# internal function to compute the global polynomial design matrix
global_poly_helper <- function (x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
}

# internal function to compute the local polynomial design matrix, handling negative and positive knots separately
local_poly_helper <- function (knots, refined_x, p = 2, neg_sign_order = 0) {
  if (min(knots) >= 0) {
    D <- get_local_poly(knots, refined_x, p)
  }
  else if (max(knots) <= 0) {
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    D <- get_local_poly(knots_neg, refined_x_neg, p)
    D <- D * ((-1)^neg_sign_order)
  }
  else {
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    D1 <- get_local_poly(knots_neg, refined_x_neg, p)
    D1 <- D1 * ((-1)^neg_sign_order)
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    D2 <- get_local_poly(knots_pos, refined_x_pos, p)
    D <- cbind(D1, D2)
  }
  D
}

# internal function to compute the local polynomial design matrix
get_local_poly <- function (knots, refined_x, p) {
  dif <- diff(knots)
  nn <- length(refined_x)
  n <- length(knots)
  D <- matrix(0, nrow = nn, ncol = n - 1)
  for (j in 1:nn) {
    for (i in 1:(n - 1)) {
      if (refined_x[j] <= knots[i]) {
        D[j, i] <- 0
      }
      else if (refined_x[j] <= knots[i + 1] & refined_x[j] >=
               knots[i]) {
        D[j, i] <- (1/factorial(p)) * (refined_x[j] -
                                         knots[i])^p
      }
      else {
        k <- 1:p
        D[j, i] <- sum((dif[i]^k) * ((refined_x[j] -
                                        knots[i + 1])^(p - k))/(factorial(k) * factorial(p -
                                                                                           k)))
      }
    }
  }
  D
}


#' Define the L-GP Object
#'
#' @description Creates an object from the one-parameter L-GP family.
#'
#' @param scale A numeric value representing the scale parameter of the L-GP.
#' @param beta Optional intercept state stored with the L-GP object.
#' @param beta_prec Optional non-negative prior precision stored with the L-GP
#'   object.
#'
#' @return An object of class `"LGP"`, a data frame containing the scale parameter.
#'
#' @examples
#' gp <- LGP(1.5)
#' print(gp)
#'
#' @export
LGP <- function(scale = 0, beta = NULL, beta_prec = NULL) {
  scale <- .check_single_numeric(scale, "scale")
  beta <- .check_optional_beta_vector(beta, "beta")
  beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")

  structure(
    list(
      scale = scale,
      beta = beta,
      beta_prec = beta_prec
    ),
    class = "LGP"
  )
}



#' Build design/penalty matrices for L-GP smoothing (TMB backend)
#'
#' @description
#' Precomputes the design matrices and penalty components used by the
#' L-GP smoother. The returned list can be passed directly to
#' \code{\link{ebnm_LGP_generator}} for known-standard-error workflows, or
#' reused through \code{\link{eb_smoother}} with
#' \code{family = "lgp"}.
#'
#' The model uses a global polynomial basis \code{X} (degree \code{p-1})
#' and a local spline-like basis \code{B}. The local coefficients are
#' assigned an L-GP prior implemented in TMB through a precision matrix
#' \code{P} and a single scale parameter \code{theta}.
#'
#' @param t Numeric vector of grid locations (length \code{n}).
#' @param p Integer polynomial degree for the global trend basis
#'   (default \code{2}, corresponding to an intercept + linear term).
#' @param num_knots Integer number of knots used to build the local basis
#'   (default \code{30}).
#' @param betaprec Numeric scalar controlling the legacy fallback precision on the
#'   global coefficients \code{beta}:
#'   \itemize{
#'     \item \code{betaprec > 0}: proper Gaussian prior with precision \code{betaprec}.
#'     \item \code{betaprec = 0}: diffuse (improper/flat) prior on \code{beta}.
#'     \item \code{betaprec < 0}: empirical-Bayes mode for \code{beta}
#'       (optimized rather than treated as random in Step A; see
#'       \code{\link{ebnm_LGP_generator}}).
#'   }
#'   This value is retained for backward compatibility. New code should prefer
#'   the public \code{beta_prec} argument on \code{\link{ebnm_LGP_generator}}
#'   and \code{\link{eb_smoother}}.
#' @param link One of \code{"identity"}, \code{"log"}, or
#'   \code{"softplus"}. This value is stored
#'   as \code{link_id} in the returned object for convenience; the final
#'   link used by the fitter is determined by the \code{link} argument in
#'   \code{\link{ebnm_LGP_generator}}.
#'
#' @return A list containing the sparse matrices \code{X}, \code{B}, \code{P},
#'   and additional scalars used by the TMB objective:
#'   \itemize{
#'     \item \code{logPdet}: \eqn{\log|P|} for the penalty/precision component.
#'     \item \code{betaprec}: the value supplied via \code{betaprec}.
#'     \item \code{link_id}: integer encoding of the requested link
#'       (\code{0} identity, \code{1} log, \code{2} softplus).
#'     \item \code{model_id}: internal TMB objective selector; \code{0} for L-GP.
#'   }
#'
#' @export
LGP_setup <- function(t, p = 2, num_knots = 30, betaprec = 0, link = "identity") {
  link <- match.arg(link, choices = c("identity", "log", "softplus"))
  if (!is.numeric(betaprec) || length(betaprec) != 1) stop("betaprec must be a single numeric.")
  t <- as.numeric(t); if (anyNA(t)) stop("t contains NA.")

  knots <- seq(min(t), max(t), length.out = num_knots)
  X <- global_poly_helper(x = t, p = p)
  P <- compute_weights_precision_helper(knots)
  B <- local_poly_helper(knots = knots, refined_x = t, p = p)

  tmbdat <- list(
    X = as(Matrix::Matrix(X, sparse = TRUE), "TsparseMatrix"),
    B = as(Matrix::Matrix(B, sparse = TRUE), "TsparseMatrix"),
    P = as(Matrix::Matrix(P, sparse = TRUE), "TsparseMatrix"),
    logPdet = sum(log(diff(knots))),  # OK for your diag(diff(knots)) P
    betaprec = as.numeric(betaprec),
    link_id = if (identical(link, "identity")) 0L else if (identical(link, "log")) 1L else 2L,
    model_id = 0L
  )
  tmbdat
}

.lgp_observation_terms <- function(eta,
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
      stop("The log-link linear predictor overflowed during L-GP optimization.")
    }
    mu <- exp(eta)
    grad <- mu * (mu - x) / (s^2)
    hess_diag <- if (identical(laplace_curvature, "fisher")) {
      mu^2 / (s^2)
    } else {
      mu * (2 * mu - x) / (s^2)
    }
  } else {
    softplus <- .softplus_stable(eta)
    sigmoid <- .sigmoid_stable(eta)
    mu <- softplus
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
    stop("Non-finite observation terms in L-GP Laplace objective.")
  }
  list(nll = as.numeric(nll), grad = as.numeric(grad), hess_diag = as.numeric(hess_diag), mu = mu)
}

.lgp_response_moments_from_eta <- function(eta_mean, eta_var, link = c("identity", "log", "softplus")) {
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

.lgp_prior_precision <- function(theta, P, pX, beta_mode, beta_prec = NULL) {
  P <- Matrix::Matrix(P, sparse = TRUE)
  Q_u <- Matrix::forceSymmetric(Matrix::Matrix(exp(theta) * P, sparse = TRUE))
  if (!(beta_mode %in% c("prior_flat", "prior_proper"))) {
    return(Q_u)
  }
  beta_prec0 <- if (identical(beta_mode, "prior_proper")) {
    beta_prec <- .check_optional_beta_prec(beta_prec, "beta_prec")
    if (is.null(beta_prec) || beta_prec <= 0) {
      stop("`beta_prec` must be positive for proper-prior beta.")
    }
    beta_prec
  } else {
    0
  }
  Matrix::forceSymmetric(Matrix::bdiag(
    Q_u,
    Matrix::Diagonal(n = pX, x = beta_prec0)
  ))
}

.lgp_laplace_inner_objective <- function(x,
                                         s,
                                         B,
                                         X,
                                         P,
                                         logPdet,
                                         theta,
                                         beta_mode,
                                         beta_value = NULL,
                                         beta_prec = NULL,
                                         beta_init = NULL,
                                         initial_mode = NULL,
                                         link = c("identity", "log", "softplus"),
                                         laplace_curvature = c("observed", "fisher")) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  x <- as.numeric(x)
  s <- as.numeric(s)
  B <- Matrix::Matrix(B, sparse = TRUE)
  X <- Matrix::Matrix(X, sparse = TRUE)
  P <- Matrix::Matrix(P, sparse = TRUE)
  pB <- ncol(B)
  pX <- ncol(X)
  integrate_beta <- beta_mode %in% c("prior_flat", "prior_proper")
  prior_precision <- .lgp_prior_precision(theta, P, pX, beta_mode, beta_prec)
  if (integrate_beta) {
    A_eta <- cbind(B, X)
    default_z0 <- c(rep(0, pB), as.numeric(beta_init))
  } else {
    beta_value <- .check_optional_beta_vector(beta_value, "beta_value", expected_length = pX, allow_null = FALSE)
    A_eta <- B
    default_z0 <- rep(0, pB)
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
    if (!integrate_beta) eta <- eta + as.numeric(X %*% beta_value)
    obs <- .lgp_observation_terms(eta = eta, x = x, s = s, link = link, laplace_curvature = laplace_curvature)
    prior_quad <- as.numeric(0.5 * sum(z * as.numeric(prior_precision %*% z)))
    obs$nll + prior_quad
  }
  gradient <- function(z) {
    z <- as.numeric(z)
    eta <- as.numeric(A_eta %*% z)
    if (!integrate_beta) eta <- eta + as.numeric(X %*% beta_value)
    obs <- .lgp_observation_terms(eta = eta, x = x, s = s, link = link, laplace_curvature = laplace_curvature)
    as.numeric(Matrix::t(A_eta) %*% obs$grad + prior_precision %*% z)
  }
  opt <- stats::nlminb(
    start = z0,
    objective = objective,
    gradient = gradient,
    control = list(eval.max = 1000, iter.max = 1000)
  )
  if (!is.finite(opt$objective)) {
    stop("L-GP Laplace inner optimization failed.")
  }

  z_mode <- as.numeric(opt$par)
  eta_mode <- as.numeric(A_eta %*% z_mode)
  if (!integrate_beta) eta_mode <- eta_mode + as.numeric(X %*% beta_value)
  obs <- .lgp_observation_terms(eta = eta_mode, x = x, s = s, link = link, laplace_curvature = laplace_curvature)
  W_eta <- Matrix::Diagonal(x = obs$hess_diag)
  H <- prior_precision + Matrix::t(A_eta) %*% W_eta %*% A_eta
  H <- Matrix::forceSymmetric(Matrix::Matrix(H, sparse = TRUE))
  H_factor <- .factorize_spd(H)

  U_mode <- z_mode[seq_len(pB)]
  log_prior <- 0.5 * (pB * theta + logPdet) - 0.5 * pB * log(2 * pi) -
    0.5 * exp(theta) * sum(U_mode * as.numeric(P %*% U_mode))
  if (integrate_beta) {
    beta_hat <- z_mode[pB + seq_len(pX)]
    if (identical(beta_mode, "prior_proper")) {
      beta_prec0 <- .check_optional_beta_prec(beta_prec, "beta_prec")
      log_prior <- log_prior + 0.5 * pX * log(beta_prec0) -
        0.5 * pX * log(2 * pi) -
        0.5 * beta_prec0 * sum(beta_hat^2)
    }
    eta_design <- A_eta
  } else {
    beta_hat <- as.numeric(beta_value)
    eta_design <- B
  }
  log_joint <- -obs$nll + log_prior
  log_marginal <- log_joint + 0.5 * length(z_mode) * log(2 * pi) - 0.5 * H_factor$logdet
  eta_var <- .compute_diag_A_Qinv_At(eta_design, H_factor)
  response <- .lgp_response_moments_from_eta(eta_mode, eta_var, link = link)

  list(
    log_marginal = as.numeric(log_marginal),
    log_joint = as.numeric(log_joint),
    mode = z_mode,
    precision = H,
    precision_factor = H_factor,
    fitted_beta = as.numeric(beta_hat),
    U_hat = as.numeric(U_mode),
    eta_mean = eta_mode,
    eta_var = eta_var,
    posterior = data.frame(
      mean = response$mean,
      var = response$var,
      second_moment = response$mean^2 + response$var
    ),
    integrated_beta = integrate_beta,
    eta_design = eta_design,
    inner_convergence = opt$convergence,
    inner_message = opt$message
  )
}

.optimize_lgp_objective <- function(par0, fixed_names, safe_objective, eval_objective, failure_label) {
  fixed_names <- intersect(unique(fixed_names), names(par0))
  free_names <- setdiff(names(par0), fixed_names)
  expand_par <- function(par_free) {
    out <- par0
    out[free_names] <- as.numeric(par_free)
    out
  }
  if (!length(free_names)) {
    return(list(objective = eval_objective(par0), opt = NULL, method = "fixed"))
  }
  methods <- c("BFGS", "Nelder-Mead")
  last_message <- NULL
  for (method in methods) {
    opt <- tryCatch(
      stats::optim(par = par0[free_names], fn = function(z) safe_objective(expand_par(z)), method = method),
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
    objective <- tryCatch(eval_objective(expand_par(opt$par)), error = function(e) e)
    if (inherits(objective, "error")) {
      last_message <- conditionMessage(objective)
      next
    }
    return(list(objective = objective, opt = opt, method = method, free_names = free_names, fixed_names = fixed_names))
  }
  stop(failure_label, ": ", last_message)
}





.fit_lgp_laplace_r <- function(x,
                               s,
                               LGP_setup,
                               g_init = NULL,
                               fix_g = FALSE,
                               beta_fixed = NULL,
                               beta_prec = NULL,
                               beta_prec_missing = FALSE,
                               fix_params = character(),
                               link = c("identity", "log", "softplus"),
                               learn_noise = FALSE,
                               laplace_curvature = c("observed", "fisher")) {
  link <- match.arg(link)
  laplace_curvature <- match.arg(laplace_curvature)
  if (isTRUE(fix_g)) fix_params <- unique(c(fix_params, "scale"))
  x <- as.numeric(x)
  n <- length(x)
  B <- Matrix::Matrix(LGP_setup$B, sparse = TRUE)
  X <- Matrix::Matrix(LGP_setup$X, sparse = TRUE)
  P <- Matrix::Matrix(LGP_setup$P, sparse = TRUE)
  pX <- ncol(X)
  if (is.null(s)) s <- rep(1, n)
  if (length(s) == 1L) s <- rep(s, n)
  if (length(s) != n) stop("`s` must have length 1 or length(x).")
  if (anyNA(x) || anyNA(s)) stop("`x` and `s` must not contain NA.")
  if (any(s <= 0)) stop("All entries of `s` must be > 0.")

  beta_spec <- .eb_smoother_resolve_beta_spec(
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    g_init_beta_prec = if (is.null(g_init)) NULL else g_init$beta_prec,
    legacy_beta_prec = if (isTRUE(beta_prec_missing)) LGP_setup$betaprec else NULL
  )
  beta_mode <- beta_spec$mode
  beta_fixed_use <- if (is.null(beta_spec$beta_fixed)) {
    NULL
  } else {
    .check_optional_beta_vector(beta_spec$beta_fixed, "beta_fixed", expected_length = pX, allow_null = FALSE)
  }
  beta_prec_use <- beta_spec$beta_prec
  g_init_input <- g_init
  if ("scale" %in% fix_params && (is.null(g_init_input) || is.null(g_init_input$scale))) {
    stop("`fix_params = \"scale\"` requires `g_init$scale`.")
  }
  if (is.null(g_init)) {
    g_init <- LGP(scale = 0, beta = rep(0, pX), beta_prec = beta_prec_use)
  }
  theta0 <- .check_single_numeric(g_init$scale, "g_init$scale")
  beta_init <- if (!is.null(beta_fixed_use)) {
    beta_fixed_use
  } else if (!is.null(g_init$beta)) {
    .check_optional_beta_vector(g_init$beta, "g_init$beta", expected_length = pX, allow_null = FALSE)
  } else {
    rep(0, pX)
  }
  noise_sd0 <- if (isTRUE(learn_noise)) stats::sd(x) else NA_real_
  if (isTRUE(learn_noise) && (!is.finite(noise_sd0) || noise_sd0 <= 0)) noise_sd0 <- 1

  beta_names <- paste0("beta", seq_len(pX))
  par0 <- c(theta = theta0)
  if (isTRUE(learn_noise)) par0 <- c(par0, log_noise = log(noise_sd0))
  if (identical(beta_mode, "empirical_bayes")) {
    par0 <- c(par0, stats::setNames(as.numeric(beta_init), beta_names))
  }
  fixed_names <- if ("scale" %in% fix_params) "theta" else character()
  last_inner_mode <- NULL

  eval_objective <- function(par) {
    theta <- as.numeric(par[["theta"]])
    s_use <- if (isTRUE(learn_noise)) rep(exp(as.numeric(par[["log_noise"]])), n) else as.numeric(s)
    beta_value <- if (identical(beta_mode, "fixed")) {
      beta_fixed_use
    } else if (identical(beta_mode, "empirical_bayes")) {
      as.numeric(par[beta_names])
    } else {
      beta_init
    }
    inner <- .lgp_laplace_inner_objective(
      x = x,
      s = s_use,
      B = B,
      X = X,
      P = P,
      logPdet = LGP_setup$logPdet,
      theta = theta,
      beta_mode = beta_mode,
      beta_value = beta_value,
      beta_prec = beta_prec_use,
      beta_init = beta_init,
      initial_mode = last_inner_mode,
      link = link,
      laplace_curvature = laplace_curvature
    )
    if (!is.null(inner$mode) && all(is.finite(inner$mode))) {
      last_inner_mode <<- as.numeric(inner$mode)
    }
    inner$theta <- theta
    inner$s <- s_use
    inner$fitted_noise_sd <- if (isTRUE(learn_noise)) exp(as.numeric(par[["log_noise"]])) else NULL
    inner
  }
  safe_objective <- function(par) {
    objective <- tryCatch(eval_objective(par), error = function(e) e)
    if (inherits(objective, "error")) return(1e100)
    -objective$log_marginal
  }
  opt_res <- .optimize_lgp_objective(
    par0 = par0,
    fixed_names = fixed_names,
    safe_objective = safe_objective,
    eval_objective = eval_objective,
    failure_label = "L-GP Laplace optimization failed"
  )
  objective <- opt_res$objective
  log_likelihood <- structure(objective$log_marginal, class = "logLik")
  log_likelihood_stepB_joint <- structure(objective$log_joint, class = "logLik")
  log_likelihood_stepB_laplace <- log_likelihood
  beta_hat <- as.numeric(objective$fitted_beta)
  posterior_sampler <- function(nsamp) {
    nsamp <- .check_single_numeric(nsamp, "nsamp")
    if (nsamp < 1 || nsamp != floor(nsamp)) stop("`nsamp` must be a positive integer.")
    samps <- LaplacesDemon::rmvnp(n = nsamp, mu = as.numeric(objective$mode), Omega = as.matrix(objective$precision))
    if (is.null(dim(samps))) samps <- matrix(samps, nrow = 1)
    eta_s <- as.matrix(objective$eta_design) %*% t(samps)
    if (!isTRUE(objective$integrated_beta)) {
      eta_s <- eta_s + as.numeric(X %*% beta_hat)
    }
    if (identical(link, "identity")) t(eta_s) else if (identical(link, "log")) t(exp(eta_s)) else t(.softplus_stable(eta_s))
  }
  fitted_g <- LGP(
    scale = objective$theta,
    beta = beta_hat,
    beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use
  )
  out <- list(
    posterior = objective$posterior,
    fitted_g = fitted_g,
    fitted_beta = beta_hat,
    beta_mode = beta_mode,
    beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use,
    log_likelihood = log_likelihood,
    log_likelihood_semantics = paste0(if (identical(laplace_curvature, "fisher")) "laplace_fisher_" else "laplace_", beta_mode),
    log_likelihood_stepA = log_likelihood,
    log_likelihood_stepB_joint = log_likelihood_stepB_joint,
    log_likelihood_stepB_laplace = log_likelihood_stepB_laplace,
    posterior_sampler = posterior_sampler,
    data = data.frame(x = x, s = objective$s),
    prior_family = paste0("LGP", if (isTRUE(learn_noise)) "_learned_noise" else ""),
    backend = if (identical(laplace_curvature, "fisher")) "laplace_fisher" else "laplace",
    laplace_implementation = "r",
    laplace_curvature = laplace_curvature,
    link = link,
    g_init = LGP(theta0, beta = beta_init, beta_prec = if (identical(beta_mode, "prior_flat")) 0 else beta_prec_use)
  )
  if (isTRUE(learn_noise)) {
    out$fitted_noise_sd <- as.numeric(objective$fitted_noise_sd)
  }
  out
}

#' Generate an `ebnm` function for L-GP smoothing
#'
#' @description
#' Returns a function with signature
#' \code{function(x, s, g_init = NULL, fix_g = FALSE, beta_fixed = NULL, beta_prec = NULL, fix_params = character(), output = NULL)}
#' that fits an L-GP smoother.
#'
#' This is the recommended public interface when the observation standard
#' errors \code{s} are known and you want an \code{ebnm}-compatible L-GP
#' smoother, for example inside \code{flashier} or \code{ebmf}. If the
#' standard errors are unknown and you want to learn one common noise SD
#' instead, use \code{\link{eb_smoother}} with \code{family = "lgp"} and
#' \code{s = NULL}.
#'
#' The fitted latent linear predictor is
#' \deqn{\eta = B U + X \beta,}
#' where \code{U} are local-basis coefficients (random effects) and
#' \code{beta} are global coefficients. The single L-GP scale parameter
#' \code{theta} controls the strength of smoothing through the prior/penalty
#' on \code{U}.
#'
#' The argument \code{link} controls how posterior moments are reported:
#' \itemize{
#'   \item \code{"identity"}: posterior moments are for \eqn{\eta}.
#'   \item \code{"log"}: posterior moments are for \eqn{\exp(\eta)} using
#'     a log-normal moment transform.
#'   \item \code{"softplus"}: posterior moments are for
#'     \eqn{\log(1 + \exp(\eta))}. The package evaluates deterministic
#'     Gauss-Hermite quadrature under the marginal Gaussian Laplace
#'     approximation for each \eqn{\eta_i}, so reported moments match
#'     \code{posterior_sampler()} in expectation.
#' }
#'
#' For \code{link = "log"}, \code{backend = "auto"} uses
#' \code{"laplace_fisher"}. For \code{link = "softplus"}, it uses
#' \code{"laplace"}. The Fisher backend keeps the conditional mode equal to the
#' mode of the original log posterior and replaces only the Laplace posterior
#' precision/log-determinant observation curvature by the Fisher/Gauss-Newton
#' curvature. Explicit \code{backend = "laplace"} retains observed-Hessian
#' Laplace semantics.
#'
#' @details
#' \strong{Link handling.}
#' The \code{link} argument here \emph{overrides} \code{LGP_setup$link_id} if present.
#' This is intentional so the caller can reuse a setup object while changing the link.
#'
#' \strong{Step A / Step B estimation logic.}
#' The public \code{beta} semantics are:
#' \itemize{
#'   \item \code{beta_fixed}: treat \code{beta} as known.
#'   \item \code{beta_prec = NULL}: estimate \code{beta} by empirical Bayes.
#'   \item \code{beta_prec = 0}: use a flat zero-mean Gaussian prior on \code{beta}.
#'   \item \code{beta_prec > 0}: use a proper zero-mean Gaussian prior on \code{beta}.
#' }
#'
#' When \code{beta_prec} is not supplied, the fitter falls back to the legacy
#' \code{LGP_setup$betaprec} field so older code continues to work.
#'
#' With these semantics, the internal Step A / Step B logic becomes:
#' \itemize{
#'   \item empirical-Bayes \code{beta}:
#'     \enumerate{
#'       \item Step A: integrate out \code{U} (Laplace) and optimize \code{(theta, beta)}.
#'       \item Step B: fix \code{(theta, beta)} at the Step A estimates and infer \code{U}.
#'     }
#'   \item fixed / flat-prior / proper-prior \code{beta}:
#'     \enumerate{
#'       \item Step A: integrate out \code{U} and \code{beta} and optimize \code{theta}.
#'       \item Step B: fix \code{theta} and infer \code{(U, beta)} jointly.
#'     }
#' }
#'
#' \strong{Initialization and fixed-beta handling.}
#' \code{g_init$beta}, when supplied, is treated as an initialization value for
#' the global coefficients and is stored in the returned fit state. It is not a
#' prior mean. Supplying \code{beta_fixed} fixes \code{beta} at that value
#' regardless of whether \code{fix_g} is \code{TRUE} or \code{FALSE}. When
#' \code{fix_g = TRUE}, \code{theta} is held fixed at \code{g_init$scale}.
#' The returned function also accepts \code{fix_params}. Allowed values are
#' \code{"scale"} and \code{"beta"}. Fixed scale uses
#' \code{g_init$scale}, fixed beta uses \code{beta_fixed} when supplied or
#' otherwise \code{g_init$beta}, and \code{fix_g = TRUE} is retained as a
#' shortcut for \code{fix_params = "scale"}.
#'
#' \strong{Log-likelihood outputs.}
#' The returned object includes multiple log-likelihood-style quantities:
#' \itemize{
#'   \item \code{log_likelihood}: the primary value returned, intended to be
#'     backward-compatible with \pkg{EBMFSmooth}. When Step A is run
#'     (\code{fix_g = FALSE}), this equals the integrated (marginal) objective
#'     from Step A, i.e. \code{-optA$value}.
#'   \item \code{log_likelihood_stepA}: the Step A integrated log-likelihood
#'     (same as \code{log_likelihood} when Step A is run).
#'   \item \code{log_likelihood_stepB_joint}: the joint log density at the Step B
#'     MAP estimate (no Laplace correction).
#'   \item \code{log_likelihood_stepB_laplace}: Laplace approximation around the
#'     Step B MAP. With consistent objectives, this matches \code{log_likelihood}
#'     when Step A is run; it is mainly provided as a diagnostic.
#' }
#'
#' @param LGP_setup A list returned by \code{\link{LGP_setup}}, containing
#'   \code{X}, \code{B}, \code{P}, \code{logPdet}, and \code{betaprec}.
#' @param link One of \code{"identity"}, \code{"log"}, or
#'   \code{"softplus"}. This argument overrides
#'   \code{LGP_setup$link_id} if present.
#' @param backend Backend choice. \code{"auto"} uses \code{"tmb"} for
#'   \code{link = "identity"} and \code{"laplace_fisher"} for
#'   \code{link = "log"}, and \code{"laplace"} for
#'   \code{link = "softplus"}. \code{"laplace"} uses the observed-Hessian
#'   Laplace backend. \code{"laplace_fisher"} is available for
#'   \code{link = "log"} and \code{link = "softplus"}.
#' @param dll Compiled TMB DLL name (default \code{"EBSmoothr"}).
#'
#' @return A function that returns an object of class \code{"ebnm"} (and \code{"list"}).
#'   The closure has signature
#'   \code{function(x, s, g_init = NULL, fix_g = FALSE, beta_fixed = NULL, beta_prec = NULL, fix_params = character(), output = NULL)}.
#'   The returned fit includes posterior mean/variance, fitted hyperparameters,
#'   fitted/fixed \code{beta}, the resolved \code{g_init}, and the log-likelihood
#'   diagnostics described above.
#'
#' @export
ebnm_LGP_generator <- function(LGP_setup,
                               link = c("identity", "log", "softplus"),
                               backend = c("auto", "tmb", "laplace", "laplace_fisher"),
                               dll = "EBSmoothr") {
  link <- match.arg(link)
  backend <- match.arg(backend)
  backend_use <- if (identical(backend, "auto")) {
    if (identical(link, "log")) "laplace_fisher" else if (identical(link, "softplus")) "laplace" else "tmb"
  } else {
    backend
  }
  if (identical(backend_use, "laplace_fisher") && !(identical(link, "log") || identical(link, "softplus"))) {
    stop("`backend = \"laplace_fisher\"` is only available for `link = \"log\"` or `link = \"softplus\"`.")
  }

  .check_numeric_scalar <- function(z, nm) {
    if (!is.numeric(z) || length(z) != 1L || is.na(z)) {
      stop(nm, " must be a single non-NA numeric.")
    }
    as.numeric(z)
  }

  # Force link_id from the function argument (override setup$link_id if present)
  link_id_arg <- if (identical(link, "identity")) 0L else if (identical(link, "log")) 1L else 2L
  if (!is.null(LGP_setup$link_id) && as.integer(LGP_setup$link_id) != link_id_arg) {
    warning("`link` overrides `LGP_setup$link_id` (they differ). Using link = '", link, "'.")
  }

  # Build A matrix columns in the exact order of free parameters in opt$par.
  # This is critical: the Hessian rows/cols correspond to opt$par order.
  .build_A_in_par_order <- function(B, X, par_names) {
    B <- Matrix::Matrix(B, sparse = TRUE)
    X <- Matrix::Matrix(X, sparse = TRUE)

    u_pos <- which(par_names == "U")
    b_pos <- which(par_names == "beta")

    if (length(u_pos) == ncol(B) &&
        length(b_pos) == ncol(X) &&
        (length(b_pos) == 0L || max(u_pos) < min(b_pos))) {
      # Fast path: grouped as all U then all beta
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

  ebnm_gp <- function(x,
                      s,
                      g_init = NULL,
                      fix_g = FALSE,
                      beta_fixed = NULL,
                      beta_prec = NULL,
                      fix_params = character(),
                      output = NULL) {

    if (is.null(LGP_setup$X) || is.null(LGP_setup$B) || is.null(LGP_setup$P)) {
      stop("LGP_setup must contain X, B, and P.")
    }

    n <- length(x)
    if (nrow(LGP_setup$X) != n || nrow(LGP_setup$B) != n) {
      if (length(s) == 3 && length(x) == 3) return(ebnm::ebnm_flat(x))
      stop("length(x) must match nrow(X) and nrow(B) in LGP_setup.")
    }

    if (!(length(s) == 1L || length(s) == n)) stop("s must have length 1 or length(x).")
    if (length(s) == 1L) s <- rep(s, n)
    if (anyNA(x) || anyNA(s)) stop("x and s must not contain NA.")
    if (any(s <= 0)) stop("All entries of s must be > 0.")
    if (!is.logical(fix_g) || length(fix_g) != 1L || is.na(fix_g)) {
      stop("`fix_g` must be TRUE or FALSE.")
    }

    tmbdat <- LGP_setup
    tmbdat$x <- as.numeric(x)
    tmbdat$s <- as.numeric(s)
    tmbdat$link_id <- as.integer(link_id_arg)
    tmbdat$learn_noise <- 0L
    tmbdat$model_id <- 0L

    pB <- ncol(tmbdat$B)
    pX <- ncol(tmbdat$X)
    beta_prec_supplied <- !missing(beta_prec)
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
      beta_prec_supplied = beta_prec_supplied,
      g_init = g_init,
      expected_length = pX
    )

    if (identical(backend_use, "laplace_fisher") || length(fix_params_use) > 0L) {
      fit <- .fit_lgp_laplace_r(
        x = x,
        s = s,
        LGP_setup = LGP_setup,
        g_init = g_init,
        fix_g = fix_g,
        beta_fixed = beta_fixed,
        beta_prec = beta_prec,
        beta_prec_missing = !beta_prec_supplied,
        fix_params = fix_params_use,
        link = link,
        learn_noise = FALSE,
        laplace_curvature = if (identical(backend_use, "laplace_fisher")) "fisher" else "observed"
      )
      return(structure(fit, class = c("list", "ebnm")))
    }

    beta_spec <- .eb_smoother_resolve_beta_spec(
      beta_fixed = beta_fixed,
      beta_prec = beta_prec,
      g_init_beta_prec = if (is.null(g_init)) NULL else g_init$beta_prec,
      legacy_beta_prec = if (missing(beta_prec)) tmbdat$betaprec else NULL
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
    betaprec_internal <- if (beta_mode == "empirical_bayes") -1 else beta_prec_use
    tmbdat$betaprec <- if (is.null(betaprec_internal)) -1 else betaprec_internal

    if (is.null(g_init)) {
      g_init <- LGP(scale = 0, beta = rep(0, pX), beta_prec = beta_prec_use)
    }
    theta0 <- .check_numeric_scalar(g_init$scale, "g_init$scale")
    beta_init <- if (!is.null(beta_fixed_use)) {
      beta_fixed_use
    } else if (!is.null(g_init$beta)) {
      .check_optional_beta_vector(g_init$beta, "g_init$beta", expected_length = pX, allow_null = FALSE)
    } else {
      rep(0, pX)
    }

    par0 <- list(
      theta = theta0,
      U = rep(0, pB),
      beta = as.numeric(beta_init),
      log_noise = 0
    )

    fitted_theta <- theta0
    fitted_beta <- as.numeric(beta_init)
    ll_stepA <- NA_real_

    if (!fix_g) {
      if (beta_mode == "fixed") {
        objA <- TMB::MakeADFun(
          data = tmbdat,
          parameters = within(par0, beta <- as.numeric(beta_fixed_use)),
          map = list(beta = factor(rep(NA, pX)), log_noise = factor(NA)),
          DLL = dll,
          random = "U",
          silent = TRUE
        )
        optA <- optim(par = objA$par, fn = objA$fn, gr = objA$gr, method = "BFGS")
        ll_stepA <- -as.numeric(optA$value)
        fitted_theta <- as.numeric(optA$par[["theta"]])
        fitted_beta <- as.numeric(beta_fixed_use)
      } else if (betaprec_internal < 0) {
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
        fitted_beta <- rep(0, pX)
      }
    }

    fixed_theta_and_beta <- identical(beta_mode, "fixed") || (!fix_g && betaprec_internal < 0)

    if (fixed_theta_and_beta) {
      mapB <- list(
        theta = factor(NA),
        beta = factor(rep(NA, pX)),
        log_noise = factor(NA)
      )
      parB <- list(
        theta = as.numeric(fitted_theta),
        U = rep(0, pB),
        beta = as.numeric(if (beta_mode == "fixed") beta_fixed_use else fitted_beta),
        log_noise = 0
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
      mapB <- list(theta = factor(NA), log_noise = factor(NA))
      parB <- list(
        theta = as.numeric(fitted_theta),
        U = rep(0, pB),
        beta = rep(0, pX),
        log_noise = 0
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
      U_hat <- as.numeric(optB$par[u_idx])
      beta_hat <- as.numeric(optB$par[b_idx])
      eta_mean <- as.numeric(tmbdat$B %*% U_hat + tmbdat$X %*% beta_hat)
      A <- .build_A_in_par_order(tmbdat$B, tmbdat$X, par_names)
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

    log_likelihood <- if (is.finite(ll_stepA)) ll_stepA else as.numeric(ll_stepB_laplace)
    class(log_likelihood) <- "logLik"

    posterior_sampler <- function(nsamp) {
      samps <- LaplacesDemon::rmvnp(n = nsamp, mu = as.numeric(optB$par), Omega = as.matrix(prec))
      if (fixed_theta_and_beta) {
        eta_s <- as.matrix(tmbdat$B) %*% t(samps) + as.matrix(tmbdat$X) %*% beta_hat
      } else {
        par_names <- names(optB$par)
        A <- .build_A_in_par_order(tmbdat$B, tmbdat$X, par_names)
        eta_s <- as.matrix(A) %*% t(samps)
      }
      if (tmbdat$link_id == 0L) t(eta_s) else if (tmbdat$link_id == 1L) t(exp(eta_s)) else t(.softplus_stable(eta_s))
    }

    fitted_g <- LGP(
      scale = fitted_theta,
      beta = beta_hat,
      beta_prec = if (beta_mode == "prior_flat") 0 else beta_prec_use
    )

    structure(
      list(
        posterior = posterior,
        fitted_g = fitted_g,
        fitted_beta = beta_hat,
        beta_mode = beta_mode,
        beta_prec = if (beta_mode == "prior_flat") 0 else beta_prec_use,
        log_likelihood = log_likelihood,
        log_likelihood_semantics = paste0("laplace_", beta_mode),
        log_likelihood_stepA = structure(ll_stepA, class = "logLik"),
        log_likelihood_stepB_joint = ll_stepB_joint,
        log_likelihood_stepB_laplace = ll_stepB_laplace,
        posterior_sampler = posterior_sampler,
        data = data.frame(x = x, s = s),
        prior_family = "LGP",
        backend = if (identical(backend_use, "laplace")) "laplace" else "tmb",
        laplace_implementation = "tmb",
        laplace_curvature = "observed",
        link = link,
        g_init = LGP(theta0, beta = beta_init, beta_prec = if (beta_mode == "prior_flat") 0 else beta_prec_use)
      ),
      class = c("list", "ebnm")
    )
  }

  ebnm_gp
}
