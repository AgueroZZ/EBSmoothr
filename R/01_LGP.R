# ---- internal: diag(A Q^{-1} A^T) for SPD precision Q ----
.compute_diag_A_Qinv_At <- function(A, Q) {
  Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))
  A <- Matrix::Matrix(A, sparse = TRUE)

  cholQ <- Matrix::Cholesky(Q, LDL = FALSE, perm = FALSE)
  V <- Matrix::solve(cholQ, Matrix::t(A), system = "A")
  as.numeric(Matrix::rowSums(A * Matrix::t(V)))
}

# # Old-package-compatible logdet: determinant(as.matrix(H))$modulus
.compute_logdet_spd <- function(Q) {
  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = TRUE)
  Q <- Matrix::forceSymmetric(Q)
  F <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE)  # perm=TRUE 通常更稳
  as.numeric(Matrix::determinant(F, logarithm = TRUE)$modulus)
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
#'
#' @return An object of class `"LGP"`, a data frame containing the scale parameter.
#'
#' @examples
#' gp <- LGP(1.5)
#' print(gp)
#'
#' @export
LGP <- function (scale) {
  structure(data.frame(scale), class = "LGP")
}



#' Build design/penalty matrices for L-GP smoothing (TMB backend)
#'
#' @description
#' Precomputes the design matrices and penalty components used by the
#' L-GP smoother. The returned list is intended to be passed to
#' \code{\link{ebnm_LGP_generator}}.
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
#' @param betaprec Numeric scalar controlling the prior precision on the
#'   global coefficients \code{beta}:
#'   \itemize{
#'     \item \code{betaprec > 0}: proper Gaussian prior with precision \code{betaprec}.
#'     \item \code{betaprec = 0}: diffuse (improper/flat) prior on \code{beta}.
#'     \item \code{betaprec < 0}: empirical-Bayes mode for \code{beta}
#'       (optimized rather than treated as random in Step A; see
#'       \code{\link{ebnm_LGP_generator}}).
#'   }
#' @param link Either \code{"identity"} or \code{"log"}. This value is stored
#'   as \code{link_id} in the returned object for convenience; the final
#'   link used by the fitter is determined by the \code{link} argument in
#'   \code{\link{ebnm_LGP_generator}}.
#'
#' @return A list containing the sparse matrices \code{X}, \code{B}, \code{P},
#'   and additional scalars used by the TMB objective:
#'   \itemize{
#'     \item \code{logPdet}: \eqn{\log|P|} for the penalty/precision component.
#'     \item \code{betaprec}: the value supplied via \code{betaprec}.
#'     \item \code{link_id}: integer encoding of the requested link (\code{0} identity, \code{1} log).
#'   }
#'
#' @export
LGP_setup <- function(t, p = 2, num_knots = 30, betaprec = 0, link = "identity") {
  if (!is.numeric(betaprec) || length(betaprec) != 1) stop("betaprec must be a single numeric.")
  t <- as.numeric(t); if (anyNA(t)) stop("t contains NA.")

  knots <- seq(min(t), max(t), length.out = num_knots)
  X <- global_poly_helper(x = t, p = p)
  P <- compute_weights_precision_helper(knots)
  B <- local_poly_helper(knots = knots, refined_x = t, p = p)

  tmbdat <- list(
    X = as(as(as(X, "dMatrix"), "generalMatrix"), "TsparseMatrix"),
    B = as(as(as(B, "dMatrix"), "generalMatrix"), "TsparseMatrix"),
    P = as(as(as(P, "dMatrix"), "generalMatrix"), "TsparseMatrix"),
    logPdet = sum(log(diff(knots))),  # OK for your diag(diff(knots)) P
    betaprec = as.numeric(betaprec),
    link_id = ifelse(link == "log", 1L, 0L)
  )
  tmbdat
}





#' Generate an `ebnm` function for L-GP smoothing (TMB backend)
#'
#' @description
#' Returns a function with signature
#' \code{function(x, s, g_init = NULL, fix_g = FALSE, output = NULL)}
#' that fits an L-GP smoother using a unified TMB objective.
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
#' }
#'
#' @details
#' \strong{Link handling.}
#' The \code{link} argument here \emph{overrides} \code{LGP_setup$link_id} if present.
#' This is intentional so the caller can reuse a setup object while changing the link.
#'
#' \strong{Step A / Step B estimation logic.}
#' The behavior is controlled by \code{LGP_setup$betaprec}:
#' \itemize{
#'   \item \code{betaprec < 0} (EB mode for \code{beta}):
#'     \enumerate{
#'       \item Step A: integrate out \code{U} (Laplace) and optimize \code{(theta, beta)}.
#'       \item Step B: fix \code{(theta, beta)} at the Step A estimates and infer \code{U}.
#'     }
#'   \item \code{betaprec >= 0}:
#'     \enumerate{
#'       \item Step A: integrate out \code{U} and \code{beta} and optimize \code{theta}.
#'       \item Step B: fix \code{theta} and infer \code{(U, beta)} jointly.
#'     }
#' }
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
#' @param link Either \code{"identity"} or \code{"log"}. This argument overrides
#'   \code{LGP_setup$link_id} if present.
#' @param dll Compiled TMB DLL name (default \code{"EBSmoothr"}).
#'
#' @return A function that returns an object of class \code{"ebnm"} (and \code{"list"}).
#'   The returned object includes posterior mean/variance, fitted hyperparameters,
#'   and the log-likelihood diagnostics described above.
#'
#' @export
ebnm_LGP_generator <- function(LGP_setup,
                               link = c("identity", "log"),
                               dll = "EBSmoothr") {
  link <- match.arg(link)

  .check_numeric_scalar <- function(z, nm) {
    if (!is.numeric(z) || length(z) != 1L || is.na(z)) {
      stop(nm, " must be a single non-NA numeric.")
    }
    as.numeric(z)
  }

  # Force link_id from the function argument (override setup$link_id if present)
  link_id_arg <- if (link == "log") 1L else 0L
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

  ebnm_gp <- function(x, s, g_init = NULL, fix_g = FALSE, output = NULL) {

    # ---- basic checks ----
    if (is.null(LGP_setup$X) || is.null(LGP_setup$B) || is.null(LGP_setup$P)) {
      stop("LGP_setup must contain X, B, and P.")
    }

    n <- length(x)
    if (nrow(LGP_setup$X) != n || nrow(LGP_setup$B) != n) {
      # ebnm sometimes probes with length-3 init
      if (length(s) == 3 && length(x) == 3) return(ebnm::ebnm_flat(x))
      stop("length(x) must match nrow(X) and nrow(B) in LGP_setup.")
    }

    if (!(length(s) == 1L || length(s) == n)) stop("s must have length 1 or length(x).")
    if (length(s) == 1L) s <- rep(s, n)
    if (anyNA(x) || anyNA(s)) stop("x and s must not contain NA.")
    if (any(s <= 0)) stop("All entries of s must be > 0.")

    # ---- TMB data ----
    tmbdat <- LGP_setup
    tmbdat$x <- as.numeric(x)
    tmbdat$s <- as.numeric(s)
    tmbdat$link_id <- as.integer(link_id_arg)

    betaprec <- .check_numeric_scalar(tmbdat$betaprec, "betaprec")

    pB <- ncol(tmbdat$B)
    pX <- ncol(tmbdat$X)

    if (is.null(g_init)) g_init <- LGP(0)
    theta0 <- .check_numeric_scalar(g_init$scale, "g_init$scale")

    # Parameter template required by the unified C++ objective
    par0 <- list(
      theta = theta0,
      U    = rep(0, pB),
      beta = rep(0, pX)
    )

    # ============================================================
    # Step A
    #   betaprec >= 0: integrate out (U, beta), optimize theta only
    #   betaprec <  0: integrate out U, optimize (theta, beta)
    # Primary logLik is Step A integrated objective: -optA$value
    # (matches EBMFSmooth's returned log_likelihood).
    # ============================================================
    fitted_theta <- theta0
    fitted_beta  <- rep(0, pX)
    ll_stepA <- NA_real_

    if (!fix_g) {
      if (betaprec < 0) {
        # Integrate out U, jointly optimize theta and beta
        objA <- TMB::MakeADFun(
          data = tmbdat,
          parameters = par0,
          DLL = dll,
          random = "U",
          silent = TRUE
        )
        optA <- optim(par = objA$par, fn = objA$fn, gr = objA$gr, method = "BFGS")

        ll_stepA <- -as.numeric(optA$value)

        if (!("theta" %in% names(optA$par))) stop("Step A: 'theta' not found in optA$par.")
        fitted_theta <- as.numeric(optA$par[["theta"]])

        # IMPORTANT: beta appears with repeated names ("beta","beta",...)
        beta_idx <- which(names(optA$par) == "beta")
        if (length(beta_idx) != pX) {
          stop(sprintf("Step A: expected %d beta entries, got %d.", pX, length(beta_idx)))
        }
        fitted_beta <- as.numeric(optA$par[beta_idx])

      } else {
        # Integrate out (U, beta), optimize theta only
        objA <- TMB::MakeADFun(
          data = tmbdat,
          parameters = par0,
          DLL = dll,
          random = c("U", "beta"),
          silent = TRUE
        )
        optA <- optim(par = objA$par, fn = objA$fn, gr = objA$gr, method = "BFGS")

        ll_stepA <- -as.numeric(optA$value)

        if (!("theta" %in% names(optA$par))) stop("Step A: 'theta' not found in optA$par.")
        fitted_theta <- as.numeric(optA$par[["theta"]])

        # beta is integrated out here; we will estimate it in Step B.
        fitted_beta <- rep(0, pX)
      }
    } else {
      # If g is fixed, we do not optimize in Step A.
      fitted_theta <- theta0
      fitted_beta  <- rep(0, pX)
    }

    fitted_g <- LGP(fitted_theta)

    # ============================================================
    # Step B
    #   betaprec >= 0: fix theta; infer (U, beta) jointly
    #   betaprec <  0: fix theta and beta; infer U only
    # We compute Step B joint and Laplace-corrected values for diagnostics.
    # IMPORTANT: StepB_laplace uses OLD-style logdet to match EBMFSmooth.
    # ============================================================
    if (betaprec < 0) {
      # Fix theta AND beta, infer U only
      mapB <- list(
        theta = factor(NA),
        beta  = factor(rep(NA, pX))
      )
      parB <- list(
        theta = as.numeric(fitted_theta),
        U    = rep(0, pB),
        beta = as.numeric(fitted_beta)
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

      if (length(optB$par) != pB) {
        stop(sprintf("Step B: expected %d free params (U only), got %d.", pB, length(optB$par)))
      }

      H <- numDeriv::hessian(function(w) ff$fn(w), optB$par)
      prec <- Matrix::forceSymmetric(H)

      U_hat <- as.numeric(optB$par)
      beta_hat <- as.numeric(fitted_beta)

      eta_mean <- as.numeric(tmbdat$B %*% U_hat + tmbdat$X %*% beta_hat)
      eta_var  <- .compute_diag_A_Qinv_At(tmbdat$B, prec)

      free_dim <- pB

    } else {
      # Fix theta only, infer (U, beta) jointly
      mapB <- list(theta = factor(NA))
      parB <- list(
        theta = as.numeric(fitted_theta),
        U    = rep(0, pB),
        beta = rep(0, pX)
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

      if (length(optB$par) != (pB + pX)) {
        stop(sprintf("Step B: expected %d free params (U+beta), got %d.",
                     pB + pX, length(optB$par)))
      }

      H <- numDeriv::hessian(function(w) ff$fn(w), optB$par)
      prec <- Matrix::forceSymmetric(H)

      par_names <- names(optB$par)
      u_idx <- which(par_names == "U")
      b_idx <- which(par_names == "beta")
      if (length(u_idx) != pB || length(b_idx) != pX) {
        stop("Step B: unexpected names/lengths in optB$par (expect pB 'U' and pX 'beta').")
      }

      U_hat <- as.numeric(optB$par[u_idx])
      beta_hat <- as.numeric(optB$par[b_idx])

      eta_mean <- as.numeric(tmbdat$B %*% U_hat + tmbdat$X %*% beta_hat)

      A <- .build_A_in_par_order(tmbdat$B, tmbdat$X, par_names)
      eta_var <- .compute_diag_A_Qinv_At(A, prec)

      free_dim <- pB + pX
    }

    # ---- transform moments according to link ----
    if (tmbdat$link_id == 0L) {
      post_mean <- eta_mean
      post_var  <- eta_var
    } else {
      post_mean <- exp(eta_mean + 0.5 * eta_var)
      post_var  <- exp(2 * eta_mean + eta_var) * (exp(eta_var) - 1)
    }

    posterior <- data.frame(mean = as.numeric(post_mean), var = as.numeric(post_var))
    posterior$second_moment <- posterior$mean^2 + posterior$var

    # ---- Step B diagnostics: joint and Laplace-corrected values ----
    ll_stepB_joint <- -ff$fn(optB$par)

    # OLD-style logdet here (matches EBMFSmooth behavior)
    ll_stepB_laplace <- ll_stepB_joint - 0.5 * .compute_logdet_spd(prec) +
      0.5 * free_dim * log(2 * pi)

    class(ll_stepB_joint) <- "logLik"
    class(ll_stepB_laplace) <- "logLik"

    # ---- Primary logLik returned ----
    # Match old package: return Step A integrated objective (-optA$value) when available.
    # If fix_g = TRUE (no Step A), fall back to Step B laplace.
    if (is.finite(ll_stepA)) {
      log_likelihood <- ll_stepA
    } else {
      log_likelihood <- as.numeric(ll_stepB_laplace)
    }
    class(log_likelihood) <- "logLik"

    posterior_sampler <- function(nsamp) {
      samps <- LaplacesDemon::rmvnp(n = nsamp, mu = as.numeric(optB$par), Omega = as.matrix(prec))
      if (betaprec < 0) {
        # sample U only
        U_s <- samps
        eta_s <- as.matrix(tmbdat$B) %*% t(U_s) + as.matrix(tmbdat$X) %*% beta_hat
      } else {
        # sample (U, beta) in the same order as optB$par
        par_names <- names(optB$par)
        A <- .build_A_in_par_order(tmbdat$B, tmbdat$X, par_names)
        eta_s <- as.matrix(A) %*% t(samps)
      }
      if (tmbdat$link_id == 0L) t(eta_s) else t(exp(eta_s))
    }

    result <- list(
      posterior = posterior,
      fitted_g = fitted_g,
      fitted_beta = beta_hat,
      log_likelihood = log_likelihood,                  # aligned to old package (Step A)
      log_likelihood_stepA = structure(ll_stepA, class = "logLik"),
      log_likelihood_stepB_joint = ll_stepB_joint,      # diagnostic
      log_likelihood_stepB_laplace = ll_stepB_laplace,  # diagnostic (old-style logdet)
      posterior_sampler = posterior_sampler,
      data = data.frame(x = x, s = s),
      prior_family = "LGP"
    )

    structure(result, class = c("list", "ebnm"))
  }

  ebnm_gp
}

