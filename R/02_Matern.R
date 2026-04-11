# ---- helper: normalize locations and build mesh/A for d = 1 or 2 ----
.normalize_locations <- function(locations) {
  if (is.data.frame(locations)) locations <- as.matrix(locations)

  if (is.null(dim(locations))) {
    # vector -> n x 1
    locations <- matrix(as.numeric(locations), ncol = 1)
  } else {
    locations <- as.matrix(locations)
    storage.mode(locations) <- "double"
  }

  if (!ncol(locations) %in% c(1L, 2L)) {
    stop("For Matern, `locations` must be a vector, an n×1 matrix, or an n×2 matrix.")
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
    return(list(mesh = mesh, A = A))
  }

  # d == 2
  if (is.null(max.edge)) {
    ranges <- apply(loc_mat, 2, function(z) diff(range(z, finite = TRUE)))
    max.edge <- min(ranges) / 10
    if (!is.finite(max.edge) || max.edge <= 0) max.edge <- 1
  }
  if (length(max.edge) == 1L) max.edge <- c(max.edge, max.edge)
  if (length(max.edge) != 2L) stop("For d=2, max.edge must be NULL, length-1, or length-2.")

  mesh <- INLA::inla.mesh.2d(loc = loc_mat, max.edge = max.edge)
  A <- INLA::inla.spde.make.A(mesh = mesh, loc = loc_mat)
  list(mesh = mesh, A = A)
}



#' Define the Matern GP Family
#'
#' @description Creates an object representing a Matern process, parameterized
#' by \code{theta}, the log-transformed range parameter.
#'
#' @param theta Numeric scalar. Log-transformed range parameter. Defaults to \code{NULL}.
#'
#' @return An object of class \code{"Matern"}.
#'
#' @export
Matern <- function(theta = NULL) {
  structure(list(theta = theta), class = "Matern")
}



#' Generate an `ebnm` Function for Matern Smoothing (INLA backend)
#'
#' @description
#' Returns an \code{ebnm}-compatible function that fits a Matern GP using INLA.
#' Supports \strong{d = 1} and \strong{d = 2} spatial domains.
#'
#' - If \code{locations} is a numeric vector, it is treated as 1D locations (n×1).
#' - If \code{locations} is a matrix/data.frame, it must be n×1 or n×2.
#'
#' @param locations Spatial locations (vector, n×1, or n×2).
#' @param max.edge Mesh maximum edge length.
#'   - d=1: scalar.
#'   - d=2: scalar or length-2.
#'   - NULL: use a default based on spatial extent.
#' @param alpha Smoothness parameter for SPDE construction (default 2).
#' @param suppress_warnings If TRUE, suppress INLA warnings (default TRUE).
#' @param penalty_range PC-prior range anchor. If NULL, defaults to ~1/10 spatial extent.
#' @param link One of \code{"identity"}, \code{"log"}, \code{"logit"}, \code{"probit"}.
#'
#' @return A function with signature \code{function(x, s, g_init=NULL, fix_g=FALSE, output=NULL)}.
#'
#' @export
ebnm_Matern_generator <- function(locations,
                                  max.edge = NULL,
                                  alpha = 2,
                                  suppress_warnings = TRUE,
                                  penalty_range = NULL,
                                  link = c("identity", "log", "logit", "probit")) {

  link <- match.arg(link)

  # normalize locations
  loc_info <- .normalize_locations(locations)
  loc_mat <- loc_info$loc
  d <- loc_info$d

  # build mesh & A based on d
  meshA <- .build_mesh_A(loc_mat, max.edge = max.edge)
  mesh <- meshA$mesh
  A <- meshA$A

  # penalty range (local var; do NOT mutate outer)
  penalty_range0 <- if (is.null(penalty_range)) .default_penalty_range(loc_mat) else penalty_range
  if (!is.numeric(penalty_range0) || length(penalty_range0) != 1L || is.na(penalty_range0) || penalty_range0 <= 0) {
    stop("penalty_range must be a single positive number (or NULL).")
  }

  ebnm_Matern <- function(x, s, g_init = NULL, fix_g = FALSE, output = NULL) {

    n <- nrow(loc_mat)

    # ---- checks ----
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

    if (!(length(s) == 1L || length(s) == length(x))) stop("s must have length 1 or length(x).")
    if (length(s) == 1L) s <- rep(s, length(x))
    if (anyNA(x) || anyNA(s)) stop("x and s must not contain NA.")
    if (any(s <= 0)) stop("All entries of s must be > 0.")

    # init g
    if (is.null(g_init)) g_init <- Matern(theta = log(penalty_range0))
    if (!is.null(g_init$theta) && (!is.numeric(g_init$theta) || length(g_init$theta) != 1L || is.na(g_init$theta))) {
      stop("g_init$theta must be a single non-NA numeric (or NULL).")
    }

    # Build SPDE object (PC priors)
    spde <- INLA::inla.spde2.pcmatern(
      mesh = mesh,
      alpha = alpha,
      prior.range = c(penalty_range0, 0.5),
      prior.sigma = c(1, NA)
    )

    idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

    # ---- Step A: EB optimization over theta and beta0 (unless fix_g) ----
    if (!fix_g) {
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
            control.link = list(model = link),
            hyper = list(prec = list(fixed = TRUE, initial = 0))
          ),
          control.fixed = INLA::control.fixed(prec = 0),
          data = INLA::inla.stack.data(stackA),
          control.predictor = list(A = INLA::inla.stack.A(stackA), link = 1),
          control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE),
          silent = TRUE
        )
      }

      resA <- if (suppress_warnings) suppressWarnings(runA()) else runA()

      theta_hat <- resA$mode$theta
      beta0_hat <- resA$summary.fixed$mean

      fitted_g <- Matern(theta_hat)

    } else {
      # fixed g: theta from g_init, beta0 set to 0 (consistent with your old behavior)
      fitted_g <- g_init
      beta0_hat <- 0
    }

    # ---- Step B: fit with theta fixed, beta0 as offset ----
    stackB <- INLA::inla.stack(
      data = list(Y = as.numeric(x)),
      A = list(A, matrix(beta0_hat, nrow = n, ncol = 1)),
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
          control.link = list(model = link),
          hyper = list(prec = list(fixed = TRUE, initial = 0))
        ),
        control.fixed = INLA::control.fixed(prec = 0),
        data = INLA::inla.stack.data(stackB),
        control.predictor = list(A = INLA::inla.stack.A(stackB), link = 1),
        control.mode = INLA::control.mode(theta = fitted_g$theta, fixed = TRUE),
        control.compute = list(mlik = TRUE),
        silent = TRUE
      )
    }

    resB <- if (suppress_warnings) suppressWarnings(runB()) else runB()

    ii <- INLA::inla.stack.index(stackB, tag = "est")$data

    posterior <- data.frame(
      mean = resB$summary.fitted.values$mean[ii],
      var  = resB$summary.fitted.values$sd[ii]^2
    )
    posterior$second_moment <- posterior$mean^2 + posterior$var

    posterior_spatial_field <- resB$summary.random$spatial.field
    log_likelihood <- resB$mlik[[1]]
    class(log_likelihood) <- "logLik"

    posterior_sampler <- function(n = n) {
      warning("posterior_sampler is not implemented yet for Matern.")
      matrix(NA_real_, nrow = n, ncol = length(x))
    }

    out <- list(
      posterior = posterior,
      fitted_g = fitted_g,
      g_init = g_init,
      log_likelihood = log_likelihood,
      posterior_sampler = posterior_sampler,
      data = data.frame(x = x, s = s),
      prior_family = paste0(link, "_Matern"),
      posterior_spatial_field = posterior_spatial_field,
      mesh = mesh,
      inla_result = resB
    )

    structure(out, class = c("list", "ebnm"))
  }

  if (suppress_warnings) suppressWarnings(ebnm_Matern) else ebnm_Matern
}
