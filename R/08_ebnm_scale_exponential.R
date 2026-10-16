## Scale mixture of exponentials: a non-negative ebnm solver.
##
## ebnm_scale_exponential() fits the prior
##   g = pi_0 delta_mode + sum_k pi_k (mode + Exp(scale = sigma_k))
## on a fixed grid of scales sigma_k, estimating only the weights (mixsqp).
## It is the non-negative analogue of ebnm::ebnm_normal_scale_mixture() and
## more flexible than ebnm::ebnm_point_exponential(), which has a single scale.
##
## Every quantity is written in terms of t = s / sigma - (x - mode) / s and the
## Mills ratio R(t) = (1 - Phi(t)) / phi(t), so that nothing cancels when the
## noise s is much larger than sigma (where ebnm_point_exponential() loses
## accuracy).

## Moments of Z ~ N(-t, 1) truncated to [0, Inf), and the log Mills ratio
## log R(t) = log(1 - Phi(t)) - log phi(t):
##   h  = E[Z]   = 1 / R(t) - t,
##   m2 = E[Z^2] = 1 - t h(t).
## Direct formulas for t < 5. For t >= 5 use the continued fraction
## 1 / R(t) = t + 1 / (t + c), c = 2 / (t + 3 / (t + 4 / (t + ...))), so that
## h = 1 / (t + c), m2 = c h and log R = -log(t + h), with no cancellation
## (depth 30 is exact to double precision for t >= 5).
.expmix_trunc_moments <- function(t, depth = 30L) {
  h <- m2 <- log_R <- t
  small <- t < 5
  if (any(small)) {
    ts <- t[small]
    lr <- stats::pnorm(ts, lower.tail = FALSE, log.p = TRUE) -
      stats::dnorm(ts, log = TRUE)
    hs <- exp(-lr) - ts
    log_R[small] <- lr
    h[small] <- hs
    m2[small] <- 1 - ts * hs
  }
  if (any(!small)) {
    tl <- t[!small]
    cf <- 0
    for (k in seq.int(depth, 2L)) cf <- k / (tl + cf)
    hl <- 1 / (tl + cf)
    log_R[!small] <- -log(tl + hl)
    h[!small] <- hl
    m2[!small] <- cf * hl
  }
  list(h = h, m2 = m2, log_R = log_R)
}

## Default grid: from min(s) / 10 up to 2 * max(x - mode), spaced by `mult`.
.expmix_default_scale <- function(x, s, mode, mult = sqrt(2)) {
  lo <- min(s) / 10
  hi <- max(2 * max(x - mode), lo * mult)
  n_grid <- ceiling(log(hi / lo) / log(mult)) + 1
  lo * mult^(seq_len(n_grid) - 1)
}

## Component log marginal likelihoods and posterior moments, all at once.
## ll is n x (K + 1), column 1 the point mass at mode; h and m2 are n x K
## (moments of (theta - mode) / s under each exponential component). For scale
## sigma (rate a = 1 / sigma), with z = (x - mode) / s and t = s a - z:
##   log p(x | k) = log a - log(2 pi) / 2 - z^2 / 2 + log R(t).
.expmix_components <- function(x, s, mode, scale) {
  z <- (x - mode) / s
  base <- -0.5 * log(2 * pi) - 0.5 * z^2
  tt <- outer(s, 1 / scale) - z
  mom <- .expmix_trunc_moments(tt)
  ll <- cbind(base - log(s), base + mom$log_R - rep(log(scale), each = length(x)))
  list(ll = ll, h = mom$h, m2 = mom$m2, z = z)
}

.expmix_loglik_matrix <- function(x, s, mode, scale) {
  .expmix_components(x, s, mode, scale)$ll
}

## The prior as an ebnm "gammamix" with all shapes 1 (component 1: scale 0,
## the point mass), the class ebnm::ebnm_point_exponential() also returns.
## Same object as ebnm::gammamix(), built without its data.frame() call.
.expmix_g <- function(pi, scale, mode) {
  K <- length(scale) + 1L
  structure(
    list(pi = pi, shape = rep(1, K), scale = c(0, scale), shift = rep(mode, K)),
    class = "gammamix",
    row.names = seq_len(K)
  )
}

## data.frame() is slow relative to the rest of a call inside flashier.
.expmix_df <- function(cols) {
  structure(cols, class = "data.frame", row.names = .set_row_names(length(cols[[1]])))
}

.expmix_check_g <- function(g) {
  ok <- inherits(g, "gammamix") && all(g$shape == 1) &&
    length(unique(g$shift)) == 1L && g$scale[1] == 0 && all(g$scale[-1] > 0)
  if (!ok) {
    stop("g_init must be a \"gammamix\" with all shapes 1, one common shift, ",
         "scale[1] = 0 (the point mass) and positive other scales.", call. = FALSE)
  }
  invisible(g)
}

#' Scale mixture of exponentials (non-negative) ebnm solver
#'
#' @description
#' Solves the empirical Bayes normal means problem
#' \deqn{x_j \sim N(\theta_j, s_j^2), \quad \theta_j \sim g,}
#' with the prior
#' \deqn{g = \pi_0 \delta_{m} + \sum_{k = 1}^K \pi_k \, (m + \mathrm{Exp}(\sigma_k)),}
#' a point mass at the mode \eqn{m} plus a mixture of exponentials with scales
#' (means) \eqn{\sigma_1 < \dots < \sigma_K} on a fixed grid. Only the weights
#' \eqn{\pi} are estimated, by maximizing the marginal likelihood with
#' [mixsqp::mixsqp()]. This is the non-negative analogue of
#' [ebnm::ebnm_normal_scale_mixture()]: it nests
#' [ebnm::ebnm_point_exponential()] (one scale) and can fit heavier tails and
#' multi-scale non-negative effects.
#'
#' The interface follows the `ebnm` package, so the function can be passed
#' directly as `ebnm_fn` to `flashier::flash()`, `flashier::flash_greedy()` or
#' `flashier::flash_factors_init()` (for example `ebnm_fn = c(ebnm_Matern,
#' ebnm_scale_exponential)` for a semi-non-negative fit, or on both sides for a
#' non-negative matrix factorization).
#'
#' @details
#' Each component's marginal likelihood and posterior (a normal truncated to
#' \eqn{[m, \infty)}) are computed through the Mills ratio with a continued
#' fraction in the tail, so the log-likelihood and posterior moments stay
#' accurate when \eqn{s_j \gg \sigma_k}, a regime where
#' [ebnm::ebnm_point_exponential()] loses precision.
#'
#' **Grid.** With `scale = "estimate"` the grid runs from `min(s) / 10` to
#' `2 * max(x - mode)` in multiples of `mult`. If `g_init` is supplied its grid
#' is reused (and its weights are the starting point of the optimization),
#' unless `scale` is given explicitly. A warm-start grid that no longer reaches
#' `max(x - mode)` is replaced by the default grid, so warm starts inside
#' flashier's backfit cannot truncate the prior.
#'
#' @param x Vector of observations.
#' @param s Vector of standard errors (or a scalar), all positive and finite.
#' @param mode Location of the point mass and lower end of the support. Fixed
#'   (not estimated); default 0.
#' @param scale Either `"estimate"` (default grid, see Details) or a numeric
#'   vector of positive exponential scales \eqn{\sigma_k} (the point mass is
#'   always added).
#' @param g_init Optional prior of class `"gammamix"` with all shapes 1 and
#'   `scale[1] = 0` (e.g. the `fitted_g` of a previous fit, or of
#'   [ebnm::ebnm_point_exponential()]), used as the grid and starting weights, or as the prior
#'   itself when `fix_g = TRUE`.
#' @param fix_g If `TRUE`, `g_init` is used as the prior without estimation.
#' @param output Character vector of outputs, as in [ebnm::ebnm()]; see
#'   [ebnm::ebnm_output_default()] and [ebnm::ebnm_output_all()].
#' @param mult Spacing of the default grid (ratio of consecutive scales).
#' @param control List of control parameters passed to [mixsqp::mixsqp()]
#'   (defaults here: `verbose = FALSE`, `eps = 1e-12`, `convtol.sqp = 1e-10`).
#'
#' @return An object of class `"ebnm"` with elements `data`, `posterior`
#'   (columns `mean`, `sd`, `second_moment`, `lfsr` as requested), `fitted_g`
#'   (class `"gammamix"`, as in [ebnm::ebnm_point_exponential()]: `pi`,
#'   `shape` (all 1), `scale` with `scale[1] = 0` for the point mass, and
#'   `shift` = `mode`), `log_likelihood` (class `"logLik"`) and, if requested,
#'   `posterior_sampler`. The local false sign rate is the posterior
#'   probability that \eqn{\theta_j = m}.
#'
#' @examples
#' set.seed(1)
#' theta <- c(rep(0, 200), rexp(100, 1), rexp(100, 1 / 10))
#' x <- theta + rnorm(400)
#' fit <- ebnm_scale_exponential(x, s = 1)
#' fit$fitted_g$pi[1]          # estimated null proportion
#' head(fit$posterior)
#'
#' @export
ebnm_scale_exponential <- function(x,
                                   s = 1,
                                   mode = 0,
                                   scale = "estimate",
                                   g_init = NULL,
                                   fix_g = FALSE,
                                   output = ebnm::ebnm_output_default(),
                                   mult = sqrt(2),
                                   control = NULL) {
  if (length(s) == 1L) s <- rep(s, length(x))
  if (length(s) != length(x)) {
    stop("s must be a scalar or have the same length as x.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(!is.finite(s)) || any(s <= 0)) {
    stop("x must be finite and s must be positive and finite.", call. = FALSE)
  }
  if (!is.numeric(mode) || length(mode) != 1L || !is.finite(mode)) {
    stop("mode must be a single finite number.", call. = FALSE)
  }
  if (!is.null(g_init)) .expmix_check_g(g_init)
  if (fix_g && is.null(g_init)) {
    stop("fix_g = TRUE requires g_init.", call. = FALSE)
  }

  pi_init <- NULL
  if (fix_g) {
    mode <- g_init$shift[1]
    grid <- g_init$scale[-1]
  } else if (!identical(scale, "estimate")) {
    if (!is.numeric(scale) || any(!is.finite(scale)) || any(scale <= 0)) {
      stop("scale must be \"estimate\" or a vector of positive scales.", call. = FALSE)
    }
    grid <- sort(scale)
  } else if (!is.null(g_init) && g_init$shift[1] == mode &&
             max(g_init$scale) >= max(x - mode)) {
    grid <- g_init$scale[-1]
    pi_init <- g_init$pi
  } else {
    grid <- .expmix_default_scale(x, s, mode, mult)
  }

  comp <- .expmix_components(x, s, mode, grid)
  ll <- comp$ll
  z <- comp$z
  ll_max <- ll[cbind(seq_along(x), max.col(ll, ties.method = "first"))]
  L <- exp(ll - ll_max)

  if (fix_g) {
    w <- g_init$pi
  } else {
    ## Rows of L have maximum 1, so a small eps is safe; mixsqp's default
    ## (1e-8) stops short of the optimum on fine, near-collinear grids.
    ctrl <- utils::modifyList(
      list(verbose = FALSE, eps = 1e-12, convtol.sqp = 1e-10, tol.svd = 0,
           normalize.rows = FALSE),
      as.list(control)
    )
    w <- if (is.null(pi_init)) {
      mixsqp::mixsqp(L, control = ctrl)$x
    } else {
      ## Warm start (e.g. flashier backfit): skip mixsqp's initial EM steps.
      ctrl <- utils::modifyList(list(numiter.em = 0), ctrl)
      mixsqp::mixsqp(L, x0 = pmax(pi_init, 1e-8), control = ctrl)$x
    }
    w <- pmax(w, 0)
    w <- w / sum(w)
  }

  lik_rows <- as.vector(L %*% w)
  loglik <- sum(log(lik_rows) + ll_max)

  g <- .expmix_g(w, grid, mode)
  retlist <- list()

  if ("data" %in% output) {
    retlist$data <- .expmix_df(list(x = x, s = s))
  }

  want_post <- intersect(
    output,
    c("posterior_mean", "posterior_sd", "posterior_second_moment", "lfsr")
  )
  if (length(want_post) > 0 || "posterior_sampler" %in% output) {
    ## Posterior component weights (n x (K + 1)).
    W <- (L * rep(w, each = length(x))) / lik_rows
    ## theta - mode = s * Z, Z ~ N(-t, 1) truncated to [0, Inf).
    m1 <- s * rowSums(W[, -1, drop = FALSE] * comp$h)
    m2 <- s^2 * rowSums(W[, -1, drop = FALSE] * comp$m2)
    post_var <- pmax(m2 - m1^2, 0)
    post_mean <- mode + m1
    post_m2 <- m2 + 2 * mode * m1 + mode^2

    posterior <- list()
    if ("posterior_mean" %in% output) posterior$mean <- post_mean
    if ("posterior_sd" %in% output) posterior$sd <- sqrt(post_var)
    if ("posterior_second_moment" %in% output) posterior$second_moment <- post_m2
    if ("lfsr" %in% output) posterior$lfsr <- W[, 1]
    if (length(want_post) > 0) retlist$posterior <- .expmix_df(posterior)
  }

  if ("fitted_g" %in% output) {
    retlist$fitted_g <- g
  }

  if ("log_likelihood" %in% output) {
    retlist$log_likelihood <- structure(
      loglik,
      class = "logLik",
      df = length(grid),
      nobs = length(x)
    )
  }

  if ("posterior_sampler" %in% output) {
    retlist$posterior_sampler <- function(nsamp) {
      n <- length(x)
      out <- matrix(mode, nsamp, n)
      for (j in seq_len(n)) {
        k <- sample.int(length(w), nsamp, replace = TRUE, prob = W[j, ])
        nz <- k > 1L
        if (any(nz)) {
          tj <- s[j] / grid[k[nz] - 1L] - z[j]
          ## Z' = Z + t ~ N(0, 1) truncated to [t, Inf), by inverse CDF.
          logp <- stats::pnorm(tj, lower.tail = FALSE, log.p = TRUE) +
            log(stats::runif(sum(nz)))
          zp <- stats::qnorm(logp, lower.tail = FALSE, log.p = TRUE)
          out[nz, j] <- mode + s[j] * pmax(zp - tj, 0)
        }
      }
      out
    }
  }

  structure(retlist, class = c("ebnm", "list"))
}
