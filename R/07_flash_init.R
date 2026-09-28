## Smooth rank-one initialization for flashier's greedy step.
##
## flash_greedy_init_smooth() is an `init_fn` for flashier::flash_greedy(). The
## .smooth_init_* helpers below do the matrix computation: they take the residual
## matrix and the coordinates, and know nothing about flashier, priors or noise.

#' Smooth rank-one initialization for flashier's greedy step
#'
#' @description
#' An initialization function for [flashier::flash_greedy()] (argument `init_fn`)
#' that starts each new factor from a smooth curve. It returns the rank-one
#' approximation \eqn{R \approx \ell f^\top} of the current residual matrix \eqn{R},
#' fit by alternating least squares on the observed entries only, with the vector of
#' dimension `smooth_dim` (\eqn{\ell}) restricted to a low-dimensional space of
#' functions of the coordinates `x`:
#' \deqn{\min_{\alpha, f} \sum_{(i, j) \in \Omega} \{R_{ij} - (B\alpha)_i f_j\}^2,}
#' where \eqn{\Omega} is the set of observed entries and \eqn{B} is the basis
#' evaluated at `x`.
#'
#' flashier's default, [flashier::flash_greedy_init_default()], runs the same
#' alternating least squares with \eqn{\ell} unrestricted. When each column is
#' observed at only a few rows (sparse, irregular sampling), a column observed once
#' is fit exactly by \eqn{f_j = R_{ij} / \ell_i} whatever \eqn{\ell} is, and the
#' unrestricted fit can put almost all of \eqn{f} on a few such columns. Restricting
#' \eqn{\ell} to a few smooth functions makes the starting curve one that many
#' columns share.
#'
#' @details
#' **Bases.**
#' * `"ns"` (default): natural cubic splines, `splines::ns(df = df, intercept =
#'   TRUE)`. `df` is the total number of basis functions, including the constant, so
#'   the default `df = 4` is a 4-dimensional space. The interior knots are at the
#'   quantiles of the distinct values of `x`, so the basis does not depend on the
#'   order of `x` or on how often a value is repeated.
#' * `"linear"`: intercept and slope, i.e. the best starting curve \eqn{a + b x}.
#' * `"constant"`: \eqn{\ell = 1}. This is solved in closed form: \eqn{f_j} is the
#'   mean of the observed residuals of column \eqn{j}. `nstarts`, `maxiter`, `tol` and
#'   `seed` are not used.
#'
#' The basis only restricts the starting values. Once they are returned, flashier
#' updates the factor with the priors in `ebnm_fn` as usual, and the fitted curve is
#' not confined to the basis: `df` sets the complexity of the starting curve, not of
#' the fit.
#'
#' **Algorithm.** Each start draws a random direction in the basis, then alternates
#' two exact least-squares steps on the observed entries: the scores
#' \eqn{f_j = \sum_{i \in \Omega_j} \ell_i R_{ij} / \sum_{i \in \Omega_j} \ell_i^2}
#' given \eqn{\ell}, and the basis coefficients \eqn{\alpha} given \eqn{f} (a
#' weighted least-squares problem with one row per row of \eqn{R}, solved by SVD; a
#' rank-deficient problem gets its minimum-norm solution). It stops when the relative
#' decrease of the residual sum of squares (RSS) is at most `tol`, or after `maxiter`
#' iterations, and the start with the smallest RSS is returned. No noise variance,
#' prior, ridge or roughness penalty is involved, and missing entries (`NA` in the
#' residuals) are left out, not treated as zeros.
#'
#' A column with no observed entry gets \eqn{f_j = 0}, as does a column whose observed
#' rows all have \eqn{|\ell_i|} at rounding level, at most \eqn{100\epsilon
#' \max_i |\ell_i|} (the ratio would then be rounding error). A row of the smooth
#' dimension with no observed entry gets the value of the fitted curve at its coordinate.
#'
#' **Limitations.** The restriction does not make a sparse design identifiable. A
#' column observed once is still fit exactly by \eqn{f_j = R_{ij} / \ell_i}, so one or
#' a few columns can still dominate \eqn{f}, e.g. a column whose single observation
#' falls near a zero of \eqn{\ell}; and if every column is observed once, many curves
#' reach the same RSS. For the same reason a start can settle where \eqn{\ell} is
#' nearly zero at the time of a column observed once, with that column's \eqn{f_j}
#' growing without bound. Such a start fits the other columns worse, so the
#' smallest-RSS rule avoids it as long as one start escapes; this is why the default
#' uses 10 starts (on sparse simulated data, 3 starts all got stuck in 2 of 3 data
#' sets, and 10 starts never did). There are no sign constraints: \eqn{\ell} and
#' \eqn{f} can take both signs, so this initialization does not suit non-negative
#' priors (a start that violates the prior's constraint can make greedy stop early). Only data given as an
#' ordinary numeric matrix are supported, with missing entries as `NA`: not tensors,
#' sparse `Matrix` objects or low-rank (`u`, `d`, `v`) data.
#'
#' @param flash The `flash_fit` object that flashier passes to `init_fn` (a flash
#'   object also works).
#' @param x Coordinates of dimension `smooth_dim`: a numeric vector with one finite
#'   value per row (`smooth_dim = 1`) or per column (`smooth_dim = 2`) of the data, in
#'   the same order. Values may be unsorted and repeated.
#' @param smooth_dim `1` (default) to restrict the loadings (rows of the data), `2` to
#'   restrict the factors (columns).
#' @param basis `"ns"` (default), `"linear"` or `"constant"`; see Details.
#' @param df For `basis = "ns"`: number of basis functions, including the constant.
#'   An integer of at least 2 and at most the number of distinct values of `x`
#'   (`df = 2` is the linear space). Not used by the other bases.
#' @param nstarts Number of random starts. Each costs a few ALS iterations on the
#'   observed entries, so 10 starts take a fraction of a second on matrices of a few
#'   hundred by a few hundred.
#' @param maxiter Maximum number of iterations per start.
#' @param tol Tolerance on the relative decrease of the RSS.
#' @param seed Random seed for the starts. The caller's random number stream is
#'   restored afterwards. `NULL` draws the starts from the current stream instead.
#'
#' @return A list of two vectors, as `init_fn` must return: the starting values for
#'   the rows (length `nrow` of the data) and for the columns (length `ncol`). Their
#'   outer product is the fitted rank-one approximation. The two vectors have equal
#'   norms, and the entry of the smooth vector with the largest absolute value is
#'   positive. If the observed residuals are all zero, both vectors are zero, and
#'   flashier then adds no further factor. Residuals that are only numerically
#'   orthogonal to the basis give vectors at rounding level, not zero; flashier then
#'   optimizes from them as usual.
#'
#' @seealso [flashier::flash_greedy()], [flashier::flash_greedy_init_default()].
#'
#' @examples
#' \donttest{
#' if (requireNamespace("flashier", quietly = TRUE)) {
#'   set.seed(1)
#'   tt <- seq(0, 1, length.out = 40)
#'   Y <- outer(sin(2 * pi * tt), rnorm(60)) + matrix(rnorm(40 * 60, 0, 0.5), 40, 60)
#'   for (j in 1:60) Y[-sample(40, 6), j] <- NA   # each column observed at 6 times
#'   prior_L <- ebnm_Matern_generator(setup = Matern_setup(tt, alpha = 2))
#'
#'   init <- function(f) flash_greedy_init_smooth(f, x = tt)
#'   fl <- flashier::flash_init(Y, var_type = 0)
#'   fl <- flashier::flash_greedy(fl, Kmax = 2, init_fn = init,
#'                                ebnm_fn = list(prior_L, ebnm::ebnm_point_normal))
#'
#'   # A richer starting space, or intercept + slope only:
#'   init_ns6 <- function(f) flash_greedy_init_smooth(f, x = tt, df = 6)
#'   init_lin <- function(f) flash_greedy_init_smooth(f, x = tt, basis = "linear")
#'   # A constant curve (closed form, no iterations):
#'   init_const <- function(f) flash_greedy_init_smooth(f, x = tt, basis = "constant")
#'
#'   # Smooth factors instead of loadings (time along the columns):
#'   init_t <- function(f) flash_greedy_init_smooth(f, x = tt, smooth_dim = 2)
#'   fl_t <- flashier::flash_init(t(Y), var_type = 0)
#'   fl_t <- flashier::flash_greedy(fl_t, Kmax = 2, init_fn = init_t,
#'                                  ebnm_fn = list(ebnm::ebnm_point_normal, prior_L))
#' }
#' }
#' @export
flash_greedy_init_smooth <- function(flash, x, smooth_dim = 1L,
                                     basis = c("ns", "constant", "linear"), df = 4L,
                                     nstarts = 10L, maxiter = 100L, tol = 1e-6,
                                     seed = 666L) {
  basis <- match.arg(basis)
  if (!inherits(flash, c("flash_fit", "flash"))) {
    stop("`flash` must be the flash_fit object that flashier passes to `init_fn`.",
         call. = FALSE)
  }
  # check how the data are stored before computing residuals, which would densify
  # sparse data and fails for low-rank data
  Y <- (if (inherits(flash, "flash")) flash[["flash_fit"]] else flash)[["Y"]]
  if (!is.null(Y) && !(is.matrix(Y) && is.numeric(Y))) {
    stop(sprintf(paste0("flash_greedy_init_smooth() needs the data as an ordinary numeric ",
                        "matrix (got class \"%s\"); tensors, sparse Matrix objects and ",
                        "low-rank data are not supported."),
                 paste(class(Y), collapse = "\", \"")), call. = FALSE)
  }
  out <- .smooth_init_matrix(stats::residuals(flash), x, smooth_dim = smooth_dim,
                             basis = basis, df = df, nstarts = nstarts,
                             maxiter = maxiter, tol = tol, seed = seed)
  list(out$row_values, out$column_values)
}

.smooth_init_count <- function(v, nm, lower) {
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v != round(v) || v < lower) {
    stop(sprintf("`%s` must be an integer >= %d.", nm, lower), call. = FALSE)
  }
  as.integer(v)
}

# Basis of the starting curve at the coordinates x (one row per coordinate).
.smooth_init_basis <- function(x, basis, df = NULL) {
  ux <- sort(unique(x))
  B <- switch(basis,
    constant = matrix(1, length(x), 1L),
    linear = {
      if (length(ux) < 2L) {
        stop("basis = \"linear\" needs at least 2 distinct values of `x`; use basis = \"constant\".",
             call. = FALSE)
      }
      cbind(1, (x - mean(range(ux))) / diff(range(ux)))
    },
    ns = {
      if (length(ux) < df) {
        stop(sprintf(paste0("basis = \"ns\" with df = %d needs at least %d distinct values of ",
                            "`x` (there are %d); lower `df`, or use basis = \"linear\" or ",
                            "\"constant\"."), df, df, length(ux)), call. = FALSE)
      }
      # knots from the distinct values, then evaluated at every x
      nb <- splines::ns(ux, df = df, intercept = TRUE)
      splines::ns(x, knots = attr(nb, "knots"), Boundary.knots = attr(nb, "Boundary.knots"),
                  intercept = TRUE)
    })
  B <- unname(matrix(as.numeric(B), nrow = length(x)))
  if (qr(B)$rank < ncol(B)) {
    stop(sprintf(paste0("The \"%s\" basis is rank deficient at these coordinates; lower `df`, ",
                        "or use basis = \"linear\" or \"constant\"."), basis), call. = FALSE)
  }
  B
}

# Scores given the loading: f_j = sum_i Z_ij l_i R_ij / sum_i Z_ij l_i^2. R is zero
# at the missing entries and Zn is the 0/1 observed indicator. A column whose
# observed l_i are all at rounding level (|l_i| <= 100 eps max|l|) gets f_j = 0.
# That is a numerical-zero test, not a regularization.
.smooth_init_fstep <- function(R, Zn, l) {
  den <- as.numeric(crossprod(Zn, l^2))
  num <- as.numeric(crossprod(R, l))
  ok <- as.numeric(crossprod(Zn, abs(l) > 100 * .Machine$double.eps * max(abs(l)))) > 0
  f <- numeric(length(den))
  f[ok] <- num[ok] / den[ok]
  f
}

# Basis coefficients given the scores: minimizes sum_(i,j) observed (R_ij - f_j Q_i a)^2.
# Row i enters with weight w_i = sum_j Z_ij f_j^2 and response b_i / w_i, where
# b_i = sum_j Z_ij f_j R_ij; rows with w_i = 0 carry no information. The small
# least-squares problem is solved by SVD (minimum-norm solution if rank deficient).
.smooth_init_alpha_step <- function(R, Zn, Q, f) {
  w <- as.numeric(Zn %*% f^2)
  b <- as.numeric(R %*% f)
  rows <- w > 0
  sw <- sqrt(w[rows])
  A <- sw * Q[rows, , drop = FALSE]
  y <- b[rows] / sw
  sv <- svd(A)
  keep <- sv$d > max(dim(A)) * .Machine$double.eps * sv$d[1]
  if (!any(keep)) return(rep(0, ncol(Q)))
  as.numeric(sv$v[, keep, drop = FALSE] %*% (crossprod(sv$u[, keep, drop = FALSE], y) / sv$d[keep]))
}

# Rank-one alternating least squares with l = Q a, from the starting directions in
# the columns of A0. Q has orthonormal columns, so |l| = |a|. Values of R at the
# missing entries (!Z) are placeholders and are never used.
.smooth_init_als <- function(R, Z, Q, A0, maxiter, tol) {
  R[!Z] <- 0
  Zn <- Z * 1
  tss <- sum(R^2)
  rss <- function(l, f) sum(Zn * (R - tcrossprod(l, f))^2)
  # rounding-level tolerances (not tuning parameters): an exact ALS step cannot
  # increase the RSS, and an RSS at rounding level of the total is an exact fit
  tiny <- 100 * .Machine$double.eps * tss
  starts <- vector("list", ncol(A0))
  for (s in seq_len(ncol(A0))) {
    l <- as.numeric(Q %*% (A0[, s] / sqrt(sum(A0[, s]^2))))
    f <- .smooth_init_fstep(R, Zn, l)
    r <- rss(l, f)
    best <- if (is.finite(r)) list(l = l, f = f, rss = r) else NULL
    status <- if (is.null(best)) "failed" else if (all(f == 0)) "zero" else "maxiter"
    it <- 0L
    while (identical(status, "maxiter") && it < maxiter) {
      it <- it + 1L
      a <- .smooth_init_alpha_step(R, Zn, Q, f)
      na <- sqrt(sum(a^2))
      if (!is.finite(na) || na == 0) {
        status <- "numerical"
        break
      }
      l <- as.numeric(Q %*% (a / na))
      f <- .smooth_init_fstep(R, Zn, l)
      r_new <- rss(l, f)
      if (!is.finite(r_new)) {
        status <- "numerical"
        break
      }
      if (r_new - r > sqrt(.Machine$double.eps) * r + tiny) {
        status <- "rss_increase"
        break
      }
      if (r_new < best$rss) best <- list(l = l, f = f, rss = r_new)
      if (r - r_new <= tol * r || r_new <= tiny) status <- "converged"
      r <- r_new
    }
    starts[[s]] <- list(fit = best, iter = it, status = status)
  }
  ok <- vapply(starts, function(st) !is.null(st$fit), logical(1))
  if (!any(ok)) {
    stop("The smooth initialization failed: no start gave a finite rank-one fit.", call. = FALSE)
  }
  rss_s <- vapply(starts, function(st) if (is.null(st$fit)) NA_real_ else st$fit$rss, numeric(1))
  k <- which(ok)[which.min(rss_s[ok])]
  list(l = starts[[k]]$fit$l, f = starts[[k]]$fit$f, rss = rss_s[k], tss = tss,
       iter = starts[[k]]$iter, status = starts[[k]]$status,
       starts = data.frame(start = seq_along(starts), rss = rss_s,
                           iter = vapply(starts, `[[`, integer(1), "iter"),
                           status = vapply(starts, `[[`, character(1), "status")))
}

# basis = "constant": l = 1 and f_j = mean of the observed residuals of column j.
.smooth_init_constant <- function(R, Z) {
  R[!Z] <- 0
  cnt <- colSums(Z)
  f <- ifelse(cnt > 0, colSums(R) / pmax(cnt, 1), 0)
  l <- rep(1, nrow(R))
  list(l = l, f = f, rss = sum(Z * (R - tcrossprod(l, f))^2), tss = sum(R^2),
       iter = 0L, status = if (all(f == 0)) "zero" else "exact", starts = NULL)
}

# Draw the random starting directions under a local seed; the caller's RNG state is
# restored afterwards.
.smooth_init_starts <- function(q, nstarts, seed) {
  draw <- function() matrix(stats::rnorm(q * nstarts), q, nstarts)
  if (is.null(seed)) return(draw())
  genv <- globalenv()
  old <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
    get(".Random.seed", envir = genv, inherits = FALSE)
  }
  on.exit(if (is.null(old)) rm(".Random.seed", envir = genv) else assign(".Random.seed", old, envir = genv))
  set.seed(seed)
  draw()
}

# The matrix computation behind flash_greedy_init_smooth(). R is the residual matrix
# with NA at the missing entries. Returns the starting values in the original row
# and column order, plus diagnostics: the smooth and the other vector, the observed
# RSS of the rank-one fit and the total sum of squares (tss), and, for the best
# start, its iterations and status ("converged", "maxiter", "rss_increase",
# "numerical", "zero" = every score exactly zero, e.g. all-zero residuals, "exact" =
# constant basis); starts has one row per start.
.smooth_init_matrix <- function(R, x, smooth_dim = 1L, basis = c("ns", "constant", "linear"),
                                df = 4L, nstarts = 10L, maxiter = 100L, tol = 1e-6,
                                seed = 666L) {
  basis <- match.arg(basis)
  if (!is.matrix(R) || !is.numeric(R)) {
    stop(sprintf(paste0("The residuals must be an ordinary numeric matrix (got class \"%s\"); ",
                        "tensors, sparse matrices and low-rank data are not supported."),
                 paste(class(R), collapse = "\", \"")), call. = FALSE)
  }
  if (!is.numeric(smooth_dim) || length(smooth_dim) != 1L || !isTRUE(smooth_dim %in% c(1, 2))) {
    stop("`smooth_dim` must be 1 or 2.", call. = FALSE)
  }
  if (smooth_dim == 2) R <- t(R)
  if (!is.numeric(x) || length(x) != nrow(R) || any(!is.finite(x))) {
    stop(sprintf("`x` must be a numeric vector of %d finite coordinates, one per %s of the data.",
                 nrow(R), if (smooth_dim == 1) "row" else "column"), call. = FALSE)
  }
  if (basis == "ns") df <- .smooth_init_count(df, "df", 2L)
  nstarts <- .smooth_init_count(nstarts, "nstarts", 1L)
  maxiter <- .smooth_init_count(maxiter, "maxiter", 1L)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a non-negative number.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))) {
    stop("`seed` must be NULL or a single number.", call. = FALSE)
  }
  storage.mode(R) <- "double"
  Z <- !is.na(R)
  if (!any(Z)) stop("The residual matrix has no observed entries.", call. = FALSE)
  if (any(!is.finite(R[Z]))) stop("The observed residuals must be finite.", call. = FALSE)

  B <- .smooth_init_basis(as.numeric(x), basis, df)
  fit <- if (basis == "constant") {
    .smooth_init_constant(R, Z)
  } else if (all(R[Z] == 0)) {
    list(l = rep(0, nrow(R)), f = rep(0, ncol(R)), rss = 0, tss = 0, iter = 0L,
         status = "zero", starts = NULL)
  } else {
    Q <- qr.Q(qr(B))
    .smooth_init_als(R, Z, Q, .smooth_init_starts(ncol(Q), nstarts, seed), maxiter, tol)
  }

  # equal norms and a fixed sign; the outer product is unchanged
  l <- fit$l
  f <- fit$f
  if (all(f == 0)) {
    l <- 0 * l
  } else {
    if (l[which.max(abs(l))] < 0) {
      l <- -l
      f <- -f
    }
    nl <- sqrt(sum(l^2))
    nf <- sqrt(sum(f^2))
    l <- l * sqrt(nf / nl)
    f <- f * sqrt(nl / nf)
  }
  list(row_values = if (smooth_dim == 1) l else f,
       column_values = if (smooth_dim == 1) f else l,
       smooth_values = l, other_values = f, rss = fit$rss, tss = fit$tss,
       iter = fit$iter, status = fit$status, starts = fit$starts,
       basis = basis, n_basis = ncol(B))
}
