## flash_greedy_init_smooth() and its matrix helpers (R/07_flash_init.R).

# n x p residuals: rank one with a smooth loading, plus noise; each column observed
# at m random rows (m = n: complete data).
sim_resid <- function(n = 40, p = 30, m = n, noise = 0.3, seed = 1) {
  set.seed(seed)
  tt <- seq(0, 1, length.out = n)
  R <- outer(sin(2 * pi * tt) + tt, rnorm(p, 0, 2)) + matrix(rnorm(n * p, 0, noise), n, p)
  if (m < n) for (j in seq_len(p)) R[-sample(n, m), j] <- NA
  list(R = R, t = tt)
}

in_span <- function(v, B) max(abs(qr.resid(qr(B), v))) <= 1e-8 * max(1, max(abs(v)))
cosine <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))

test_that("basis = \"constant\" gives the column means of the observed residuals", {
  d <- sim_resid(m = 7)
  R <- d$R
  R[, 5] <- NA                                     # a column with no observed entry
  out <- .smooth_init_matrix(R, d$t, basis = "constant")
  expect_identical(out$status, "exact")
  expect_equal(max(out$row_values) - min(out$row_values), 0)
  expect_gt(out$row_values[1], 0)
  fitted <- tcrossprod(out$row_values, out$column_values)
  cm <- colMeans(R, na.rm = TRUE)
  cm[is.nan(cm)] <- 0
  expect_equal(fitted[1, ], cm, tolerance = 1e-12)
  expect_identical(out$column_values[5], 0)
  expect_equal(out$rss, sum((R - fitted)^2, na.rm = TRUE), tolerance = 1e-12)
})

test_that("the starting curve lies in the requested basis, and df counts the constant", {
  d <- sim_resid(m = 8)
  x <- d$t
  for (df in c(2L, 3L, 4L, 6L)) {
    B <- .smooth_init_basis(x, "ns", df)
    expect_identical(ncol(B), df)
    expect_true(in_span(rep(1, length(x)), B))    # the constant is included
    out <- .smooth_init_matrix(d$R, x, basis = "ns", df = df)
    expect_true(in_span(out$row_values, B))
    expect_identical(out$n_basis, df)
  }
  # df = 2 is the linear space
  expect_true(in_span(x, .smooth_init_basis(x, "ns", 2L)))
  out <- .smooth_init_matrix(d$R, x, basis = "linear")
  expect_true(in_span(out$row_values, cbind(1, x)))
  expect_identical(out$n_basis, 2L)
  # the basis is splines::ns() with knots at the quantiles of the distinct x
  B4 <- splines::ns(x, df = 4, intercept = TRUE)
  expect_equal(qr.Q(qr(.smooth_init_basis(x, "ns", 4L))) %*% t(qr.Q(qr(.smooth_init_basis(x, "ns", 4L)))),
               qr.Q(qr(B4)) %*% t(qr.Q(qr(B4))), tolerance = 1e-10)
  # ... and repeating a coordinate does not move the knots
  xr <- c(x, rep(x[2], 20))
  expect_equal(.smooth_init_basis(xr, "ns", 4L)[seq_along(x), ], .smooth_init_basis(x, "ns", 4L),
               tolerance = 1e-12)
})

test_that("complete data: smooth ALS reaches the rank-one SVD of Q'R", {
  d <- sim_resid(noise = 0.5)
  for (basis in c("ns", "linear")) {
    B <- .smooth_init_basis(d$t, basis, 5L)
    Q <- qr.Q(qr(B))
    sv <- svd(crossprod(Q, d$R))
    expect_gt(sv$d[1] / sv$d[2], 3)
    best <- sv$d[1] * tcrossprod(Q %*% sv$u[, 1], sv$v[, 1])
    out <- .smooth_init_matrix(d$R, d$t, basis = basis, df = 5L, tol = 1e-13, maxiter = 1000L)
    fitted <- tcrossprod(out$row_values, out$column_values)
    expect_lt(max(abs(fitted - best)) / max(abs(best)), 1e-5)
    expect_equal(out$rss, sum((d$R - best)^2), tolerance = 1e-8)
    expect_identical(out$status, "converged")
  }
})

test_that("the aggregated loading update equals least squares on the stacked observations", {
  d <- sim_resid(n = 25, p = 12, m = 5)
  Z <- !is.na(d$R)
  R0 <- d$R
  R0[!Z] <- 0
  Q <- qr.Q(qr(.smooth_init_basis(d$t, "ns", 4L)))
  f <- rnorm(ncol(d$R))
  a <- .smooth_init_alpha_step(R0, Z * 1, Q, f)
  obs <- which(Z, arr.ind = TRUE)
  X <- Q[obs[, 1], ] * f[obs[, 2]]
  a_direct <- qr.coef(qr(X), d$R[obs])
  expect_equal(a, as.numeric(a_direct), tolerance = 1e-10)
  # the score update is least squares column by column
  l <- as.numeric(Q %*% a)
  fs <- .smooth_init_fstep(R0, Z * 1, l)
  for (j in 1:3) {
    i <- which(Z[, j])
    expect_equal(fs[j], sum(l[i] * d$R[i, j]) / sum(l[i]^2), tolerance = 1e-12)
  }
  # a column whose observed loadings are all at rounding level gets a zero score
  l <- c(1, rep(1.5e-14, 4))
  Zc <- cbind(c(0, 1, 1, 1, 1), c(1, 1, 0, 0, 0))
  fs <- .smooth_init_fstep(Zc * 2, Zc, l)
  expect_identical(fs[1], 0)
  expect_equal(fs[2], 2 * (1 + 1.5e-14) / (1 + 1.5e-14^2))
})

test_that("missing entries are left out, whatever their placeholder; observed zeros count", {
  d <- sim_resid(m = 6)
  Z <- !is.na(d$R)
  Q <- qr.Q(qr(.smooth_init_basis(d$t, "ns", 4L)))
  A0 <- .smooth_init_starts(ncol(Q), 3L, 666)
  fits <- lapply(c(0, 999, NA), function(ph) {
    Rp <- d$R
    Rp[!Z] <- ph
    .smooth_init_als(Rp, Z, Q, A0, 100L, 1e-6)
  })
  expect_identical(fits[[1]]$l, fits[[2]]$l)
  expect_identical(fits[[1]]$f, fits[[3]]$f)
  expect_identical(fits[[1]]$rss, fits[[2]]$rss)
  # turning observed values into zeros changes the fit; turning them into NA too
  R_zero <- d$R
  R_zero[which(Z[, 1])[1:3], 1] <- 0
  R_na <- d$R
  R_na[which(Z[, 1])[1:3], 1] <- NA
  base <- .smooth_init_matrix(d$R, d$t)
  expect_false(isTRUE(all.equal(.smooth_init_matrix(R_zero, d$t)$column_values, base$column_values)))
  expect_false(isTRUE(all.equal(.smooth_init_matrix(R_zero, d$t)$column_values,
                                .smooth_init_matrix(R_na, d$t)$column_values)))
})

test_that("sparse irregular designs: unobserved columns, single visits, unobserved rows", {
  d <- sim_resid(n = 60, p = 80, m = 4, noise = 0.2)
  R <- d$R
  R[, 3] <- NA                                     # no observation
  for (j in 4:10) R[which(!is.na(R[, j]))[-1], j] <- NA   # observed once
  R[17, ] <- NA                                    # an unobserved time point
  out <- .smooth_init_matrix(R, d$t)
  B <- .smooth_init_basis(d$t, "ns", 4L)
  expect_true(all(is.finite(out$row_values)) && all(is.finite(out$column_values)))
  expect_identical(out$column_values[3], 0)
  expect_true(in_span(out$row_values, B))          # the unobserved row follows the curve
  expect_gt(abs(cosine(out$row_values, sin(2 * pi * d$t) + d$t)), 0.99)
  expect_lte(out$rss, out$tss)
  expect_equal(out$rss, sum((R - tcrossprod(out$row_values, out$column_values))^2, na.rm = TRUE),
               tolerance = 1e-10)
  # here one start stalls far from the others, and the best one is kept
  expect_identical(out$rss, min(out$starts$rss))
  expect_gt(max(out$starts$rss), 10 * out$rss)
  # every column observed once: any curve that is nonzero there fits exactly;
  # only the numerics and the basis constraint are checked
  R1 <- d$R
  for (j in seq_len(ncol(R1))) R1[which(!is.na(R1[, j]))[-1], j] <- NA
  out1 <- .smooth_init_matrix(R1, d$t)
  expect_true(all(is.finite(out1$column_values)))
  expect_true(in_span(out1$row_values, B))
  expect_lt(out1$rss, 1e-8 * out1$tss)
})

test_that("rank-deficient updates get a minimum-norm solution, without a ridge", {
  d <- sim_resid(n = 30, p = 20)
  R <- d$R
  R[-c(5, 20), ] <- NA                             # only two distinct times observed, df = 4
  out <- .smooth_init_matrix(R, d$t, df = 4L, tol = 1e-12, maxiter = 1000L)
  expect_true(all(is.finite(out$row_values)))
  expect_true(in_span(out$row_values, .smooth_init_basis(d$t, "ns", 4L)))
  # the two observed rows are fit as well as any rank-one fit of that 2 x p block
  sv <- svd(R[c(5, 20), ])
  expect_equal(out$rss, sum(sv$d[-1]^2), tolerance = 1e-6)
  # minimum norm: the coefficients lie in the row space of the observed rows' design,
  # which fixes the curve at the 28 unobserved rows
  Q <- qr.Q(qr(.smooth_init_basis(d$t, "ns", 4L)))
  expect_lt(max(abs(qr.resid(qr(t(Q[c(5, 20), ])), crossprod(Q, out$row_values)))), 1e-8)
})

test_that("all-zero residuals give zero vectors", {
  R <- matrix(0, 10, 6)
  R[1:3, 2] <- NA
  for (basis in c("ns", "linear", "constant")) {
    out <- .smooth_init_matrix(R, seq_len(10), basis = basis)
    expect_identical(out$status, "zero")
    expect_true(all(out$row_values == 0) && all(out$column_values == 0))
  }
})

test_that("repeated and unsorted coordinates: equal curve values, positions preserved", {
  d <- sim_resid(n = 30, p = 25, m = 6)
  x <- d$t
  x[c(4, 9)] <- x[2]                               # repeated coordinates
  out <- .smooth_init_matrix(d$R, x, tol = 1e-12, maxiter = 1000L)
  expect_equal(out$row_values[c(4, 9)], rep(out$row_values[2], 2), tolerance = 1e-12)
  # a row permutation rotates the orthonormal basis, so the random starts differ;
  # both runs are converged to the same optimum
  perm <- sample(30)
  out_p <- .smooth_init_matrix(d$R[perm, ], x[perm], tol = 1e-12, maxiter = 1000L)
  expect_equal(abs(cosine(out_p$row_values, out$row_values[perm])), 1, tolerance = 1e-6)
  expect_equal(tcrossprod(out_p$row_values, out_p$column_values),
               tcrossprod(out$row_values, out$column_values)[perm, ], tolerance = 1e-6)
})

test_that("smooth_dim = 2 smooths the columns and returns the vectors in matrix order", {
  d <- sim_resid(m = 7)
  a <- .smooth_init_matrix(d$R, d$t, smooth_dim = 1L)
  b <- .smooth_init_matrix(t(d$R), d$t, smooth_dim = 2L)
  expect_length(b$row_values, ncol(d$R))
  expect_length(b$column_values, nrow(d$R))
  expect_identical(b$column_values, a$row_values)
  expect_identical(b$row_values, a$column_values)
  expect_error(.smooth_init_matrix(t(d$R), d$t, smooth_dim = 1L), "one per row")
})

test_that("the best start is kept, the seed is local and the result is reproducible", {
  d <- sim_resid(m = 5)
  set.seed(42)
  before <- .Random.seed
  out1 <- .smooth_init_matrix(d$R, d$t, nstarts = 5L, seed = 11)
  expect_identical(.Random.seed, before)           # caller's RNG state untouched
  out2 <- .smooth_init_matrix(d$R, d$t, nstarts = 5L, seed = 11)
  expect_identical(out1$row_values, out2$row_values)
  expect_identical(out1$column_values, out2$column_values)
  expect_identical(out1$rss, min(out1$starts$rss))
  expect_identical(nrow(out1$starts), 5L)
  # seed = NULL draws from the caller's stream
  set.seed(42)
  .smooth_init_matrix(d$R, d$t, seed = NULL)
  expect_false(identical(.Random.seed, before))
  # also when there was no RNG state before
  rm(".Random.seed", envir = globalenv())
  .smooth_init_matrix(d$R, d$t, seed = 1)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(42)
})

test_that("scale and sign are normalized without changing the fit", {
  d <- sim_resid(m = 6)
  Z <- !is.na(d$R)
  Q <- qr.Q(qr(.smooth_init_basis(d$t, "ns", 4L)))
  raw <- .smooth_init_als(d$R, Z, Q, .smooth_init_starts(ncol(Q), 10L, 666), 100L, 1e-6)
  out <- .smooth_init_matrix(d$R, d$t)
  expect_equal(tcrossprod(out$row_values, out$column_values), tcrossprod(raw$l, raw$f),
               tolerance = 1e-12)
  expect_equal(sum(out$row_values^2), sum(out$column_values^2), tolerance = 1e-12)
  expect_gt(out$row_values[which.max(abs(out$row_values))], 0)
})

test_that("invalid inputs are rejected with clear errors", {
  d <- sim_resid(m = 6)
  expect_error(.smooth_init_matrix(matrix(NA_real_, 5, 4), 1:5), "no observed entries")
  expect_error(.smooth_init_matrix(d$R, d$t[-1]), "one per row")
  expect_error(.smooth_init_matrix(d$R, replace(d$t, 3, NA)), "finite coordinates")
  expect_error(.smooth_init_matrix(d$R, d$t, smooth_dim = 3), "`smooth_dim` must be 1 or 2")
  expect_error(.smooth_init_matrix(d$R, d$t, df = 1), "`df` must be an integer >= 2")
  expect_error(.smooth_init_matrix(d$R, d$t, df = 4.5), "`df` must be an integer")
  expect_error(.smooth_init_matrix(d$R, rep(1:3, length.out = 40), df = 4), "at least 4 distinct")
  expect_error(.smooth_init_matrix(d$R, rep(1, 40), basis = "linear"), "at least 2 distinct")
  expect_error(.smooth_init_matrix(d$R, d$t, nstarts = 0), "`nstarts`")
  expect_error(.smooth_init_matrix(d$R, d$t, tol = -1), "`tol`")
  expect_error(.smooth_init_matrix(d$R, d$t, seed = "a"), "`seed`")
  expect_error(.smooth_init_matrix(array(0, c(4, 3, 2)), 1:4), "ordinary numeric matrix")
  expect_error(.smooth_init_matrix(Matrix::Matrix(1, 4, 3), 1:4), "ordinary numeric matrix")
  R <- d$R
  R[2, 2] <- Inf
  expect_error(.smooth_init_matrix(R, d$t), "must be finite")
  expect_error(flash_greedy_init_smooth(d$R, d$t), "flash_fit object")
})

test_that("low-rank and sparse Matrix data are rejected before any residuals are formed", {
  skip_if_not_installed("flashier")
  set.seed(2)
  M <- matrix(rnorm(200), 20, 10)
  fl_lr <- flashier::flash_init(svd(M), var_type = 0)
  expect_error(flash_greedy_init_smooth(fl_lr$flash_fit, x = 1:20), "low-rank data are not supported")
  fl_sp <- flashier::flash_init(Matrix::rsparsematrix(20, 10, 0.3), var_type = 0)
  expect_error(flash_greedy_init_smooth(fl_sp, x = 1:20), "\"dgCMatrix\"")
})

test_that("flash_greedy() runs with the smooth initialization and recovers the loading", {
  skip_if_not_installed("flashier")
  set.seed(3)
  n <- 40; p <- 60
  tt <- seq(0, 1, length.out = n)
  l1 <- sin(2 * pi * tt)
  Y <- outer(l1, rnorm(p, 0, 3)) + matrix(rnorm(n * p, 0, 0.5), n, p)
  for (j in seq_len(p)) Y[-sample(n, 5), j] <- NA
  fl0 <- flashier::flash_init(Y, var_type = 0)
  init <- function(f) flash_greedy_init_smooth(f, x = tt)
  # flashier's residuals have NA exactly at the missing entries
  expect_identical(is.na(stats::residuals(fl0$flash_fit)), is.na(Y))
  # the adapter works on the flash_fit that flashier passes, and on a flash object,
  # and passes every setting on
  expect_identical(init(fl0$flash_fit), init(fl0))
  ref <- .smooth_init_matrix(Y, tt)
  expect_identical(init(fl0$flash_fit), list(ref$row_values, ref$column_values))
  for (a in list(list(basis = "linear"), list(basis = "constant"), list(df = 6L),
                 list(nstarts = 1L, seed = 7L, maxiter = 5L, tol = 1e-3))) {
    r <- do.call(.smooth_init_matrix, c(list(Y, tt), a))
    expect_identical(do.call(flash_greedy_init_smooth, c(list(fl0$flash_fit, x = tt), a)),
                     list(r$row_values, r$column_values))
  }
  # ... also at later greedy steps, after factors have been added
  na_ok <- logical(0)
  init_chk <- function(f) {
    na_ok <<- c(na_ok, identical(is.na(stats::residuals(f)), is.na(Y)))
    init(f)
  }
  priors <- list(Matern = ebnm_Matern_generator(setup = Matern_setup(tt, alpha = 2)),
                 LGP = ebnm_LGP_generator(LGP_setup(tt, num_knots = 20, betaprec = -1)))
  for (pr in names(priors)) {
    fl <- suppressWarnings(flashier::flash_greedy(fl0, Kmax = 2, init_fn = init_chk, verbose = 0,
                                                  ebnm_fn = list(priors[[pr]], ebnm::ebnm_point_normal)))
    expect_gte(fl$n_factors, 1)
    fl <- suppressWarnings(flashier::flash_backfit(fl, maxiter = 50, verbose = 0))
    expect_gt(max(abs(apply(fl$L_pm, 2, cosine, l1))), 0.98)
  }
  expect_gte(length(na_ok), 4)
  expect_true(all(na_ok))
  # smoothing prior on the columns: time runs along the columns of t(Y)
  init_t <- function(f) flash_greedy_init_smooth(f, x = tt, smooth_dim = 2L)
  flt <- suppressWarnings(flashier::flash_init(t(Y), var_type = 0) |>
    flashier::flash_greedy(Kmax = 1, init_fn = init_t, verbose = 0,
                           ebnm_fn = list(ebnm::ebnm_point_normal, priors$Matern)))
  expect_equal(flt$n_factors, 1)
  expect_gt(abs(cosine(flt$F_pm[, 1], l1)), 0.98)
  # residuals with nothing left give a zero start, and greedy stops without error
  zero_fit <- suppressWarnings(flashier::flash_init(matrix(0, n, p), var_type = 0) |>
    flashier::flash_greedy(Kmax = 1, init_fn = init, verbose = 0))
  expect_equal(zero_fit$n_factors, 0)
})
