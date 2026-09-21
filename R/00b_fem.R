# ---- Native FEM assembly for the Matern SPDE -------------------------------
#
# The Matern SPDE (kappa^2 - Laplacian)^(alpha/2) (tau u) = W has a FEM
# (piecewise-linear basis) precision that, for integer alpha, is
#
#   Q = tau^2 * sum_{k=0}^{alpha} choose(alpha, k) kappa^(2 (alpha - k)) G_k
#
# with G_0 = C (the lumped mass matrix), G_1 = G (stiffness), and
# G_k = G C^{-1} G_{k-1}. This is the binomial expansion of the
# Lindgren-Rue-Lindstrom recursion Q_alpha = K C^{-1} Q_{alpha-2} C^{-1} K with
# K = kappa^2 C + G, and it reproduces INLA's `inla.spde2.matern` precision to
# machine precision for the alpha in {1, 2} that INLA's spde2 model supports.
#
# `fmesher` supplies C and the G_k directly (`fm_fem`), so nothing here needs
# INLA. `fm_fem` implements order <= 2 for 1-d meshes, so higher orders are
# extended with the recursion above; C is diagonal, which makes that cheap and
# exact.

.matern_max_alpha <- function() 8L

.check_matern_alpha <- function(alpha, d) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || !is.finite(alpha)) {
    stop("`alpha` must be a single finite numeric.")
  }
  if (abs(alpha - round(alpha)) > sqrt(.Machine$double.eps)) {
    stop(
      "`alpha` must be an integer. Fractional smoothness needs a rational SPDE ",
      "approximation, which EBSmoothr does not implement."
    )
  }
  alpha <- as.integer(round(alpha))
  if (alpha < 1L) stop("`alpha` must be a positive integer.")
  if (alpha > .matern_max_alpha()) {
    stop("`alpha` must be at most ", .matern_max_alpha(), ".")
  }
  if (alpha <= d / 2) {
    stop("`alpha` must satisfy alpha > d / 2.")
  }
  alpha
}

.as_csparse_sym <- function(M) {
  M <- Matrix::forceSymmetric(as(Matrix::Matrix(M, sparse = TRUE), "CsparseMatrix"))
  as(M, "CsparseMatrix")
}

# C, G_1, ..., G_alpha for a mesh. `fm_fem` returns c0 (lumped mass) plus
# g1..g_order; for 1-d meshes it caps `order` at 2, so anything above that is
# built here from G_k = G C^{-1} G_{k-1}.
.matern_fem_matrices <- function(mesh, alpha) {
  alpha <- as.integer(alpha)
  native_order <- if (inherits(mesh, "fm_mesh_1d") || inherits(mesh, "inla.mesh.1d")) {
    min(alpha, 2L)
  } else {
    alpha
  }

  fem <- suppressWarnings(fmesher::fm_fem(mesh, order = max(1L, native_order)))
  if (is.null(fem$c0) || is.null(fem$g1)) {
    stop("`fmesher::fm_fem()` did not return the expected `c0` and `g1` matrices.")
  }

  C <- .as_csparse_sym(fem$c0)
  c_diag <- Matrix::diag(C)
  if (any(!is.finite(c_diag)) || any(c_diag <= 0)) {
    stop("The FEM mass matrix has non-positive diagonal entries; the mesh is degenerate.")
  }
  Cinv <- Matrix::Diagonal(x = 1 / c_diag)
  G1 <- .as_csparse_sym(fem$g1)

  out <- vector("list", alpha + 1L)
  out[[1L]] <- C
  if (alpha >= 1L) out[[2L]] <- G1
  for (k in seq_len(alpha)[-1L]) {
    gk <- fem[[paste0("g", k)]]
    out[[k + 1L]] <- if (!is.null(gk)) {
      .as_csparse_sym(gk)
    } else {
      .as_csparse_sym(G1 %*% Cinv %*% out[[k]])
    }
  }
  names(out) <- paste0("G", seq_along(out) - 1L)
  out
}

# Align every G_k onto one shared sparsity pattern so that building a precision
# at new hyperparameters is a scalar combination of cached x-slots rather than
# a fresh sparse matrix sum.
.matern_fem_build <- function(mesh, alpha, d) {
  alpha <- .check_matern_alpha(alpha, d)
  G <- .matern_fem_matrices(mesh, alpha)

  # Union pattern via absolute values so no entry can cancel away.
  U <- G[[1L]]
  for (k in seq_along(G)[-1L]) U <- U + abs(G[[k]])
  U <- .as_csparse_sym(abs(U))

  Ut <- as(U, "TsparseMatrix")
  n <- nrow(U)
  key_u <- as.numeric(Ut@i) + as.numeric(Ut@j) * n

  align <- function(M) {
    Mt <- as(M, "TsparseMatrix")
    key_m <- as.numeric(Mt@i) + as.numeric(Mt@j) * n
    idx <- match(key_m, key_u)
    if (anyNA(idx)) stop("FEM component entry outside the union sparsity pattern.")
    x <- numeric(length(key_u))
    x[idx] <- Mt@x
    x
  }

  structure(
    list(
      alpha = alpha,
      d = as.integer(d),
      n.spde = as.integer(n),
      U = U,
      x = lapply(G, align),
      binom = choose(alpha, 0:alpha),
      matrices = G
    ),
    class = c("matern_fem", "list")
  )
}

.is_matern_fem <- function(x) inherits(x, "matern_fem")

# Q = tau^2 * sum_k choose(alpha, k) kappa^(2 (alpha - k)) G_k
.matern_fem_precision <- function(fem, kappa, tau) {
  if (!.is_matern_fem(fem)) stop("`fem` must be a `matern_fem` object.")
  alpha <- fem$alpha
  kappa2 <- kappa^2
  w <- fem$binom * kappa2^(alpha - (0:alpha))

  xs <- fem$x
  acc <- w[[1L]] * xs[[1L]]
  for (k in seq_along(xs)[-1L]) acc <- acc + w[[k]] * xs[[k]]

  Q <- fem$U
  Q@x <- tau^2 * acc
  # Matrix caches factorizations inside @factors by reference; drop anything
  # inherited from the shared pattern template so no stale factor can be
  # picked up for the new numeric values.
  Q@factors <- list()
  Q
}
