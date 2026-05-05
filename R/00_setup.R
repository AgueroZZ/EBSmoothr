#' Roxygen commands
#'
#' @useDynLib EBSmoothr
#' @importFrom methods as
#' @importFrom stats nlminb optim
#'
dummy <- function() {
  return(NULL)
}

.softplus_stable <- function(x) {
  x_dim <- dim(x)
  x <- as.numeric(x)
  out <- numeric(length(x))
  pos <- x > 0
  out[pos] <- x[pos] + log1p(exp(-x[pos]))
  out[!pos] <- log1p(exp(x[!pos]))
  if (!is.null(x_dim)) dim(out) <- x_dim
  out
}

.sigmoid_stable <- function(x) {
  x_dim <- dim(x)
  x <- as.numeric(x)
  out <- numeric(length(x))
  pos <- x >= 0
  out[pos] <- 1 / (1 + exp(-x[pos]))
  exp_x <- exp(x[!pos])
  out[!pos] <- exp_x / (1 + exp_x)
  if (!is.null(x_dim)) dim(out) <- x_dim
  out
}

.gauss_hermite_rule_cache <- new.env(parent = emptyenv())

.gauss_hermite_rule <- function(n = 35L) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("`n` must be an integer greater than 1.")
  }
  key <- as.character(n)
  cached <- .gauss_hermite_rule_cache[[key]]
  if (!is.null(cached)) return(cached)

  J <- matrix(0, nrow = n, ncol = n)
  offdiag <- sqrt(seq_len(n - 1L) / 2)
  J[cbind(seq_len(n - 1L), seq.int(2L, n))] <- offdiag
  J[cbind(seq.int(2L, n), seq_len(n - 1L))] <- offdiag
  eig <- eigen(J, symmetric = TRUE)
  ord <- order(eig$values)
  rule <- list(
    nodes = as.numeric(eig$values[ord]),
    weights = as.numeric(sqrt(pi) * eig$vectors[1L, ord]^2)
  )
  .gauss_hermite_rule_cache[[key]] <- rule
  rule
}

.softplus_gaussian_moments <- function(eta_mean, eta_var, quadrature_n = 35L) {
  eta_mean <- as.numeric(eta_mean)
  eta_var <- pmax(as.numeric(eta_var), 0)
  if (length(eta_mean) != length(eta_var)) {
    stop("`eta_mean` and `eta_var` must have the same length.")
  }

  mean <- .softplus_stable(eta_mean)
  var <- numeric(length(eta_mean))
  positive_var <- eta_var > 0
  if (!any(positive_var)) {
    return(list(mean = mean, var = var))
  }

  rule <- .gauss_hermite_rule(quadrature_n)
  normalizer <- sqrt(pi)
  for (i in which(positive_var)) {
    eta_draws <- eta_mean[i] + sqrt(2 * eta_var[i]) * rule$nodes
    values <- .softplus_stable(eta_draws)
    m1 <- sum(rule$weights * values) / normalizer
    m2 <- sum(rule$weights * values^2) / normalizer
    mean[i] <- m1
    var[i] <- max(m2 - m1^2, 0)
  }

  list(mean = as.numeric(mean), var = as.numeric(var))
}
