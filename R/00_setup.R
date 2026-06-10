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

.inverse_softplus_stable <- function(y) {
  y_dim <- dim(y)
  y <- as.numeric(y)
  if (anyNA(y) || any(!is.finite(y)) || any(y <= 0)) {
    stop("`y` must contain positive finite values.")
  }
  out <- numeric(length(y))
  large <- y > 40
  out[large] <- y[large]
  out[!large] <- log(expm1(y[!large]))
  if (!is.null(y_dim)) dim(out) <- y_dim
  out
}

.positive_response_floor <- function(x, s = NULL) {
  x <- as.numeric(x)
  positive_x <- x[is.finite(x) & x > 0]
  floor0 <- if (length(positive_x)) {
    min(positive_x) / 2
  } else if (!is.null(s)) {
    positive_s <- as.numeric(s)[is.finite(as.numeric(s)) & as.numeric(s) > 0]
    if (length(positive_s)) min(positive_s) / 10 else 1e-6
  } else {
    1e-6
  }
  if (!is.finite(floor0) || floor0 <= 0) floor0 <- 1e-6
  as.numeric(floor0)
}

.response_mean_for_init <- function(x, s = NULL) {
  x <- as.numeric(x)
  if (is.null(s)) {
    return(as.numeric(mean(x)))
  }
  w <- 1 / (as.numeric(s)^2)
  as.numeric(stats::weighted.mean(x, w = w))
}

.response_eta_from_mean <- function(mu, x, s = NULL, link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  mu <- as.numeric(mu)[1]
  if (identical(link, "identity")) {
    return(mu)
  }
  mu <- max(mu, .positive_response_floor(x, s = s))
  if (identical(link, "log")) {
    return(as.numeric(log(mu)))
  }
  as.numeric(.inverse_softplus_stable(mu))
}

.response_beta_init <- function(x, s = NULL, link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  .response_eta_from_mean(
    mu = .response_mean_for_init(x, s = s),
    x = x,
    s = s,
    link = link
  )
}

.response_mean_from_eta <- function(eta, link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  eta <- as.numeric(eta)
  if (identical(link, "identity")) {
    return(eta)
  }
  if (identical(link, "log")) {
    return(as.numeric(exp(eta)))
  }
  as.numeric(.softplus_stable(eta))
}

.response_raw_residual_scale <- function(x, eta, link = c("identity", "log", "softplus")) {
  link <- match.arg(link)
  x <- as.numeric(x)
  eta <- as.numeric(eta)
  if (length(eta) == 1L) eta <- rep(eta, length(x))
  if (length(eta) != length(x)) {
    stop("`eta` must have length 1 or length(x).")
  }
  resid <- x - .response_mean_from_eta(eta, link = link)
  out <- stats::sd(resid)
  if (!is.finite(out) || out <= 0) {
    out <- sqrt(mean(resid^2))
  }
  if (!is.finite(out) || out <= 0) {
    out <- stats::sd(x)
  }
  if (!is.finite(out) || out <= 0) {
    out <- 1
  }
  as.numeric(out)
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
