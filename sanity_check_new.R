### Sanity check for EBSmoothr::LGP smoothing
library(EBSmoothr)

set.seed(1)

# --- Identity link example ---
f <- function(x) sin(x)
t <- seq(0, 10, length.out = 100)
sd_vec <- sample(c(0.1, 0.5, 1), length(t), replace = TRUE)
y <- f(t) + rnorm(length(t), sd = sd_vec)

plot(t, y, main = "Identity link", xlab = "t", ylab = "y")

setup <- EBSmoothr::LGP_setup(
  t = t,
  p = 2,
  beta_mode = "fixed"   # new syntax (equivalent to betaprec = 0)
)

fit_fun <- EBSmoothr::ebnm_LGP_generator(setup, link = "identity")
fit <- fit_fun(y, sd_vec)

lines(t, fit$posterior$mean, col = "blue", lwd = 2)
polygon(
  c(t, rev(t)),
  c(fit$posterior$mean + 2 * sqrt(fit$posterior$var),
    rev(fit$posterior$mean - 2 * sqrt(fit$posterior$var))),
  col = rgb(0, 0, 1, 0.2), border = NA
)
lines(t, f(t), col = "red", lty = 2, lwd = 2)

# --- Log link example ---
y_log <- exp(f(t)) + rnorm(length(t), sd = sd_vec)

plot(t, y_log, main = "Log link", xlab = "t", ylab = "y")

setup_log <- EBSmoothr::LGP_setup(
  t = t,
  p = 2,
  beta_mode = "fixed"
)

fit_fun_log <- EBSmoothr::ebnm_LGP_generator(setup_log, link = "log")
fit_log <- fit_fun_log(y_log, sd_vec)

lines(t, fit_log$posterior$mean, col = "blue", lwd = 2)
polygon(
  c(t, rev(t)),
  c(fit_log$posterior$mean + 2 * sqrt(fit_log$posterior$var),
    rev(fit_log$posterior$mean - 2 * sqrt(fit_log$posterior$var))),
  col = rgb(0, 0, 1, 0.2), border = NA
)
lines(t, exp(f(t)), col = "red", lty = 2, lwd = 2)
