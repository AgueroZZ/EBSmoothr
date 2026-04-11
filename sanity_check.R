### A simple example to check LGP works
library(EBSmoothr)

# identity link case:
f <- function(x) sin(x)
x <- seq(0, 10, length.out = 100)
sd_vec <- sample(c(0.1, 0.5, 1), length(x), replace = TRUE)
y <- f(x) + rnorm(length(x), sd = sd_vec)
LGP_setup <- EBSmoothr::LGP_setup(t = x, p = 2, betaprec = -1)
LGP_fun <- EBSmoothr::ebnm_LGP_generator(LGP_setup)
LGP_fun_old <- EBMFSmooth::ebnm_LGP_generator(LGP_setup)
LGP_res_old <- LGP_fun_old(y, sd_vec)
LGP_res <- LGP_fun(y, sd_vec)
plot(LGP_res$posterior$mean ~ x, type = "l", col = "blue")
points(x, y, col = "black", pch = 1)
polygon(c(x, rev(x)), c(LGP_res$posterior$mean + 2 * sqrt(LGP_res$posterior$var),
                      rev(LGP_res$posterior$mean - 2 * sqrt(LGP_res$posterior$var))),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(x, f(x), col = "red", lty = 2)
LGP_res$log_likelihood - LGP_res_old$log_likelihood
LGP_res$log_likelihood - LGP_res$log_likelihood_stepB_laplace
max(abs(LGP_res$posterior$mean - LGP_res_old$posterior$mean))
max(abs(sqrt(LGP_res$posterior$var) - sqrt(LGP_res_old$posterior$var)))

# log link case:
y <- exp(f(x)) + rnorm(length(x), sd = sd_vec)
LGP_setup <- EBSmoothr::LGP_setup(t = x, p = 2, betaprec = -1, link = "log")
LGP_fun <- EBSmoothr::ebnm_LGP_generator(LGP_setup, link = "log")
LGP_fun_old <- EBMFSmooth::ebnm_LGP_generator(LGP_setup, link = "log")
LGP_res_old <- LGP_fun_old(y, sd_vec)
LGP_res <- LGP_fun(y, sd_vec)
plot(LGP_res$posterior$mean ~ x, type = "l", col = "blue")
points(x, y, col = "black", pch = 1)
polygon(c(x, rev(x)), c(LGP_res$posterior$mean + 2 * sqrt(LGP_res$posterior$var),
                      rev(LGP_res$posterior$mean - 2 * sqrt(LGP_res$posterior$var))),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(x, exp(f(x)), col = "red", lty = 2)

LGP_res$log_likelihood - LGP_res_old$log_likelihood
LGP_res$log_likelihood - LGP_res$log_likelihood_stepB_laplace
max(abs(LGP_res$posterior$mean - LGP_res_old$posterior$mean))
max(abs(sqrt(LGP_res$posterior$var) - sqrt(LGP_res_old$posterior$var)))

