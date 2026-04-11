# identity link case:
f <- function(x) sin(x)
x <- seq(0, 10, length.out = 100)
sd_vec <- sample(c(0.1, 0.5, 1), length(x), replace = TRUE)
y <- f(x) + rnorm(length(x), sd = sd_vec)

Matern_fun <- EBSmoothr::ebnm_Matern_generator(locations = x)
Matern_fun_old <- EBMFSmooth::ebnm_Matern_generator(locations = cbind(x))
Matern_res_old <- Matern_fun_old(y, sd_vec)
Matern_res <- Matern_fun(y, sd_vec)
plot(Matern_res$posterior$mean ~ x, type = "l", col = "blue")
points(x, y, col = "black", pch = 1)
polygon(c(x, rev(x)), c(Matern_res$posterior$mean + 2 * sqrt(Matern_res$posterior$var),
                        rev(Matern_res$posterior$mean - 2 * sqrt(Matern_res$posterior$var))),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(x, f(x), col = "red", lty = 2)
Matern_res$log_likelihood - Matern_res_old$log_likelihood
max(abs(Matern_res$posterior$mean - Matern_res_old$posterior$mean))
max(abs(sqrt(Matern_res$posterior$var) - sqrt(Matern_res_old$posterior$var)))

# log link case:
y <- exp(f(x)) + rnorm(length(x), sd = sd_vec)
Matern_fun <- EBSmoothr::ebnm_Matern_generator(locations = x, link = "log")
Matern_fun_old <- EBMFSmooth::ebnm_Matern_generator(locations = cbind(x), link = "log")
Matern_res_old <- Matern_fun_old(y, sd_vec)
Matern_res <- Matern_fun(y, sd_vec)
plot(Matern_res$posterior$mean ~ x, type = "l", col = "blue")
points(x, y, col = "black", pch = 1)
polygon(c(x, rev(x)), c(Matern_res$posterior$mean + 2 * sqrt(Matern_res$posterior$var),
                        rev(Matern_res$posterior$mean - 2 * sqrt(Matern_res$posterior$var))),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(x, exp(f(x)), col = "red", lty = 2)
Matern_res$log_likelihood - Matern_res_old$log_likelihood
max(abs(Matern_res$posterior$mean - Matern_res_old$posterior$mean))
max(abs(sqrt(Matern_res$posterior$var) - sqrt(Matern_res_old$posterior$var)))
