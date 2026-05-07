if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all("EBSmoothr", quiet = TRUE, export_all = FALSE)
} else {
  library(EBSmoothr)
}

# Exact identity-link Matern smoother in one dimension
set.seed(1)

truth_1d <- function(x) 0.8 + sin(1.5 * x)
x_1d <- seq(0, 10, length.out = 120)
s_1d <- sample(c(0.1, 0.2, 0.35), length(x_1d), replace = TRUE)
y_1d <- truth_1d(x_1d) + rnorm(length(x_1d), sd = s_1d)

fit_fun_1d <- ebnm_Matern_generator(locations = x_1d)
fit_1d <- fit_fun_1d(y_1d, s_1d)

plot(fit_1d$posterior$mean ~ x_1d, type = "l", col = "blue")
points(x_1d, y_1d, col = "black", pch = 1)
polygon(
  c(x_1d, rev(x_1d)),
  c(
    fit_1d$posterior$mean + 2 * sqrt(fit_1d$posterior$var),
    rev(fit_1d$posterior$mean - 2 * sqrt(fit_1d$posterior$var))
  ),
  col = rgb(0, 0, 1, 0.2),
  border = NA
)
lines(x_1d, truth_1d(x_1d), col = "red", lty = 2)
fit_1d$log_likelihood
fit_1d$fitted_g
fit_1d$fitted_beta
head(fit_1d$posterior_spatial_field)


# Fixed hyperparameters and fixed intercept
fit_1d_fixed <- fit_fun_1d(
  y_1d,
  s_1d,
  g_init = Matern(theta = log(1.5), sigma = 0.8),
  fix_g = TRUE,
  beta_fixed = 0.8
)

fit_1d_fixed$log_likelihood
fit_1d_fixed$fitted_g
fit_1d_fixed$fitted_beta


# Exact identity-link Matern smoother in two dimensions
set.seed(2)

grid_2d <- expand.grid(
  x = seq(0, 1, length.out = 12),
  y = seq(0, 1, length.out = 10)
)
locations_2d <- as.matrix(grid_2d)
truth_2d <- with(grid_2d, sin(2 * pi * x) * cos(2 * pi * y))
s_2d <- rep(0.15, nrow(locations_2d))
y_2d <- truth_2d + rnorm(nrow(locations_2d), sd = s_2d)

fit_fun_2d <- ebnm_Matern_generator(
  locations = locations_2d,
  max.edge = c(0.15, 0.3)
)
fit_2d <- fit_fun_2d(y_2d, s_2d)

image(
  x = sort(unique(grid_2d$x)),
  y = sort(unique(grid_2d$y)),
  z = matrix(fit_2d$posterior$mean, nrow = 12, ncol = 10),
  main = "Exact Matern posterior mean",
  xlab = "x",
  ylab = "y",
  col = hcl.colors(20, "BluYl")
)
fit_2d$log_likelihood


# ---- inlabru cross-check on softplus link ----
if (requireNamespace("inlabru", quietly = TRUE) &&
    requireNamespace("INLA", quietly = TRUE)) {
  set.seed(42)
  loc_sp <- seq(0, 1, length.out = 100)
  s_sp <- rep(0.08, length(loc_sp))
  eta_sp <- -0.4 + 1.2 * sin(2 * pi * loc_sp)
  x_sp <- log1p(exp(eta_sp)) + rnorm(length(loc_sp), sd = s_sp)

  fit_lap <- ebnm_Matern_generator(
    locations = loc_sp, link = "softplus", backend = "laplace"
  )(x_sp, s_sp)
  fit_bru <- ebnm_Matern_generator(
    locations = loc_sp, link = "softplus", backend = "inlabru"
  )(x_sp, s_sp)

  cat("\nSoftplus posterior mean: laplace vs inlabru (head)\n")
  print(cbind(
    laplace = head(fit_lap$posterior$mean),
    inlabru = head(fit_bru$posterior$mean)
  ))
  cat("\nRange of (laplace - inlabru) posterior mean:",
      paste(format(range(fit_lap$posterior$mean - fit_bru$posterior$mean),
                   digits = 3), collapse = ", "),
      "\n")
  cat("Range of (laplace - inlabru) posterior var:",
      paste(format(range(fit_lap$posterior$var - fit_bru$posterior$var),
                   digits = 3), collapse = ", "),
      "\n")
  cat("Hyperparameters (theta, sigma, beta):\n")
  print(rbind(
    laplace = c(theta = fit_lap$fitted_g$theta,
                sigma = fit_lap$fitted_g$sigma,
                beta  = fit_lap$fitted_beta),
    inlabru = c(theta = fit_bru$fitted_g$theta,
                sigma = fit_bru$fitted_g$sigma,
                beta  = fit_bru$fitted_beta)
  ))

  # Cross-check eb_smoother (s = NULL) path
  pc_sp <- list(range = c(0.1, 0.5), sigma = c(0.35, 0.5), noise = c(0.2, 0.5))
  fit_lap_full <- eb_smoother(
    x_sp, s = NULL, family = "matern", locations = loc_sp,
    link = "softplus", backend = "laplace", pc.penalty = pc_sp
  )
  fit_bru_full <- eb_smoother(
    x_sp, s = NULL, family = "matern", locations = loc_sp,
    link = "softplus", backend = "inlabru", pc.penalty = pc_sp
  )
  cat("\neb_smoother fitted_noise_sd: laplace =",
      format(fit_lap_full$fitted_noise_sd, digits = 4),
      "  inlabru =",
      format(fit_bru_full$fitted_noise_sd, digits = 4), "\n")

  if (nzchar(Sys.getenv("EBSMOOTHR_BENCH")) &&
      requireNamespace("microbenchmark", quietly = TRUE)) {
    cat("\nBenchmark (5 reps, n =", length(loc_sp), "):\n")
    print(microbenchmark::microbenchmark(
      laplace = ebnm_Matern_generator(
        locations = loc_sp, link = "softplus", backend = "laplace"
      )(x_sp, s_sp),
      inlabru = ebnm_Matern_generator(
        locations = loc_sp, link = "softplus", backend = "inlabru"
      )(x_sp, s_sp),
      times = 5L
    ))
  }
}
