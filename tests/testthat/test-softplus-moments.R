test_that("softplus Gaussian moments match high-accuracy numerical integration", {
  reference_moments <- function(mu, eta_var) {
    if (eta_var <= 0) {
      value <- .softplus_stable(mu)
      return(c(mean = value, var = 0))
    }
    sd <- sqrt(eta_var)
    mean <- stats::integrate(
      function(z) .softplus_stable(mu + sd * z) * stats::dnorm(z),
      lower = -Inf,
      upper = Inf,
      rel.tol = 1e-12
    )$value
    second <- stats::integrate(
      function(z) .softplus_stable(mu + sd * z)^2 * stats::dnorm(z),
      lower = -Inf,
      upper = Inf,
      rel.tol = 1e-12
    )$value
    c(mean = mean, var = max(second - mean^2, 0))
  }

  grid <- expand.grid(
    eta_mean = c(-4, -2, 0, 2, 4),
    eta_var = c(0, 0.01, 0.1, 0.5, 1)
  )

  for (i in seq_len(nrow(grid))) {
    mu <- grid$eta_mean[i]
    eta_var <- grid$eta_var[i]
    reference <- reference_moments(mu, eta_var)
    actual <- .softplus_gaussian_moments(mu, eta_var)

    expect_equal(actual$mean, unname(reference["mean"]), tolerance = 1e-8)
    expect_equal(actual$var, unname(reference["var"]), tolerance = 1e-8)
  }
})
