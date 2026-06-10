# Benchmark + reference-output harness for the 2D non-identity Matern path.
# Usage: Rscript bench/bench_2d_nonidentity.R <tag>
# Writes bench/results_<tag>.rds with timings and reference outputs.

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[[1]] else "baseline"

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
})

make_problem <- function(nx, ny, link, seed = 7) {
  set.seed(seed)
  grid <- expand.grid(x = seq(0, 1, length.out = nx), y = seq(0, 1, length.out = ny))
  loc <- as.matrix(grid)
  eta <- -0.3 + 1.1 * sin(2 * pi * grid$x) * cos(2 * pi * grid$y)
  mu <- if (link == "log") exp(eta) else log1p(exp(eta))
  s <- rep(0.1, nrow(loc))
  x <- mu + rnorm(nrow(loc), sd = s)
  list(loc = loc, x = x, s = s, eta = eta)
}

extract_ref <- function(fit) {
  list(
    log_likelihood = as.numeric(fit$log_likelihood),
    theta = as.numeric(fit$fitted_g$theta),
    sigma = as.numeric(fit$fitted_g$sigma),
    beta = as.numeric(fit$fitted_beta),
    fitted_noise_sd = if (is.null(fit$fitted_noise_sd)) NA_real_ else as.numeric(fit$fitted_noise_sd),
    post_mean = as.numeric(fit$posterior$mean),
    post_var = as.numeric(fit$posterior$var),
    backend = fit$backend
  )
}

run_case <- function(label, expr_fun) {
  cat(sprintf("== %-44s ", label))
  t <- system.time(fit <- expr_fun())
  cat(sprintf("%.2fs\n", t[["elapsed"]]))
  list(elapsed = t[["elapsed"]], ref = extract_ref(fit))
}

results <- list()

## ---- 2D problems ----
for (size in list(c(15L, 15L), c(22L, 22L), c(35L, 35L))) {
  nx <- size[1]; ny <- size[2]
  for (link in c("log", "softplus")) {
    prob <- make_problem(nx, ny, link)
    setup <- Matern_setup(prob$loc)
    key <- sprintf("2d_%dx%d_%s", nx, ny, link)
    cat(sprintf("-- problem %s: n_obs=%d n_spde=%d\n", key, nrow(prob$loc), ncol(setup$A)))

    fn_pql <- ebnm_Matern_generator(setup = setup, link = link, backend = "fisher_pql")
    results[[paste0(key, "_fisher_pql")]] <- run_case(
      paste0(key, " fisher_pql"),
      function() fn_pql(prob$x, prob$s)
    )

    fn_lf <- ebnm_Matern_generator(setup = setup, link = link, backend = "laplace_fisher")
    results[[paste0(key, "_laplace_fisher")]] <- run_case(
      paste0(key, " laplace_fisher"),
      function() fn_lf(prob$x, prob$s)
    )
  }
}

## ---- learned-noise 2D softplus (eb_smoother path) ----
prob <- make_problem(15L, 15L, "softplus")
results[["2d_15x15_softplus_learnnoise_pql"]] <- run_case(
  "2d_15x15 softplus learn-noise fisher_pql",
  function() eb_smoother(prob$x, s = NULL, family = "matern", locations = prob$loc,
                         link = "softplus", backend = "fisher_pql")
)

## ---- 1D identity exact (equivalence anchor for exact-path changes) ----
set.seed(11)
x1 <- seq(0, 10, length.out = 200)
s1 <- rep(0.15, length(x1))
y1 <- 0.8 + sin(1.5 * x1) + rnorm(length(x1), sd = s1)
fn_1d <- ebnm_Matern_generator(locations = x1)
results[["1d_identity_exact"]] <- run_case(
  "1d identity exact",
  function() fn_1d(y1, s1)
)

saveRDS(results, file.path("bench", paste0("results_", tag, ".rds")))
cat(sprintf("\nSaved bench/results_%s.rds\n", tag))
cat("\nTimings summary:\n")
for (nm in names(results)) cat(sprintf("  %-46s %8.2fs\n", nm, results[[nm]]$elapsed))
