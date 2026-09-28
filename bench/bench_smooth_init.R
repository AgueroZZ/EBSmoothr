# Small comparison of greedy initializations for flashier with an EBSmoothr Matern prior:
# flashier's default (unrestricted rank-one ALS) vs flash_greedy_init_smooth() with the
# "ns" (df = 4), "linear" and "constant" bases, on dense and on sparse irregular data.
# Usage: Rscript bench/bench_smooth_init.R [n_seeds]   (from the package root)
# Writes bench/results_smooth_init.rds and prints a summary.
#
# Per design and seed it records, for the first greedy start on the raw data: the observed
# RSS of the start, the time the init takes, |cos| of the starting curve with the closest true
# loading, and the share of sum f0^2 on the largest column. It then runs greedy (Kmax = 6) ->
# backfit (200) -> nullcheck and records the final number of factors, the ELBO and the
# loading recovery (mean over true loadings of the best |cos| with a fitted loading).
# These are diagnostics. The unrestricted problem has the lower optimal start RSS, since its
# parameter space is larger, but on sparse data the default's ALS does not always reach it.

args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1) as.integer(args[1]) else 3L

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
  library(flashier)
})

n <- 60
tt <- seq(0, 1, length.out = n)
loadings <- cbind(sin(2 * pi * tt), exp(-(tt - 0.5)^2 / (2 * 0.1^2)), cos(2 * pi * tt) * tt)
loadings <- sweep(loadings, 2, sqrt(colSums(loadings^2)), "/")

# K true factors with point-normal scores; "sparse": column j observed at m_j random times,
# with m_j drawn from 1..10 (about 10% of the columns observed once).
simulate <- function(K, sparse, seed, p = 150, sd_f = c(12, 9, 7)) {
  set.seed(seed)
  L <- loadings[, seq_len(K), drop = FALSE]
  F <- sapply(seq_len(K), function(k) rbinom(p, 1, 0.6) * rnorm(p, 0, sd_f[k]))
  Y <- L %*% t(F) + matrix(rnorm(n * p), n, p)
  if (sparse) {
    for (j in seq_len(p)) Y[-sample(n, sample.int(10, 1)), j] <- NA
    stopifnot(all(rowSums(!is.na(Y)) > 0))   # flashier cannot fit an unobserved row
  }
  list(Y = Y, L = L, F = F)
}

inits <- list(
  default  = function(f) flash_greedy_init_default(f),
  ns4      = function(f) flash_greedy_init_smooth(f, x = tt),
  linear   = function(f) flash_greedy_init_smooth(f, x = tt, basis = "linear"),
  constant = function(f) flash_greedy_init_smooth(f, x = tt, basis = "constant")
)
prior_L <- ebnm_Matern_generator(setup = Matern_setup(tt, alpha = 2))
abscos <- function(a, b) abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2))

rows <- list()
for (design in c("dense", "sparse")) for (K in c(1, 3)) for (seed in seq_len(n_seeds)) {
  sim <- simulate(K, design == "sparse", seed)
  fl0 <- flash_init(sim$Y, var_type = 0)
  R <- stats::residuals(fl0)
  for (nm in names(inits)) {
    t0 <- proc.time()[["elapsed"]]
    e <- inits[[nm]](fl0$flash_fit)
    t_init <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    fl <- suppressWarnings(flash_init(sim$Y, var_type = 0) |>
      flash_greedy(Kmax = 6, ebnm_fn = list(prior_L, ebnm::ebnm_point_normal), init_fn = inits[[nm]],
                   verbose = 0) |>
      flash_backfit(maxiter = 200, verbose = 0) |>
      flash_nullcheck(verbose = 0))
    t_fit <- proc.time()[["elapsed"]] - t0
    rec <- if (fl$n_factors > 0) {
      mean(apply(sim$L, 2, function(l) max(apply(fl$L_pm, 2, abscos, l))))
    } else 0
    rows[[length(rows) + 1]] <- data.frame(
      design = design, K = K, seed = seed, init = nm,
      start_rss = sum((R - tcrossprod(e[[1]], e[[2]]))^2, na.rm = TRUE),
      start_secs = t_init,
      start_cos = max(apply(sim$L, 2, abscos, e[[1]])),
      start_top_share = max(e[[2]]^2) / sum(e[[2]]^2),
      K_final = fl$n_factors, elbo = fl$elbo, recovery = rec, fit_secs = t_fit)
  }
}
res <- do.call(rbind, rows)
res$init <- factor(res$init, levels = names(inits))
saveRDS(res, "bench/results_smooth_init.rds")

# ELBO relative to the default init of the same data set
res$elbo_gain <- res$elbo - ave(ifelse(res$init == "default", res$elbo, NA), res$design, res$K, res$seed,
                                FUN = function(v) v[!is.na(v)][1])
agg <- aggregate(cbind(start_rss, start_secs, start_cos, start_top_share, K_final, elbo_gain, recovery,
                       fit_secs) ~ init + K + design, res, mean)
num <- vapply(agg, is.numeric, logical(1))
agg[num] <- lapply(agg[num], signif, 4)
cat(sprintf("means over %d seeds; true K = 1 or 3; elbo_gain = ELBO - ELBO(default init)\n", n_seeds))
print(agg, row.names = FALSE)
