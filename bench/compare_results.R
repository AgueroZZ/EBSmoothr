# Compare two bench result RDS files (reference vs candidate).
# Usage: Rscript bench/compare_results.R baseline after

args <- commandArgs(trailingOnly = TRUE)
ref_tag <- if (length(args) >= 1) args[[1]] else "baseline"
new_tag <- if (length(args) >= 2) args[[2]] else "after"

ref <- readRDS(file.path("bench", paste0("results_", ref_tag, ".rds")))
new <- readRDS(file.path("bench", paste0("results_", new_tag, ".rds")))

cat(sprintf("%-46s %9s %9s %8s  %10s %10s %10s\n",
            "case", ref_tag, new_tag, "speedup", "d.loglik", "d.postmean", "d.postvar"))
for (nm in intersect(names(ref), names(new))) {
  r <- ref[[nm]]$ref; n <- new[[nm]]$ref
  d_ll <- abs(r$log_likelihood - n$log_likelihood)
  d_pm <- max(abs(r$post_mean - n$post_mean))
  d_pv <- max(abs(r$post_var - n$post_var))
  cat(sprintf("%-46s %8.2fs %8.2fs %7.2fx  %10.3g %10.3g %10.3g\n",
              nm, ref[[nm]]$elapsed, new[[nm]]$elapsed,
              ref[[nm]]$elapsed / new[[nm]]$elapsed, d_ll, d_pm, d_pv))
}
