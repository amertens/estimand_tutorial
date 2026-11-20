# 69_informative_ltfu.R -- supplementary scenario d2e_informative_ltfu (review 2026-11-19, item A12):
# loss to follow-up 5 times as fast once measured eGFR falls below 60, so censoring depends on the
# post-baseline kidney function that also drives the renal event. Every other parameter is the
# base case. Runs the one-feature estimator set ("deviation": the Cox analyses in both windows and
# SDR for treatment policy and sustained use) on 200 cohorts, as for the other one-feature
# scenarios. The truths do not depend on censoring, so they reproduce the base case's (a check).
#
#   Rscript --vanilla sandbox_rwe/run/69_informative_ltfu.R pilot   # truth, then 1 cohort, timed
#   Rscript --vanilla sandbox_rwe/run/69_informative_ltfu.R         # 200 cohorts, RWE_CORES workers
setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
for (f in list.files("sandbox_rwe/R", pattern = "[.]R$", full.names = TRUE)) source(f)
S <- "d2e_informative_ltfu"
cores <- as.integer(Sys.getenv("RWE_CORES", "4"))
pilot <- length(commandArgs(trailingOnly = TRUE)) && commandArgs(trailingOnly = TRUE)[1] == "pilot"
t0 <- Sys.time()
p <- get_scenario_params(S)                                  # writes results/params_<S>.rds
tf <- sprintf("sandbox_rwe/results/truth_%s.rds", S)
if (!file.exists(tf) || !identical(readRDS(tf)$fingerprint, params_fingerprint(p))) {
  saveRDS(scenario_truth(S, p), tf)
  cat(sprintf("truth in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
}
tr <- get_truth(S)
print(cbind(tr$strategies[, c("strategy", "rd", "rr")],
            base_rd = get_truth("base")$strategies$rd), digits = 4)

if (pilot) {
  t1 <- Sys.time()
  x <- run_replicate(S, 1, p, "deviation")
  print(x[, c("family", "strategy", "method", "window", "est", "lo", "hi", "ok", "error", "secs")], digits = 3)
  cat(sprintf("one cohort in %.1f min\n", as.numeric(Sys.time() - t1, units = "mins")))
} else {
  rows <- run_scenario(S, 1:200, "deviation", cores = cores)
  perf <- summarise_runs(rows, tr)
  saveRDS(list(rows = rows, perf = perf, code_version = code_version()), "sandbox_rwe/results/informative_ltfu_rows.rds")
  cat("rows:", nrow(rows), " failed:", sum(!rows$ok %in% TRUE), "\n")
  cat(sprintf("done in %.1f h\n", as.numeric(Sys.time() - t0, units = "hours")))
}
