# 72_wot_comp_scenarios.R -- SDR for while on treatment and the composite in the five
# one-feature scenarios, which ran SDR for treatment policy and sustained use only (Andrew
# 2026-11-20: complete Figure 1, Table S4 and Table S14). Same 200 cohorts per scenario
# (seed_for(scenario, rep)) and the same per-estimator seeds as the base case ("wot", "comp").
#
#   Rscript --vanilla sandbox_rwe/run/72_wot_comp_scenarios.R pilot   # 1 cohort, timed
#   Rscript --vanilla sandbox_rwe/run/72_wot_comp_scenarios.R         # 5 x 200 cohorts, RWE_CORES workers
setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
for (f in list.files("sandbox_rwe/R", pattern = "[.]R$", full.names = TRUE)) source(f)
S5 <- c("c1_ph_rare_switch", "d2a_ph", "d2b_heavy_ltfu", "d2c_renal_switch", "d2d_common_outcome")
cores <- as.integer(Sys.getenv("RWE_CORES", "4"))
t0 <- Sys.time()
if (length(commandArgs(trailingOnly = TRUE)) && commandArgs(trailingOnly = TRUE)[1] == "pilot") {
  x <- run_replicate("d2d_common_outcome", 1, get_scenario_params("d2d_common_outcome"), "wotcomp")
  print(x[, c("scenario", "strategy", "method", "est", "lo", "hi", "rr", "max_w", "ok", "error", "secs")], digits = 3)
  print(get_truth("d2d_common_outcome")$strategies[, c("strategy", "rd", "rr")], digits = 3)
  cat(sprintf("one cohort in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
} else {
  rows <- do.call(rbind, lapply(S5, function(s) {
    cat(format(Sys.time(), "%m-%d %H:%M"), "-", s, "\n")
    run_scenario(s, 1:200, "wotcomp", cores = cores)
  }))
  perf <- do.call(rbind, lapply(S5, function(s) summarise_runs(rows[rows$scenario == s, ], get_truth(s))))
  saveRDS(list(rows = rows, perf = perf, code_version = code_version()), "sandbox_rwe/results/wotcomp_rows.rds")
  cat("rows:", nrow(rows), " failed:", sum(!rows$ok %in% TRUE), "\n")
  cat(sprintf("done in %.1f h\n", as.numeric(Sys.time() - t0, units = "hours")))
}
