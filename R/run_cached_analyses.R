# run_cached_analyses.R
# Run all long-running analyses (LMTP, ground truth, simulations, support)
# and save results to results/ as .rds files.
# The QMD document loads these cached files instead of running inline.
#
# Usage: Rscript R/run_cached_analyses.R
# Or source from RStudio: source("R/run_cached_analyses.R")
#
# Speed settings: bin_width=7 uses weekly time bins (~26 columns instead
# of 180) for LMTP, cutting runtime ~7x. Increase to bin_width=1 for
# daily resolution in production.

library(here)
library(dplyr)
library(survival)

source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("R", "support_diagnostics.R"))
source(here("R", "run_simulations.R"))

dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "sim_results"), showWarnings = FALSE, recursive = TRUE)

tau <- 180
# Weekly bins: ~26 time points instead of 180. Major speedup for LMTP.
# Set to 1 for daily resolution in production analyses.
BIN_WIDTH <- 7

# ── 1. Ground truth (two estimands) ─────────────────────────────────────────
# This is fast (~1 min) -- pure DGP simulation, no LMTP.
source(here("calc_truth.R"))

truth_cache <- here("results", "sim_results", "ground_truth.rds")
if (!file.exists(truth_cache)) {
  message("=== Computing ground truth (2 x 200K subjects) ===")
  set.seed(9999)
  truth_tp <- truth_treatment_policy(
    N = 200000, tau = tau, seed = 9999,
    np_hazard = TRUE, complexity = TRUE, switch_on = TRUE
  )
  truth_ns <- truth_no_switch(
    N = 200000, tau = tau, seed = 9999,
    np_hazard = TRUE, complexity = TRUE
  )
  truth_all <- list(treatment_policy = truth_tp, no_switch = truth_ns)
  saveRDS(truth_all, truth_cache)
  message("Saved ground truth to ", truth_cache)
} else {
  message("Ground truth already cached.")
  truth_all <- readRDS(truth_cache)
}

message("  TP RD = ", round(truth_all$treatment_policy$true_rd, 6))
message("  NS RD = ", round(truth_all$no_switch$true_rd, 6))

# ── 2. Main LMTP analysis (single dataset) ──────────────────────────────────
# N=2000 with weekly bins: ~2-5 min
lmtp_cache <- here("results", "lmtp_main.rds")
if (!file.exists(lmtp_cache)) {
  message("\n=== Running main LMTP analysis (N=2000, bin_width=", BIN_WIDTH, ") ===")
  requireNamespace("lmtp", quietly = TRUE)
  requireNamespace("SuperLearner", quietly = TRUE)

  set.seed(2026)
  dat <- generate_hep_data(
    N = 2000, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, policy = "treatment_policy", seed = 2026
  )

  lmtp_prep <- prepare_lmtp_data(
    dat, tau = tau, bin_width = BIN_WIDTH,
    baseline = c("age", "sex_male", "ckd", "diabetes",
                 "hypertension", "heart_failure")
  )
  message("  LMTP columns: ", lmtp_prep$n_bins, " time bins")

  lmtp_res <- run_lmtp_analysis(lmtp_prep, folds = 2,
                                learners = c("SL.glm"))
  saveRDS(lmtp_res, lmtp_cache)
  message("Saved LMTP results to ", lmtp_cache)
} else {
  message("\nMain LMTP already cached.")
}

# ── 3. Support estimation across scenarios ───────────────────────────────────
# 3 scenarios x LMTP at N=1000 with weekly bins: ~5-10 min total
support_cache <- here("results", "support_estimation.rds")
if (!file.exists(support_cache)) {
  message("\n=== Running support scenario estimation ===")
  requireNamespace("lmtp", quietly = TRUE)
  requireNamespace("SuperLearner", quietly = TRUE)

  set.seed(101)
  dat_good <- generate_hep_data(
    N = 1000, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, policy = "treatment_policy", seed = 101
  )
  set.seed(102)
  dat_strained <- generate_hep_data(
    N = 1000, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, policy = "treatment_policy",
    gamma_A = 1.5, gamma_ckd = 1.2, lambda_sw0 = 5e-5, seed = 102
  )
  set.seed(103)
  dat_poor <- generate_hep_data(
    N = 1000, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, policy = "treatment_policy",
    gamma_A = 2.5, gamma_ckd = 2.0, lambda_sw0 = 1e-4, seed = 103
  )

  scenarios <- list(dat_good, dat_strained, dat_poor)
  scenario_names <- c("Good support", "Strained support", "Poor support")

  pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
  support_results <- lapply(seq_along(scenarios), function(i) {
    d <- scenarios[[i]]
    nm <- scenario_names[i]
    message("\n  Scenario: ", nm)

    cox_n <- tryCatch(
      fit_cox_naive(d, tau = tau),
      error = function(e) list(hr = NA, ci_low = NA, ci_high = NA)
    )
    cox_c <- tryCatch(
      fit_cox_censor_switch(d, tau = tau),
      error = function(e) list(hr = NA, ci_low = NA, ci_high = NA)
    )
    lmtp_rd <- tryCatch({
      prep <- prepare_lmtp_data(d, tau = tau, bin_width = BIN_WIDTH)
      res <- run_lmtp_analysis(prep, folds = 2)
      res$contrast_rd$estimates$estimate
    }, error = function(e) NA_real_)

    setTxtProgressBar(pb, i)
    tibble::tibble(
      Scenario        = nm,
      `Cox naive HR`  = as.numeric(cox_n$hr),
      `Cox censor HR` = as.numeric(cox_c$hr),
      `LMTP RD`       = as.numeric(lmtp_rd)
    )
  }) %>% bind_rows()
  close(pb)

  saveRDS(support_results, support_cache)
  message("Saved support results to ", support_cache)
} else {
  message("\nSupport estimation already cached.")
}

# ── 4. Repeated simulation study ─────────────────────────────────────────────
# 10 iters x 2 estimands x N=500 with weekly bins: ~10-20 min
sim_cache <- here("results", "sim_study_main.rds")
if (!file.exists(sim_cache)) {
  message("\n=== Running simulation study (10 iters x 2 estimands, N=500) ===")
  message("    (Increase n_iter and sample_size for production)")
  sim_results <- run_simulation_study(
    n_iter      = 10,
    sample_size = 500,
    tau         = tau,
    estimands   = c("treatment_policy", "no_switch"),
    run_lmtp    = TRUE,
    run_cox_td  = FALSE,
    dgp_args    = list(np_hazard = TRUE, dep_censor = TRUE,
                       complexity = TRUE),
    cache_file  = sim_cache
  )
  message("Saved simulation results to ", sim_cache)
} else {
  message("\nSimulation study already cached.")
}

message("\n=== All analyses complete. ===")
message("Cached files:")
message("  ", truth_cache)
message("  ", lmtp_cache)
message("  ", support_cache)
message("  ", sim_cache)
message("\nYou can now render the QMD document -- it will load these cached results.")
