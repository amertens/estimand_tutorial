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
BIN_WIDTH <- 30

# ── 1. Ground truth (five estimands) ─────────────────────────────────────────
# Fast (~3-5 min) -- pure DGP simulation, no LMTP.
source(here("calc_truth.R"))

truth_cache <- here("results", "sim_results", "ground_truth.rds")
if (!file.exists(truth_cache)) {
  message("=== Computing ground truth (5 estimands x 500K subjects) ===")
  dgp_args <- list(np_hazard = TRUE, complexity = TRUE, switch_on = TRUE)

  message("  Treatment-policy...")
  truth_tp <- do.call(truth_treatment_policy,
    c(list(N = 500000, tau = tau, seed = 9999), dgp_args))
  message("  No-switch...")
  truth_ns <- do.call(truth_no_switch,
    c(list(N = 500000, tau = tau, seed = 9999), dgp_args))
  message("  While-on-treatment...")
  truth_wot <- do.call(truth_while_on_treatment,
    c(list(N = 500000, tau = tau, seed = 9999), dgp_args))
  message("  Composite...")
  truth_comp <- do.call(truth_composite,
    c(list(N = 500000, tau = tau, seed = 9999), dgp_args))
  message("  Principal stratum...")
  truth_ps <- do.call(truth_principal_stratum,
    c(list(N = 500000, tau = tau, seed = 9999), dgp_args))

  truth_all <- list(
    treatment_policy   = truth_tp,
    no_switch          = truth_ns,
    while_on_treatment = truth_wot,
    composite          = truth_comp,
    principal_stratum  = truth_ps
  )
  saveRDS(truth_all, truth_cache)
  message("Saved ground truth to ", truth_cache)
} else {
  message("Ground truth already cached.")
  truth_all <- readRDS(truth_cache)
}

message("  TP   RD = ", round(truth_all$treatment_policy$true_rd, 6))
message("  NS   RD = ", round(truth_all$no_switch$true_rd, 6))
message("  WOT  RD = ", round(truth_all$while_on_treatment$true_rd, 6))
message("  COMP RD = ", round(truth_all$composite$true_rd, 6))
message("  PS   RD = ", round(truth_all$principal_stratum$true_rd, 6))
if (!is.null(truth_all$principal_stratum$pct_never_switch))
  message("  PS never-switcher %: ",
          round(truth_all$principal_stratum$pct_never_switch * 100, 1))

# ── 2. Main LMTP analysis (single dataset, N=10000) ─────────────────────────
lmtp_cache <- here("results", "lmtp_main.rds")
if (!file.exists(lmtp_cache)) {
  message("\n=== Running main LMTP analysis (N=10000, bin_width=", BIN_WIDTH, ") ===")
  requireNamespace("lmtp", quietly = TRUE)
  requireNamespace("SuperLearner", quietly = TRUE)
  requireNamespace("arm", quietly = TRUE)

  SL_LIBRARY <- c("SL.mean", "SL.bayesglm")

  set.seed(2026)
  dat <- generate_hep_data(
    N = 10000, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, policy = "treatment_policy", seed = 2026
  )

  message("  Event rate: ", round(mean(dat$event), 4))
  message("  Switch rate: ", round(mean(dat$switched), 4))

  lmtp_prep <- prepare_lmtp_data(
    dat, tau = tau, bin_width = BIN_WIDTH,
    baseline = c("age", "sex_male", "ckd", "diabetes",
                 "hypertension", "heart_failure")
  )
  message("  LMTP columns: ", lmtp_prep$n_bins, " time bins")
  message("  SL library: ", paste(SL_LIBRARY, collapse = ", "))

  lmtp_res <- run_lmtp_analysis(lmtp_prep, folds = 2,
                                learners = SL_LIBRARY)
  saveRDS(lmtp_res, lmtp_cache)
  message("Saved LMTP results to ", lmtp_cache)
} else {
  message("\nMain LMTP already cached.")
}

# ── 3. Support estimation across scenarios ───────────────────────────────────
# 3 scenarios x LMTP at N=5000 with weekly bins
support_cache <- here("results", "support_estimation.rds")
if (!file.exists(support_cache)) {
  message("\n=== Running support scenario estimation (N=5000 each) ===")
  requireNamespace("lmtp", quietly = TRUE)
  requireNamespace("SuperLearner", quietly = TRUE)

  set.seed(101)
  dat_good <- generate_hep_data(
    N = 5000, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, policy = "treatment_policy", seed = 101
  )
  set.seed(102)
  dat_strained <- generate_hep_data(
    N = 5000, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, policy = "treatment_policy",
    gamma_A = 1.5, gamma_ckd = 1.2, lambda_sw0 = 5e-3, seed = 102
  )
  set.seed(103)
  dat_poor <- generate_hep_data(
    N = 5000, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, policy = "treatment_policy",
    gamma_A = 2.5, gamma_ckd = 2.0, lambda_sw0 = 1e-2, seed = 103
  )

  scenarios <- list(dat_good, dat_strained, dat_poor)
  scenario_names <- c("Good support", "Strained support", "Poor support")

  pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
  support_results <- lapply(seq_along(scenarios), function(i) {
    d <- scenarios[[i]]
    nm <- scenario_names[i]
    message("\n  Scenario: ", nm)
    message("    Event rate: ", round(mean(d$event), 4),
            "  Switch rate: ", round(mean(d$switched), 4))

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
      res <- run_lmtp_analysis(prep, folds = 2,
                               learners = c("SL.mean", "SL.bayesglm"))
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
# 10 iters x 2 estimands x N=10000 with weekly bins
sim_cache <- here("results", "sim_study_main.rds")
if (!file.exists(sim_cache)) {
  message("\n=== Running simulation study (10 iters x 5 estimands, N=10000) ===")
  message("    (Increase n_iter for production; 10 is for fast iteration)")
  sim_results <- run_simulation_study(
    n_iter      = 10,
    sample_size = 10000,
    tau         = tau,
    bin_width   = BIN_WIDTH,
    estimands   = c("treatment_policy", "no_switch", "while_on_treatment",
                    "composite", "principal_stratum"),
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
