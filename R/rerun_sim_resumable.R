# rerun_sim_resumable.R -- run the main sim study estimand-by-estimand,
# caching each intermediate result so we can recover from a crash without
# re-running the estimands that already completed.
#
# After all five estimands cache successfully, results are concatenated and
# written to results/sim_study_main.rds. If a particular estimand keeps
# failing, lower N_CORES below or run that estimand in serial mode.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(survival); library(broom)
  library(parallel); library(lmtp); library(SuperLearner); library(arm)
})
source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("R", "run_simulations.R"))

tau <- 180
BIN_WIDTH <- 14
N_ITER <- 500    # bumped from 200 to stabilise coverage estimates
N_CORES <- 4

# Expanded Super Learner library: adds elastic net (SL.glmnet) and random
# forest (SL.ranger) on top of the original mean/glm/bayesglm to mitigate
# overfitting in the longitudinal nuisance models. SL.ranger is optional --
# drop it if it materially slows production runs.
LMTP_LEARNERS <- c("SL.mean", "SL.glm", "SL.bayesglm", "SL.glmnet", "SL.ranger")

dgp_args <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, dep_censor = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60,
  return_potential_switching = TRUE
)

ESTIMANDS <- c("treatment_policy", "no_switch", "while_on_treatment",
               "composite", "principal_stratum")

partial_dir <- here("results", "sim_partial")
dir.create(partial_dir, showWarnings = FALSE, recursive = TRUE)

for (est in ESTIMANDS) {
  cache <- file.path(partial_dir, paste0(est, ".rds"))
  if (file.exists(cache)) {
    cat(sprintf("[skip] %s already cached at %s\n", est, cache))
    next
  }
  cat(sprintf("\n=== %s (n_iter=%d, n_cores=%d) ===\n", est, N_ITER, N_CORES))
  t0 <- Sys.time()
  res <- run_simulation_study(
    n_iter        = N_ITER,
    sample_size   = 10000,
    tau           = tau,
    bin_width     = BIN_WIDTH,
    estimands     = est,
    run_lmtp      = TRUE,
    run_cox_td    = FALSE,
    n_cores       = N_CORES,
    dgp_args      = dgp_args,
    lmtp_learners = LMTP_LEARNERS,
    cache_file    = cache
  )
  cat(sprintf("  %s done: %d rows in %.1f min\n",
              est, nrow(res),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# Concatenate and save
all_res <- bind_rows(lapply(ESTIMANDS, function(e) {
  readRDS(file.path(partial_dir, paste0(e, ".rds")))
}))
saveRDS(all_res, here("results", "sim_study_main.rds"))
cat(sprintf("\nWrote %d rows to results/sim_study_main.rds\n", nrow(all_res)))
