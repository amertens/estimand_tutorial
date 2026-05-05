# run_depcensor_sensitivity.R
#
# Sensitivity to informative censoring.
#
# Baseline simulation has dep_censor = TRUE: administrative censoring depends
# on baseline outcome covariates, so loss to follow-up is informative. This
# script runs a contrasting arm with dep_censor = FALSE (random administrative
# censoring) so the contribution of informative censoring to bias and
# coverage degradation can be isolated from the rest of the DGP.
#
# Output: results/sim_depcensor_sensitivity.rds (rows from both arms with
#   `censor_arm` column = "informative" or "non-informative")
#
# Usage: Rscript R/run_depcensor_sensitivity.R

suppressPackageStartupMessages({
  library(here); library(dplyr); library(parallel); library(lmtp)
  library(SuperLearner); library(arm)
})
source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("R", "run_simulations.R"))

tau       <- 180
BIN_WIDTH <- 14
N_ITER    <- 200
N_CORES   <- 4
LMTP_LEARNERS <- c("SL.mean", "SL.glm", "SL.bayesglm", "SL.glmnet")

base_dgp <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60,
  return_potential_switching = TRUE
)

ESTIMANDS <- c("treatment_policy", "no_switch", "while_on_treatment")

partial_dir <- here("results", "sim_partial_depcensor")
dir.create(partial_dir, showWarnings = FALSE, recursive = TRUE)

for (censor_arm in c("informative", "non_informative")) {
  dgp_args <- c(base_dgp,
                list(dep_censor = (censor_arm == "informative")))
  for (est in ESTIMANDS) {
    cache <- file.path(partial_dir, sprintf("%s_%s.rds", est, censor_arm))
    if (file.exists(cache)) {
      cat(sprintf("[skip] %s [%s]\n", est, censor_arm))
      next
    }
    cat(sprintf("\n=== %s [%s] (n_iter=%d) ===\n", est, censor_arm, N_ITER))
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
    res$censor_arm <- censor_arm
    saveRDS(res, cache)
    cat(sprintf("  done: %d rows\n", nrow(res)))
  }
}

files <- list.files(partial_dir, pattern = "\\.rds$", full.names = TRUE)
all_res <- bind_rows(lapply(files, function(f) {
  d <- readRDS(f)
  if (!"censor_arm" %in% names(d)) {
    arm <- if (grepl("informative", basename(f)) &&
               !grepl("non_informative", basename(f))) {
      "informative"
    } else "non_informative"
    d$censor_arm <- arm
  }
  d
}))
saveRDS(all_res, here("results", "sim_depcensor_sensitivity.rds"))
cat(sprintf("\nWrote %d rows to results/sim_depcensor_sensitivity.rds\n",
            nrow(all_res)))
