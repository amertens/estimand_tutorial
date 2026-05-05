# run_markov_sensitivity.R
#
# Sensitivity analysis: compare LMTP SDR fits using
#   (a) full-history nuisance models (k = Inf, default)
#   (b) Markov restriction (k = 0; only the most recent time-varying covariate
#       enters each nuisance model)
#
# Run on the no-switch and while-on-treatment estimands, where time-varying
# nuisance modelling is most consequential. 100 iterations per arm is enough
# to detect a meaningful coverage shift.
#
# Output: results/sim_markov_sensitivity.rds
#
# Usage: Rscript R/run_markov_sensitivity.R

suppressPackageStartupMessages({
  library(here); library(dplyr); library(parallel); library(lmtp)
  library(SuperLearner); library(arm)
})
source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("R", "run_simulations.R"))

tau       <- 180
BIN_WIDTH <- 14
N_ITER    <- 100
N_CORES   <- 4
LMTP_LEARNERS <- c("SL.mean", "SL.glm", "SL.bayesglm", "SL.glmnet")

dgp_args <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, dep_censor = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60,
  return_potential_switching = TRUE
)

ESTIMANDS <- c("no_switch", "while_on_treatment")

out <- list()
for (k_val in c(Inf, 0)) {
  k_label <- if (is.infinite(k_val)) "full_history" else "markov_k0"
  for (est in ESTIMANDS) {
    cache <- here("results", "sim_partial_markov",
                  sprintf("%s_%s.rds", est, k_label))
    dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
    if (file.exists(cache)) {
      cat(sprintf("[skip] %s [%s]\n", est, k_label))
      next
    }
    cat(sprintf("\n=== %s [%s] (n_iter=%d) ===\n", est, k_label, N_ITER))
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
      lmtp_k        = k_val,
      cache_file    = cache
    )
    cat(sprintf("  done: %d rows\n", nrow(res)))
  }
}

# Concatenate
files <- list.files(here("results", "sim_partial_markov"),
                    pattern = "\\.rds$", full.names = TRUE)
all_res <- bind_rows(lapply(files, function(f) {
  d <- readRDS(f)
  parts <- strsplit(basename(tools::file_path_sans_ext(f)), "_")[[1]]
  d$markov_arm <- if (tail(parts, 1) == "history") "full_history" else "markov_k0"
  d
}))
saveRDS(all_res, here("results", "sim_markov_sensitivity.rds"))
cat(sprintf("\nWrote %d rows to results/sim_markov_sensitivity.rds\n",
            nrow(all_res)))
