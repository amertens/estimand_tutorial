# smoke_one_iter.R -- run a single serial iteration of treatment_policy
# to isolate where the parallel run is failing.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(survival); library(broom)
  library(lmtp); library(SuperLearner); library(arm)
})

source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("R", "run_simulations.R"))

dgp_args <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, dep_censor = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60,
  return_potential_switching = TRUE
)

cat("Running one TP iteration serially...\n")
t0 <- Sys.time()
res <- run_one_iter(
  i = 1, estimand = "treatment_policy",
  sample_size = 10000, tau = 180, bin_width = 14,
  lmtp_learners = c("SL.mean", "SL.glm", "SL.bayesglm"),
  run_lmtp = TRUE, run_cox_td = FALSE,
  dgp_args = dgp_args
)
cat(sprintf("Done in %.1f sec; %d rows.\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(res)))
print(res[, c("method", "risk_diff", "ci_low", "ci_high")])
