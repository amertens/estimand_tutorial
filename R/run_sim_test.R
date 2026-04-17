# run_sim_test.R
# Small test (5 iterations x 5 estimands) to verify bootstrap CI columns
# work correctly in the parallel simulation pipeline.

library(here)
library(dplyr)
library(survival)

source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("R", "run_simulations.R"))

tau <- 180
BIN_WIDTH <- 14

test_cache <- here("results", "sim_study_test.rds")
if (file.exists(test_cache)) file.remove(test_cache)  # always fresh

sim_results <- run_simulation_study(
  n_iter      = 5,
  sample_size = 10000,
  tau         = tau,
  bin_width   = BIN_WIDTH,
  estimands   = c("treatment_policy", "no_switch", "while_on_treatment",
                  "composite", "principal_stratum"),
  run_lmtp    = TRUE,
  run_cox_td  = FALSE,
  dgp_args    = list(h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
                     np_hazard = TRUE, dep_censor = TRUE,
                     complexity = FALSE,
                     lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60,
                     return_potential_switching = TRUE),
  cache_file  = test_cache
)

# Quick sanity check on columns
cat("\n=== Column names in results ===\n")
print(names(sim_results))

cat("\n=== Method × estimand counts ===\n")
print(table(sim_results$estimand, sim_results$method))

cat("\n=== Sample bootstrap CIs for Cox methods (first iter of TP) ===\n")
tp1 <- sim_results[sim_results$estimand == "treatment_policy" & sim_results$iter == 1, ]
print(tp1[, c("method", "hr", "hr_ci_low", "hr_ci_high", "risk_diff", "ci_low", "ci_high")])

cat("\nDone. Results saved to ", test_cache, "\n")
