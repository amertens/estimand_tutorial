# run_production_check.R
# Production-level validation: richer SL library, 5-fold CV, 100 iterations.
# Targets the treatment-policy estimand under good support (default DGP)
# to validate whether LMTP recovers the truth with adequate specification.
#
# Rationale for choices:
#   - Treatment-policy: the primary estimand, uses baseline-only treatment
#     (simplest LMTP specification, isolates SL/CV performance from
#     time-varying treatment complexity)
#   - Good support: default switching parameters (~12%), avoids positivity
#     issues that would confound the SL performance assessment
#   - 100 iterations: enough for stable bias/coverage estimates
#     (95% CI for 95% coverage with n=100 is roughly 89-98%)
#
# Expected runtime: ~4-8 hours depending on hardware.
# Run: Rscript R/run_production_check.R

library(here)
library(dplyr)
library(survival)

source(here("DGP.R"))
source(here("R", "helpers.R"))

tau <- 180
BIN_WIDTH <- 14
N_SIM <- 10000
N_ITER <- 100

# Production SL library
SL_PROD <- c("SL.mean", "SL.glm", "SL.bayesglm", "SL.glm.interaction")
CV_FOLDS <- 5

# Load truth
truth_cache <- here("results", "sim_results", "ground_truth.rds")
if (!file.exists(truth_cache)) stop("Run R/run_cached_analyses.R first for ground truth.")
truth_all <- readRDS(truth_cache)
true_rd <- truth_all$treatment_policy$true_rd
true_hr <- truth_all$treatment_policy$true_hr

message("=== Production-Level Validation Check ===")
message("  Estimand: treatment_policy")
message("  True RD: ", round(true_rd, 5))
message("  True HR: ", round(true_hr, 4))
message("  N per iter: ", N_SIM)
message("  Iterations: ", N_ITER)
message("  SL library: ", paste(SL_PROD, collapse = ", "))
message("  CV folds: ", CV_FOLDS)
message("  Bin width: ", BIN_WIDTH, " days")
message("")

cache_file <- here("results", "production_check.rds")
if (file.exists(cache_file)) {
  message("Production check already cached. Loading...")
  results <- readRDS(cache_file)
} else {
  requireNamespace("lmtp", quietly = TRUE)
  requireNamespace("SuperLearner", quietly = TRUE)
  requireNamespace("arm", quietly = TRUE)

  results <- vector("list", N_ITER)
  pb <- txtProgressBar(min = 0, max = N_ITER, style = 3)

  for (i in seq_len(N_ITER)) {
    set.seed(2000 + i)

    dat <- generate_hep_data(
      N = N_SIM, np_hazard = TRUE, dep_censor = TRUE,
      complexity = TRUE, seed = 2000 + i
    )

    covars <- c("age", "ckd", "cirrhosis", "diabetes", "heart_failure")

    # Cox naive
    cox_res <- tryCatch({
      fit <- fit_cox_naive(dat, tau = tau, covars = covars)
      km <- km_risk_at(fit$survfit, t = tau)
      list(cox_hr = fit$hr, cox_rd = km$risk_diff, cox_converged = TRUE)
    }, error = function(e) {
      list(cox_hr = NA, cox_rd = NA, cox_converged = FALSE)
    })

    # LMTP — development spec (current: SL.mean + SL.glm + SL.bayesglm, 2-fold)
    lmtp_dev <- tryCatch({
      prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH,
                                time_varying_trt = FALSE)
      res <- run_lmtp_analysis(prep, folds = 2,
                               learners = c("SL.mean", "SL.glm", "SL.bayesglm"))
      rd_obj <- res$contrast_rd
      rd_est <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
                else rd_obj$vals$theta
      rd_se <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
               else rd_obj$vals$std.error
      rd_ci <- rd_est + c(-1, 1) * 1.96 * rd_se
      list(dev_rd = rd_est, dev_ci_low = rd_ci[1], dev_ci_high = rd_ci[2],
           dev_converged = TRUE)
    }, error = function(e) {
      list(dev_rd = NA, dev_ci_low = NA, dev_ci_high = NA,
           dev_converged = FALSE)
    })

    # LMTP — production spec (richer SL, 5-fold)
    lmtp_prod <- tryCatch({
      prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH,
                                time_varying_trt = FALSE)
      res <- run_lmtp_analysis(prep, folds = CV_FOLDS,
                               learners = SL_PROD)
      rd_obj <- res$contrast_rd
      rd_est <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
                else rd_obj$vals$theta
      rd_se <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
               else rd_obj$vals$std.error
      rd_ci <- rd_est + c(-1, 1) * 1.96 * rd_se
      list(prod_rd = rd_est, prod_ci_low = rd_ci[1], prod_ci_high = rd_ci[2],
           prod_converged = TRUE)
    }, error = function(e) {
      list(prod_rd = NA, prod_ci_low = NA, prod_ci_high = NA,
           prod_converged = FALSE)
    })

    results[[i]] <- tibble(
      iter = i,
      cox_hr = cox_res$cox_hr,
      cox_rd = cox_res$cox_rd,
      cox_converged = cox_res$cox_converged,
      dev_rd = lmtp_dev$dev_rd,
      dev_ci_low = lmtp_dev$dev_ci_low,
      dev_ci_high = lmtp_dev$dev_ci_high,
      dev_converged = lmtp_dev$dev_converged,
      prod_rd = lmtp_prod$prod_rd,
      prod_ci_low = lmtp_prod$prod_ci_low,
      prod_ci_high = lmtp_prod$prod_ci_high,
      prod_converged = lmtp_prod$prod_converged
    )

    setTxtProgressBar(pb, i)
  }
  close(pb)

  results <- bind_rows(results)
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, cache_file)
  message("Saved to ", cache_file)
}

# ── Summarize ──
message("\n=== RESULTS ===")
message("True TP RD: ", round(true_rd, 5))

summarize_method <- function(rd, ci_low, ci_high, converged, label) {
  ok <- converged & !is.na(rd)
  n <- sum(ok)
  if (n == 0) { message("  ", label, ": 0 converged"); return(invisible()) }
  rd_vals <- rd[ok]
  bias <- mean(rd_vals) - true_rd
  rmse <- sqrt(mean((rd_vals - true_rd)^2))
  emp_se <- sd(rd_vals)
  cov <- mean(ci_low[ok] <= true_rd & ci_high[ok] >= true_rd, na.rm = TRUE)
  message(sprintf("  %-30s n=%3d  mean_RD=%+.5f  bias=%+.5f  RMSE=%.5f  emp_SE=%.5f  coverage=%.0f%%",
                  label, n, mean(rd_vals), bias, rmse, emp_se, cov * 100))
}

message("\nCox naive:")
summarize_method(results$cox_rd, NA, NA, results$cox_converged, "Cox KM RD")

message("\nLMTP development spec (SL.mean+glm+bayesglm, 2-fold):")
summarize_method(results$dev_rd, results$dev_ci_low, results$dev_ci_high,
                 results$dev_converged, "LMTP dev")

message("\nLMTP production spec (", paste(SL_PROD, collapse="+"), ", ", CV_FOLDS, "-fold):")
summarize_method(results$prod_rd, results$prod_ci_low, results$prod_ci_high,
                 results$prod_converged, "LMTP prod")

message("\nDone.")
