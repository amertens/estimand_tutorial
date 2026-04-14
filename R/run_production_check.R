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
BIN_WIDTH <- 28      # monthly bins (~7 time points) for speed
N_SIM <- 5000        # smaller N for faster SL fits
N_ITER <- 100

# Production SL library — adding xgboost to test flexible learner effect
SL_PROD <- c("SL.mean", "SL.glm", "SL.bayesglm", "SL.xgboost")
CV_FOLDS <- 2  # keep 2-fold for speed (5-fold showed no improvement)

# How often to save incremental results
SAVE_EVERY <- 5

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

  # Load any existing results (preserves dev columns from prior runs)
  if (file.exists(cache_file)) {
    results_so_far <- readRDS(cache_file)
    message("Loaded ", nrow(results_so_far), " cached rows (",
            max(results_so_far$iter), " iterations)")
  } else {
    results_so_far <- NULL
  }

  # Determine which iterations need prod rerun
  # (dev results are preserved from cache; only prod is recomputed)
  iters_need_prod <- if (is.null(results_so_far)) {
    seq_len(N_ITER)
  } else {
    # Rerun prod for all iters; keep dev from cache where available
    seq_len(N_ITER)
  }

  results_list <- list()

  for (i in iters_need_prod) {
    t0 <- Sys.time()
    set.seed(2000 + i)

    dat <- generate_hep_data(
      N = N_SIM, np_hazard = TRUE, dep_censor = TRUE,
      complexity = TRUE, seed = 2000 + i
    )

    covars <- c("age", "ckd", "cirrhosis", "diabetes", "heart_failure")

    # Check if dev results exist in cache for this iteration
    cached_row <- if (!is.null(results_so_far) && i %in% results_so_far$iter) {
      results_so_far[results_so_far$iter == i, ]
    } else NULL

    has_dev <- !is.null(cached_row) && !is.na(cached_row$dev_rd[1])

    # Cox naive — reuse from cache if available
    if (!is.null(cached_row) && !is.na(cached_row$cox_hr[1])) {
      cox_res <- list(cox_hr = cached_row$cox_hr[1],
                      cox_rd = cached_row$cox_rd[1],
                      cox_converged = cached_row$cox_converged[1])
    } else {
      cox_res <- tryCatch({
        fit <- fit_cox_naive(dat, tau = tau, covars = covars)
        km <- km_risk_at(fit$survfit, t = tau)
        list(cox_hr = fit$hr, cox_rd = km$risk_diff, cox_converged = TRUE)
      }, error = function(e) {
        list(cox_hr = NA, cox_rd = NA, cox_converged = FALSE)
      })
    }

    # LMTP dev — reuse from cache if available, skip recomputation
    if (has_dev) {
      lmtp_dev <- list(dev_rd = cached_row$dev_rd[1],
                       dev_ci_low = cached_row$dev_ci_low[1],
                       dev_ci_high = cached_row$dev_ci_high[1],
                       dev_converged = cached_row$dev_converged[1])
    } else {
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
    }

    # LMTP prod — always recompute (this is the new spec being tested)
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

    results_list[[as.character(i)]] <- tibble(
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

      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

      # Save and print running summary every SAVE_EVERY iterations
      if (i %% SAVE_EVERY == 0 || i == N_ITER) {
        results <- bind_rows(results_list)
        dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
        saveRDS(results, cache_file)

        n_done <- sum(!is.na(results$dev_rd) & results$dev_converged)
        n_done_p <- sum(!is.na(results$prod_rd) & results$prod_converged)
        dev_bias <- if (n_done > 0) mean(results$dev_rd[results$dev_converged], na.rm = TRUE) - true_rd else NA
        prod_bias <- if (n_done_p > 0) mean(results$prod_rd[results$prod_converged], na.rm = TRUE) - true_rd else NA
        dev_cov <- if (n_done > 0) mean(results$dev_ci_low[results$dev_converged] <= true_rd &
                                          results$dev_ci_high[results$dev_converged] >= true_rd, na.rm = TRUE) else NA
        prod_cov <- if (n_done_p > 0) mean(results$prod_ci_low[results$prod_converged] <= true_rd &
                                              results$prod_ci_high[results$prod_converged] >= true_rd, na.rm = TRUE) else NA

        message(sprintf("\n  [iter %d/%d, %.1f min] Saved. Running summary (n=%d/%d):",
                        i, N_ITER, elapsed, n_done, n_done_p))
        message(sprintf("    Dev  (2-fold):  bias=%+.5f  coverage=%.0f%%",
                        dev_bias, dev_cov * 100))
        message(sprintf("    Prod (%d-fold): bias=%+.5f  coverage=%.0f%%",
                        CV_FOLDS, prod_bias, prod_cov * 100))
      } else {
        message(sprintf("  [iter %d, %.1f min] dev_rd=%+.4f  prod_rd=%+.4f",
                        i, elapsed,
                        ifelse(lmtp_dev$dev_converged, lmtp_dev$dev_rd, NA),
                        ifelse(lmtp_prod$prod_converged, lmtp_prod$prod_rd, NA)))
      }
    }
    results <- bind_rows(results_list)
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
