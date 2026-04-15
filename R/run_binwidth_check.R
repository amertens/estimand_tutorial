# run_binwidth_check.R
# Diagnose whether LMTP bias is driven by time discretisation.
# Compares bin_width = 7 (weekly), 14 (biweekly), 28 (monthly), 60 (bimonthly)
# with minimal SL library (SL.mean + SL.glm) and 2-fold CV.
#
# Treatment-policy estimand, good support, N=5000, 50 iterations.
# Saves incrementally every 5 iterations; resume-safe.
#
# Run: Rscript R/run_binwidth_check.R

library(here)
library(dplyr)
library(survival)

source(here("DGP.R"))
source(here("R", "helpers.R"))

tau <- 180
N_SIM <- 5000
N_ITER <- 50
SL_LIB <- c("SL.mean", "SL.glm")
CV_FOLDS <- 2
SAVE_EVERY <- 5

BIN_WIDTHS <- c(7, 14, 28, 60)

# Load truth
truth_cache <- here("results", "sim_results", "ground_truth.rds")
if (!file.exists(truth_cache)) stop("Run R/run_cached_analyses.R first.")
truth_all <- readRDS(truth_cache)
true_rd <- truth_all$treatment_policy$true_rd

message("=== Bin Width Diagnostic ===")
message("  True TP RD: ", round(true_rd, 5))
message("  Bin widths: ", paste(BIN_WIDTHS, collapse = ", "), " days")
message("  N: ", N_SIM, " | Iters: ", N_ITER)
message("  SL: ", paste(SL_LIB, collapse = "+"), " | CV: ", CV_FOLDS, "-fold")
message("")

cache_file <- here("results", "binwidth_check.rds")

# Load partial results if they exist
if (file.exists(cache_file)) {
  results_so_far <- readRDS(cache_file)
  message("Loaded ", nrow(results_so_far), " cached rows")
} else {
  results_so_far <- NULL
}

requireNamespace("lmtp", quietly = TRUE)
requireNamespace("SuperLearner", quietly = TRUE)

results_list <- list()

for (i in seq_len(N_ITER)) {
  t0 <- Sys.time()
  set.seed(3000 + i)

  dat <- generate_hep_data(
    N = N_SIM, np_hazard = TRUE, dep_censor = TRUE,
    complexity = TRUE, seed = 3000 + i
  )

  row <- tibble(iter = i)

  for (bw in BIN_WIDTHS) {
    col_rd  <- paste0("rd_bw", bw)
    col_low <- paste0("ci_low_bw", bw)
    col_hi  <- paste0("ci_high_bw", bw)
    col_ok  <- paste0("conv_bw", bw)

    # Check if this iter+bw already cached
    cached <- if (!is.null(results_so_far) && i %in% results_so_far$iter) {
      r <- results_so_far[results_so_far$iter == i, ]
      if (col_rd %in% names(r) && !is.na(r[[col_rd]])) r else NULL
    } else NULL

    if (!is.null(cached)) {
      row[[col_rd]]  <- cached[[col_rd]]
      row[[col_low]] <- cached[[col_low]]
      row[[col_hi]]  <- cached[[col_hi]]
      row[[col_ok]]  <- cached[[col_ok]]
    } else {
      res <- tryCatch(
        suppressWarnings({
          prep <- prepare_lmtp_data(dat, tau = tau, bin_width = bw,
                                    time_varying_trt = FALSE)
          fit <- run_lmtp_analysis(prep, folds = CV_FOLDS, learners = SL_LIB)
          rd_obj <- fit$contrast_rd
          rd_est <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
                    else rd_obj$vals$theta
          rd_se  <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
                    else rd_obj$vals$std.error
          rd_ci  <- rd_est + c(-1, 1) * 1.96 * rd_se
          list(rd = rd_est, ci_low = rd_ci[1], ci_high = rd_ci[2], ok = TRUE)
        }),
        error = function(e) {
          list(rd = NA, ci_low = NA, ci_high = NA, ok = FALSE)
        }
      )
      row[[col_rd]]  <- res$rd
      row[[col_low]] <- res$ci_low
      row[[col_hi]]  <- res$ci_high
      row[[col_ok]]  <- res$ok
    }
  }

  results_list[[i]] <- row
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  if (i %% SAVE_EVERY == 0 || i == N_ITER) {
    results <- bind_rows(results_list)
    saveRDS(results, cache_file)

    message(sprintf("\n  [iter %d/%d, %.1f min] Saved. Running summary:", i, N_ITER, elapsed))
    for (bw in BIN_WIDTHS) {
      col_rd <- paste0("rd_bw", bw)
      col_ok <- paste0("conv_bw", bw)
      col_low <- paste0("ci_low_bw", bw)
      col_hi  <- paste0("ci_high_bw", bw)
      vals <- results[[col_rd]][results[[col_ok]] == TRUE]
      lows <- results[[col_low]][results[[col_ok]] == TRUE]
      his  <- results[[col_hi]][results[[col_ok]] == TRUE]
      n <- length(vals[!is.na(vals)])
      if (n > 0) {
        bias <- mean(vals, na.rm = TRUE) - true_rd
        cov <- mean(lows <= true_rd & his >= true_rd, na.rm = TRUE)
        message(sprintf("    bw=%2d (%2d bins): n=%2d  bias=%+.5f  coverage=%.0f%%",
                        bw, ceiling(tau / bw), n, bias, cov * 100))
      }
    }
  } else {
    rd_str <- paste(sapply(BIN_WIDTHS, function(bw) {
      v <- row[[paste0("rd_bw", bw)]]
      sprintf("bw%d=%s", bw, ifelse(is.na(v), "NA", sprintf("%+.4f", v)))
    }), collapse = "  ")
    message(sprintf("  [iter %d, %.1f min] %s", i, elapsed, rd_str))
  }
}

results <- bind_rows(results_list)

# Final summary
message("\n=== FINAL RESULTS ===")
message("True TP RD: ", round(true_rd, 5), "\n")
for (bw in BIN_WIDTHS) {
  col_rd <- paste0("rd_bw", bw)
  col_ok <- paste0("conv_bw", bw)
  col_low <- paste0("ci_low_bw", bw)
  col_hi  <- paste0("ci_high_bw", bw)
  vals <- results[[col_rd]][results[[col_ok]] == TRUE]
  lows <- results[[col_low]][results[[col_ok]] == TRUE]
  his  <- results[[col_hi]][results[[col_ok]] == TRUE]
  n <- length(vals[!is.na(vals)])
  if (n > 0) {
    bias <- mean(vals, na.rm = TRUE) - true_rd
    rmse <- sqrt(mean((vals - true_rd)^2, na.rm = TRUE))
    emp_se <- sd(vals, na.rm = TRUE)
    cov <- mean(lows <= true_rd & his >= true_rd, na.rm = TRUE)
    message(sprintf("  bw=%2d (%2d bins): n=%2d  mean_RD=%+.5f  bias=%+.5f  RMSE=%.5f  emp_SE=%.5f  coverage=%.0f%%",
                    bw, ceiling(tau / bw), n, mean(vals, na.rm = TRUE), bias, rmse, emp_se, cov * 100))
  } else {
    message(sprintf("  bw=%2d: 0 converged", bw))
  }
}

message("\nDone.")
