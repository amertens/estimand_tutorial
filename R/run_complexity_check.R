# run_complexity_check.R
# Diagnose whether LMTP bias is driven by nonlinear confounding.
# Compares complexity=TRUE (nonlinear PS & outcome) vs complexity=FALSE
# (linear) with SL.mean+SL.glm, 2-fold CV, bin_width=14.
#
# Treatment-policy estimand, N=5000, 30 iterations.
# Saves incrementally every 5 iterations.

library(here)
library(dplyr)
library(survival)

source(here("DGP.R"))
source(here("R", "helpers.R"))

tau <- 180
BIN_WIDTH <- 14
N_SIM <- 5000
N_ITER <- 30
SL_LIB <- c("SL.mean", "SL.glm")
CV_FOLDS <- 2
SAVE_EVERY <- 5

cache_file <- here("results", "complexity_check.rds")

# Compute truth for BOTH complexity settings
message("=== Computing truth for complexity=TRUE and FALSE ===")
requireNamespace("lmtp", quietly = TRUE)
requireNamespace("SuperLearner", quietly = TRUE)

# Truth: complexity=TRUE (default)
set.seed(9999)
df_a1_c <- generate_hep_data(N = 500000, np_hazard = TRUE, dep_censor = FALSE,
                              switch_on = TRUE, complexity = TRUE,
                              treat_override = "all_treated", seed = 9999)
df_a0_c <- generate_hep_data(N = 500000, np_hazard = TRUE, dep_censor = FALSE,
                              switch_on = TRUE, complexity = TRUE,
                              treat_override = "all_control", seed = 9999)
true_rd_complex <- mean(df_a1_c$event == 1 & df_a1_c$follow_time <= tau) -
                   mean(df_a0_c$event == 1 & df_a0_c$follow_time <= tau)

# Truth: complexity=FALSE
set.seed(9999)
df_a1_s <- generate_hep_data(N = 500000, np_hazard = TRUE, dep_censor = FALSE,
                              switch_on = TRUE, complexity = FALSE,
                              treat_override = "all_treated", seed = 9999)
df_a0_s <- generate_hep_data(N = 500000, np_hazard = TRUE, dep_censor = FALSE,
                              switch_on = TRUE, complexity = FALSE,
                              treat_override = "all_control", seed = 9999)
true_rd_simple <- mean(df_a1_s$event == 1 & df_a1_s$follow_time <= tau) -
                  mean(df_a0_s$event == 1 & df_a0_s$follow_time <= tau)

message("  True RD (complex): ", round(true_rd_complex, 5))
message("  True RD (simple):  ", round(true_rd_simple, 5))

# Load partial results if they exist
if (file.exists(cache_file)) {
  results_so_far <- readRDS(cache_file)
  message("Loaded ", nrow(results_so_far), " cached rows")
} else {
  results_so_far <- NULL
}

results_list <- list()

for (i in seq_len(N_ITER)) {
  t0 <- Sys.time()

  row <- tibble(iter = i)

  for (cmplx in c(TRUE, FALSE)) {
    tag <- if (cmplx) "complex" else "simple"
    col_rd  <- paste0("rd_", tag)
    col_low <- paste0("ci_low_", tag)
    col_hi  <- paste0("ci_high_", tag)
    col_ok  <- paste0("conv_", tag)

    # Check cache
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
      set.seed(4000 + i)
      dat <- generate_hep_data(N = N_SIM, np_hazard = TRUE, dep_censor = TRUE,
                                complexity = cmplx, seed = 4000 + i)

      res <- tryCatch(suppressWarnings({
        prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH,
                                  time_varying_trt = FALSE)
        fit <- run_lmtp_analysis(prep, folds = CV_FOLDS, learners = SL_LIB)
        rd_obj <- fit$contrast_rd
        rd_est <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
                  else rd_obj$vals$theta
        rd_se  <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
                  else rd_obj$vals$std.error
        rd_ci  <- rd_est + c(-1, 1) * 1.96 * rd_se
        list(rd = rd_est, ci_low = rd_ci[1], ci_high = rd_ci[2], ok = TRUE)
      }), error = function(e) {
        list(rd = NA, ci_low = NA, ci_high = NA, ok = FALSE)
      })

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

    # Complex summary
    v_c <- results$rd_complex[results$conv_complex == TRUE & !is.na(results$rd_complex)]
    l_c <- results$ci_low_complex[results$conv_complex == TRUE & !is.na(results$rd_complex)]
    h_c <- results$ci_high_complex[results$conv_complex == TRUE & !is.na(results$rd_complex)]
    # Simple summary
    v_s <- results$rd_simple[results$conv_simple == TRUE & !is.na(results$rd_simple)]
    l_s <- results$ci_low_simple[results$conv_simple == TRUE & !is.na(results$rd_simple)]
    h_s <- results$ci_high_simple[results$conv_simple == TRUE & !is.na(results$rd_simple)]

    message(sprintf("\n  [iter %d/%d, %.1f min] Saved.", i, N_ITER, elapsed))
    if (length(v_c) > 0)
      message(sprintf("    Complex (true=%+.5f): n=%d  bias=%+.5f  coverage=%.0f%%",
                      true_rd_complex, length(v_c),
                      mean(v_c) - true_rd_complex,
                      mean(l_c <= true_rd_complex & h_c >= true_rd_complex) * 100))
    if (length(v_s) > 0)
      message(sprintf("    Simple  (true=%+.5f): n=%d  bias=%+.5f  coverage=%.0f%%",
                      true_rd_simple, length(v_s),
                      mean(v_s) - true_rd_simple,
                      mean(l_s <= true_rd_simple & h_s >= true_rd_simple) * 100))
  } else {
    message(sprintf("  [iter %d, %.1f min] complex=%s  simple=%s", i, elapsed,
                    ifelse(is.na(row$rd_complex), "NA", sprintf("%+.4f", row$rd_complex)),
                    ifelse(is.na(row$rd_simple), "NA", sprintf("%+.4f", row$rd_simple))))
  }
}

message("\nDone.")
