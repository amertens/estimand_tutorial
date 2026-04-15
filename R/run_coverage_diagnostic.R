# run_coverage_diagnostic.R
# Iterative diagnostic to achieve 95% LMTP coverage under the simple DGP.
#
# Strategy: start with the simplest possible setup where LMTP should work
# well, verify coverage, then add complexity one factor at a time.
#
# Each "scenario" is a named configuration. We run all scenarios in each
# iteration so they share the same random seed for baseline covariates.
# 30 iterations per scenario, N=5000, treatment-policy estimand.
#
# Scenarios tested (cumulative complexity):
#   A. "oracle"       — no confounding (treat_override, no dep_censor), bw=14
#   B. "linear_noC"   — linear confounding, no dep_censor, bw=14
#   C. "linear_depC"  — linear confounding, dep_censor=TRUE, bw=14
#   D. "linear_bw7"   — linear confounding, dep_censor=TRUE, bw=7
#   E. "linear_N10k"  — linear confounding, dep_censor=TRUE, bw=14, N=10000
#
# If scenario A has 95% coverage, the LMTP framework is correct.
# If B drops coverage, confounding adjustment is the issue.
# If C drops further, dependent censoring is the issue.
# D and E test whether finer bins or larger N recover coverage.

library(here)
library(dplyr)
library(survival)

source(here("DGP.R"))
source(here("R", "helpers.R"))

tau <- 180
SL_LIB <- c("SL.mean", "SL.glm")
CV_FOLDS <- 2
N_ITER <- 30
SAVE_EVERY <- 5

requireNamespace("lmtp", quietly = TRUE)
requireNamespace("SuperLearner", quietly = TRUE)

cache_file <- here("results", "coverage_diagnostic.rds")

# ── Define scenarios ──
scenarios <- list(
  A_oracle = list(
    label = "A: No confounding (oracle)",
    N = 5000, complexity = FALSE, dep_censor = FALSE,
    treat_override = "simulate",  # still simulate trt but no confounding effect
    bin_width = 14,
    note = "Simplest case. If coverage is bad here, the LMTP framework itself has a problem."
  ),
  B_linear_noC = list(
    label = "B: Linear confounding, no dep censor",
    N = 5000, complexity = FALSE, dep_censor = FALSE,
    treat_override = "simulate",
    bin_width = 14,
    note = "Adds confounded treatment assignment. Tests confounding adjustment."
  ),
  C_linear_depC = list(
    label = "C: Linear confounding + dep censor",
    N = 5000, complexity = FALSE, dep_censor = TRUE,
    treat_override = "simulate",
    bin_width = 14,
    note = "Adds dependent censoring. Tests censoring model."
  ),
  D_linear_bw7 = list(
    label = "D: Same as C but bw=7 (weekly)",
    N = 5000, complexity = FALSE, dep_censor = TRUE,
    treat_override = "simulate",
    bin_width = 7,
    note = "Finer time bins. Tests discretisation effect."
  ),
  E_linear_N10k = list(
    label = "E: Same as C but N=10000",
    N = 10000, complexity = FALSE, dep_censor = TRUE,
    treat_override = "simulate",
    bin_width = 14,
    note = "Larger sample. Tests finite-sample bias."
  )
)

# ── Compute truth for each scenario ──
message("=== Computing truth for each scenario ===")
truths <- list()
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  set.seed(9999)
  df1 <- generate_hep_data(N = 500000, np_hazard = TRUE,
                            dep_censor = FALSE, switch_on = TRUE,
                            complexity = sc$complexity,
                            treat_override = "all_treated", seed = 9999)
  df0 <- generate_hep_data(N = 500000, np_hazard = TRUE,
                            dep_censor = FALSE, switch_on = TRUE,
                            complexity = sc$complexity,
                            treat_override = "all_control", seed = 9999)
  r1 <- mean(df1$event == 1 & df1$follow_time <= tau)
  r0 <- mean(df0$event == 1 & df0$follow_time <= tau)
  truths[[sc_name]] <- r1 - r0
  message(sprintf("  %-25s true_rd = %+.5f  event_rate = %.3f",
                  sc$label, truths[[sc_name]],
                  mean(df1$event == 1 & df1$follow_time <= tau)))
}

# ── Load partial results ──
if (file.exists(cache_file)) {
  results_so_far <- readRDS(cache_file)
  message("Loaded ", nrow(results_so_far), " cached rows")
} else {
  results_so_far <- NULL
}

# ── Run iterations ──
message("\n=== Running ", N_ITER, " iterations x ", length(scenarios), " scenarios ===")
results_list <- list()

for (i in seq_len(N_ITER)) {
  t0 <- Sys.time()
  row <- tibble(iter = i)

  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    col_rd  <- paste0("rd_", sc_name)
    col_low <- paste0("lo_", sc_name)
    col_hi  <- paste0("hi_", sc_name)
    col_ok  <- paste0("ok_", sc_name)

    # Check cache
    cached <- if (!is.null(results_so_far) && i %in% results_so_far$iter &&
                  col_rd %in% names(results_so_far)) {
      r <- results_so_far[results_so_far$iter == i, ]
      if (!is.na(r[[col_rd]])) r else NULL
    } else NULL

    if (!is.null(cached)) {
      row[[col_rd]]  <- cached[[col_rd]]
      row[[col_low]] <- cached[[col_low]]
      row[[col_hi]]  <- cached[[col_hi]]
      row[[col_ok]]  <- cached[[col_ok]]
    } else {
      set.seed(5000 + i)
      dat <- generate_hep_data(N = sc$N, np_hazard = TRUE,
                                dep_censor = sc$dep_censor,
                                complexity = sc$complexity,
                                seed = 5000 + i)

      res <- tryCatch(suppressWarnings({
        prep <- prepare_lmtp_data(dat, tau = tau, bin_width = sc$bin_width,
                                  time_varying_trt = FALSE)
        fit <- run_lmtp_analysis(prep, folds = CV_FOLDS, learners = SL_LIB)
        rd_obj <- fit$contrast_rd
        rd_est <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
                  else rd_obj$vals$theta
        rd_se  <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
                  else rd_obj$vals$std.error
        rd_ci  <- rd_est + c(-1, 1) * 1.96 * rd_se
        list(rd = rd_est, lo = rd_ci[1], hi = rd_ci[2], ok = TRUE)
      }), error = function(e) {
        list(rd = NA, lo = NA, hi = NA, ok = FALSE)
      })

      row[[col_rd]]  <- res$rd
      row[[col_low]] <- res$lo
      row[[col_hi]]  <- res$hi
      row[[col_ok]]  <- res$ok
    }
  }

  results_list[[i]] <- row
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  if (i %% SAVE_EVERY == 0 || i == N_ITER) {
    results <- bind_rows(results_list)
    saveRDS(results, cache_file)

    message(sprintf("\n  [iter %d/%d, %.1f min] Saved.", i, N_ITER, elapsed))
    for (sc_name in names(scenarios)) {
      sc <- scenarios[[sc_name]]
      col_rd <- paste0("rd_", sc_name)
      col_ok <- paste0("ok_", sc_name)
      col_lo <- paste0("lo_", sc_name)
      col_hi <- paste0("hi_", sc_name)
      vals <- results[[col_rd]][results[[col_ok]] == TRUE & !is.na(results[[col_rd]])]
      lows <- results[[col_lo]][results[[col_ok]] == TRUE & !is.na(results[[col_rd]])]
      his  <- results[[col_hi]][results[[col_ok]] == TRUE & !is.na(results[[col_rd]])]
      n <- length(vals)
      tr <- truths[[sc_name]]
      if (n > 0) {
        bias <- mean(vals) - tr
        cov <- mean(lows <= tr & his >= tr)
        message(sprintf("    %-30s n=%2d  bias=%+.5f  coverage=%.0f%%",
                        sc$label, n, bias, cov * 100))
      }
    }
  } else {
    parts <- sapply(names(scenarios), function(sc_name) {
      v <- row[[paste0("rd_", sc_name)]]
      sprintf("%s=%s", sc_name, ifelse(is.na(v), "NA", sprintf("%+.4f", v)))
    })
    message(sprintf("  [iter %d, %.1f min] %s", i, elapsed, paste(parts, collapse = "  ")))
  }
}

message("\n=== FINAL RESULTS ===")
results <- bind_rows(results_list)
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  col_rd <- paste0("rd_", sc_name)
  col_ok <- paste0("ok_", sc_name)
  col_lo <- paste0("lo_", sc_name)
  col_hi <- paste0("hi_", sc_name)
  vals <- results[[col_rd]][results[[col_ok]] == TRUE & !is.na(results[[col_rd]])]
  lows <- results[[col_lo]][results[[col_ok]] == TRUE & !is.na(results[[col_rd]])]
  his  <- results[[col_hi]][results[[col_ok]] == TRUE & !is.na(results[[col_rd]])]
  n <- length(vals)
  tr <- truths[[sc_name]]
  if (n > 0) {
    bias <- mean(vals) - tr
    rmse <- sqrt(mean((vals - tr)^2))
    se <- sd(vals)
    cov <- mean(lows <= tr & his >= tr)
    message(sprintf("  %-30s truth=%+.5f  n=%2d  bias=%+.5f  RMSE=%.5f  SE=%.5f  cov=%.0f%%",
                    sc$label, tr, n, bias, rmse, se, cov * 100))
  }
}
message("\nDone.")
