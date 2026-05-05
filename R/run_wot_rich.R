# run_wot_rich.R
# Test richer SuperLearner libraries for the WOT estimand.
#
# Production LMTP WOT achieved 48% coverage with the minimal library
# (SL.mean, SL.glm, SL.bayesglm). The censoring model for WOT must capture
# the treatment-dependent switching process, which simple parametric
# learners cannot fit adequately. This script tests richer libraries.
#
# Learner descriptions:
#   minimal    : Baseline (matches production); SL.mean + SL.glm.
#                Reference point for the ~48% coverage we observed.
#   ranger_2f  : Adds SL.ranger (random forest) with 2-fold CV.
#                RF captures nonlinear interactions in the switching model.
#   rich_2f    : Adds SL.glmnet (penalised regression) + SL.ranger.
#                The SL ensemble can adaptively weight parametric and
#                nonparametric learners based on CV risk.
#   rich_5f    : Same learners as rich_2f but with 5-fold CV.
#                More stable CV risk estimates, at ~2.5x the cost.
#
# Usage: Rscript R/run_wot_rich.R (runs in parallel via clusterApply)

library(here)
library(dplyr)
library(survival)
library(parallel)

source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("calc_truth.R"))

tau <- 180
BIN_WIDTH <- 14
N_ITER <- 20
SAMPLE_SIZE <- 10000
N_CORES <- 2  # Reduced from 4: SL.ranger + SL.glmnet on N=10k blew RAM with 4 workers

requireNamespace("lmtp", quietly = TRUE)
requireNamespace("SuperLearner", quietly = TRUE)
requireNamespace("ranger", quietly = TRUE)
requireNamespace("glmnet", quietly = TRUE)

dgp_args <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60
)

# Compute WOT truth (cached)
truth_cache <- here("results", "sim_results", "ground_truth.rds")
truth_all <- readRDS(truth_cache)
truth_wot_rd <- truth_all$while_on_treatment$true_rd
message(sprintf("WOT truth RD = %+.4f", truth_wot_rd))

# SL library configurations
sl_configs <- list(
  minimal   = list(learners = c("SL.mean", "SL.glm"), folds = 2),
  ranger_2f = list(learners = c("SL.mean", "SL.glm", "SL.ranger"), folds = 2),
  rich_2f   = list(learners = c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger"), folds = 2),
  rich_5f   = list(learners = c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger"), folds = 5)
)

# Single-iteration worker: runs all SL configs on one dataset
run_one_iter <- function(i, sl_configs, dgp_args, tau, bin_width, sample_size) {
  dat <- do.call(generate_hep_data,
    c(list(N = sample_size, switch_on = TRUE, seed = 4000 + i), dgp_args))

  out <- list()
  for (sl_name in names(sl_configs)) {
    cfg <- sl_configs[[sl_name]]
    rd <- tryCatch(suppressWarnings({
      prep <- prepare_lmtp_data(dat, tau = tau, bin_width = bin_width,
                                time_varying_trt = FALSE,
                                competing_risk_at_switch = TRUE)
      fit <- run_lmtp_analysis(prep, folds = cfg$folds, learners = cfg$learners)
      rd_obj <- fit$contrast_rd
      rd_s <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
              else rd_obj$vals$theta
      se <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
            else rd_obj$vals$std.error
      # Negate for risk scale; swap CI bounds
      ci <- -(rd_s + c(1, -1) * 1.96 * se)
      list(rd = -rd_s, lo = ci[1], hi = ci[2], ok = TRUE)
    }), error = function(e) list(rd = NA, lo = NA, hi = NA, ok = FALSE))

    out[[length(out) + 1]] <- data.frame(
      iter = i, sl = sl_name, rd = rd$rd, lo = rd$lo, hi = rd$hi,
      ok = rd$ok, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

# Run in parallel
message(sprintf("Running %d iterations x %d SL configs on %d cores...",
                N_ITER, length(sl_configs), N_CORES))

cl <- makeCluster(N_CORES)
tryCatch({
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(dplyr); library(survival); library(here)
      library(lmtp); library(SuperLearner); library(ranger); library(glmnet); library(arm)
    })
    source(here("DGP.R")); source(here("R", "helpers.R"))
  })
  clusterExport(cl, c("run_one_iter", "sl_configs", "dgp_args", "tau",
                      "BIN_WIDTH", "SAMPLE_SIZE"), envir = environment())
  t0 <- Sys.time()
  iter_results <- clusterApply(cl, seq_len(N_ITER), function(i) {
    run_one_iter(i, sl_configs, dgp_args, tau, BIN_WIDTH, SAMPLE_SIZE)
  })
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  message(sprintf("  Completed in %.1f min", elapsed))
}, finally = stopCluster(cl))

res <- bind_rows(iter_results)
saveRDS(res, here("results", "wot_rich_sl.rds"))

# Summary
message("\n=== WOT BIAS/COVERAGE BY SL LIBRARY (truth RD = ",
        sprintf("%+.4f", truth_wot_rd), ") ===\n")
message(sprintf("  %-12s  n_ok  mean_RD    bias       coverage", "SL config"))
message(paste(rep("-", 55), collapse = ""))
for (sl_name in names(sl_configs)) {
  sub <- res[res$sl == sl_name & res$ok == TRUE, ]
  n <- nrow(sub)
  if (n == 0) { message(sprintf("  %-12s  0     ALL FAILED", sl_name)); next }
  mean_rd <- mean(sub$rd)
  bias <- mean_rd - truth_wot_rd
  cov <- mean(sub$lo <= truth_wot_rd & sub$hi >= truth_wot_rd, na.rm = TRUE)
  message(sprintf("  %-12s  %2d    %+.4f    %+.4f     %.0f%%",
    sl_name, n, mean_rd, bias, cov * 100))
}
message("\nSaved to results/wot_rich_sl.rds")
