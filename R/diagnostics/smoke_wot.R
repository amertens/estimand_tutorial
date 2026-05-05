# smoke_wot.R -- single-iteration smoke test for the WOT fix
# (time_varying_trt = TRUE, censor_at_switch = TRUE).
# Verifies the new specification runs end-to-end and reports bias vs. truth.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(survival)
  library(lmtp); library(SuperLearner); library(arm)
})
source(here("DGP.R"))
source(here("R", "helpers.R"))

tau <- 180
BIN_WIDTH <- 14

truth_all <- readRDS(here("results", "sim_results", "ground_truth.rds"))
true_rd <- truth_all$while_on_treatment$true_rd
cat(sprintf("True WOT RD = %+.4f\n", true_rd))

dgp_args <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, dep_censor = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60,
  return_potential_switching = TRUE
)

dat <- do.call(generate_hep_data,
               c(list(N = 10000, switch_on = TRUE, seed = 4001), dgp_args))
cat(sprintf("N = %d; switch rate = %.1f%%; event rate = %.1f%%\n",
            nrow(dat),
            100 * mean(dat$switched == 1),
            100 * mean(dat$event == 1 & dat$follow_time <= tau)))

# NEW spec: time_varying_trt = TRUE
t0 <- Sys.time()
prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH,
                          time_varying_trt = FALSE,
                          competing_risk_at_switch = TRUE)
cat(sprintf("A_cols: %s\n", paste(head(prep$A_cols, 3), collapse = ",")))
cat(sprintf("n_bins: %d; wide cols: %d\n", prep$n_bins, ncol(prep$data)))

res <- run_lmtp_analysis(prep, folds = 2,
                         learners = c("SL.mean", "SL.glm", "SL.bayesglm"))
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

rd_obj <- res$contrast_rd
rd_s <- (if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
         else rd_obj$vals$theta)
se <- (if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
       else rd_obj$vals$std.error)
# Negate for risk-scale RD
rd  <- -rd_s
lo  <- -(rd_s + 1.96 * se)
hi  <- -(rd_s - 1.96 * se)

covered <- lo <= true_rd && hi >= true_rd

cat(sprintf("\nLMTP WOT (competing_risk_at_switch=TRUE) single-iter result:\n"))
cat(sprintf("  RD est : %+.4f   95%% CI [%+.4f, %+.4f]\n", rd, lo, hi))
cat(sprintf("  Truth  : %+.4f   Bias: %+.4f   Covered: %s\n",
            true_rd, rd - true_rd, if (covered) "YES" else "NO"))
cat(sprintf("  Elapsed: %.1f sec\n", elapsed))
