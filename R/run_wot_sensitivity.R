# run_wot_sensitivity.R
# Sensitivity analysis for the WOT estimand: isolate the bottleneck.
#
# WOT LMTP achieves only ~48% coverage in the 60-iteration production run
# (bias +0.036, true RD = +0.064). Richer SL libraries barely move the
# needle. This script tests six hypotheses about the source of the bias:
#
#   A  LMTP baseline           Current production spec (reference)
#   B  concrete (continuous)   Continuous-time TMLE via `concrete` package
#                              — avoids discretization; models switch-
#                              censoring as continuous-time intensity
#                              (same estimand as SUSTAIN-6 OT-RC)
#   C  LMTP no admin censor    dep_censor = FALSE in sim data. Tests
#                              whether pooling admin + switch censoring
#                              into one C column is the bottleneck.
#   D  LMTP admin as covar.    dep_censor = TRUE, but the admin-censor
#                              linear predictor is added as a baseline
#                              covariate so the censoring model doesn't
#                              have to learn it. Isolates switching as
#                              the remaining informative component.
#   E  LMTP bin_width = 7      Finer time bins (26 vs 13). Tests whether
#                              coarse discretization smears the switching-
#                              time distribution.
#   F  LMTP N = 20,000         Larger sample. Tests finite-sample bias in
#                              the censoring model vs structural misspec.
#
# Each condition uses the same DGP (h0=5e-3, HR_early=2.0, HR_late=0.90,
# lambda_sw0=1e-3, gamma_A=0.80, gamma_ckd=0.60, complexity=FALSE) except
# where noted.
#
# Usage:
#   Test mode: `Rscript R/run_wot_sensitivity.R` (N_ITER = 5)
#   Edit N_ITER below to scale up once verified.

library(here)
library(dplyr)
library(data.table)
library(survival)
library(parallel)

source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("calc_truth.R"))

tau <- 180
N_ITER <- 20          # Focused rerun: A vs E only, bumped from 5 -> 20
N_CORES <- 4
BIN_WIDTH <- 14
SL_LIB <- c("SL.mean", "SL.glm", "SL.bayesglm")
# Focused rerun: only conditions A (baseline) and E (bin_width=7), since the
# 5-iter pilot showed B/C/D/F not moving the bias and only E was promising.
RUN_CONDITIONS <- c("A", "E")

requireNamespace("lmtp", quietly = TRUE)
requireNamespace("concrete", quietly = TRUE)
requireNamespace("SuperLearner", quietly = TRUE)

# Baseline DGP args (shared across most conditions)
dgp_base <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60
)

# WOT truth (already cached)
truth_all <- readRDS(here("results", "sim_results", "ground_truth.rds"))
truth_wot_rd <- truth_all$while_on_treatment$true_rd
cat(sprintf("WOT truth RD = %+.4f\n", truth_wot_rd))

# ── Helper: LMTP WOT estimator for a given data and options ────────────────
fit_lmtp_wot <- function(dat, bin_width = BIN_WIDTH, baseline_extra = NULL,
                         learners = SL_LIB, folds = 2) {
  baseline <- c("age", "sex_male", "ckd", "cirrhosis", "nsaid",
                "diabetes", "hypertension", "heart_failure", baseline_extra)
  baseline <- baseline[baseline %in% names(dat)]
  prep <- prepare_lmtp_data(dat, tau = tau, bin_width = bin_width,
                            time_varying_trt = FALSE,
                            competing_risk_at_switch = TRUE,
                            baseline = baseline)
  fit <- run_lmtp_analysis(prep, folds = folds, learners = learners)
  rd_obj <- fit$contrast_rd
  rd_s <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
          else rd_obj$vals$theta
  se <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
        else rd_obj$vals$std.error
  # Negate for risk scale; swap CI bounds
  list(rd = -rd_s,
       lo = -(rd_s + 1.96 * se),
       hi = -(rd_s - 1.96 * se),
       ok = TRUE)
}

# ── Helper: concrete WOT estimator ──────────────────────────────────────────
# Model switching as right-censoring (same as SUSTAIN-6 OT-RC).
#   status = 1 if renal failure occurred before switching
#   status = 0 otherwise (switched or admin-censored)
# Do NOT cap wot_time at tau — concrete needs follow-up past the target time
# to know that subjects are still at risk at tau. Let wot_time run to the
# actual min of (event, switch, admin censor), which may exceed tau.
fit_concrete_wot <- function(dat) {
  baseline <- c("age", "sex_male", "ckd", "cirrhosis", "nsaid",
                "diabetes", "hypertension", "heart_failure")
  baseline <- baseline[baseline %in% names(dat)]

  # Construct WOT-derived follow-up: end at min(event, switch, admin)
  # Without the tau cap, subjects who reach admin censor past day 180
  # contribute person-time to the day-180 cumulative incidence.
  dat_wot <- dat %>%
    mutate(
      wot_time = pmin(event_time, switch_time, censor_admin),
      wot_status = as.integer(event_time <= switch_time &
                              event_time <= censor_admin)
    )

  dt <- data.table::as.data.table(
    dat_wot[, c("wot_time", "wot_status", "treatment", baseline)]
  )
  setnames(dt, c("wot_time", "wot_status", "treatment"),
           c("time", "status", "trt"))

  args <- concrete::formatArguments(
    DataTable = dt,
    EventTime = "time", EventType = "status", Treatment = "trt",
    TargetTime = tau, Intervention = 0:1
  )
  # Use concrete's default Model: MainTerms + TrtOnly Cox for cause-specific
  # hazard and censoring; SL.xgboost + SL.glmnet for treatment propensity.
  # Overriding with the string "coxnet" breaks internal parsing.

  est <- concrete::doConcrete(ConcreteArgs = args)
  out <- concrete::getOutput(est, GComp = TRUE,
                             Estimand = c("RD", "Risk"),
                             Simultaneous = FALSE)

  # Pull the TMLE risk-difference estimate at tau for event type 1.
  # Note: concrete labels the RD estimand "Risk Diff" (not "RD") and the
  # estimator "tmle" (lowercase).
  rd_row <- out[out$Estimand == "Risk Diff" & out$Estimator == "tmle" &
                out$Event == 1 & out$Time == tau, ]
  if (nrow(rd_row) == 0) return(list(rd = NA, lo = NA, hi = NA, ok = FALSE))

  list(rd = as.numeric(rd_row$`Pt Est`[1]),
       lo = as.numeric(rd_row$`CI Low`[1]),
       hi = as.numeric(rd_row$`CI Hi`[1]),
       ok = TRUE)
}

# ── Condition-specific data generators ──────────────────────────────────────
gen_data <- function(condition, seed, n_override = NULL) {
  n <- if (!is.null(n_override)) n_override else 10000
  args <- modifyList(dgp_base,
                     list(N = n, switch_on = TRUE, seed = seed))
  if (condition == "C") args$dep_censor <- FALSE
  else                  args$dep_censor <- TRUE
  dat <- do.call(generate_hep_data, args)

  # For condition D: precompute admin-censoring linear predictor
  if (condition == "D") {
    # Admin censor rate: c_rate = censor_base * exp(0.04*lp_out + 0.03*trt)
    # lp_out for complexity=FALSE: -2.8 + 0.03*age + 0.7*ckd + 0.5*cirrhosis
    #                              + 0.3*heart_failure + 0.25*nsaid
    lp_out <- with(dat,
                   -2.8 + 0.03 * age + 0.7 * ckd + 0.5 * cirrhosis +
                     0.3 * heart_failure + 0.25 * nsaid)
    dat$admin_cens_pred <- 0.04 * lp_out + 0.03 * dat$treatment
  }
  dat
}

# ── Worker: run one iteration across the requested conditions ──────────────
# RUN_CONDITIONS controls which sub-experiments fire on each iteration.
# Default ("A","E") focuses on baseline vs finer-bin discretization.
run_one <- function(i, tau, truth_wot_rd, dgp_base, BIN_WIDTH, SL_LIB,
                    RUN_CONDITIONS = c("A", "B", "C", "D", "E", "F")) {
  out <- list()

  if ("A" %in% RUN_CONDITIONS) {
    out[["A_LMTP_base"]] <- tryCatch(
      fit_lmtp_wot(gen_data("A", seed = 5000 + i)),
      error = function(e) list(rd = NA, lo = NA, hi = NA, ok = FALSE))
  }
  if ("B" %in% RUN_CONDITIONS) {
    out[["B_concrete"]] <- tryCatch(
      fit_concrete_wot(gen_data("B", seed = 5000 + i)),
      error = function(e) list(rd = NA, lo = NA, hi = NA, ok = FALSE))
  }
  if ("C" %in% RUN_CONDITIONS) {
    out[["C_LMTP_no_admin"]] <- tryCatch(
      fit_lmtp_wot(gen_data("C", seed = 5000 + i)),
      error = function(e) list(rd = NA, lo = NA, hi = NA, ok = FALSE))
  }
  if ("D" %in% RUN_CONDITIONS) {
    out[["D_LMTP_admin_covar"]] <- tryCatch(
      fit_lmtp_wot(gen_data("D", seed = 5000 + i),
                   baseline_extra = "admin_cens_pred"),
      error = function(e) list(rd = NA, lo = NA, hi = NA, ok = FALSE))
  }
  if ("E" %in% RUN_CONDITIONS) {
    out[["E_LMTP_bw7"]] <- tryCatch(
      fit_lmtp_wot(gen_data("E", seed = 5000 + i), bin_width = 7),
      error = function(e) list(rd = NA, lo = NA, hi = NA, ok = FALSE))
  }
  if ("F" %in% RUN_CONDITIONS) {
    out[["F_LMTP_N20k"]] <- tryCatch(
      fit_lmtp_wot(gen_data("F", seed = 5000 + i, n_override = 20000)),
      error = function(e) list(rd = NA, lo = NA, hi = NA, ok = FALSE))
  }

  do.call(rbind, lapply(names(out), function(nm) {
    r <- out[[nm]]
    data.frame(iter = i, condition = nm, rd = r$rd, lo = r$lo, hi = r$hi,
               ok = r$ok, stringsAsFactors = FALSE)
  }))
}

# ── Parallel execution ─────────────────────────────────────────────────────
cat(sprintf("\nRunning %d iterations x 6 conditions on %d cores...\n",
            N_ITER, N_CORES))

cl <- makeCluster(N_CORES)
tryCatch({
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(dplyr); library(data.table); library(survival); library(here)
      library(lmtp); library(SuperLearner); library(arm); library(concrete)
    })
    source(here("DGP.R")); source(here("R", "helpers.R"))
  })
  clusterExport(cl, c("run_one", "fit_lmtp_wot", "fit_concrete_wot",
                      "gen_data", "tau", "truth_wot_rd", "dgp_base",
                      "BIN_WIDTH", "SL_LIB", "RUN_CONDITIONS"),
                envir = environment())

  t0 <- Sys.time()
  iter_results <- clusterApply(cl, seq_len(N_ITER), function(i) {
    run_one(i, tau, truth_wot_rd, dgp_base, BIN_WIDTH, SL_LIB,
            RUN_CONDITIONS = RUN_CONDITIONS)
  })
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  cat(sprintf("  Completed in %.1f min\n", elapsed))
}, finally = stopCluster(cl))

res <- dplyr::bind_rows(iter_results)
saveRDS(res, here("results", "wot_sensitivity.rds"))

# ── Summary ────────────────────────────────────────────────────────────────
cat(sprintf("\n=== WOT SENSITIVITY (truth RD = %+.4f) ===\n\n", truth_wot_rd))
cat(sprintf("  %-22s  %4s  %8s  %8s  %7s\n",
            "Condition", "n_ok", "Mean RD", "Bias", "Coverage"))
cat(paste(rep("-", 60), collapse = ""), "\n")

cond_full <- c("A_LMTP_base", "B_concrete", "C_LMTP_no_admin",
               "D_LMTP_admin_covar", "E_LMTP_bw7", "F_LMTP_N20k")
cond_order <- intersect(cond_full, unique(res$condition))
for (cn in cond_order) {
  s <- res[res$condition == cn & res$ok, ]
  n <- nrow(s)
  if (n == 0) {
    cat(sprintf("  %-22s  %4d  %8s  %8s  %7s\n", cn, 0, "—", "—", "ALL FAIL"))
    next
  }
  mean_rd <- mean(s$rd, na.rm = TRUE)
  bias <- mean_rd - truth_wot_rd
  cov <- mean(s$lo <= truth_wot_rd & s$hi >= truth_wot_rd, na.rm = TRUE)
  cat(sprintf("  %-22s  %4d  %+.4f    %+.4f    %5.0f%%\n",
              cn, n, mean_rd, bias, cov * 100))
}
cat(sprintf("\nSaved to results/wot_sensitivity.rds\n"))
