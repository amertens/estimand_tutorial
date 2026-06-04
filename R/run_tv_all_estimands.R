# run_tv_all_estimands.R
#
# Path 2 sensitivity: under the time-varying-confounder DGP, run all 5 ICH
# E9(R1) estimands x 3 estimator specifications (estimand-aligned Cox,
# L-TMLE baseline-only, L-TMLE +L(t)).
#
# Two DGP variants:
#   - "pure":    delta_trt = 0,   delta_recover = 0   (L(t) is exogenous)
#   - "fb":      delta_trt = 0.10, delta_recover = 0.06 (A->L feedback)
#
# Output:
#   results/sim_tv_all_estimands_pure.rds
#   results/sim_tv_all_estimands_fb.rds
# Per-iter checkpoints:
#   results/sim_partial_tv_pure/iter_<i>.rds
#   results/sim_partial_tv_fb/iter_<i>.rds
# Truth caches:
#   results/sim_results/ground_truth_tv_pure.rds
#   results/sim_results/ground_truth_tv_fb.rds
#
# Usage: Rscript R/run_tv_all_estimands.R

suppressPackageStartupMessages({
  library(here); library(dplyr); library(parallel); library(survival)
  library(lmtp); library(SuperLearner); library(arm)
})
source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("R", "run_simulations.R"))
source(here("R", "dgp_tv_confounder.R"))

tau           <- 180
BIN_WIDTH     <- 14
N_ITER        <- as.integer(Sys.getenv("TV_N_ITER", "15"))   # env-configurable; default 15 pilot
N_CORES       <- as.integer(Sys.getenv("TV_N_CORES", "4"))
CHUNK         <- as.integer(Sys.getenv("TV_CHUNK", "5"))
N_PER_ITER    <- 3000    # reduced from 5000 for tractability
N_TRUTH       <- 200000
LMTP_LEARNERS <- c("SL.mean", "SL.glm", "SL.bayesglm")
# Dropped principal_stratum: subset-based, less stable at N=3000, less
# central to the "Cox biased / lmtp recovers under TV confounding" story.
ESTIMANDS     <- c("treatment_policy", "no_switch", "while_on_treatment",
                   "composite")

# DGP parameter sets
DGP_VARIANTS <- list(
  pure = list(delta_trt = 0.00, delta_recover = 0.00),
  fb   = list(delta_trt = 0.10, delta_recover = 0.06)
)

# ── 1. Truth computation per DGP variant ─────────────────────────────────
compute_truths <- function(variant_args, cache_path) {
  if (file.exists(cache_path)) return(readRDS(cache_path))
  cat("Computing truths on N =", N_TRUTH, "for variant", basename(cache_path), "...\n")

  gen <- function(arm, switching_enabled, seed) {
    do.call(generate_hep_data_tv,
            c(list(N = N_TRUTH, seed = seed, treat_override = arm,
                   switching_enabled = switching_enabled),
              variant_args))
  }

  # Switching-enabled cohorts under each arm (for TP, WOT, composite, PS)
  d_sw_1 <- gen("all_treated", TRUE,  1001)
  d_sw_0 <- gen("all_control", TRUE,  1001)
  # Switching-disabled cohorts (for no-switch truth)
  d_ns_1 <- gen("all_treated", FALSE, 1002)
  d_ns_0 <- gen("all_control", FALSE, 1002)

  # Helper: tau-day risk under different definitions
  risk_tp <- function(d) mean(d$event == 1 & d$follow_time <= tau)
  risk_wot <- function(d) mean(d$event == 1 & d$event_time <= d$switch_time &
                                 d$event_time <= tau)
  risk_composite <- function(d) mean(pmin(d$event_time, d$switch_time) <= tau)
  risk_ps <- function(d) {  # observed-non-switcher approximation
    nev <- d[d$switch_time > tau, ]
    if (nrow(nev) == 0) return(NA_real_)
    mean(nev$event == 1 & nev$follow_time <= tau)
  }

  mk <- function(r1, r0) list(true_risk_1 = r1, true_risk_0 = r0,
                              true_rd = r1 - r0, true_rr = r1 / r0)

  truths <- list(
    treatment_policy   = mk(risk_tp(d_sw_1),       risk_tp(d_sw_0)),
    no_switch          = mk(risk_tp(d_ns_1),       risk_tp(d_ns_0)),
    while_on_treatment = mk(risk_wot(d_sw_1),      risk_wot(d_sw_0)),
    composite          = mk(risk_composite(d_sw_1), risk_composite(d_sw_0)),
    principal_stratum  = mk(risk_ps(d_sw_1),       risk_ps(d_sw_0))
  )
  dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(truths, cache_path)
  cat("  Truths written to", cache_path, "\n")
  for (nm in names(truths)) {
    cat(sprintf("    %-20s RD = %+0.4f  RR = %0.3f\n",
                nm, truths[[nm]]$true_rd, truths[[nm]]$true_rr))
  }
  truths
}

# ── 2. Per-iter analysis: 1 cohort, all 5 estimands, all 3 specs ─────────
tv_all_one_iter <- function(i, variant_args) {
  set.seed(2000 + i)
  dat <- do.call(generate_hep_data_tv,
                 c(list(N = N_PER_ITER, seed = 2000 + i,
                        treat_override = "simulate", switching_enabled = TRUE),
                   variant_args))
  L_cols <- grep("^L[0-9]+$", names(dat), value = TRUE)
  covars <- intersect(c("age", "ckd", "cirrhosis", "diabetes", "heart_failure"),
                      names(dat))

  rows <- list()

  add <- function(estimand, method, hr = NA, hr_lo = NA, hr_hi = NA,
                  rd = NA, rd_lo = NA, rd_hi = NA, conv = TRUE, err = NA) {
    rows[[length(rows) + 1L]] <<- data.frame(
      estimand = estimand, method = method,
      hr = hr, hr_lo = hr_lo, hr_hi = hr_hi,
      rd = rd, rd_lo = rd_lo, rd_hi = rd_hi,
      converged = conv, error_msg = err, stringsAsFactors = FALSE
    )
  }

  # Helper: run lmtp with optional time-varying L
  run_lmtp_one <- function(prep_dat, A_spec, label_suffix, time_vary_L = FALSE) {
    tv_arg <- if (time_vary_L) lapply(seq_along(L_cols), function(k) L_cols[k]) else NULL
    res1 <- lmtp::lmtp_sdr(
      data = prep_dat$data, trt = A_spec, outcome = prep_dat$Y_cols,
      cens = prep_dat$C_cols, baseline = prep_dat$baseline,
      time_vary = tv_arg, outcome_type = "survival", folds = 2,
      learners_trt = LMTP_LEARNERS, learners_outcome = LMTP_LEARNERS,
      shift = lmtp::static_binary_on
    )
    res0 <- lmtp::lmtp_sdr(
      data = prep_dat$data, trt = A_spec, outcome = prep_dat$Y_cols,
      cens = prep_dat$C_cols, baseline = prep_dat$baseline,
      time_vary = tv_arg, outcome_type = "survival", folds = 2,
      learners_trt = LMTP_LEARNERS, learners_outcome = LMTP_LEARNERS,
      shift = lmtp::static_binary_off
    )
    ctr <- lmtp::lmtp_contrast(res1, ref = res0, type = "additive")
    rd_surv <- ctr$estimates$estimate
    rd_se   <- ctr$estimates$std.error
    list(rd = -rd_surv,
         rd_lo = -(rd_surv + 1.96 * rd_se),
         rd_hi = -(rd_surv - 1.96 * rd_se))
  }

  # ── treatment_policy ──
  tryCatch({
    f <- fit_cox_naive(dat, tau = tau, covars = covars)
    add("treatment_policy", "Cox naive", hr = f$hr, hr_lo = f$ci_low, hr_hi = f$ci_high)
  }, error = function(e) add("treatment_policy", "Cox naive", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH, time_varying_trt = FALSE)
    r <- run_lmtp_one(prep, "treatment", "tp_base", FALSE)
    add("treatment_policy", "LMTP (baseline-only)", rd = r$rd, rd_lo = r$rd_lo, rd_hi = r$rd_hi)
  }, error = function(e) add("treatment_policy", "LMTP (baseline-only)", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH, time_varying_trt = FALSE)
    prep$data <- left_join(prep$data, dat[, c("id", L_cols)], by = "id")
    r <- run_lmtp_one(prep, "treatment", "tp_tv", TRUE)
    add("treatment_policy", "LMTP (+L(t))", rd = r$rd, rd_lo = r$rd_lo, rd_hi = r$rd_hi)
  }, error = function(e) add("treatment_policy", "LMTP (+L(t))", conv = FALSE, err = conditionMessage(e)))

  # ── no_switch ──
  tryCatch({
    f <- fit_cox_censor_switch(dat, tau = tau, covars = covars)
    add("no_switch", "Cox censor-at-switch", hr = f$hr, hr_lo = f$ci_low, hr_hi = f$ci_high)
  }, error = function(e) add("no_switch", "Cox censor-at-switch", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    f <- fit_cox_ipcw(dat, tau = tau, covars = covars)
    add("no_switch", "Cox IPCW", hr = f$hr, hr_lo = f$ci_low, hr_hi = f$ci_high)
  }, error = function(e) add("no_switch", "Cox IPCW", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH, time_varying_trt = TRUE)
    r <- run_lmtp_one(prep, prep$A_cols, "ns_base", FALSE)
    add("no_switch", "LMTP (baseline-only)", rd = r$rd, rd_lo = r$rd_lo, rd_hi = r$rd_hi)
  }, error = function(e) add("no_switch", "LMTP (baseline-only)", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH, time_varying_trt = TRUE)
    prep$data <- left_join(prep$data, dat[, c("id", L_cols)], by = "id")
    r <- run_lmtp_one(prep, prep$A_cols, "ns_tv", TRUE)
    add("no_switch", "LMTP (+L(t))", rd = r$rd, rd_lo = r$rd_lo, rd_hi = r$rd_hi)
  }, error = function(e) add("no_switch", "LMTP (+L(t))", conv = FALSE, err = conditionMessage(e)))

  # ── while_on_treatment (cause-specific) ──
  dat_wot <- derive_estimand(dat, "while_on_treatment")
  tryCatch({
    f <- fit_cox_naive(dat_wot, tau = tau, covars = covars)
    add("while_on_treatment", "Cox cause-specific", hr = f$hr, hr_lo = f$ci_low, hr_hi = f$ci_high)
  }, error = function(e) add("while_on_treatment", "Cox cause-specific", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH,
                              time_varying_trt = FALSE, competing_risk_at_switch = TRUE)
    r <- run_lmtp_one(prep, "treatment", "wot_base", FALSE)
    add("while_on_treatment", "LMTP (baseline-only)", rd = r$rd, rd_lo = r$rd_lo, rd_hi = r$rd_hi)
  }, error = function(e) add("while_on_treatment", "LMTP (baseline-only)", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    prep <- prepare_lmtp_data(dat, tau = tau, bin_width = BIN_WIDTH,
                              time_varying_trt = FALSE, competing_risk_at_switch = TRUE)
    prep$data <- left_join(prep$data, dat[, c("id", L_cols)], by = "id")
    r <- run_lmtp_one(prep, "treatment", "wot_tv", TRUE)
    add("while_on_treatment", "LMTP (+L(t))", rd = r$rd, rd_lo = r$rd_lo, rd_hi = r$rd_hi)
  }, error = function(e) add("while_on_treatment", "LMTP (+L(t))", conv = FALSE, err = conditionMessage(e)))

  # ── composite ──
  dat_comp <- derive_estimand(dat, "composite")
  tryCatch({
    f <- fit_cox_naive(dat_comp, tau = tau, covars = covars)
    add("composite", "Cox naive (composite)", hr = f$hr, hr_lo = f$ci_low, hr_hi = f$ci_high)
  }, error = function(e) add("composite", "Cox naive (composite)", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    prep <- prepare_lmtp_data(dat_comp, tau = tau, bin_width = BIN_WIDTH, time_varying_trt = FALSE)
    r <- run_lmtp_one(prep, "treatment", "comp_base", FALSE)
    add("composite", "LMTP (baseline-only)", rd = r$rd, rd_lo = r$rd_lo, rd_hi = r$rd_hi)
  }, error = function(e) add("composite", "LMTP (baseline-only)", conv = FALSE, err = conditionMessage(e)))

  tryCatch({
    prep <- prepare_lmtp_data(dat_comp, tau = tau, bin_width = BIN_WIDTH, time_varying_trt = FALSE)
    prep$data <- left_join(prep$data, dat[, c("id", L_cols)], by = "id")
    r <- run_lmtp_one(prep, "treatment", "comp_tv", TRUE)
    add("composite", "LMTP (+L(t))", rd = r$rd, rd_lo = r$rd_lo, rd_hi = r$rd_hi)
  }, error = function(e) add("composite", "LMTP (+L(t))", conv = FALSE, err = conditionMessage(e)))

  # ── principal_stratum: dropped for tractability (see ESTIMANDS comment) ──

  out <- do.call(rbind, rows)
  out$iter <- i
  out$n    <- N_PER_ITER
  out
}

# ── 3. Driver: chunked checkpointing per DGP variant ─────────────────────
run_variant <- function(variant_name) {
  variant_args <- DGP_VARIANTS[[variant_name]]
  truth_cache <- here("results", "sim_results",
                      sprintf("ground_truth_tv_%s.rds", variant_name))
  partial_dir <- here("results",
                      sprintf("sim_partial_tv_all_%s", variant_name))
  final_cache <- here("results",
                      sprintf("sim_tv_all_estimands_%s.rds", variant_name))
  dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)

  compute_truths(variant_args, truth_cache)

  done_files <- list.files(partial_dir, pattern = "^iter_\\d+\\.rds$",
                           full.names = TRUE)
  done_iters <- as.integer(gsub("[^0-9]", "",
                                basename(tools::file_path_sans_ext(done_files))))
  remaining <- setdiff(seq_len(N_ITER), done_iters)
  cat(sprintf("\n[variant=%s] %d of %d iters remaining\n",
              variant_name, length(remaining), N_ITER))
  if (length(remaining) == 0) {
    all_files <- list.files(partial_dir, pattern = "^iter_\\d+\\.rds$",
                            full.names = TRUE)
    saveRDS(bind_rows(lapply(all_files, readRDS)), final_cache)
    cat(sprintf("[variant=%s] cached %s\n", variant_name, final_cache))
    return(invisible())
  }

  chunks <- split(remaining, ceiling(seq_along(remaining) / CHUNK))
  for (chunk_idx in seq_along(chunks)) {
    iters <- chunks[[chunk_idx]]
    t0 <- Sys.time()
    cat(sprintf("[variant=%s] chunk %d/%d: iters %d..%d (n=%d)\n",
                variant_name, chunk_idx, length(chunks),
                min(iters), max(iters), length(iters)))
    n_workers <- min(N_CORES, length(iters))
    Sys.sleep(2)
    cl <- makeCluster(n_workers)
    chunk_results <- NULL
    tryCatch({
      clusterEvalQ(cl, {
        suppressPackageStartupMessages({
          library(dplyr); library(tibble); library(tidyr)
          library(survival); library(broom); library(here)
          library(lmtp); library(SuperLearner); library(arm)
        })
        source(here::here("DGP.R"))
        source(here::here("R", "helpers.R"))
        source(here::here("R", "run_simulations.R"))
        source(here::here("R", "dgp_tv_confounder.R"))
      })
      clusterExport(cl, c("tv_all_one_iter", "tau", "BIN_WIDTH",
                          "LMTP_LEARNERS", "N_PER_ITER", "variant_args"),
                    envir = environment())
      chunk_results <- clusterApply(cl, iters, function(i)
        tv_all_one_iter(i, variant_args))
    }, error = function(e) {
      message("    [chunk error] ", conditionMessage(e))
    }, finally = {
      try(stopCluster(cl), silent = TRUE)
      gc(verbose = FALSE)
    })
    if (!is.null(chunk_results)) {
      for (k in seq_along(iters)) {
        if (k > length(chunk_results)) next
        out <- chunk_results[[k]]
        if (is.null(out)) next
        tryCatch(saveRDS(out, file.path(partial_dir,
                                        sprintf("iter_%d.rds", iters[[k]]))),
                 error = function(e) message("[save error] iter ", iters[[k]]))
      }
    }
    cat(sprintf("  done in %.1fs\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  all_files <- list.files(partial_dir, pattern = "^iter_\\d+\\.rds$",
                          full.names = TRUE)
  saveRDS(bind_rows(lapply(all_files, readRDS)), final_cache)
  cat(sprintf("[variant=%s] wrote %d rows to %s\n",
              variant_name, length(all_files), final_cache))
}

# ── 4. Run variants sequentially (order/selection via TV_VARIANT env) ─────
# Default order runs the feedback calibration first so its result is available
# before the exogenous calibration finishes.
.tv_variant <- Sys.getenv("TV_VARIANT", "fb,pure")
for (.v in strsplit(.tv_variant, ",")[[1]]) run_variant(trimws(.v))
cat("\nAll variants done.\n")
