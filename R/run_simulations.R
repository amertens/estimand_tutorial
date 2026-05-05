# run_simulations.R
# Repeated simulation study comparing estimators across the five
# ICH E9(R1) estimand strategies.
#
# Per iteration and estimand, the following methods are run:
#   - Cox naive (KM RD + Cox g-comp RD)
#   - Cox censor-at-switch (KM RD + Cox g-comp RD)
#   - Cox IPCW (KM RD + Cox g-comp RD)
#   - LMTP SDR (one fit; principal_stratum runs two fits, oracle + approx)
#
# Estimand-specific LMTP specifications:
#   treatment_policy:   baseline-only treatment, no switch handling
#   no_switch:          time-varying A_j, static intervention holds A constant
#   while_on_treatment: baseline-only treatment, switching as a competing
#                       event in the outcome (competing_risk_at_switch=TRUE);
#                       targets the M1 / Rufibach cause-specific incidence.
#   composite:          baseline-only treatment, composite outcome
#   principal_stratum:  baseline-only treatment on (a) the oracle
#                       never-switcher subset and (b) observed non-switchers.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(broom)
  library(here)
  library(parallel)
})

source(here("DGP.R"))
source(here("R", "helpers.R"))

# Number of parallel workers. Set to 1 for serial (with progress bar).
N_CORES <- 8L

# ── Single iteration for one estimand ────────────────────────────────────────
#' Run one iteration of the simulation study for a specific estimand.
#' @param i integer; iteration index.
#' @param estimand character; "treatment_policy" or "no_switch".
#' @param sample_size integer; per-iteration sample size.
#' @param tau integer; follow-up horizon.
#' @param run_lmtp logical; whether to run LMTP (slow).
#' @param run_cox_td logical; whether to run Cox time-dependent (slow for large N).
#' @param dgp_args list of additional arguments passed to generate_hep_data().
#' @return data.frame with one row per method.
run_one_iter <- function(i, estimand = "treatment_policy",
                         sample_size = 1000, tau = 180,
                         bin_width = 14,
                         lmtp_learners = c("SL.mean", "SL.glm", "SL.bayesglm", "SL.glmnet", "SL.ranger"),
                         lmtp_k = Inf,
                         run_lmtp = TRUE, run_cox_td = FALSE,
                         dgp_args = list()) {
  set.seed(1000 + i)

  # Generate dataset(s) for this iteration.
  #
  # For no_switch: we need TWO datasets with the same seed:
  #   - dat_noswitch (switch_on=FALSE): used by Cox, where truth = outcomes
  #     under sustained treatment with no switching at all.
  #   - dat_raw (switch_on=TRUE): used by LMTP, which observes the switching
  #     process and intervenes to hold treatment constant via time-varying A_j.
  #     LMTP needs to see the observed switching to learn the switching model,
  #     then overrides it with the static intervention.
  #
  # For all other estimands: one dataset with switch_on=TRUE.
  args_sw <- modifyList(
    list(N = sample_size, switch_on = TRUE, seed = 1000 + i),
    dgp_args
  )
  dat_raw <- do.call(generate_hep_data, args_sw)

  if (estimand == "no_switch") {
    args_nosw <- modifyList(
      list(N = sample_size, switch_on = FALSE, seed = 1000 + i),
      dgp_args
    )
    dat_noswitch <- do.call(generate_hep_data, args_nosw)
    dat_cox <- dat_noswitch
  } else if (estimand %in% c("while_on_treatment", "composite")) {
    dat_cox <- derive_estimand(dat_raw, estimand)
  } else {
    dat_cox <- dat_raw
  }

  covars <- intersect(
    c("age", "ckd", "cirrhosis", "diabetes", "heart_failure"),
    names(dat_raw)
  )

  results <- list()

  # Helper: build Cox result row with HR + HR CI only.
  # We deliberately do NOT report a Cox-derived risk difference. The
  # unadjusted KM RD does not use the Cox fit at all (the Cox model
  # contributes nothing beyond the HR), and the covariate-standardized
  # Cox g-comp RD inherits the proportional-hazards assumption that
  # the simulation DGP intentionally violates. Cox is judged on its
  # native scale (HR vs true marginal HR); LMTP is judged on the
  # risk-scale parameters (RD and RR) it natively targets.
  make_cox_row <- function(fit, label) {
    data.frame(
      method = label,
      hr = fit$hr, hr_ci_low = fit$ci_low, hr_ci_high = fit$ci_high,
      risk_diff = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      risk_ratio = NA_real_, rr_ci_low = NA_real_, rr_ci_high = NA_real_,
      converged = TRUE
    )
  }

  cox_fail <- function(label) {
    data.frame(
      method = label,
      hr = NA_real_, hr_ci_low = NA_real_, hr_ci_high = NA_real_,
      risk_diff = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      risk_ratio = NA_real_, rr_ci_low = NA_real_, rr_ci_high = NA_real_,
      converged = FALSE
    )
  }

  # ── A. Naive Cox (ignores switching) ──
  # Cox is judged on its native scale: the HR. We do not report a
  # Cox-derived risk difference; see make_cox_row() rationale above.
  cox_naive_fit <- tryCatch(
    fit_cox_naive(dat_cox, tau = tau, covars = covars),
    error = function(e) NULL
  )
  results[["cox_naive"]] <- if (!is.null(cox_naive_fit)) {
    make_cox_row(cox_naive_fit, "Cox naive")
  } else cox_fail("Cox naive")

  # ── B. Cox censor at switch ──
  cox_cs_fit <- tryCatch(
    fit_cox_censor_switch(dat_cox, tau = tau, covars = covars),
    error = function(e) NULL
  )
  results[["cox_censor"]] <- if (!is.null(cox_cs_fit)) {
    make_cox_row(cox_cs_fit, "Cox censor-at-switch")
  } else cox_fail("Cox censor-at-switch")

  # ── C. Cox IPCW (censor at switch with inverse probability weights) ──
  cox_ipcw_fit <- tryCatch(
    fit_cox_ipcw(dat_cox, tau = tau, covars = covars),
    error = function(e) NULL
  )
  results[["cox_ipcw"]] <- if (!is.null(cox_ipcw_fit)) {
    make_cox_row(cox_ipcw_fit, "Cox IPCW")
  } else cox_fail("Cox IPCW")

  # ── D. Cox time-dependent treatment (optional) ──
  if (run_cox_td) {
    results[["cox_td"]] <- tryCatch({
      fit <- fit_cox_td(dat_cox, tau = tau, covars = covars)
      data.frame(
        method = "Cox time-dependent",
        hr = fit$hr, hr_ci_low = fit$ci_low, hr_ci_high = fit$ci_high,
        risk_diff = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        risk_ratio = NA_real_, rr_ci_low = NA_real_, rr_ci_high = NA_real_,
        converged = TRUE
      )
    }, error = function(e) cox_fail("Cox time-dependent"))
  }

  # ── D. LMTP SDR ──
  # Estimand-specific dispatch (see estimand_truth_lmtp_audit.md):
  #
  #   treatment_policy: BASELINE-ONLY treatment (time_varying_trt=FALSE).
  #     Assigns A at baseline; LMTP censoring model handles informative
  #     follow-up loss from switching. static_binary_on/off sets A=1/0.
  #
  #   no_switch: TIME-VARYING A_j (time_varying_trt=TRUE) on dat_raw.
  #     A_j reflects observed switching trajectory. static_binary_on/off
  #     OVERRIDES all A_j to hold treatment constant = no-switch intervention.
  #
  #   while_on_treatment: BASELINE-ONLY treatment, censor_at_switch=TRUE.
  #     Switching enters through C columns (censoring nodes), not A columns.
  #     LMTP censoring model adjusts for informative switch-censoring.
  #
  #   composite: BASELINE-ONLY treatment (time_varying_trt=FALSE) on dat_cox.
  #     Switching IS part of the outcome; should not intervene on A_j.
  #
  #   principal_stratum: BASELINE-ONLY on subsets of dat_raw.
  #     Never-switchers have constant A_j anyway.

  lmtp_fail <- function(label) {
    data.frame(
      method = label,
      hr = NA_real_, hr_ci_low = NA_real_, hr_ci_high = NA_real_,
      risk_diff = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      risk_ratio = NA_real_, rr_ci_low = NA_real_, rr_ci_high = NA_real_,
      converged = FALSE
    )
  }

  # Helper to extract LMTP contrasts.
  # Note: lmtp_contrast(type = "additive") returns S(1) - S(0) on the
  # survival scale. The risk-scale RD is R(1) - R(0) = -(S(1) - S(0)),
  # so we negate the point estimate and swap the CI bounds.
  # The risk-scale RR (R(1)/R(0)) and its CI are computed in
  # run_lmtp_analysis via an EIF-based delta method (see helpers.R).
  extract_lmtp <- function(res, label) {
    rd_obj <- res$contrast_rd
    rd_surv <- (if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
                else rd_obj$vals$theta)
    rd_se   <- (if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
                else rd_obj$vals$std.error)
    rd_est <- -rd_surv
    rd_ci  <- -(rd_surv + c(1, -1) * 1.96 * rd_se)

    rr_est     <- res$rr
    rr_ci_low  <- res$rr_ci_low
    rr_ci_high <- res$rr_ci_high

    data.frame(
      method = label,
      hr = NA_real_, hr_ci_low = NA_real_, hr_ci_high = NA_real_,
      risk_diff = rd_est, ci_low = rd_ci[1], ci_high = rd_ci[2],
      risk_ratio = rr_est, rr_ci_low = rr_ci_low, rr_ci_high = rr_ci_high,
      converged = TRUE
    )
  }

  if (run_lmtp) {
    requireNamespace("lmtp", quietly = TRUE)

    if (estimand == "treatment_policy") {
      # ── Baseline-only treatment: let switching happen naturally ──
      results[["lmtp"]] <- tryCatch({
        n_events <- sum(dat_raw$event == 1 & dat_raw$follow_time <= tau)
        if (n_events < 5) stop("Too few events: ", n_events)
        prep <- prepare_lmtp_data(dat_raw, tau = tau, bin_width = bin_width,
                                  time_varying_trt = FALSE)
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners, k = lmtp_k)
        extract_lmtp(res, "LMTP SDR")
      }, error = function(e) lmtp_fail("LMTP SDR"))

    } else if (estimand == "no_switch") {
      # ── Time-varying A_j on switching data: static intervention holds constant ──
      results[["lmtp"]] <- tryCatch({
        n_events <- sum(dat_raw$event == 1 & dat_raw$follow_time <= tau)
        if (n_events < 5) stop("Too few events: ", n_events)
        prep <- prepare_lmtp_data(dat_raw, tau = tau, bin_width = bin_width,
                                  time_varying_trt = TRUE)
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners, k = lmtp_k)
        extract_lmtp(res, "LMTP SDR (no-switch)")
      }, error = function(e) lmtp_fail("LMTP SDR (no-switch)"))

    } else if (estimand == "while_on_treatment") {
      # ── Baseline-only treatment, switching as a competing event ──
      # WOT here is the M1 / Rufibach interpretation: the cause-specific
      # cumulative incidence of the outcome before switching by tau, under
      # an intervention on baseline treatment. Switching is absorbed into
      # the outcome definition (Y locks at 0 once a subject switches
      # without a prior event); only admin censoring before any switch
      # contributes to the C nodes. This aligns the LMTP estimator with
      # the truth in calc_truth.R::truth_while_on_treatment.
      #
      # Earlier iterations of this block used censor_at_switch = TRUE
      # (Latimer / IPCW interpretation), which targets the no-switch
      # counterfactual under non-informative censoring -- a different
      # estimand from the truth. That mismatch produced the apparent
      # +0.035 "bias" reported in earlier drafts; see Discussion.
      results[["lmtp"]] <- tryCatch({
        prep <- prepare_lmtp_data(dat_raw, tau = tau, bin_width = bin_width,
                                  time_varying_trt = FALSE,
                                  competing_risk_at_switch = TRUE)
        n_events <- sum(prep$data$Y1 == 1, na.rm = TRUE)  # events in first bin
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners, k = lmtp_k)
        extract_lmtp(res, "LMTP SDR (WOT)")
      }, error = function(e) lmtp_fail("LMTP SDR (WOT)"))

    } else if (estimand == "composite") {
      # ── Baseline-only treatment on composite-derived data ──
      results[["lmtp"]] <- tryCatch({
        n_events <- sum(dat_cox$event == 1 & dat_cox$follow_time <= tau)
        if (n_events < 5) stop("Too few events: ", n_events)
        prep <- prepare_lmtp_data(dat_cox, tau = tau, bin_width = bin_width,
                                  time_varying_trt = FALSE)
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners, k = lmtp_k)
        extract_lmtp(res, "LMTP SDR (composite)")
      }, error = function(e) lmtp_fail("LMTP SDR (composite)"))

    } else if (estimand == "principal_stratum") {
      # ── Oracle principal stratum (primary) ──
      # Restrict to true never-switchers (never_switcher == 1) identified
      # via paired potential switching draws. This is an oracle analysis
      # using simulation-only information not available in practice.
      if ("never_switcher" %in% names(dat_raw)) {
        results[["lmtp"]] <- tryCatch(
          suppressWarnings({
            sub <- dat_raw[dat_raw$never_switcher == 1, ]
            n_events <- sum(sub$event == 1 & sub$follow_time <= tau)
            if (n_events < 5) stop("Too few events: ", n_events)
            prep <- prepare_lmtp_data(sub, tau = tau, bin_width = bin_width,
                                      time_varying_trt = FALSE)
            res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners, k = lmtp_k)
            extract_lmtp(res, "LMTP SDR (oracle PS)")
          }),
          error = function(e) lmtp_fail("LMTP SDR (oracle PS)")
        )
      }

      # ── Observed non-switcher approximation (supplementary) ──
      results[["lmtp_approx"]] <- tryCatch({
        sub <- dat_raw[dat_raw$switched == 0, ]
        n_events <- sum(sub$event == 1 & sub$follow_time <= tau)
        if (n_events < 5) stop("Too few events: ", n_events)
        prep <- prepare_lmtp_data(sub, tau = tau, bin_width = bin_width,
                                  time_varying_trt = FALSE)
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners, k = lmtp_k)
        extract_lmtp(res, "LMTP SDR (obs. non-sw, approx)")
      }, error = function(e) lmtp_fail("LMTP SDR (obs. non-sw, approx)"))
    }
  }

  out <- bind_rows(results)
  out$estimand <- estimand
  out$iter     <- i
  out$n        <- sample_size
  out
}


# ── Full simulation runner ───────────────────────────────────────────────────
#' Run the repeated simulation study across estimands.
#' @param n_iter integer; number of iterations (default 200).
#' @param sample_size integer; sample size per iteration.
#' @param tau integer; follow-up horizon.
#' @param estimands character vector; which estimands to run.
#' @param run_lmtp logical.
#' @param run_cox_td logical.
#' @param dgp_args list; extra args for generate_hep_data().
#' @param cache_file character or NULL; path to cache results.
#' @return data.frame of aggregated results.
run_simulation_study <- function(n_iter = 200, sample_size = 1000,
                                 tau = 180, bin_width = 14,
                                 lmtp_learners = c("SL.mean", "SL.glm", "SL.bayesglm", "SL.glmnet", "SL.ranger"),
                                 lmtp_k = Inf,
                                 estimands = c("treatment_policy", "no_switch"),
                                 run_lmtp = TRUE, run_cox_td = FALSE,
                                 n_cores = N_CORES,
                                 dgp_args = list(),
                                 cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file)) {
    message("Loading cached results from ", cache_file)
    return(readRDS(cache_file))
  }

  total_iters <- n_iter * length(estimands)
  message("Running ", n_iter, " iterations x ", length(estimands), " estimands (",
          total_iters, " total) on ", n_cores, " cores...")

  all_res <- NULL

  for (est in estimands) {
    message("  Estimand: ", est)

    if (n_cores > 1) {
      # ── Parallel via PSOCK cluster (Windows-compatible) ──
      # Each worker sources DGP.R and helpers.R for data generation and
      # analysis functions. run_one_iter is defined as a standalone
      # wrapper function on each worker to avoid closure serialization
      # issues. Arguments are passed explicitly via clusterExport.
      cl <- makeCluster(n_cores)

      tryCatch({
        # Set up workers: packages + source files
        clusterEvalQ(cl, {
          suppressPackageStartupMessages({
            library(dplyr); library(survival); library(broom)
            library(here); library(lmtp); library(SuperLearner); library(arm)
          })
          source(here("DGP.R"))
          source(here("R", "helpers.R"))
        })

        # Export run_one_iter and all arguments as plain objects
        # to workers' global environments. No closures.
        clusterExport(cl, "run_one_iter", envir = environment())
        worker_args <- list(
          est = est, sample_size = sample_size, tau = tau,
          bin_width = bin_width, lmtp_learners = lmtp_learners,
          lmtp_k = lmtp_k,
          run_lmtp = run_lmtp, run_cox_td = run_cox_td,
          dgp_args = dgp_args
        )
        clusterExport(cl, "worker_args", envir = environment())

        # Use clusterApply with a minimal wrapper — no closure captures
        iter_results <- clusterApply(cl, seq_len(n_iter), function(i) {
          run_one_iter(
            i, estimand = worker_args$est,
            sample_size = worker_args$sample_size,
            tau = worker_args$tau,
            bin_width = worker_args$bin_width,
            lmtp_learners = worker_args$lmtp_learners,
            lmtp_k = worker_args$lmtp_k,
            run_lmtp = worker_args$run_lmtp,
            run_cox_td = worker_args$run_cox_td,
            dgp_args = worker_args$dgp_args
          )
        })

        est_res <- dplyr::bind_rows(iter_results)
      }, error = function(e) {
        message("    Parallel failed: ", conditionMessage(e))
        message("    Falling back to serial execution.")
        est_res <<- NULL
      }, finally = {
        stopCluster(cl)
      })

      # If parallel failed, fall back to serial
      if (is.null(est_res)) {
        pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
        est_res <- NULL
        for (i in seq_len(n_iter)) {
          iter_res <- run_one_iter(
            i, estimand = est, sample_size = sample_size,
            tau = tau, bin_width = bin_width,
            lmtp_learners = lmtp_learners,
            lmtp_k = lmtp_k,
            run_lmtp = run_lmtp, run_cox_td = run_cox_td,
            dgp_args = dgp_args
          )
          est_res <- bind_rows(est_res, iter_res)
          setTxtProgressBar(pb, i)
        }
        close(pb)
      }
    } else {
      # ── Serial fallback ──
      pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
      est_res <- NULL
      for (i in seq_len(n_iter)) {
        iter_res <- run_one_iter(
          i, estimand = est, sample_size = sample_size,
          tau = tau, bin_width = bin_width,
          lmtp_learners = lmtp_learners,
          run_lmtp = run_lmtp, run_cox_td = run_cox_td,
          dgp_args = dgp_args
        )
        est_res <- bind_rows(est_res, iter_res)
        setTxtProgressBar(pb, i)
      }
      close(pb)
    }

    message("    Done: ", nrow(est_res), " rows")
    all_res <- bind_rows(all_res, est_res)
  }

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(all_res, cache_file)
    message("Results cached to ", cache_file)
  }

  all_res
}


# ── Support scenario runner ──────────────────────────────────────────────────
#' Run estimators under three support scenarios for a given estimand.
#' @param estimand character; "treatment_policy" or "no_switch".
#' @param n_iter integer; iterations per scenario.
#' @param sample_size integer.
#' @param tau integer.
#' @param cache_file character or NULL.
#' @return data.frame with scenario column.
run_support_scenarios <- function(estimand = "treatment_policy",
                                  n_iter = 50, sample_size = 2000,
                                  tau = 180, cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file)) {
    message("Loading cached support results from ", cache_file)
    return(readRDS(cache_file))
  }

  scenarios <- list(
    "Good support" = list(),
    "Strained support" = list(gamma_A = 1.5, gamma_ckd = 1.2, lambda_sw0 = 5e-5),
    "Poor support" = list(gamma_A = 2.5, gamma_ckd = 2.0, lambda_sw0 = 1e-4)
  )

  all_res <- NULL
  for (sc_name in names(scenarios)) {
    message("  Scenario: ", sc_name)
    sc_args <- scenarios[[sc_name]]
    for (i in seq_len(n_iter)) {
      if (i %% 10 == 0) message("    iteration ", i, " / ", n_iter)
      iter_res <- run_one_iter(
        i, estimand = estimand, sample_size = sample_size,
        tau = tau, run_lmtp = TRUE, run_cox_td = FALSE,
        dgp_args = sc_args
      )
      iter_res$scenario <- sc_name
      all_res <- bind_rows(all_res, iter_res)
    }
  }

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(all_res, cache_file)
    message("Cached to ", cache_file)
  }

  all_res
}


# ── Summarize simulation results ─────────────────────────────────────────────
#' Summarize simulation results: bias, empirical SE, RMSE, coverage.
#' @param sim_results data.frame from run_simulation_study().
#' @param truth_rd numeric; true risk difference (for LMTP and KM-based RD).
#' @param truth_hr numeric; true hazard ratio (for Cox, if known).
#' @return tibble with summary statistics grouped by estimand and method.
summarize_simulation <- function(sim_results, truth_rd = NA, truth_rr = NA,
                                 truth_hr = NA) {
  has_rr_ci <- all(c("rr_ci_low", "rr_ci_high") %in% names(sim_results))
  has_hr_ci <- all(c("hr_ci_low", "hr_ci_high") %in% names(sim_results))
  sim_results %>%
    group_by(estimand, method) %>%
    summarise(
      n_iter       = n(),
      n_converged  = sum(converged, na.rm = TRUE),
      # HR summaries (Cox methods only)
      mean_hr      = mean(hr, na.rm = TRUE),
      sd_hr        = sd(hr, na.rm = TRUE),
      bias_hr      = if (!is.na(truth_hr)) mean(hr, na.rm = TRUE) - truth_hr
                     else NA_real_,
      coverage_hr  = if (!is.na(truth_hr) && has_hr_ci)
                       mean(hr_ci_low <= truth_hr & hr_ci_high >= truth_hr,
                            na.rm = TRUE)
                     else NA_real_,
      # RD summaries (LMTP only)
      mean_rd      = mean(risk_diff, na.rm = TRUE),
      emp_se_rd    = sd(risk_diff, na.rm = TRUE),
      bias_rd      = if (!is.na(truth_rd)) mean(risk_diff, na.rm = TRUE) - truth_rd
                     else NA_real_,
      rmse_rd      = if (!is.na(truth_rd))
                       sqrt(mean((risk_diff - truth_rd)^2, na.rm = TRUE))
                     else NA_real_,
      coverage_rd  = if (!is.na(truth_rd))
                       mean(ci_low <= truth_rd & ci_high >= truth_rd,
                            na.rm = TRUE)
                     else NA_real_,
      # RR summaries (LMTP only; risk-scale RR with EIF-based delta-method CI)
      mean_rr      = mean(risk_ratio, na.rm = TRUE),
      sd_rr        = sd(risk_ratio, na.rm = TRUE),
      bias_rr      = if (!is.na(truth_rr)) mean(risk_ratio, na.rm = TRUE) - truth_rr
                     else NA_real_,
      coverage_rr  = if (!is.na(truth_rr) && has_rr_ci)
                       mean(rr_ci_low <= truth_rr & rr_ci_high >= truth_rr,
                            na.rm = TRUE)
                     else NA_real_,
      .groups = "drop"
    )
}
