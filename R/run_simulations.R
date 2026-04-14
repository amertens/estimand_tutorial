# run_simulations.R
# Repeated simulation study comparing estimators under each DGP policy.
#
# For each DGP policy (treatment_policy, no_switch, while_on_treatment,
# composite, principal_stratum), data are generated and three estimators
# are applied:
#   1. Baseline-treatment Cox HR (ignores switching)
#   2. Censor-at-switch Cox HR
#   3. LMTP SDR (static baseline intervention)
#
# The LMTP analysis uses the same static intervention specification for
# most policies. For the principal stratum, it restricts to observed
# non-switchers (an approximation). The composite LMTP is known to be
# misspecified (the censoring model adjusts away switching, whereas the
# composite counts it as an event).

# TODO(Joy): fill in exact runtime notes after benchmarking.

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
                         lmtp_learners = c("SL.mean", "SL.glm", "SL.bayesglm"),
                         run_lmtp = TRUE, run_cox_td = FALSE,
                         dgp_args = list()) {
  set.seed(1000 + i)

  # Generate ONE dataset under treatment-policy (switching occurs, does not censor).
  # Derive estimand-specific outcomes via derive_estimand().
  # For no_switch: generate separately with switch_on=FALSE.
  if (estimand == "no_switch") {
    args <- modifyList(
      list(N = sample_size, switch_on = FALSE, seed = 1000 + i),
      dgp_args
    )
  } else {
    args <- modifyList(
      list(N = sample_size, switch_on = TRUE, seed = 1000 + i),
      dgp_args
    )
  }
  dat_raw <- do.call(generate_hep_data, args)

  # Apply estimand-specific outcome/censoring definition for Cox analyses.
  # LMTP uses dat_raw (or dat_raw subsets) directly — it handles switching
  # through its own treatment/censoring model, not through pre-derived
  # follow_time/event columns.
  if (estimand %in% c("while_on_treatment", "composite")) {
    dat_cox <- derive_estimand(dat_raw, estimand)
  } else {
    dat_cox <- dat_raw
  }

  covars <- intersect(
    c("age", "ckd", "cirrhosis", "diabetes", "heart_failure"),
    names(dat_raw)
  )

  results <- list()

  # ── A. Naive Cox (ignores switching) ──
  results[["cox_naive"]] <- tryCatch({
    fit <- fit_cox_naive(dat_cox, tau = tau, covars = covars)
    km  <- km_risk_at(fit$survfit, t = tau)
    km_rr <- if (km$risk_ctrl > 0) km$risk_trt / km$risk_ctrl else NA_real_
    data.frame(
      method = "Cox naive", hr = fit$hr,
      ci_low = fit$ci_low, ci_high = fit$ci_high,
      risk_diff = km$risk_diff, risk_ratio = km_rr, converged = TRUE
    )
  }, error = function(e) {
    data.frame(method = "Cox naive", hr = NA, ci_low = NA, ci_high = NA,
               risk_diff = NA, risk_ratio = NA, converged = FALSE)
  })

  # ── B. Cox censor at switch ──
  results[["cox_censor"]] <- tryCatch({
    fit <- fit_cox_censor_switch(dat_cox, tau = tau, covars = covars)
    km  <- km_risk_at(fit$survfit, t = tau)
    km_rr <- if (km$risk_ctrl > 0) km$risk_trt / km$risk_ctrl else NA_real_
    data.frame(
      method = "Cox censor-at-switch", hr = fit$hr,
      ci_low = fit$ci_low, ci_high = fit$ci_high,
      risk_diff = km$risk_diff, risk_ratio = km_rr, converged = TRUE
    )
  }, error = function(e) {
    data.frame(method = "Cox censor-at-switch", hr = NA, ci_low = NA,
               ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
  })

  # ── C. Cox time-dependent treatment (optional) ──
  if (run_cox_td) {
    results[["cox_td"]] <- tryCatch({
      fit <- fit_cox_td(dat_cox, tau = tau, covars = covars)
      data.frame(
        method = "Cox time-dependent", hr = fit$hr,
        ci_low = fit$ci_low, ci_high = fit$ci_high,
        risk_diff = NA_real_, risk_ratio = NA_real_, converged = TRUE
      )
    }, error = function(e) {
      data.frame(method = "Cox time-dependent", hr = NA, ci_low = NA,
                 ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
    })
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
  #   while_on_treatment: TIME-VARYING A_j on dat_cox (WOT-derived).
  #     Y/C columns reflect censoring at switch. static_binary holds constant.
  #     Distinct label: "LMTP SDR (WOT)".
  #
  #   composite: BASELINE-ONLY treatment (time_varying_trt=FALSE) on dat_cox.
  #     Switching IS part of the outcome; should not intervene on A_j.
  #
  #   principal_stratum: BASELINE-ONLY on subsets of dat_raw.
  #     Never-switchers have constant A_j anyway.

  # Helper to extract LMTP contrasts (handles old and new API)
  extract_lmtp <- function(res, label) {
    rd_obj <- res$contrast_rd
    rd_est <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
              else rd_obj$vals$theta
    rd_se  <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
              else rd_obj$vals$std.error
    rd_ci  <- rd_est + c(-1, 1) * 1.96 * rd_se
    rr_obj <- res$contrast_rr
    rr_est <- if (!is.null(rr_obj$estimates)) rr_obj$estimates$estimate
              else rr_obj$vals$theta
    data.frame(
      method = label, hr = NA_real_,
      ci_low = rd_ci[1], ci_high = rd_ci[2],
      risk_diff = rd_est, risk_ratio = rr_est, converged = TRUE
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
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners)
        extract_lmtp(res, "LMTP SDR")
      }, error = function(e) {
        data.frame(method = "LMTP SDR", hr = NA, ci_low = NA,
                   ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
      })

    } else if (estimand == "no_switch") {
      # ── Time-varying A_j on no-switch data: hold treatment constant ──
      results[["lmtp"]] <- tryCatch({
        n_events <- sum(dat_raw$event == 1 & dat_raw$follow_time <= tau)
        if (n_events < 5) stop("Too few events: ", n_events)
        prep <- prepare_lmtp_data(dat_raw, tau = tau, bin_width = bin_width,
                                  time_varying_trt = TRUE)
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners)
        extract_lmtp(res, "LMTP SDR (no-switch)")
      }, error = function(e) {
        data.frame(method = "LMTP SDR (no-switch)", hr = NA, ci_low = NA,
                   ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
      })

    } else if (estimand == "while_on_treatment") {
      # ── Time-varying A_j on dat_raw: hold treatment constant ──
      # WOT is a conditioning estimand, but LMTP implements an intervention.
      # Using dat_raw (not dat_cox) with time-varying A_j and static_binary
      # holds treatment constant = no-switch intervention. This approximates
      # WOT when the no-switch intervention is close to the natural course
      # for non-switchers. See Limitations for the conceptual distinction.
      results[["lmtp"]] <- tryCatch({
        n_events <- sum(dat_raw$event == 1 & dat_raw$follow_time <= tau)
        if (n_events < 5) stop("Too few events: ", n_events)
        prep <- prepare_lmtp_data(dat_raw, tau = tau, bin_width = bin_width,
                                  time_varying_trt = TRUE)
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners)
        extract_lmtp(res, "LMTP SDR (WOT)")
      }, error = function(e) {
        data.frame(method = "LMTP SDR (WOT)", hr = NA, ci_low = NA,
                   ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
      })

    } else if (estimand == "composite") {
      # ── Baseline-only treatment on composite-derived data ──
      results[["lmtp"]] <- tryCatch({
        n_events <- sum(dat_cox$event == 1 & dat_cox$follow_time <= tau)
        if (n_events < 5) stop("Too few events: ", n_events)
        prep <- prepare_lmtp_data(dat_cox, tau = tau, bin_width = bin_width,
                                  time_varying_trt = FALSE)
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners)
        extract_lmtp(res, "LMTP SDR (composite)")
      }, error = function(e) {
        data.frame(method = "LMTP SDR (composite)", hr = NA, ci_low = NA,
                   ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
      })

    } else if (estimand == "principal_stratum") {
      # ── Baseline-only on subsets: obs non-switchers + true PS ──
      results[["lmtp_approx"]] <- tryCatch({
        sub <- dat_raw[dat_raw$switched == 0, ]
        n_events <- sum(sub$event == 1 & sub$follow_time <= tau)
        if (n_events < 5) stop("Too few events: ", n_events)
        prep <- prepare_lmtp_data(sub, tau = tau, bin_width = bin_width,
                                  time_varying_trt = FALSE)
        res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners)
        extract_lmtp(res, "LMTP SDR (obs. non-switchers)")
      }, error = function(e) {
        data.frame(method = "LMTP SDR (obs. non-switchers)", hr = NA, ci_low = NA,
                   ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
      })

      if ("never_switcher" %in% names(dat_raw)) {
        results[["lmtp_true_ps"]] <- tryCatch(
          suppressWarnings({
            sub <- dat_raw[dat_raw$never_switcher == 1, ]
            n_events <- sum(sub$event == 1 & sub$follow_time <= tau)
            if (n_events < 5) stop("Too few events: ", n_events)
            prep <- prepare_lmtp_data(sub, tau = tau, bin_width = bin_width,
                                      time_varying_trt = FALSE)
            res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners)
            extract_lmtp(res, "LMTP SDR (true PS)")
          }),
          error = function(e) {
            data.frame(method = "LMTP SDR (true PS)", hr = NA, ci_low = NA,
                       ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
          }
        )
      }
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
                                 lmtp_learners = c("SL.mean", "SL.glm", "SL.bayesglm"),
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
  sim_results %>%
    group_by(estimand, method) %>%
    summarise(
      n_iter       = n(),
      n_converged  = sum(converged, na.rm = TRUE),
      # HR summaries (Cox methods)
      mean_hr      = mean(hr, na.rm = TRUE),
      sd_hr        = sd(hr, na.rm = TRUE),
      bias_hr      = if (!is.na(truth_hr)) mean(hr, na.rm = TRUE) - truth_hr
                     else NA_real_,
      # RD summaries (LMTP and KM-based)
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
      # RR summaries
      mean_rr      = mean(risk_ratio, na.rm = TRUE),
      bias_rr      = if (!is.na(truth_rr)) mean(risk_ratio, na.rm = TRUE) - truth_rr
                     else NA_real_,
      .groups = "drop"
    )
}
