# run_simulations.R
# Repeated simulation study comparing estimators WITHIN each estimand.
# For each estimand (treatment-policy and hypothetical no-switch), we fit:
#   1. Naive Cox HR (baseline treatment only) -- misaligned
#   2. Cox censoring at switch -- misaligned for treatment-policy; naive for hypothetical
#   3. Cox with time-dependent treatment -- misaligned
#   4. LMTP SDR targeting the correct intervention -- aligned
# We also run three support scenarios (good, strained, poor) to show how
# estimation degrades with weak positivity.

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

# Number of parallel workers (set to number of cores)
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

  # Generate data under the appropriate policy
  args <- modifyList(
    list(N = sample_size, policy = estimand, seed = 1000 + i),
    dgp_args
  )
  dat <- do.call(generate_hep_data, args)

  covars <- intersect(
    c("age", "ckd", "cirrhosis", "diabetes", "heart_failure"),
    names(dat)
  )

  results <- list()

  # ── A. Naive Cox (ignores switching) ──
  results[["cox_naive"]] <- tryCatch({
    fit <- fit_cox_naive(dat, tau = tau, covars = covars)
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
    fit <- fit_cox_censor_switch(dat, tau = tau, covars = covars)
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
      fit <- fit_cox_td(dat, tau = tau, covars = covars)
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
  # Estimand-specific LMTP handling:
  #   treatment_policy:   static on/off on full data (default)
  #   no_switch:          static on/off on data generated without switching
  #   while_on_treatment: static on/off; DGP censors at switch, LMTP
  #                       censoring model accounts for this
  #   composite:          static on/off; outcome = AKI or switch (already
  #                       encoded in DGP event column under composite policy)
  #   principal_stratum:  subset to observed non-switchers, then static on/off.
  #                       This is an approximation -- the true principal stratum
  #                       requires both potential switching indicators.
  if (run_lmtp) {
    results[["lmtp"]] <- tryCatch({
      requireNamespace("lmtp", quietly = TRUE)

      # For principal stratum: restrict to observed non-switchers
      lmtp_dat <- if (estimand == "principal_stratum") {
        dat[dat$switched == 0, ]
      } else {
        dat
      }

      # Check we have enough data and events
      n_events <- sum(lmtp_dat$event == 1 & lmtp_dat$follow_time <= tau)
      if (n_events < 5) stop("Too few events for LMTP: ", n_events)

      prep <- prepare_lmtp_data(lmtp_dat, tau = tau, bin_width = bin_width)
      res  <- run_lmtp_analysis(prep, folds = 2, learners = lmtp_learners)

      # Handle both old ($vals$theta) and new ($estimates$estimate) lmtp API
      rd_obj <- res$contrast_rd
      rd_est <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$estimate
                else rd_obj$vals$theta
      rd_se  <- if (!is.null(rd_obj$estimates)) rd_obj$estimates$std.error
                else rd_obj$vals$std.error
      rd_ci  <- rd_est + c(-1, 1) * 1.96 * rd_se
      rr_obj <- res$contrast_rr
      rr_est <- if (!is.null(rr_obj$estimates)) rr_obj$estimates$estimate
                else rr_obj$vals$theta

      method_label <- if (estimand == "principal_stratum") {
        "LMTP SDR (non-switchers)"
      } else {
        "LMTP SDR"
      }

      data.frame(
        method = method_label, hr = NA_real_,
        ci_low = rd_ci[1], ci_high = rd_ci[2],
        risk_diff = rd_est, risk_ratio = rr_est, converged = TRUE
      )
    }, error = function(e) {
      data.frame(method = "LMTP SDR", hr = NA, ci_low = NA,
                 ci_high = NA, risk_diff = NA, risk_ratio = NA, converged = FALSE)
    })
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
      # ── Parallel execution via PSOCK cluster (Windows-compatible) ──
      cl <- makeCluster(n_cores)
      on.exit(stopCluster(cl), add = TRUE)

      # Export required functions and packages to workers
      clusterEvalQ(cl, {
        suppressPackageStartupMessages({
          library(dplyr)
          library(survival)
          library(here)
        })
        source(here("DGP.R"))
        source(here("R", "helpers.R"))
      })

      # Export the iteration function and parameters
      clusterExport(cl, c("run_one_iter", "fit_cox_naive",
                          "fit_cox_censor_switch", "fit_cox_td",
                          "prepare_lmtp_data", "run_lmtp_analysis",
                          "km_risk_at", "generate_hep_data"),
                    envir = environment())

      iter_results <- parLapply(cl, seq_len(n_iter), function(i) {
        run_one_iter(
          i, estimand = est, sample_size = sample_size,
          tau = tau, bin_width = bin_width,
          lmtp_learners = lmtp_learners,
          run_lmtp = run_lmtp, run_cox_td = run_cox_td,
          dgp_args = dgp_args
        )
      })

      stopCluster(cl)
      on.exit(NULL)  # clear the on.exit since we stopped manually

      est_res <- bind_rows(iter_results)
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
