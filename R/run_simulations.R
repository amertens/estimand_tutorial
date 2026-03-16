# run_simulations.R
# Repeated simulation study comparing Cox and LMTP estimators within estimands.
# Adapted from archive/ps_aki_analysis_par.R.

# TODO(Joy): fill in exact runtime notes after benchmarking.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(broom)
  library(here)
})

source(here("DGP.R"))
source(here("R", "helpers.R"))

# ── Single iteration ─────────────────────────────────────────────────────────
#' Run one iteration of the simulation study.
#' @param i integer; iteration index.
#' @param sample_size integer; per-iteration sample size.
#' @param tau integer; follow-up horizon.
#' @param run_lmtp logical; whether to run LMTP (slow).
#' @param dgp_args list of additional arguments passed to generate_hcv_data().
#' @return data.frame with one row of results.
run_one_iter <- function(i, sample_size = 1000, tau = 180,
                         run_lmtp = TRUE, dgp_args = list()) {
  set.seed(1000 + i)

  # Generate fresh data each iteration
  args <- modifyList(
    list(N = sample_size, seed = 1000 + i),
    dgp_args
  )
  dat <- do.call(generate_hcv_data, args)

  covars <- intersect(
    c("age", "ckd", "cirrhosis", "diabetes", "heart_failure"),
    names(dat)
  )

  # ── A. Naive Cox (ignores switching) ──
  cox_naive <- tryCatch({
    fit <- fit_cox_naive(dat, tau = tau, covars = covars)
    data.frame(
      method = "Cox naive (ignore switch)",
      hr = fit$hr, ci_low = fit$ci_low, ci_high = fit$ci_high,
      risk_diff = NA_real_, converged = TRUE
    )
  }, error = function(e) {
    data.frame(method = "Cox naive (ignore switch)",
               hr = NA, ci_low = NA, ci_high = NA,
               risk_diff = NA, converged = FALSE)
  })

  # ── B. Cox censor at switch ──
  cox_censor <- tryCatch({
    fit <- fit_cox_censor_switch(dat, tau = tau, covars = covars)
    data.frame(
      method = "Cox censor at switch",
      hr = fit$hr, ci_low = fit$ci_low, ci_high = fit$ci_high,
      risk_diff = NA_real_, converged = TRUE
    )
  }, error = function(e) {
    data.frame(method = "Cox censor at switch",
               hr = NA, ci_low = NA, ci_high = NA,
               risk_diff = NA, converged = FALSE)
  })

  # ── C. LMTP ──
  lmtp_row <- data.frame(
    method = "LMTP SDR", hr = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_,
    risk_diff = NA_real_, converged = FALSE
  )

  if (run_lmtp) {
    lmtp_row <- tryCatch({
      requireNamespace("lmtp", quietly = TRUE)
      prep <- prepare_lmtp_data(dat, tau = tau)
      res  <- run_lmtp_analysis(prep, folds = 2, learners = c("SL.glm"))
      rd_est <- res$contrast_rd$vals$theta
      rd_ci  <- rd_est + c(-1, 1) * 1.96 * res$contrast_rd$vals$std.error
      data.frame(
        method = "LMTP SDR",
        hr = NA_real_,
        ci_low = rd_ci[1], ci_high = rd_ci[2],
        risk_diff = rd_est, converged = TRUE
      )
    }, error = function(e) {
      data.frame(method = "LMTP SDR", hr = NA,
                 ci_low = NA, ci_high = NA,
                 risk_diff = NA, converged = FALSE)
    })
  }

  results <- bind_rows(cox_naive, cox_censor, lmtp_row)
  results$iter <- i
  results$n    <- sample_size
  results
}


# ── Full simulation runner ────────────────────────────────────────────────────
#' Run the repeated simulation study.
#' @param n_iter integer; number of iterations (default 50).
#' @param sample_size integer; sample size per iteration.
#' @param tau integer; follow-up horizon.
#' @param run_lmtp logical.
#' @param dgp_args list; extra args for generate_hcv_data().
#' @param cache_file character or NULL; path to cache results.
#' @return data.frame of aggregated results.
run_simulation_study <- function(n_iter = 50, sample_size = 1000,
                                 tau = 180, run_lmtp = TRUE,
                                 dgp_args = list(),
                                 cache_file = NULL) {
  # Check cache

if (!is.null(cache_file) && file.exists(cache_file)) {
    message("Loading cached results from ", cache_file)
    return(readRDS(cache_file))
  }

  message("Running ", n_iter, " simulation iterations...")
  all_res <- NULL
  for (i in seq_len(n_iter)) {
    if (i %% 10 == 0) message("  iteration ", i, " / ", n_iter)
    iter_res <- run_one_iter(i, sample_size = sample_size,
                             tau = tau, run_lmtp = run_lmtp,
                             dgp_args = dgp_args)
    all_res <- bind_rows(all_res, iter_res)
  }

  # Cache
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(all_res, cache_file)
    message("Results cached to ", cache_file)
  }

  all_res
}


# ── Summarize simulation results ─────────────────────────────────────────────
#' Summarize simulation results: bias, empirical SE, RMSE, coverage.
#' @param sim_results data.frame from run_simulation_study().
#' @param truth_rd numeric; true risk difference (for LMTP).
#' @param truth_hr numeric; true hazard ratio (for Cox, if known).
#' @return tibble with summary statistics.
summarize_simulation <- function(sim_results, truth_rd = NA, truth_hr = NA) {
  sim_results %>%
    group_by(method) %>%
    summarise(
      n_iter       = n(),
      n_converged  = sum(converged, na.rm = TRUE),
      # HR-based summaries (Cox methods)
      mean_hr      = mean(hr, na.rm = TRUE),
      bias_hr      = if (!is.na(truth_hr)) mean(hr, na.rm = TRUE) - truth_hr
                     else NA_real_,
      # RD-based summaries (LMTP)
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
      .groups = "drop"
    )
}
