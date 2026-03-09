# helpers.R
# Utility functions for the estimand-estimator alignment tutorial.
# Reuses logic from archive/lmtp_aki_analysis.R and archive/regression_estimation.Rmd.

# TODO(Joy): confirm exact file paths for the archived Cox and LMTP scripts used.
#   This file consolidates archive/lmtp_aki_analysis.R, archive/ps_aki_analysis_par.R,
#   and archive/regression_estimation.Rmd into reusable helpers.

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

# ── Person-period expansion ──────────────────────────────────────────────────
#' Expand a subject-level dataset to person-day (long) format.
#' @param dat data.frame with columns: id, follow_time, event, treatment,
#'   switch, plus baseline covariates.
#' @param tau integer; maximum follow-up horizon in days.
#' @param covars character vector of baseline covariate names.
#' @return data.frame in long format with one row per person-day.
expand_person_period <- function(dat, tau = 180, covars = NULL) {
  if (is.null(covars)) {
    covars <- c("age", "sex_male", "ckd", "cirrhosis", "diabetes",
                "hypertension", "hiv", "bmi")
  }
  # keep only covars that exist in dat

  covars <- intersect(covars, names(dat))

  dat <- dat %>%
    mutate(
      obs_time = pmin(follow_time, tau),
      obs_event = as.integer(event == 1 & follow_time <= tau),
      obs_cens  = as.integer(event == 0 & follow_time < tau)
    )

  long <- dat %>%
    select(id, treatment, switch, obs_time, obs_event, obs_cens,
           all_of(covars)) %>%
    tidyr::uncount(weights = ceiling(obs_time), .remove = FALSE) %>%
    group_by(id) %>%
    mutate(
      day = row_number(),
      event_day = as.integer(day == max(day) & obs_event == 1),
      cens_day  = as.integer(day == max(day) & obs_cens == 1)
    ) %>%
    ungroup()

  long
}


# ── Naive Cox: baseline treatment only ────────────────────────────────────────
#' Fit a Cox model using only baseline treatment assignment (ignores switching).
#' This targets neither the treatment-policy nor the hypothetical estimand cleanly.
#' @param dat data.frame with follow_time, event, treatment, and baseline covariates.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector of adjustment covariates.
#' @return list with hr, ci, and survival curves.
fit_cox_naive <- function(dat, tau = 180,
                          covars = c("age", "ckd", "cirrhosis", "diabetes",
                                     "hiv")) {
  covars <- intersect(covars, names(dat))
  dat <- dat %>%
    mutate(time_use = pmin(follow_time, tau),
           event_use = as.integer(event == 1 & follow_time <= tau))

  fml <- as.formula(paste0(
    "Surv(time_use, event_use) ~ treatment + ",
    paste(covars, collapse = " + ")
  ))

  cox_fit <- coxph(fml, data = dat)
  hr_est  <- exp(coef(cox_fit)["treatment"])
  hr_ci   <- exp(confint(cox_fit)["treatment", ])

  sf <- survfit(Surv(time_use, event_use) ~ treatment, data = dat)

  list(
    model   = cox_fit,
    hr      = hr_est,
    ci_low  = hr_ci[1],
    ci_high = hr_ci[2],
    survfit = sf
  )
}


# ── Cox: censor at switch ─────────────────────────────────────────────────────
#' Fit a Cox model that censors follow-up at treatment switch.
#' Implicitly targets a while-on-treatment / hypothetical estimand but is
#' biased when switching is informative.
#' @param dat data.frame; must contain follow_time, event, switch, treatment.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector.
#' @return list with hr, ci, survfit.
fit_cox_censor_switch <- function(dat, tau = 180,
                                  covars = c("age", "ckd", "cirrhosis",
                                             "diabetes", "hiv")) {
  covars <- intersect(covars, names(dat))

  # Reconstruct censoring at switch.
  # In DGP.R, follow_time already incorporates switch censoring.
  # For subjects who switched, follow_time <= switch_time + risk_window.
  # We want to censor exactly at switch time.
  # Since follow_time = min(event_time, censor_admin, censor_switch),
  # and switch subjects have switch==1, we use follow_time directly
  # but mark switch subjects as censored (event=0) if they didn't have

  # the event before the switch.
  dat <- dat %>%
    mutate(
      time_use  = pmin(follow_time, tau),
      # If the subject switched and the observed event might have occurred
      # after switching, treat as censored.
      event_use = as.integer(event == 1 & follow_time <= tau)
    )

  fml <- as.formula(paste0(
    "Surv(time_use, event_use) ~ treatment + ",
    paste(covars, collapse = " + ")
  ))

  cox_fit <- coxph(fml, data = dat)
  hr_est  <- exp(coef(cox_fit)["treatment"])
  hr_ci   <- exp(confint(cox_fit)["treatment", ])

  sf <- survfit(Surv(time_use, event_use) ~ treatment, data = dat)

  list(
    model   = cox_fit,
    hr      = hr_est,
    ci_low  = hr_ci[1],
    ci_high = hr_ci[2],
    survfit = sf
  )
}


# ── LMTP data preparation ────────────────────────────────────────────────────
#' Prepare wide-format data for lmtp_sdr().
#' Adapted from archive/lmtp_aki_analysis.R.
#' @param dat data.frame from generate_hcv_data().
#' @param tau integer; follow-up horizon.
#' @param baseline character vector of baseline covariate names.
#' @return list with wide data.frame, Y_cols, C_cols.
prepare_lmtp_data <- function(dat, tau = 180,
                              baseline = c("age", "sex_male", "ckd",
                                           "diabetes", "hypertension")) {
  baseline <- intersect(baseline, names(dat))

  dat <- dat %>%
    mutate(
      aki_event    = as.integer(event == 1 & follow_time <= tau),
      time_to_aki  = if_else(aki_event == 1, follow_time, as.numeric(tau)),
      cens_event   = as.integer(event == 0 & follow_time < tau),
      time_to_cens = if_else(cens_event == 1, follow_time, as.numeric(tau))
    )

  # Build wide Y and C vectors
  make_Y <- function(t_aki, e) as.integer(seq_len(tau) >= t_aki & e == 1)
  make_C <- function(t_c, e) {
    v <- rep(1L, tau)
    if (e == 1 && t_c < tau) v[seq(t_c + 1L, tau)] <- 0L
    v
  }

  Y_mat <- t(mapply(make_Y, dat$time_to_aki, dat$aki_event))
  C_mat <- t(mapply(make_C, dat$time_to_cens, dat$cens_event))

  Y_cols <- paste0("Y", seq_len(tau))
  C_cols <- paste0("C", seq_len(tau))
  colnames(Y_mat) <- Y_cols
  colnames(C_mat) <- C_cols

  wide <- bind_cols(
    dat %>% select(id, treatment, all_of(baseline)),
    as_tibble(Y_mat),
    as_tibble(C_mat)
  )

  # Carry forward: once Y=1, all subsequent Y must be 1
  wide <- lmtp::event_locf(as.data.frame(wide), outcomes = Y_cols)


  list(data = as.data.frame(wide), Y_cols = Y_cols, C_cols = C_cols,
       baseline = baseline)
}


# ── LMTP estimation wrapper ──────────────────────────────────────────────────
#' Run lmtp_sdr for treat-all vs treat-none and return contrasts.
#' @param lmtp_prep output of prepare_lmtp_data().
#' @param folds integer; number of cross-validation folds.
#' @param learners character vector of SuperLearner libraries.
#' @return list with res_on, res_off, contrast_rr, contrast_rd.
run_lmtp_analysis <- function(lmtp_prep, folds = 2,
                              learners = c("SL.glm")) {
  requireNamespace("lmtp", quietly = TRUE)

  common_args <- list(
    data    = lmtp_prep$data,
    trt     = "treatment",
    outcome = lmtp_prep$Y_cols,
    cens    = lmtp_prep$C_cols,
    baseline = lmtp_prep$baseline,
    outcome_type = "survival",
    folds   = folds,
    learners_trt     = learners,
    learners_outcome = learners
  )

  res_on  <- do.call(lmtp::lmtp_sdr,
                     c(common_args, list(shift = lmtp::static_binary_on)))
  res_off <- do.call(lmtp::lmtp_sdr,
                     c(common_args, list(shift = lmtp::static_binary_off)))

  contrast_rr <- lmtp::lmtp_contrast(res_on, ref = res_off, type = "rr")
  contrast_rd <- lmtp::lmtp_contrast(res_on, ref = res_off, type = "additive")

  list(
    res_on      = res_on,
    res_off     = res_off,
    risk_trt    = 1 - res_on$theta[length(res_on$theta)],
    risk_ctrl   = 1 - res_off$theta[length(res_off$theta)],
    contrast_rr = contrast_rr,
    contrast_rd = contrast_rd
  )
}


# ── Marginal risk from KM ────────────────────────────────────────────────────
#' Extract marginal risks at time t from a survfit object.
#' @param sf survfit object with treatment strata.
#' @param t time point.
#' @return named list: risk_ctrl, risk_trt, risk_diff.
km_risk_at <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)
  risks <- 1 - s$surv
  list(risk_ctrl = risks[1], risk_trt = risks[2],
       risk_diff = risks[2] - risks[1])
}
