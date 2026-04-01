# helpers.R
# Utility functions for the estimand-estimator alignment tutorial.
# Provides Cox regression variants, LMTP data preparation, and risk utilities.

# TODO(Joy): confirm exact file paths for the archived Cox and LMTP scripts used.

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

# ── Person-period expansion ──────────────────────────────────────────────────
#' Expand a subject-level dataset to person-day (long) format.
#' @param dat data.frame with columns: id, follow_time, event, treatment,
#'   switched, plus baseline covariates.
#' @param tau integer; maximum follow-up horizon in days.
#' @param covars character vector of baseline covariate names.
#' @return data.frame in long format with one row per person-day.
expand_person_period <- function(dat, tau = 180, covars = NULL) {
  if (is.null(covars)) {
    covars <- c("age", "sex_male", "ckd", "cirrhosis", "diabetes",
                "hypertension", "heart_failure", "nsaid", "acearb", "statin")
  }
  covars <- intersect(covars, names(dat))

  dat <- dat %>%
    mutate(
      obs_time  = pmin(follow_time, tau),
      obs_event = as.integer(event == 1 & follow_time <= tau),
      obs_cens  = as.integer(event == 0 & follow_time < tau)
    )

  long <- dat %>%
    select(id, treatment, switched, obs_time, obs_event, obs_cens,
           any_of("switch_time"), all_of(covars)) %>%
    tidyr::uncount(weights = ceiling(obs_time), .remove = FALSE) %>%
    group_by(id) %>%
    mutate(
      day       = row_number(),
      event_day = as.integer(day == max(day) & obs_event == 1),
      cens_day  = as.integer(day == max(day) & obs_cens == 1)
    ) %>%
    ungroup()

  long
}


# ── Naive Cox: baseline treatment only ───────────────────────────────────────
#' Fit a Cox model using only baseline treatment assignment (ignores switching).
#' This does not cleanly target either the treatment-policy or the hypothetical
#' estimand: the Cox HR is a conditional parameter and the risk set may be
#' altered by switching-induced censoring in the original data.
#' @param dat data.frame with follow_time, event, treatment, and baseline covariates.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector of adjustment covariates.
#' @return list with hr, ci, and survival curves.
fit_cox_naive <- function(dat, tau = 180,
                          covars = c("age", "ckd", "cirrhosis", "diabetes",
                                     "heart_failure")) {
  covars <- intersect(covars, names(dat))

  dat <- dat %>%
    mutate(time_use  = pmin(follow_time, tau),
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


# ── Cox: censor at switch ────────────────────────────────────────────────────
#' Fit a Cox model that censors follow-up at the observed treatment switch time.
#' Targets a while-on-treatment / hypothetical estimand but is biased when
#' switching is informative.
#' @param dat data.frame; must contain follow_time, event, switch_time,
#'   switched, treatment.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector.
#' @return list with hr, ci, survfit.
fit_cox_censor_switch <- function(dat, tau = 180,
                                  covars = c("age", "ckd", "cirrhosis",
                                             "diabetes", "heart_failure")) {
  covars <- intersect(covars, names(dat))

  # Censor exactly at the observed switch time.
  # If a subject switched before event AND before tau, treat as censored.
  dat <- dat %>%
    mutate(
      # effective time: event_time, admin censor, switch_time, or tau
      time_use = pmin(follow_time, tau),
      event_use = as.integer(event == 1 & follow_time <= tau)
    )

  # For subjects who switched before their event: censor at switch_time
  if ("switch_time" %in% names(dat)) {
    dat <- dat %>%
      mutate(
        sw_before_event = switched == 1 & switch_time < time_use,
        time_use  = ifelse(sw_before_event, pmin(switch_time, tau), time_use),
        event_use = ifelse(sw_before_event, 0L, event_use)
      ) %>%
      select(-sw_before_event)
  }

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


# ── Cox: time-dependent treatment covariate ─────────────────────────────────
#' Fit a Cox model with treatment as a time-dependent covariate that changes
#' at the switch time. This is sometimes used to "adjust" for switching, but
#' it does not target the treatment-policy estimand without strong assumptions
#' (no unmeasured time-varying confounding of the switch decision).
#' @param dat data.frame from generate_hep_data() with switch_time and switched.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector of baseline covariates.
#' @return list with hr, ci.
fit_cox_td <- function(dat, tau = 180,
                       covars = c("age", "ckd", "cirrhosis", "diabetes",
                                  "heart_failure")) {
  covars <- intersect(covars, names(dat))

  dat <- dat %>%
    mutate(
      time_use  = pmin(follow_time, tau),
      event_use = as.integer(event == 1 & follow_time <= tau)
    )

  # Build time-dependent dataset: split each person at their switch time
  # Before switch: trt_current = treatment (baseline)
  # After switch: trt_current = 1 - treatment (switched)
  rows <- list()
  for (i in seq_len(nrow(dat))) {
    r <- dat[i, ]
    sw_t <- if ("switch_time" %in% names(r)) r$switch_time else Inf
    if (r$switched == 1 && sw_t < r$time_use && sw_t > 0) {
      # interval 1: [0, switch_time) -- on original treatment
      rows[[length(rows) + 1]] <- data.frame(
        id = r$id, tstart = 0, tstop = sw_t,
        event_td = 0L, trt_current = r$treatment,
        r[, covars, drop = FALSE],
        stringsAsFactors = FALSE
      )
      # interval 2: [switch_time, time_use] -- on switched treatment
      rows[[length(rows) + 1]] <- data.frame(
        id = r$id, tstart = sw_t, tstop = r$time_use,
        event_td = r$event_use, trt_current = 1L - r$treatment,
        r[, covars, drop = FALSE],
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        id = r$id, tstart = 0, tstop = r$time_use,
        event_td = r$event_use, trt_current = r$treatment,
        r[, covars, drop = FALSE],
        stringsAsFactors = FALSE
      )
    }
  }
  td_dat <- bind_rows(rows)

  fml <- as.formula(paste0(
    "Surv(tstart, tstop, event_td) ~ trt_current + ",
    paste(covars, collapse = " + ")
  ))

  cox_fit <- coxph(fml, data = td_dat)
  hr_est  <- exp(coef(cox_fit)["trt_current"])
  hr_ci   <- exp(confint(cox_fit)["trt_current", ])

  list(
    model   = cox_fit,
    hr      = hr_est,
    ci_low  = hr_ci[1],
    ci_high = hr_ci[2]
  )
}


# ── LMTP intervention functions ─────────────────────────────────────────────
# These define the target interventions for lmtp_sdr().

#' Static intervention: set treatment to 1 (initiate active treatment).
#' Equivalent to lmtp::static_binary_on for a single time-point treatment.
static_binary_on <- function(data, trt) {
  rep(1L, nrow(data))
}

#' Static intervention: set treatment to 0 (initiate control).
static_binary_off <- function(data, trt) {
  rep(0L, nrow(data))
}

# TODO(Joy): finalise the no_switch intervention function for the
#   hypothetical estimand in the longitudinal (person-period) LMTP setup.
#   The function below is a placeholder for the case where treatment is
#   time-varying and the intervention sets post-baseline treatment to the
#   baseline value (i.e., no switching allowed).

#' No-switch intervention: set post-baseline treatment to baseline value.
#' For use with a time-varying treatment column in LMTP, this ensures
#' each subject remains on their initially assigned treatment throughout.
#' @param data data.frame with treatment column.
#' @param trt character; name of the treatment variable at the current time.
#' @return vector of treatment values equal to baseline treatment.
intervention_no_switch <- function(data, trt) {
  # In the current setup, treatment is time-fixed at baseline.
  # This function preserves that value across all time points.
  data[["treatment"]]
}


# ── LMTP data preparation ───────────────────────────────────────────────────
#' Prepare wide-format data for lmtp_sdr().
#' @param dat data.frame from generate_hep_data().
#' @param tau integer; follow-up horizon.
#' @param baseline character vector of baseline covariate names.
#' @return list with wide data.frame, Y_cols, C_cols, baseline.
prepare_lmtp_data <- function(dat, tau = 180,
                              baseline = c("age", "sex_male", "ckd",
                                           "diabetes", "hypertension",
                                           "heart_failure")) {
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


# ── LMTP estimation wrapper ─────────────────────────────────────────────────
#' Run lmtp_sdr for two interventions and return contrasts.
#' @param lmtp_prep output of prepare_lmtp_data().
#' @param shift_on function; intervention for "treated" arm.
#'   Defaults to lmtp::static_binary_on.
#' @param shift_off function; intervention for "control" arm.
#'   Defaults to lmtp::static_binary_off.
#' @param folds integer; number of cross-validation folds.
#'   For production: V >= 10 when n_eff < 5000 (Gruber et al. 2022).
#' @param learners character vector of SuperLearner libraries.
#'   For production: c("SL.glm", "SL.glmnet", "SL.xgboost") at minimum.
#'   Tutorial default is SL.glm only for speed.
#' @return list with res_on, res_off, risk_trt, risk_ctrl, contrast_rr, contrast_rd.
run_lmtp_analysis <- function(lmtp_prep,
                              shift_on  = NULL,
                              shift_off = NULL,
                              folds = 2,
                              learners = c("SL.glm")) {
  requireNamespace("lmtp", quietly = TRUE)

  if (is.null(shift_on))  shift_on  <- lmtp::static_binary_on
  if (is.null(shift_off)) shift_off <- lmtp::static_binary_off

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
                     c(common_args, list(shift = shift_on)))
  res_off <- do.call(lmtp::lmtp_sdr,
                     c(common_args, list(shift = shift_off)))

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


# ── Risk extraction utilities ────────────────────────────────────────────────

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

#' Extract hazard ratio and CI from a coxph fit.
#' @param cox_fit coxph object.
#' @param var character; name of the treatment variable.
#' @return named list: hr, ci_low, ci_high.
extract_hr <- function(cox_fit, var = "treatment") {
  hr  <- exp(coef(cox_fit)[var])
  ci  <- exp(confint(cox_fit)[var, ])
  list(hr = hr, ci_low = ci[1], ci_high = ci[2])
}

#' Compute risk difference from two scalar risks.
#' @param risk_trt numeric; risk under treatment.
#' @param risk_ctrl numeric; risk under control.
#' @return list with rd and rr.
risk_contrast <- function(risk_trt, risk_ctrl) {
  list(
    rd = risk_trt - risk_ctrl,
    rr = risk_trt / risk_ctrl
  )
}
