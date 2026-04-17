# helpers.R
# Utility functions for the estimand-estimator simulation study.
# Provides Cox regression variants, LMTP data preparation, and risk utilities.
#
# LMTP specification per estimand:
#   treatment_policy:   time_varying_trt=FALSE, censor_at_switch=FALSE.
#                       Baseline treatment only; censoring model handles admin censoring.
#   no_switch:          time_varying_trt=TRUE, censor_at_switch=FALSE.
#                       A_j columns with static_binary = hold treatment constant.
#   while_on_treatment: time_varying_trt=FALSE, censor_at_switch=TRUE.
#                       Baseline treatment only; switching enters through C columns.
#   composite:          time_varying_trt=FALSE, censor_at_switch=FALSE.
#                       Composite outcome in Y columns.
#   principal_stratum:  time_varying_trt=FALSE, censor_at_switch=FALSE.
#                       Baseline treatment on oracle never-switcher subset.


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
    dplyr::select(id, treatment, switched, obs_time, obs_event, obs_cens,
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
#' Estimates a conditional hazard ratio, not a marginal risk contrast. Often
#' interpreted as a treatment-policy analysis, but the Cox HR does not
#' correspond to the marginal risk difference defined by that estimand.
#'
#' Intended estimand: treatment-policy (Cox version). Targets the conditional
#' HR for baseline treatment assignment over a [0, tau] follow-up window that
#' ignores any post-baseline switching.
#'
#' @param dat data.frame with follow_time, event, treatment, and baseline covariates.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector of adjustment covariates.
#' @return list with:
#'   - `model`: fitted coxph object
#'   - `hr`, `ci_low`, `ci_high`: HR point estimate and 95% HR CI
#'     (from the Cox model; on the HR scale, not the RD scale)
#'   - `survfit`: stratified survfit for KM risk extraction via km_risk_at()
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
#' Often used to approximate a hypothetical no-switch or while-on-treatment
#' analysis. Estimates a conditional HR among non-censored person-time; biased
#' when switching is informative (i.e., when the censoring-at-random assumption
#' is violated because CKD or other confounders drive both switching and the
#' outcome).
#'
#' Intended estimand: WOT or no-switch (Cox version). Both targets are
#' approximated by this model under the strong assumption of non-informative
#' censoring; the resulting HR is a conditional within-model summary and
#' does not target the marginal RD of any ICH E9(R1) estimand.
#'
#' @param dat data.frame; must contain follow_time, event, switch_time,
#'   switched, treatment.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector of adjustment covariates.
#' @return list with:
#'   - `model`: fitted coxph object
#'   - `hr`, `ci_low`, `ci_high`: HR point estimate and 95% HR CI
#'   - `survfit`: stratified survfit under the censor-at-switch follow-up
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
      dplyr::select(-sw_before_event)
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
#' at the switch time. This is sometimes used to "adjust" for switching by
#' attributing post-switch person-time to the switched treatment, but it does
#' not target any ICH E9(R1) marginal estimand without strong assumptions
#' (no unmeasured time-varying confounding of the switch decision, and a
#' proportional-hazards assumption on the current-treatment indicator).
#'
#' Intended estimand: none of the ICH E9(R1) strategies map cleanly to this
#' model. Reported for completeness as a commonly-used switching adjustment.
#'
#' @param dat data.frame from generate_hep_data() with switch_time and switched.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector of baseline covariates.
#' @return list with:
#'   - `model`: fitted coxph object with (tstart, tstop) intervals
#'   - `hr`, `ci_low`, `ci_high`: HR for the current-treatment indicator
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


# ── Cox: IPCW-weighted (censor at switch, adjust for informative censoring) ─
#' Fit an IPCW-weighted Cox model that censors at switch but reweights
#' non-switchers to account for informative censoring. A commonly-used
#' pharmacoepi approach for WOT / hypothetical no-switch estimands.
#'
#' The switching model is a subject-level logistic regression of
#' `did_switch ~ treatment + weight_covars`. Stabilised IPCW weights are
#' computed as `(1 - p_marginal) / (1 - p_conditional)` for non-switchers
#' and truncated at the 99th percentile; switchers receive weight 1 since
#' they are censored (their weight is irrelevant to the outcome model).
#' If no subject switches, IPCW reduces to naive Cox and the function
#' returns `fit_cox_naive()` on the input data.
#'
#' Simplification: this implementation uses a subject-level binary switching
#' model rather than a time-varying pooled logistic model. A production
#' IPCW analysis would specify the time-varying censoring process with
#' hazard-based weights; see Limitations in the manuscript.
#'
#' Intended estimand: WOT / hypothetical no-switch (Cox version). Targets
#' the conditional HR under follow-up censored at switch, reweighted to
#' recover the marginal among covariate-matched non-switchers.
#'
#' @param dat data.frame from generate_hep_data() with switch_time, switched.
#' @param tau numeric; follow-up horizon.
#' @param covars character vector of baseline covariates for the outcome model.
#' @param weight_covars character vector of covariates for the switching model.
#'   Defaults to covars (treatment is always added to the switching model).
#' @return list with:
#'   - `model`: fitted coxph object with robust SE
#'   - `hr`, `ci_low`, `ci_high`: HR point estimate and 95% HR CI
#'     (CI uses the robust sandwich variance to account for weighting)
#'   - `survfit`: weighted stratified survfit under censor-at-switch follow-up
fit_cox_ipcw <- function(dat, tau = 180,
                         covars = c("age", "ckd", "cirrhosis", "diabetes",
                                    "heart_failure"),
                         weight_covars = NULL) {
  covars <- intersect(covars, names(dat))
  if (is.null(weight_covars)) weight_covars <- covars

  dat <- dat %>%
    mutate(
      time_use  = pmin(follow_time, tau),
      event_use = as.integer(event == 1 & follow_time <= tau),
      sw_time   = if ("switch_time" %in% names(.)) switch_time else Inf,
      did_switch = as.integer(switched == 1 & sw_time < time_use)
    )

  # If nobody switches, IPCW reduces to naive Cox (no reweighting needed)
  if (sum(dat$did_switch) == 0) {
    return(fit_cox_naive(dat, tau = tau, covars = covars))
  }

  # Censor at switch (same as fit_cox_censor_switch)
  dat <- dat %>%
    mutate(
      time_ipcw  = ifelse(did_switch == 1, pmin(sw_time, tau), time_use),
      event_ipcw = ifelse(did_switch == 1, 0L, event_use)
    )

  # Fit switching model: P(switch by min(time_use, tau) | treatment, W)
  sw_fml <- as.formula(paste0(
    "did_switch ~ treatment + ", paste(weight_covars, collapse = " + ")
  ))
  sw_fit <- glm(sw_fml, data = dat, family = binomial())
  p_switch <- predict(sw_fit, type = "response")

  # Stabilised IPCW weights:
  #   For non-switchers: w = P(no switch | marginal) / P(no switch | A, W)
  #   For switchers: censored, so their weight doesn't matter (event_ipcw = 0)
  p_switch_marginal <- mean(dat$did_switch)
  w_unstab <- 1 / pmax(1 - p_switch, 0.01)
  w_stab   <- (1 - p_switch_marginal) * w_unstab
  # Switchers get weight 1 (they're censored, weight is irrelevant)
  dat$ipcw <- ifelse(dat$did_switch == 1, 1, w_stab)
  # Truncate extreme weights at 99th percentile
  w99 <- quantile(dat$ipcw, 0.99)
  dat$ipcw <- pmin(dat$ipcw, w99)

  fml <- as.formula(paste0(
    "Surv(time_ipcw, event_ipcw) ~ treatment + ",
    paste(covars, collapse = " + ")
  ))

  cox_fit <- coxph(fml, data = dat, weights = ipcw, robust = TRUE)
  hr_est  <- exp(coef(cox_fit)["treatment"])
  # Use robust SE for CI (vcov extracts the correct variance matrix)
  robust_se <- sqrt(vcov(cox_fit)["treatment", "treatment"])
  hr_ci <- exp(coef(cox_fit)["treatment"] + c(-1, 1) * 1.96 * robust_se)

  sf <- survfit(Surv(time_ipcw, event_ipcw) ~ treatment, data = dat,
                weights = dat$ipcw)

  list(
    model   = cox_fit,
    hr      = hr_est,
    ci_low  = hr_ci[1],
    ci_high = hr_ci[2],
    survfit = sf
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

#' Build a shifted dataset that intervenes on baseline treatment only.
#' Used by the optional `add_baseline_only` branch of `run_lmtp_analysis()`
#' (not invoked by the main simulation pipeline). Sets the first A column to
#' a fixed value and leaves later A columns at their observed values, with
#' all C columns set to 1 (required by `lmtp_sdr()` when `shifted=` is used
#' in a survival setting).
#'
#' @param data data.frame (wide format, with A and C columns).
#' @param baseline_trt character; name of the baseline treatment column
#'   (e.g., "A1" or "treatment").
#' @param value numeric 0/1; treatment value to assign at baseline.
#' @param C_cols character vector of censoring column names (e.g., C1..CK).
#' @return data.frame with baseline_trt set to `value` and C_cols set to 1.
make_shifted_baseline <- function(data, baseline_trt, value, C_cols = NULL) {
  shifted <- data

  # Intervene only on baseline treatment
  shifted[[baseline_trt]] <- value

  # Required for lmtp when using shifted= in survival settings
  if (!is.null(C_cols)) {
    for (cc in C_cols) {
      shifted[[cc]] <- 1
    }
  }

  shifted
}


#' No-switch intervention: preserve baseline treatment at every time point.
#' Not currently called by the main simulation pipeline — the no-switch
#' estimand is implemented via `static_binary_on/off` on time-varying A_j
#' columns (which override all A_j to a constant, equivalent to holding
#' baseline treatment). Retained for reference and possible future use.
#'
#' @param data data.frame with treatment column.
#' @param trt character; name of the treatment variable at the current time.
#' @return vector of treatment values equal to baseline treatment.
intervention_no_switch <- function(data, trt) {
  data[["treatment"]]
}


# ── LMTP data preparation ───────────────────────────────────────────────────
#' Prepare wide-format data for lmtp_sdr().
#'
#' Treatment specification:
#'   time_varying_trt = TRUE: A_j columns reflecting observed switching
#'     trajectory. Used for no-switch (static_binary holds constant).
#'   time_varying_trt = FALSE: baseline "treatment" only. Used for
#'     treatment-policy, WOT, composite, principal stratum.
#'
#' Censoring specification:
#'   censor_at_switch = FALSE (default): C columns encode only
#'     administrative censoring (event == 0 & follow_time < tau).
#'   censor_at_switch = TRUE: C columns encode both administrative
#'     censoring AND switching-induced censoring. Switching is treated
#'     as a censoring event: after switch_time, C_j = 0. Used for the
#'     while-on-treatment (WOT) estimand.
#'
#' @param dat data.frame from generate_hep_data().
#' @param tau integer; follow-up horizon in days.
#' @param bin_width integer; days per time bin.
#' @param time_varying_trt logical; create A_j columns if TRUE.
#' @param censor_at_switch logical; include switching in C columns if TRUE.
#' @param trt_var character; baseline treatment column name.
#' @param switch_var character; switch time column name.
#' @param baseline character vector of baseline covariate names.
#' @return list with data, Y_cols, C_cols, A_cols (or NULL), baseline, n_bins.
prepare_lmtp_data <- function(dat, tau = 180, bin_width = 1,
                              time_varying_trt = TRUE,
                              censor_at_switch = FALSE,
                              trt_var = "treatment",
                              switch_var = "switch_time",
                              baseline = c("age", "sex_male", "ckd",
                                           "cirrhosis", "nsaid",
                                           "diabetes", "hypertension",
                                           "heart_failure")) {
  baseline <- intersect(baseline, names(dat))

  bin_edges <- seq(0, tau, by = bin_width)
  if (tail(bin_edges, 1) < tau) bin_edges <- c(bin_edges, tau)
  n_bins <- length(bin_edges) - 1

  # When censor_at_switch = TRUE, redefine follow-up to end at switch
  # so that Y and C columns reflect WOT censoring.
  if (censor_at_switch && switch_var %in% names(dat)) {
    dat <- dat %>%
      mutate(
        wot_follow = pmin(event_time, .data[[switch_var]], follow_time, tau),
        wot_event  = as.integer(event_time <= wot_follow &
                                  event_time <= .data[[switch_var]]),
        aki_event    = wot_event,
        time_to_aki  = if_else(aki_event == 1, wot_follow, as.numeric(tau)),
        # Censoring = not having the event before tau (includes both
        # admin censoring AND switch-induced censoring)
        cens_event   = as.integer(aki_event == 0 & wot_follow < tau),
        time_to_cens = if_else(cens_event == 1, wot_follow, as.numeric(tau))
      )
  } else {
    dat <- dat %>%
      mutate(
        aki_event    = as.integer(event == 1 & follow_time <= tau),
        time_to_aki  = if_else(aki_event == 1, follow_time, as.numeric(tau)),
        cens_event   = as.integer(event == 0 & follow_time < tau),
        time_to_cens = if_else(cens_event == 1, follow_time, as.numeric(tau))
      )
  }

  make_Y <- function(t_aki, e) {
    as.integer(bin_edges[-1] >= t_aki & e == 1)
  }
  make_C <- function(t_c, e) {
    v <- rep(1L, n_bins)
    if (e == 1 && t_c < tau) {
      cens_bin <- which(bin_edges[-1] > t_c)[1]
      if (!is.na(cens_bin) && cens_bin <= n_bins) v[seq(cens_bin, n_bins)] <- 0L
    }
    v
  }

  Y_mat <- t(mapply(make_Y, dat$time_to_aki, dat$aki_event))
  C_mat <- t(mapply(make_C, dat$time_to_cens, dat$cens_event))

  Y_cols <- paste0("Y", seq_len(n_bins))
  C_cols <- paste0("C", seq_len(n_bins))
  colnames(Y_mat) <- Y_cols
  colnames(C_mat) <- C_cols

  if (time_varying_trt) {
    # Time-varying A_j: reflects observed treatment at each bin.
    # Before switch_time: A_j = baseline treatment.
    # After switch_time: A_j = 1 - baseline treatment.
    # Used for no-switch/WOT: static_binary intervention overrides all A_j
    # to hold treatment constant = implements "no switching" counterfactual.
    make_A <- function(a0, sw_time) {
      A <- rep(as.integer(a0), n_bins)
      if (is.finite(sw_time) && sw_time < tau) {
        sw_bin <- which(bin_edges[-1] > sw_time)[1]
        if (!is.na(sw_bin) && sw_bin <= n_bins)
          A[seq(sw_bin, n_bins)] <- 1L - as.integer(a0)
      }
      A
    }
    A_mat <- t(mapply(make_A, dat[[trt_var]], dat[[switch_var]]))
    A_cols <- paste0("A", seq_len(n_bins))
    colnames(A_mat) <- A_cols

    wide <- bind_cols(
      dat %>% dplyr::select(id, treatment, all_of(baseline)),
      as_tibble(A_mat),
      as_tibble(Y_mat),
      as_tibble(C_mat)
    )

    # When A_cols are used as trt, include baseline treatment in adjustment set
    baseline_out <- if (!"treatment" %in% baseline) c("treatment", baseline)
                    else baseline
  } else {
    # Baseline-only treatment: single "treatment" column.
    # Used for treatment-policy (let switching happen naturally; LMTP
    # censoring model handles informative follow-up loss) and composite
    # (switching is part of the outcome, not a treatment node).
    A_cols <- NULL

    wide <- bind_cols(
      dat %>% dplyr::select(id, treatment, all_of(baseline)),
      as_tibble(Y_mat),
      as_tibble(C_mat)
    )
    baseline_out <- baseline
  }

  wide <- lmtp::event_locf(as.data.frame(wide), outcomes = Y_cols)

  list(data = as.data.frame(wide), A_cols = A_cols, Y_cols = Y_cols,
       C_cols = C_cols, baseline = baseline_out, n_bins = n_bins,
       bin_width = bin_width)
}


# ── LMTP estimation wrapper ─────────────────────────────────────────────────
#' Run `lmtp_sdr()` for two static interventions and return both interventions'
#' fits plus the risk-ratio and risk-difference contrasts.
#'
#' The intervention dispatch is driven by what `prepare_lmtp_data()` produced:
#'   - If `lmtp_prep$A_cols` is non-NULL (time-varying A_j), `static_binary_on/off`
#'     overrides every A_j to 1 or 0 — implementing the no-switch counterfactual
#'     "hold treatment constant at every time point."
#'   - If `lmtp_prep$A_cols` is NULL (baseline-only), `static_binary_on/off`
#'     sets the single baseline treatment column to 1 or 0 — implementing the
#'     treatment-policy / WOT / composite / principal-stratum targets.
#'
#' Return-value sign convention:
#'   - `risk_trt`, `risk_ctrl` are on the risk scale: `1 - theta[K]`.
#'   - `contrast_rd` returns the survival-scale contrast `S(1) - S(0)` directly
#'     from `lmtp_contrast(type = "additive")`. Callers that want the
#'     risk-scale RD `R(1) - R(0)` must negate the point estimate and swap the
#'     CI bounds (done in `extract_lmtp()` in R/run_simulations.R).
#'
#' @param lmtp_prep output of `prepare_lmtp_data()`.
#' @param shift_on function; intervention for the "treated" arm.
#'   Defaults to `lmtp::static_binary_on`.
#' @param shift_off function; intervention for the "control" arm.
#'   Defaults to `lmtp::static_binary_off`.
#' @param add_baseline_only logical; if TRUE and `lmtp_prep$A_cols` is non-NULL,
#'   also runs a baseline-only intervention using `shifted=` (treats the first
#'   A column as the intervention point and leaves later A columns at observed
#'   values). Optional diagnostic; not used by the main simulation pipeline.
#' @param folds integer; number of cross-validation folds. For production:
#'   V >= 10 when n_eff < 5000 (Gruber et al. 2022). Simulation default is 2.
#' @param learners character vector of SuperLearner libraries. For production:
#'   `c("SL.glm", "SL.glmnet", "SL.xgboost")` at minimum. Tutorial default
#'   is `c("SL.mean", "SL.glm", "SL.bayesglm")` for speed.
#' @return list with:
#'   - `res_on`, `res_off`: full `lmtp_sdr` fits for the two interventions
#'   - `risk_trt`, `risk_ctrl`: marginal risks on the risk scale
#'   - `contrast_rr`: `lmtp_contrast(type = "rr")` output
#'   - `contrast_rd`: `lmtp_contrast(type = "additive")` output, on the
#'     survival scale (negate to get the risk-scale RD)
#'   - (if `add_baseline_only = TRUE`) additional `*_baseline_only` fields
run_lmtp_analysis <- function(lmtp_prep,
                              shift_on  = NULL,
                              shift_off = NULL,
                              add_baseline_only = FALSE,
                              folds = 2,
                              learners = c("SL.glm")) {
  requireNamespace("lmtp", quietly = TRUE)
  
  if (is.null(shift_on))  shift_on  <- lmtp::static_binary_on
  if (is.null(shift_off)) shift_off <- lmtp::static_binary_off
  
  # Use time-varying treatment columns if available (from prepare_lmtp_data_tv)
  trt_spec <- if (!is.null(lmtp_prep$A_cols)) lmtp_prep$A_cols else "treatment"
  
   
  
  common_args <- list(
    data    = lmtp_prep$data,
    trt     = trt_spec,
    outcome = lmtp_prep$Y_cols,
    cens    = lmtp_prep$C_cols,
    baseline = lmtp_prep$baseline,
    outcome_type = "survival",
    folds   = folds,
    learners_trt     = learners,
    learners_outcome = learners
  )
  
  # Full always_on vs always_off
  res_on  <- do.call(lmtp::lmtp_sdr,
                     c(common_args, list(shift = shift_on)))
  res_off <- do.call(lmtp::lmtp_sdr,
                     c(common_args, list(shift = shift_off)))
  
  contrast_rr <- lmtp::lmtp_contrast(res_on, ref = res_off, type = "rr")
  contrast_rd <- lmtp::lmtp_contrast(res_on, ref = res_off, type = "additive")
  
 out<- list(
    res_on      = res_on,
    res_off     = res_off,
    risk_trt    = 1 - res_on$theta[length(res_on$theta)],
    risk_ctrl   = 1 - res_off$theta[length(res_off$theta)],
    contrast_rr = contrast_rr,
    contrast_rd = contrast_rd
  )
  
 # Optional: baseline-only intervention on wide-format time-varying data.
 # Intervenes on A1 only, leaving A2..An at observed values.
 # Uses lmtp's `shifted=` argument rather than a shift function.
 # Not used by the main simulation pipeline, where the treatment-policy
 # estimand is implemented via `time_varying_trt = FALSE` in prepare_lmtp_data().
  if (add_baseline_only && !is.null(lmtp_prep$A_cols)) {
    baseline_trt <- lmtp_prep$A_cols[1]

    shifted_baseline_on <- make_shifted_baseline(
      data = lmtp_prep$data, baseline_trt = baseline_trt,
      value = 1, C_cols = lmtp_prep$C_cols)
    shifted_baseline_off <- make_shifted_baseline(
      data = lmtp_prep$data, baseline_trt = baseline_trt,
      value = 0, C_cols = lmtp_prep$C_cols)

    res_baseline_on <- do.call(lmtp::lmtp_sdr,
      c(common_args, list(shifted = shifted_baseline_on, mtp = FALSE)))
    res_baseline_off <- do.call(lmtp::lmtp_sdr,
      c(common_args, list(shifted = shifted_baseline_off, mtp = FALSE)))

    out$res_baseline_on  <- res_baseline_on
    out$res_baseline_off <- res_baseline_off
    out$risk_baseline_on  <- 1 - res_baseline_on$theta[length(res_baseline_on$theta)]
    out$risk_baseline_off <- 1 - res_baseline_off$theta[length(res_baseline_off$theta)]
    out$contrast_rr_baseline_only <-
      lmtp::lmtp_contrast(res_baseline_on, ref = res_baseline_off, type = "rr")
    out$contrast_rd_baseline_only <-
      lmtp::lmtp_contrast(res_baseline_on, ref = res_baseline_off, type = "additive")
  }
 
 
 out
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

#' Bootstrap confidence interval for the KM risk difference.
#' @param dat data.frame with time_use, event_use, treatment columns.
#' @param t time point for risk extraction.
#' @param n_boot integer; number of bootstrap replicates.
#' @param weights optional numeric vector of observation weights.
#' @return named list: rd, rd_ci_low, rd_ci_high.
boot_km_rd_ci <- function(dat, t, n_boot = 200, weights = NULL) {
  rd_boot <- numeric(n_boot)
  n <- nrow(dat)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    d_b <- dat[idx, ]
    if (!is.null(weights)) {
      w_b <- weights[idx]
      sf_b <- tryCatch(
        survfit(Surv(time_use, event_use) ~ treatment, data = d_b, weights = w_b),
        error = function(e) NULL)
    } else {
      sf_b <- tryCatch(
        survfit(Surv(time_use, event_use) ~ treatment, data = d_b),
        error = function(e) NULL)
    }
    if (is.null(sf_b)) { rd_boot[b] <- NA; next }
    km_b <- km_risk_at(sf_b, t)
    rd_boot[b] <- km_b$risk_diff
  }
  rd_boot <- rd_boot[!is.na(rd_boot)]
  if (length(rd_boot) < 10) return(list(rd_ci_low = NA, rd_ci_high = NA))
  list(rd_ci_low = quantile(rd_boot, 0.025),
       rd_ci_high = quantile(rd_boot, 0.975))
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