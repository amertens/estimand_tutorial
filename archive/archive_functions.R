# archive_functions.R
# ============================================================================
# Useful functions salvaged from archive/ before deletion.
# Each section is annotated with the original source file.
#
# Contents:
#   1. Ground-truth computation    (from archive/DGP_test.R)
#   2. Landmark analysis           (from archive/simstudy3.R)
#   3. Multi-estimator comparison  (from archive/simstudy3.R)
#   4. Simulation performance      (from archive/simstudy3.R)
#   5. PS matching + Cox           (from archive/ps_aki_analysis_par.R,
#                                        archive/PS_analysis.Rmd)
#   6. Cox diagnostics             (from archive/PS_analysis.Rmd)
#   7. IPTW Cox                    (from archive/PS_analysis.Rmd,
#                                        archive/simstudy3.R)
#   8. IPCW no-switch estimand     (from archive/regression_estimation.Rmd)
#   9. survtmle + cumulative HR    (from archive/survTMLE_example.R)
#  10. concrete (competing risks)  (from archive/concrete_example.R)
#  11. survivalSL ensemble         (from archive/survivalSL.R)
#  12. tmle3 / tlverse workflow    (from archive/tmle_diagnostics.R)
#  13. E-value sensitivity         (from archive/tmle_diagnostics.R)
#  14. LMTP shift function factory (from archive/lmtp_tutorial.Rmd)
#  15. Longitudinal ART DGP        (from archive/lmtp_tutorial.Rmd)
#  16. Incidence rate table         (from archive/gen_data_v2.R)
# ============================================================================


# ============================================================================
# 1. GROUND-TRUTH COMPUTATION
#    Source: archive/DGP_test.R
#    Computes true causal risk difference and risk ratio from large
#    counterfactual datasets (all-treated vs all-control, no censoring).
# ============================================================================

#' Compute true risk at a given time horizon from counterfactual datasets.
#' @param df0 data.frame from generate_hcv_data(treat_override="all_control").
#' @param df1 data.frame from generate_hcv_data(treat_override="all_treated").
#' @param t_star numeric; time horizon in days.
#' @return data.frame with columns t, RD, RR.
calc_true_risk <- function(df0, df1, t_star) {
  risk_untreated <- mean(df0$follow_time <= t_star & df0$event == 1)
  risk_treated   <- mean(df1$follow_time <= t_star & df1$event == 1)
  risk_diff      <- risk_treated - risk_untreated
  risk_ratio     <- risk_treated / risk_untreated

  cat(t_star, "-day risk (untreated):", round(risk_untreated, 4), "\n")
  cat(t_star, "-day risk (treated):  ", round(risk_treated, 4), "\n")
  cat("Risk difference:        ", round(risk_diff, 4), "\n")
  cat("Risk ratio:             ", round(risk_ratio, 4), "\n")

  data.frame(t = t_star, RD = risk_diff, RR = risk_ratio)
}


#' Compute true period-specific hazard ratio via counting-process Cox
#' on large counterfactual data.
#' @param Nbig integer; large sample size (default 1e6).
#' @param start numeric; start of time window (days).
#' @param end numeric; end of time window (days).
#' @param ... additional args passed to generate_hcv_data().
#' @return named numeric; the true HR in [start, end].
true_hr_period <- function(Nbig = 1000000, start = 0, end = 90, ...) {
  requireNamespace("survival", quietly = TRUE)

  ctrl <- generate_hcv_data(Nbig, treat_override = "all_control",
                            max_follow = 360, ...) |>
    dplyr::select(follow_time, event) |>
    dplyr::mutate(A = 0)

  trt <- generate_hcv_data(Nbig, treat_override = "all_treated",
                           max_follow = 360, ...) |>
    dplyr::select(follow_time, event) |>
    dplyr::mutate(A = 1)

  dat <- dplyr::bind_rows(ctrl, trt)

  dat_cp <- dat |>
    dplyr::mutate(
      start  = pmax(0, pmin(follow_time, start)),
      stop   = pmin(follow_time, end),
      eventW = as.integer(event == 1 & follow_time <= end &
                            follow_time > start)
    ) |>
    dplyr::filter(stop > start)

  fit <- survival::coxph(survival::Surv(start, stop, eventW) ~ A,
                         data = dat_cp)
  exp(coef(fit))
}


# ============================================================================
# 2. LANDMARK ANALYSIS
#    Source: archive/simstudy3.R
#    Converts survival data to a binary landmark outcome at time t0.
# ============================================================================

#' Create a landmark binary outcome from survival data.
#' @param dat data.frame with follow_time and event columns.
#' @param t0 numeric; landmark time.
#' @return dat with Y (1=event by t0, 0=survived past t0, NA=censored before
#'   t0) and Delta (1=observed, 0=censored) columns added.
make_landmark <- function(dat, t0) {
  dat |> dplyr::mutate(
    Y     = dplyr::case_when(
      event == 1 & follow_time <= t0 ~ 1,
      follow_time >= t0              ~ 0,
      TRUE ~ NA_real_
    ),
    Delta = ifelse(is.na(Y), 0, 1)
  )
}


# ============================================================================
# 3. MULTI-ESTIMATOR COMPARISON (SINGLE ITERATION)
#    Source: archive/simstudy3.R
#    Runs TMLE, crude Cox, IPTW Cox, and PS-matched Cox in a single
#    iteration of the DGP.
# ============================================================================

#' Run TMLE + crude Cox + IPTW Cox + PS-matched Cox for one dataset.
#' Requires: tmle, survival, survey, MatchIt packages.
#' @param N integer; sample size.
#' @param horizons numeric vector; landmark time horizons.
#' @return tibble with R1, R0, RD, RR, hr_crude, hr_iptw, hr_psm.
analyze_once <- function(N = 20000, horizons = c(30, 90, 180)) {
  dat   <- generate_hcv_data(N = N, np_hazard = TRUE,
                             dep_censor = TRUE, complexity = TRUE)
  Wvars <- c("age", "sex_male", "ckd", "cirrhosis", "diabetes",
             "hypertension", "bmi", "race", "region")

  lapply(horizons, function(t0) {
    d0 <- make_landmark(dat, t0)
    W  <- d0[, Wvars]

    # TMLE
    fit <- tmle::tmle(Y = d0$Y, A = d0$treatment, W = W, Delta = d0$Delta,
                      family = "binomial",
                      Q.SL.library = c("SL.glm", "SL.mean"),
                      g.SL.library = "SL.glm",
                      g.Delta.SL.library = "SL.glm")
    R1 <- fit$estimates$EY1$psi
    R0 <- fit$estimates$EY0$psi
    RD <- R1 - R0
    RR <- R1 / R0

    # Crude Cox HR
    hr_crude <- exp(coef(survival::coxph(
      survival::Surv(follow_time, event) ~ treatment, data = dat)))

    # IPTW Cox via survey::svycoxph
    ps <- glm(treatment ~ ., data = dat[, c("treatment", Wvars)],
              family = binomial)$fitted
    w  <- ifelse(dat$treatment == 1, 1 / ps, 1 / (1 - ps))
    des <- survey::svydesign(ids = ~1, weights = ~w, data = dat)
    hr_iptw <- exp(coef(survey::svycoxph(
      survival::Surv(follow_time, event) ~ treatment, design = des)))

    # PS-matched Cox
    ps_form <- as.formula(paste("treatment ~", paste(Wvars, collapse = "+")))
    m.out <- MatchIt::matchit(ps_form,
                              data = dat[, c("treatment", "follow_time",
                                             "event", Wvars)],
                              method = "nearest", ratio = 1, caliper = 0.20)
    d.m <- MatchIt::match.data(m.out, data = "all")
    hr_psm <- exp(coef(survival::coxph(
      survival::Surv(follow_time, event) ~ treatment,
      data = d.m, cluster = subclass)))

    tibble::tibble(t0, R1, R0, RD, RR, hr_crude, hr_iptw, hr_psm)
  }) |> dplyr::bind_rows()
}


# ============================================================================
# 4. SIMULATION PERFORMANCE SUMMARY
#    Source: archive/simstudy3.R
#    Computes bias, oracle SE, and coverage for multiple estimands.
# ============================================================================

#' Summarize simulation performance across estimands.
#' @param res data.frame of replicate results.
#' @param truth named list/data.frame with true values for each estimand.
#' @param est_cols character vector of estimand column names.
#' @param lwr_suffix, upr_suffix suffixes for CI columns.
#' @return tibble with mean, sd, bias, coverage per estimand.
perf_stats <- function(res, truth,
                       est_cols   = c("RD", "RR", "hr_crude", "hr_psm"),
                       lwr_suffix = "_lwr", upr_suffix = "_upr") {
  one_stat <- function(col) {
    x <- res[[col]]
    tibble::tibble(
      estimand = col,
      mean     = mean(x),
      sd       = sd(x),
      bias     = mean(x) - truth[[col]],
      cover    = {
        lwr <- res[[paste0(col, lwr_suffix)]]
        upr <- res[[paste0(col, upr_suffix)]]
        if (!is.null(lwr) && !is.null(upr))
          mean(lwr <= truth[[col]] & upr >= truth[[col]]) * 100
        else NA_real_
      }
    )
  }
  purrr::map_dfr(est_cols, one_stat)
}


# ============================================================================
# 5. PS MATCHING + MATCHED COX
#    Source: archive/ps_aki_analysis_par.R, archive/PS_analysis.Rmd
#    Propensity score matching via MatchIt and Cox on matched data.
# ============================================================================

#' Fit a PS-matched Cox model.
#' @param dat data.frame with treatment, follow_time, event, and covariates.
#' @param ps_formula formula for propensity score model.
#' @param caliper numeric; matching caliper (default 0.2).
#' @return list with hr, ci_low, ci_high, match_object, matched_data.
fit_cox_ps_matched <- function(dat,
                               ps_formula = treatment ~ age + sex_male +
                                 ckd + cirrhosis + hiv + diabetes +
                                 hypertension + bmi,
                               caliper = 0.2) {
  ps_mod <- glm(ps_formula, data = dat, family = binomial)
  dat$ps <- predict(ps_mod, type = "response")

  m_out <- MatchIt::matchit(ps_formula, data = dat, method = "nearest",
                            caliper = caliper, ratio = 1)
  m_dat <- MatchIt::match.data(m_out)

  cox_fit <- survival::coxph(survival::Surv(follow_time, event) ~ treatment,
                             data = m_dat)
  hr_tab <- broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE)
  hr_row <- hr_tab[hr_tab$term == "treatment", ]

  list(
    hr           = hr_row$estimate,
    ci_low       = hr_row$conf.low,
    ci_high      = hr_row$conf.high,
    match_object = m_out,
    matched_data = m_dat,
    cox_model    = cox_fit
  )
}


# ============================================================================
# 6. COX DIAGNOSTICS
#    Source: archive/PS_analysis.Rmd
#    Schoenfeld residual (PH) test and Cox-Snell residual plot.
# ============================================================================

#' Run proportional hazards diagnostic (Schoenfeld test + plot).
#' @param cox_fit coxph object.
#' @param plot logical; whether to plot Schoenfeld residuals.
#' @return cox.zph test result.
cox_ph_diagnostic <- function(cox_fit, plot = TRUE) {
  scho <- survival::cox.zph(cox_fit)
  if (plot) plot(scho)
  scho
}

#' Plot Cox-Snell residuals for model fit assessment.
#' @param cox_fit coxph object.
cox_snell_plot <- function(cox_fit) {
  M <- residuals(cox_fit, type = "martingale")
  event_m <- cox_fit$y[, "status"]
  cs_resid <- event_m - M
  km_cs <- survival::survfit(survival::Surv(cs_resid, event_m) ~ 1)
  plot(km_cs$time, -log(km_cs$surv), type = "l",
       xlab = "Cox-Snell residual", ylab = "Cumulative hazard",
       main = "Cox-Snell residual plot")
  abline(0, 1, col = "red", lty = 2)
}


# ============================================================================
# 7. IPTW COX
#    Source: archive/PS_analysis.Rmd, archive/simstudy3.R
#    Inverse probability of treatment weighted Cox model.
# ============================================================================

#' Fit an IPTW-weighted Cox model using survey::svycoxph.
#' @param dat data.frame with treatment, follow_time, event, and covariates.
#' @param covars character vector of covariate names for PS model.
#' @return list with hr, model, weights.
fit_cox_iptw <- function(dat,
                         covars = c("age", "ckd", "cirrhosis", "diabetes",
                                    "hiv", "hypertension", "bmi")) {
  covars <- intersect(covars, names(dat))
  fml <- as.formula(paste("treatment ~", paste(covars, collapse = " + ")))
  ps <- glm(fml, data = dat, family = binomial)$fitted.values
  w  <- ifelse(dat$treatment == 1, 1 / ps, 1 / (1 - ps))
  dat$w_iptw <- w

  des <- survey::svydesign(ids = ~1, weights = ~w_iptw, data = dat)
  fit <- survey::svycoxph(
    survival::Surv(follow_time, event) ~ treatment, design = des)

  list(
    hr    = exp(coef(fit))[["treatment"]],
    model = fit,
    weights = w
  )
}


# ============================================================================
# 8. IPCW NO-SWITCH ESTIMAND
#    Source: archive/regression_estimation.Rmd
#    Inverse probability of censoring weighting for the hypothetical
#    no-switch estimand via a pooled logistic censoring model.
# ============================================================================

#' Expand person-day data and compute IPCW weights for switch censoring.
#' @param dat data.frame with id, follow_time, event, switch, treatment,
#'   and baseline covariates.
#' @param cutoff numeric; analysis time horizon.
#' @param covars character vector of baseline covariates for censoring model.
#' @return data.frame with IPCW weights (column w) merged back to subject
#'   level.
compute_ipcw_weights <- function(dat, cutoff = 90,
                                 covars = c("age", "ckd", "cirrhosis",
                                            "diabetes", "bmi",
                                            "sex_male")) {
  covars <- intersect(covars, names(dat))

  dat <- dat %>%
    dplyr::mutate(
      trt = treatment,
      Tstar = pmin(follow_time, cutoff),
      admin_cens   = as.integer(follow_time > cutoff),
      switch_cens  = as.integer(switch == 1 & follow_time > 0),
      C = pmax(admin_cens, switch_cens),
      event_ipcw = as.integer(event == 1 & follow_time <= cutoff & C == 0)
    )

  # Person-day expansion
  long <- purrr::map_dfr(seq_len(nrow(dat)), function(i) {
    tibble::tibble(
      id       = dat$id[i],
      trt      = dat$trt[i],
      day      = seq_len(dat$Tstar[i]),
      censored = as.integer(day == dat$Tstar[i] & dat$C[i] == 1),
      event    = as.integer(day == dat$Tstar[i] & dat$event_ipcw[i] == 1)
    ) %>%
      dplyr::bind_cols(dat[i, covars, drop = FALSE])
  })

  # Pooled logistic model for censoring with natural splines for time
  fitC <- glm(censored ~ trt + splines::ns(day, 4) + .,
              family = binomial(), data = long)
  long <- long %>%
    dplyr::mutate(pC = predict(fitC, type = "response")) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(Surv_uncens = cumprod(1 - pC)) %>%
    dplyr::ungroup()

  weights <- long %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(w = dplyr::last(Surv_uncens), .groups = "drop")

  dat <- dat %>% dplyr::left_join(weights, by = "id")
  dat$w <- pmin(dat$w, quantile(dat$w, 0.995, na.rm = TRUE))
  dat
}


# ============================================================================
# 9. survtmle + CUMULATIVE HAZARD RATIO FROM INFLUENCE CURVES
#    Source: archive/survTMLE_example.R
#    Fits survtmle and derives a cumulative-hazard ratio with IC-based SE.
# ============================================================================

#' Fit survtmle and compute cumulative-hazard ratio with IC-based CI.
#' Requires: survtmle, SuperLearner.
#' @param ftime integer vector of event/censoring times.
#' @param ftype integer vector (0=censored, 1=event).
#' @param trt integer vector (0/1 treatment).
#' @param W data.frame of adjustment covariates.
#' @param t0 numeric; time horizon.
#' @param sl_lib_f, sl_lib_c, sl_lib_g character vectors of SL learners.
#' @return list with cumHR, ci, survtmle_fit.
fit_survtmle_cumhr <- function(ftime, ftype, trt, W, t0 = 5,
                               sl_lib_f = c("SL.glm", "SL.mean"),
                               sl_lib_c = "SL.glm",
                               sl_lib_g = "SL.glm") {
  fit <- survtmle::survtmle(
    ftime      = ftime,
    ftype      = ftype,
    trt        = trt,
    adjustVars = W,
    t0         = t0,
    SL.ftime   = sl_lib_f,
    SL.ctime   = sl_lib_c,
    SL.trt     = sl_lib_g,
    method     = "hazard"
  )

  S1 <- fit$est[1, 1]
  S0 <- fit$est[2, 1]
  Lambda1 <- -log(S1)
  Lambda0 <- -log(S0)
  cumHR   <- Lambda1 / Lambda0

  # IC-based SE for log(cumHR)
  n  <- length(ftime)
  ic <- (fit$ic[, 2] - fit$ic[, 1]) / (S1 * log(S0 / S1))
  se <- sd(ic) / sqrt(n)
  ci <- exp(log(cumHR) + c(-1, 1) * 1.96 * se)

  list(cumHR = cumHR, ci = ci, survtmle_fit = fit)
}


# ============================================================================
# 10. COMPETING-RISKS TMLE VIA concrete PACKAGE
#     Source: archive/concrete_example.R
#     One-step TMLE for competing risks using the PBC dataset as a template.
# ============================================================================

#' Run concrete TMLE for competing risks.
#' Requires: concrete, survival.
#' @param data data.frame with time, status, trt columns.
#' @param target_times numeric vector of evaluation times.
#' @param target_events integer vector of competing event types.
#' @param cv_folds integer; cross-validation folds.
#' @return list with concrete_est, output (RD with simultaneous CIs).
run_concrete_tmle <- function(data, target_times, target_events = 1:2,
                              cv_folds = 10) {
  args <- concrete::formatArguments(
    DataTable    = data,
    EventTime    = "time",
    EventType    = "status",
    Treatment    = "trt",
    Intervention = 0:1,
    TargetTime   = target_times,
    TargetEvent  = target_events,
    CVArg        = list(V = cv_folds),
    Verbose      = FALSE
  )
  est <- concrete::doConcrete(args)
  out <- concrete::getOutput(est, Estimand = "RD", Simultaneous = TRUE)
  list(concrete_est = est, output = out)
}


# ============================================================================
# 11. SURVIVAL SUPER LEARNER VIA survivalSL
#     Source: archive/survivalSL.R
#     Ensemble survival prediction using multiple survival model libraries.
# ============================================================================

#' Fit a survival Super Learner ensemble.
#' Requires: survivalSL.
#' @param train_data data.frame with follow_time, event, and covariates.
#' @param cov_quanti character vector of continuous covariate names.
#' @param cov_quali character vector of binary/factor covariate names.
#' @param methods character vector of survivalSL methods (default:
#'   Cox elastic-net, AFT gamma, PH exponential).
#' @return survivalSL model object.
fit_survival_sl <- function(train_data,
                            cov_quanti = c("age"),
                            cov_quali = c("treatment", "sex_male", "ckd"),
                            methods = c("LIB_COXen", "LIB_AFTgamma",
                                        "LIB_PHexponential")) {
  survivalSL::survivalSL(
    methods    = methods,
    metric     = "ci",
    data       = as.data.frame(train_data),
    times      = "follow_time",
    failures   = "event",
    cov.quanti = cov_quanti,
    cov.quali  = cov_quali,
    progress   = TRUE
  )
}


# ============================================================================
# 12. tmle3 / tlverse WORKFLOW
#     Source: archive/tmle_diagnostics.R
#     Full tmle3 specification for survival TMLE with sl3 learner stack.
#     NOTE: requires the tlverse ecosystem (sl3, tmle3).
# ============================================================================

#' Set up and fit tmle3 survival TMLE.
#' Requires: tlverse (tmle3, sl3).
#' @param dat data.frame with treatment, event, follow_time, and covariates.
#' @param t_max numeric; target time for survival estimand.
#' @param sl_learners list of sl3 learner objects.
#' @return tmle3 fit object.
fit_tmle3_survival <- function(dat, t_max = 180, sl_learners = NULL) {
  node_list <- list(
    W = setdiff(names(dat), c("treatment", "event", "follow_time")),
    A = "treatment",
    Y = "event",
    Ttilde = "follow_time"
  )

  if (is.null(sl_learners)) {
    sl_learners <- list(
      sl3::Lrnr_glm$new(),
      sl3::Lrnr_glmnet$new(alpha = 0),
      sl3::Lrnr_glmnet$new(alpha = 1)
    )
  }

  tmle_sp <- tmle3::tmle3_Spec_survival(
    time          = t_max,
    contrast      = list(A = 1, B = 0),
    censoring_node = NULL,
    outcome_type  = "binary",
    marginal      = "RD"
  )

  task <- tmle_sp$make_tmle_task(dat, node_list)
  learner_stack <- sl3::make_learner(sl3::Stack, sl_learners)
  learner_list  <- list(Y = learner_stack, A = learner_stack)

  tmle3::tmle3(tmle_sp, task, learner_list)
}


# ============================================================================
# 13. E-VALUE SENSITIVITY ANALYSIS
#     Source: archive/tmle_diagnostics.R
#     Computes the E-value for unmeasured confounding sensitivity.
# ============================================================================

#' Compute E-value for a risk difference estimate.
#' Requires: EValue.
#' @param risk_rd numeric; estimated risk difference.
#' @param se_rd numeric; standard error of the risk difference.
#' @param baseline_risk numeric; baseline event rate for RR approximation.
#' @return E-value result from EValue::evalue().
compute_evalue <- function(risk_rd, se_rd, baseline_risk) {
  rr_est <- 1 + risk_rd / baseline_risk
  EValue::evalue(
    est  = rr_est,
    lo   = rr_est - 1.96 * se_rd,
    hi   = rr_est + 1.96 * se_rd,
    true = 1,
    type = "RR"
  )
}


# ============================================================================
# 14. LMTP SHIFT FUNCTION FACTORIES (DYNAMIC & STOCHASTIC)
#     Source: archive/lmtp_tutorial.Rmd
#     Creates shift functions for dynamic and stochastic modified treatment
#     policies for use with lmtp_sdr().
# ============================================================================

#' Factory for a dynamic "boost to threshold" shift function.
#' Forces regimen 1 (A=1) when adherence drops below thr.
#' @param thr numeric; adherence threshold (e.g., 0.8).
#' @return function(data, trt) suitable for lmtp shift argument.
make_dynamic_shift <- function(thr = 0.8) {
  function(data, trt) {
    pdc <- data[[gsub("^A", "PDC_prev", trt)]]
    at  <- data[[trt]]
    ifelse(pdc < thr, 1, at)
  }
}

#' Factory for a stochastic shift function.
#' When adherence < thr, switches to regimen 1 with probability p.
#' @param thr numeric; adherence threshold.
#' @param p numeric; probability of switching when below threshold.
#' @return function(data, trt) suitable for lmtp shift argument.
make_stochastic_shift <- function(thr = 0.8, p = 0.7) {
  function(data, trt) {
    pdc  <- data[[gsub("^A", "PDC_prev", trt)]]
    at   <- data[[trt]]
    need <- (pdc < thr)
    flip <- rbinom(length(at), 1, p)
    ifelse(need & flip == 1, 1, at)
  }
}


# ============================================================================
# 15. LONGITUDINAL ART COHORT DGP
#     Source: archive/lmtp_tutorial.Rmd
#     Simulates a longitudinal antiretroviral therapy cohort with
#     time-varying treatment, adherence (PDC), and censoring.
# ============================================================================

#' Simulate longitudinal ART cohort data.
#' @param N integer; number of subjects.
#' @param K integer; number of time blocks (e.g., 5 = 15 months).
#' @param censor_prob numeric; probability of being censored.
#' @param seed integer; random seed.
#' @return data.frame in long format with id, t, A, PDC, VL, C, Y, and
#'   baseline covariates.
sim_art_cohort <- function(N = 4000, K = 5, censor_prob = 0.15,
                           seed = 123) {
  set.seed(seed)
  id  <- 1:N
  age <- rnorm(N, 45, 10)
  sex <- rbinom(N, 1, 0.4)
  cd4 <- pmax(50, rnorm(N, 500, 150))

  DT <- data.table::CJ(id, t = 0:(K - 1))
  DT <- merge(DT, data.table::data.table(id, age, sex, cd4), by = "id")

  DT[t == 0, `:=`(
    A   = rbinom(.N, 1, plogis(0.5 - 0.01 * age + 0.5 * sex)),
    PDC = runif(.N, 0.5, 1)
  )]

  for (tt in 1:(K - 1)) {
    lag <- DT[t == tt - 1, .(id, A_prev = A, PDC_prev = PDC)]
    DT <- merge(DT, lag, by = "id", all.x = TRUE)
    DT[t == tt, A := rbinom(.N, 1,
                            plogis(1.5 * A_prev - 3 * (PDC_prev < 0.6) -
                                     0.01 * age))]
    DT[t == tt, PDC := pmin(1, pmax(0,
                                    rnorm(.N, 0.85 * A + 0.6 * (1 - A),
                                          0.2)))]
    DT[, c("A_prev", "PDC_prev") := NULL]
  }

  DT[, VL := rnorm(.N, 4 - 1.2 * A - 2 * PDC, 0.6)]

  Y <- DT[t == K - 1,
           .(id, Y = rbinom(.N, 1,
                            plogis(-3 - 0.4 * A - 2 * PDC +
                                     0.3 * (VL > 4))))]

  drop <- data.table::data.table(
    id, t_cens = rbinom(N, 1, censor_prob) * sample(1:K, N, TRUE))
  DT <- merge(DT, drop, by = "id")
  DT[, C := as.integer(t >= t_cens & t_cens > 0)]
  DT[, t_cens := NULL]

  long <- merge(DT, Y, by = "id")
  data.table::setorder(long, id, t)
  as.data.frame(long)
}


# ============================================================================
# 16. INCIDENCE RATE TABLE (PUBLICATION STYLE)
#     Source: archive/gen_data_v2.R
#     Computes arm-wise person-years, event counts, incidence rates,
#     rate ratios with Wald CIs, and rate differences.
# ============================================================================

#' Create a publication-style incidence rate table by treatment arm.
#' @param dat data.frame with treatment, follow_time, event.
#' @param tmax numeric; analysis window in days.
#' @return tibble formatted as a rate table with RR and RD rows.
make_rate_table <- function(dat, tmax = 180) {
  dat <- dat %>%
    dplyr::mutate(
      fup_d    = pmin(follow_time, tmax),
      event_t  = as.integer(event == 1 & follow_time <= tmax),
      fup_yrs  = fup_d / 365
    )

  tab <- dat %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarise(
      n_patients = dplyr::n(),
      person_yrs = sum(fup_yrs),
      n_events   = sum(event_t),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(rate_1000py = n_events / person_yrs * 1e3)

  rate0 <- tab$rate_1000py[tab$treatment == 0]
  rate1 <- tab$rate_1000py[tab$treatment == 1]
  py0   <- tab$person_yrs[tab$treatment == 0]
  py1   <- tab$person_yrs[tab$treatment == 1]
  ev0   <- tab$n_events[tab$treatment == 0]
  ev1   <- tab$n_events[tab$treatment == 1]

  RR       <- rate1 / rate0
  se_logRR <- sqrt(1 / ev1 + 1 / ev0)
  ciRR     <- exp(log(RR) + c(-1, 1) * 1.96 * se_logRR)

  RD    <- rate1 - rate0
  se_RD <- sqrt(ev1 / py1^2 + ev0 / py0^2) * 1e3
  ciRD  <- RD + c(-1, 1) * 1.96 * se_RD

  list(
    table     = tab,
    rate_ratio = list(RR = RR, ci = ciRR),
    rate_diff  = list(RD = RD, ci = ciRD)
  )
}
