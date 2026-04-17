library(tidyverse)


generate_hep_data <- function(
    N               = 125000,
    ## treatment assignment
    p_trt_target    = 0.36,
    ## baseline hazard parameters
    h0              = 5e-4,
    HR_early        = 1.30,        # TDF vs ETV early HR (first tau days)
    HR_late         = 0.90,        # TDF vs ETV late HR (attenuates; modest reversal)
    tau             = 90,          # change-point for non-PH (days)
    max_follow      = 720,         # admin cut-off (days)
    ## behaviour toggles
    np_hazard       = TRUE,       # TRUE -> HR varies over time (non-PH)
    dep_censor      = TRUE,       # TRUE -> censoring depends on risk
    complexity      = TRUE,       # TRUE -> nonlinear PS & outcome models
    ## informative switching
    switch_on       = TRUE,       # TRUE -> informative switching
    lambda_sw0      = 1.0e-3,     # baseline switch hazard (~15-20% by 180d)
    gamma_A         = 0.80,       # log-HR for treatment on switch
    gamma_ckd       = 0.60,       # log-HR per CKD on switch
    ## misc
    censor_base     = 1/100,      # admin censoring rate
    treat_override  = c("simulate", "random", "all_treated", "all_control"),
    return_potential_switching = FALSE,
    add_missing     = FALSE,
    impute          = FALSE,
    seed            = NULL)
{
  if (!is.null(seed)) set.seed(seed)
  treat_override <- match.arg(treat_override)
  if (impute) requireNamespace("missForest")
  requireNamespace("tidyverse")

  ##--------------------------------------------------------------------------
  ## 1. Demographics ---------------------------------------------------------
  raw <- tibble::tibble(
    id           = seq_len(N),
    age          = pmax(rnorm(N, 48, 13), 18),
    sex_male     = rbinom(N, 1, 0.58),
    race         = sample(c("white", "black", "hispanic", "asian", "other"),
                          N, TRUE, prob = c(.48, .14, .06, .02, .30)))

  ## 2. Comorbidities & co-medications ---------------------------------------
  add_bin <- function(p) rbinom(N, 1, p)
  raw <- raw %>%
    mutate(
      ckd            = add_bin(.08),
      diabetes       = add_bin(.20),
      hypertension   = add_bin(.45),
      cirrhosis      = add_bin(.18),
      heart_failure  = add_bin(.07),
      nsaid          = add_bin(.25),
      acearb         = add_bin(.30),
      statin         = add_bin(.15))

  ##--------------------------------------------------------------------------
  ## 3. Treatment assignment -------------------------------------------------
  cohort <- raw
  if (treat_override == "simulate") {
    lp0 <- with(cohort,
                0.015 * age + 0.30 * cirrhosis + 0.25 * ckd +
                  0.10 * diabetes + 0.05 * heart_failure +
                  rnorm(nrow(cohort), 0, 0.6))
    if (complexity) {
      lp0 <- with(cohort, lp0 +
                    0.5 * (age / 50)^3 + 1.5 * ckd * cirrhosis)
    }
    alpha0  <- qlogis(p_trt_target) - mean(lp0)
    p_trt   <- plogis(alpha0 + lp0) |> pmin(0.95) |> pmax(0.05)
    cohort$treatment <- rbinom(nrow(cohort), 1, p_trt)
  } else if (treat_override == "random") {
    cohort$treatment <- rbinom(nrow(cohort), 1, p_trt_target)
  } else {
    cohort$treatment <- ifelse(treat_override == "all_treated", 1L, 0L)
  }

  ##--------------------------------------------------------------------------
  ## 4. Individual baseline hazard -------------------------------------------
  if (!complexity) {
    lp_out <- with(cohort,
                   -2.8 + 0.03 * age + 0.7 * ckd + 0.5 * cirrhosis +
                     0.3 * heart_failure + 0.25 * nsaid)
  } else {
    lp_out <- with(cohort,
                   -2.8 + 0.03 * age + 0.0005 * age^2 + 0.7 * ckd +
                     0.5 * cirrhosis + 0.4 * heart_failure * acearb +
                     0.6 * nsaid * cohort$treatment)
  }
  base_rate <- h0 * exp(lp_out)

  ##--------------------------------------------------------------------------
  ## 5. Treatment switching --------------------------------------------------
  if (!switch_on) {
    cohort$switch_time <- Inf
    if (return_potential_switching) cohort$never_switcher <- 1L
  } else {
    if (return_potential_switching) {
      rng_state <- .Random.seed
      lambda_sw_a1 <- lambda_sw0 * exp(gamma_A * 1 + gamma_ckd * cohort$ckd)
      sw_a1 <- rexp(nrow(cohort), rate = lambda_sw_a1)
      .Random.seed <<- rng_state
      lambda_sw_a0 <- lambda_sw0 * exp(gamma_A * 0 + gamma_ckd * cohort$ckd)
      sw_a0 <- rexp(nrow(cohort), rate = lambda_sw_a0)
      cohort$switch_time <- ifelse(cohort$treatment == 1, sw_a1, sw_a0)
      cohort$never_switcher <- as.integer(sw_a1 > max_follow &
                                            sw_a0 > max_follow)
    } else {
      lambda_sw <- lambda_sw0 *
        exp(gamma_A * cohort$treatment + gamma_ckd * cohort$ckd)
      cohort$switch_time <- rexp(nrow(cohort), rate = lambda_sw)
    }
  }

  ##--------------------------------------------------------------------------
  ## 6. Event times (renal failure) ------------------------------------------
  # When switch_on = TRUE, patients who switch have their post-switch
  # hazard modified to reflect the new treatment. The event time is
  # generated in two pieces:
  #   [0, switch_time): original treatment hazard
  #   [switch_time, Inf): switched treatment hazard

  if (!np_hazard) {
    ## ── Proportional hazards ──
    rate_orig   <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    rate_switch <- base_rate * ifelse(cohort$treatment == 1, 1, HR_early)
    event_time_orig <- rexp(nrow(cohort), rate = rate_orig)

    if (switch_on) {
      did_sw <- cohort$switch_time < event_time_orig
      residual <- rexp(nrow(cohort), rate = rate_switch)
      cohort$event_time <- ifelse(did_sw,
                                  cohort$switch_time + residual,
                                  event_time_orig)
    } else {
      cohort$event_time <- event_time_orig
    }

  } else {
    ## ── Non-PH: change-point in HR at tau days ──
    rpexp_piece <- function(n, r1, r2, cp) {
      u  <- runif(n); p1 <- 1 - exp(-r1 * cp); t <- numeric(n)
      e  <- u <= p1
      t[e]  <- -log(1 - u[e]) / r1[e]
      t[!e] <- cp - log((1 - u[!e]) / (1 - p1[!e])) / r2[!e]
      t
    }

    r1_orig <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    r2_orig <- base_rate * ifelse(cohort$treatment == 1, HR_late,  1)
    event_time_orig <- rpexp_piece(nrow(cohort), r1_orig, r2_orig, tau)

    if (switch_on) {
      r1_switch <- base_rate * ifelse(cohort$treatment == 1, 1, HR_early)
      r2_switch <- base_rate * ifelse(cohort$treatment == 1, 1, HR_late)
      did_sw <- cohort$switch_time < event_time_orig
      sw_t <- cohort$switch_time
      post_rate <- ifelse(sw_t <= tau, r1_switch, r2_switch)
      residual <- rexp(nrow(cohort), rate = post_rate)
      cohort$event_time <- ifelse(did_sw, sw_t + residual, event_time_orig)
    } else {
      cohort$event_time <- event_time_orig
    }
  }

  ##--------------------------------------------------------------------------
  ## 7. Administrative censoring ---------------------------------------------
  if (!dep_censor) {
    censor_admin <- rep(max_follow, nrow(cohort))
  } else {
    c_rate <- censor_base * exp(0.04 * lp_out + 0.03 * cohort$treatment)
    censor_admin <- rexp(nrow(cohort), rate = c_rate)
  }
  cohort$censor_admin <- pmin(censor_admin, max_follow)

  ##--------------------------------------------------------------------------
  ## 8. Construct analysis variables -----------------------------------------
  # The dataset is generated under the treatment-policy regime:
  # switching occurs, modifies the hazard, but does NOT censor.
  # All raw timing variables are retained so that the analysis code
  # can derive any estimand by redefining outcome/censoring:
  #
  # Treatment-policy:   follow_time = min(event_time, censor_admin)
  #                     event = I(event_time <= follow_time)
  # While-on-treatment: follow_time = min(event_time, switch_time, censor_admin)
  #                     event = I(event_time <= follow_time)
  # Composite:          follow_time = min(event_time, switch_time, censor_admin)
  #                     event = I(min(event_time, switch_time) <= censor_admin)
  # No-switch:          requires a separate DGP call with switch_on=FALSE

  # Default follow-up: treatment-policy (no censoring at switch)
  cohort$follow_time <- pmin(cohort$event_time, cohort$censor_admin)
  cohort$event       <- as.integer(cohort$event_time <= cohort$follow_time)
  cohort$switched    <- as.integer(cohort$switch_time < cohort$follow_time)
  cohort$would_switch <- as.integer(cohort$switch_time <= max_follow)

  ##--------------------------------------------------------------------------
  ## 9. Analysis dataset -----------------------------------------------------
  ana <- cohort %>%
    dplyr::select(
      id, age, sex_male, race,
      ckd, diabetes, hypertension, cirrhosis, heart_failure,
      nsaid, acearb, statin,
      treatment, event_time, switch_time, follow_time, event, switched,
      censor_admin, would_switch,
      any_of("never_switcher")
    )

  ##--------------------------------------------------------------------------
  ## 10. Optional missingness / imputation -----------------------------------
  if (add_missing) {
    ana$race[sample(nrow(ana), 0.05 * nrow(ana))] <- NA
    ana$ckd[sample(nrow(ana), 0.10 * nrow(ana))]  <- NA
    if (impute) {
      imp_vars <- c("age", "race", "ckd", "cirrhosis", "diabetes",
                     "hypertension")
      imp_in   <- ana %>% dplyr::select(dplyr::all_of(imp_vars)) %>%
        dplyr::mutate(across(c(race), as.factor))
      ana[, imp_vars] <- missForest::missForest(as.data.frame(imp_in),
                                                verbose = FALSE)$ximp
    }
  }

  return(ana)
}


# ── Utility: derive estimand-specific outcome from a treatment-policy dataset ──
#' Given a dataset generated under treatment-policy, construct the
#' follow_time and event columns for a specific estimand strategy.
#' This allows one dataset to support all five estimand analyses.
#' @param dat data.frame from generate_hep_data() (treatment-policy).
#' @param strategy character; one of "treatment_policy", "while_on_treatment",
#'   "composite".
#' @return data.frame with follow_time and event redefined.
derive_estimand <- function(dat, strategy = c("treatment_policy",
                                               "while_on_treatment",
                                               "composite")) {
  strategy <- match.arg(strategy)
  out <- dat

  if (strategy == "treatment_policy") {
    # Already the default from generate_hep_data()
    out$follow_time <- pmin(dat$event_time, dat$censor_admin)
    out$event       <- as.integer(dat$event_time <= out$follow_time)

  } else if (strategy == "while_on_treatment") {
    # Censor at switch time
    out$follow_time <- pmin(dat$event_time, dat$switch_time, dat$censor_admin)
    out$event       <- as.integer(dat$event_time <= out$follow_time)

  } else if (strategy == "composite") {
    # Outcome = renal failure OR switching (whichever first)
    out$follow_time <- pmin(dat$event_time, dat$switch_time, dat$censor_admin)
    out$event       <- as.integer(
      pmin(dat$event_time, dat$switch_time) <= dat$censor_admin
    )
  }

  out$switched <- as.integer(dat$switch_time < out$follow_time)
  out
}


# ── Generate and cache example datasets ──────────────────────────────────────
if (sys.nframe() == 0) {
  set.seed(1234)
  df <- generate_hep_data(
    h0 = 3e-4, np_hazard = FALSE, dep_censor = FALSE,
    complexity = FALSE, seed = 1234
  )
  c(mean(df$follow_time), mean(df$event) * 100,
    mean(df$switched) * 100)
  write.csv(df, here::here("data/sim_hep_renal.csv"), row.names = FALSE)

  df_complex <- generate_hep_data(
    np_hazard = TRUE, dep_censor = TRUE, complexity = TRUE, seed = 1234
  )
  c(mean(df_complex$follow_time), mean(df_complex$event) * 100,
    mean(df_complex$switched) * 100)
  write.csv(df_complex, here::here("data/sim_hep_renal_complex.csv"),
            row.names = FALSE)
}
