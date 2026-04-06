library(tidyverse)


generate_hep_data <- function(
    N               = 125000,
    ## treatment assignment
    p_trt_target    = 0.36,
    ## baseline hazard parameters
    h0              = 5e-4,
    HR_early        = 1.25,
    HR_late         = 0.70,
    tau             = 45,          # change-point for non-PH (days)
    max_follow      = 720,         # admin cut-off (days)
    ## behaviour toggles
    np_hazard       = TRUE,       # TRUE -> HR varies over time (non-PH)
    dep_censor      = TRUE,       # TRUE -> censoring depends on risk; FALSE -> no censoring (only max_follow)
    complexity      = TRUE,       # TRUE -> nonlinear PS & outcome models
    ## informative switching
    switch_on       = TRUE,       # TRUE -> informative switching; FALSE -> no switching at all
    lambda_sw0      = 1.0e-3,     # baseline switch hazard (targets ~15-20% switching by 180d)
    gamma_A         = 0.80,       # log-HR for treatment on switch
    gamma_ckd       = 0.60,       # log-HR per CKD on switch
    ## estimand policy
    policy          = c("treatment_policy", "no_switch", "while_on_treatment",
                        "composite", "principal_stratum"),
    ## misc
    censor_base     = 1/100,      # admin censoring rate (1/days), only used when dep_censor=TRUE
    treat_override  = c("simulate", "all_treated", "all_control"),
    return_potential_switching = FALSE,  # TRUE: add never_switcher column for PS
    add_missing     = FALSE,
    impute          = FALSE,
    seed            = NULL)
{
  if (!is.null(seed)) set.seed(seed)
  treat_override <- match.arg(treat_override)
  policy         <- match.arg(policy)
  if (impute) requireNamespace("missForest")
  requireNamespace("tidyverse")

  # Under the no_switch policy, disable switching entirely so that

  # potential outcomes are generated in a world where nobody switches.
  if (policy == "no_switch") {
    switch_on <- FALSE
  }

  ##--------------------------------------------------------------------------
  ## 1. Demographics ---------------------------------------------------------
  raw <- tibble::tibble(
    id           = seq_len(N),
    age          = pmax(rnorm(N, 48, 13), 18),
    sex_male     = rbinom(N, 1, 0.58),
    race         = sample(c("white", "black", "hispanic", "asian", "other"), N, TRUE,
                          prob = c(.48, .14, .06, .02, .30)))

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

    ## baseline PS linear predictor
    lp0 <- with(cohort,
                0.015 * age + 0.30 * cirrhosis + 0.25 * ckd + 0.10 * diabetes +
                  0.05 * heart_failure + rnorm(nrow(cohort), 0, 0.6))

    if (complexity) {
      lp0 <- with(cohort, lp0 +
                    0.5 * (age / 50)^3 + 1.5 * ckd * cirrhosis)
    }
    alpha0  <- qlogis(p_trt_target) - mean(lp0)
    p_trt   <- plogis(alpha0 + lp0) |> pmin(0.95) |> pmax(0.05)
    cohort$treatment <- rbinom(nrow(cohort), 1, p_trt)
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
                   -2.8 + 0.03 * age + 0.0005 * age^2 + 0.7 * ckd + 0.5 * cirrhosis +
                     0.4 * heart_failure * acearb + 0.6 * nsaid * cohort$treatment)
  }
  base_rate <- h0 * exp(lp_out)

  ##--------------------------------------------------------------------------
  ## 5. Treatment switching --------------------------------------------------
  # Generate switch times BEFORE event times so that switching can modify
  # the post-switch hazard (step 6).
  # switch_on = FALSE -> no switching at all
  # switch_on = TRUE  -> informative switching modeled with its own hazard
  if (!switch_on) {
    cohort$switch_time   <- Inf
    cohort$switched      <- 0L
    if (return_potential_switching) {
      cohort$never_switcher <- 1L
    }
  } else {
    if (return_potential_switching) {
      # Draw potential switching times under both treatment assignments
      # using the same uniform random draws (save/restore RNG state).
      rng_state <- .Random.seed
      lambda_sw_a1 <- lambda_sw0 * exp(gamma_A * 1 + gamma_ckd * cohort$ckd)
      sw_a1 <- rexp(nrow(cohort), rate = lambda_sw_a1)
      .Random.seed <<- rng_state
      lambda_sw_a0 <- lambda_sw0 * exp(gamma_A * 0 + gamma_ckd * cohort$ckd)
      sw_a0 <- rexp(nrow(cohort), rate = lambda_sw_a0)
      # Use actual treatment to determine observed switch time
      cohort$switch_time <- ifelse(cohort$treatment == 1, sw_a1, sw_a0)
      cohort$never_switcher <- as.integer(sw_a1 > max_follow & sw_a0 > max_follow)
    } else {
      lambda_sw <- lambda_sw0 *
        exp(gamma_A * cohort$treatment + gamma_ckd * cohort$ckd)
      cohort$switch_time <- rexp(nrow(cohort), rate = lambda_sw)
    }
    cohort$switched <- 0L  # will be refined after event times
  }

  ##--------------------------------------------------------------------------
  ## 6. Event times (renal failure) ------------------------------------------
  # When switch_on = TRUE and policy = "treatment_policy", a patient who
  # switches at time S has:
  #   - hazard based on ORIGINAL treatment for t in [0, S)
  #   - hazard based on OPPOSITE treatment for t >= S
  # This means the treatment-policy and no-switch truths will differ.
  #
  # Implementation: generate event time in two pieces.
  #   Piece 1: [0, switch_time) under original treatment hazard
  #   Piece 2: [switch_time, Inf) under switched treatment hazard
  # If no event in piece 1, the residual survival is memoryless (exponential)
  # so we draw a new residual time under the post-switch rate.

  if (!np_hazard) {
    ## ── Proportional hazards ──
    rate_orig <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    rate_switch <- base_rate * ifelse(cohort$treatment == 1, 1, HR_early)

    # Draw event time under original treatment
    event_time_orig <- rexp(nrow(cohort), rate = rate_orig)

    if (switch_on) {
      # For patients who switch before their original event time,
      # redraw the residual event time under the post-switch hazard.
      did_switch_before_event <- cohort$switch_time < event_time_orig
      residual <- rexp(nrow(cohort), rate = rate_switch)
      event_time_new <- ifelse(
        did_switch_before_event,
        cohort$switch_time + residual,
        event_time_orig
      )
      cohort$event_time <- event_time_new
    } else {
      cohort$event_time <- event_time_orig
    }

  } else {
    ## ── Non-PH: change-point in HR at tau days ──
    rpexp_piece <- function(n, r1, r2, tau) {
      u  <- runif(n); p1 <- 1 - exp(-r1 * tau); t <- numeric(n)
      e  <- u <= p1
      t[e]  <- -log(1 - u[e]) / r1[e]
      t[!e] <- tau - log((1 - u[!e]) / (1 - p1[!e])) / r2[!e]
      t
    }

    # Rates under original treatment
    r1_orig <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    r2_orig <- base_rate * ifelse(cohort$treatment == 1, HR_late,  1)
    event_time_orig <- rpexp_piece(nrow(cohort), r1_orig, r2_orig, tau)

    if (switch_on) {
      # Rates under switched treatment (opposite arm)
      r1_switch <- base_rate * ifelse(cohort$treatment == 1, 1, HR_early)
      r2_switch <- base_rate * ifelse(cohort$treatment == 1, 1, HR_late)

      did_switch_before_event <- cohort$switch_time < event_time_orig

      # For switchers: draw residual time under post-switch rate.
      # Use the rate appropriate for the calendar time after switch.
      # Simplification: use the late-period rate if switch_time > tau,
      # early-period rate otherwise, then add a piecewise residual.
      # For tractability, we use a constant-rate residual at the
      # time-appropriate post-switch rate.
      sw_t <- cohort$switch_time
      post_switch_rate <- ifelse(sw_t <= tau, r1_switch, r2_switch)
      residual <- rexp(nrow(cohort), rate = post_switch_rate)

      event_time_new <- ifelse(
        did_switch_before_event,
        sw_t + residual,
        event_time_orig
      )
      cohort$event_time <- event_time_new
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
  ## 8. Observed follow-up & event indicator ---------------------------------
  # Under treatment_policy: switching does NOT censor. Subjects are
  #   followed until the event or administrative censoring regardless of
  #   whether they switched treatments.
  # Under no_switch: switch_on is already FALSE (set above), so
  #   switch_time = Inf and has no effect.
  # Under while_on_treatment: subjects are censored at switch_time.
  #   Only events occurring before switching contribute.

  if (policy == "while_on_treatment") {
    # Censor at switch time
    cohort$follow_time <- pmin(cohort$event_time, cohort$censor_admin,
                               cohort$switch_time)
  } else if (policy == "composite") {
    # Composite: switching is absorbed into the outcome.
    # follow_time = min(event_time, switch_time, censor_admin)
    # event = 1 if either AKI or switch came first (before admin censor)
    cohort$follow_time <- pmin(cohort$event_time, cohort$switch_time,
                               cohort$censor_admin)
  } else {
    # treatment_policy, no_switch, principal_stratum:
    # follow until event or admin censor (switching does not censor)
    cohort$follow_time <- pmin(cohort$event_time, cohort$censor_admin)
  }

  if (policy == "composite") {
    # Event = 1 if either the original event or switching occurred before admin censor
    cohort$event <- as.integer(
      cohort$follow_time < cohort$censor_admin |
        (cohort$event_time <= cohort$censor_admin) |
        (cohort$switch_time <= cohort$censor_admin)
    )
    # More precisely: event if the earliest of event/switch is <= admin censor
    cohort$event <- as.integer(
      pmin(cohort$event_time, cohort$switch_time) <= cohort$censor_admin
    )
  } else {
    cohort$event <- as.integer(cohort$event_time <= cohort$follow_time)
  }

  # Refine switched indicator: only count switching before follow-up end
  cohort$switched <- as.integer(
    cohort$switch_time < cohort$follow_time
  )

  # For principal stratum: flag never-switchers.
  # A subject is a "never-switcher" if they would not switch under EITHER
  # treatment assignment. Since we only observe one potential switching
  # outcome, we approximate by flagging subjects who did not switch under
  # their observed treatment. The true principal stratum requires both
  # potential switching indicators, which we store for truth calculations
  # when treat_override is used.
  cohort$would_switch <- as.integer(cohort$switch_time <= max_follow)

  ##--------------------------------------------------------------------------
  ## 9. Analysis dataset -----------------------------------------------------
  # Keep switch_time for downstream analyses (e.g., censor-at-switch Cox)
  ana <- cohort %>%
    dplyr::select(
      id, age, sex_male, race,
      ckd, diabetes, hypertension, cirrhosis, heart_failure,
      nsaid, acearb, statin,
      treatment, event_time, switch_time, follow_time, event, switched,
      any_of("never_switcher"),
      would_switch
    )

  ##--------------------------------------------------------------------------
  ## 10. Optional missingness / imputation -----------------------------------
  if (add_missing) {
    ana$race[sample(nrow(ana), 0.05 * nrow(ana))] <- NA
    ana$ckd[sample(nrow(ana), 0.10 * nrow(ana))]  <- NA
    if (impute) {
      imp_vars <- c("age", "race", "ckd", "cirrhosis", "diabetes", "hypertension")
      imp_in   <- ana %>% dplyr::select(dplyr::all_of(imp_vars)) %>%
        dplyr::mutate(across(c(race), as.factor))
      ana[, imp_vars] <- missForest::missForest(as.data.frame(imp_in),
                                                verbose = FALSE)$ximp
    }
  }

  return(ana)
}


# TODO(Joy): refine hepatitis B DGP parameters (h0, HR_early, HR_late)
#   to match realistic renal-failure incidence in tenofovir vs entecavir.

# ── Generate and cache example datasets ──────────────────────────────────────
# Only run when this script is executed directly, not when sourced by other files.

if (sys.nframe() == 0) {
  set.seed(1234)
  df <- generate_hep_data(
    h0         = 3e-4,
    np_hazard  = FALSE,
    dep_censor = FALSE,
    complexity = FALSE,
    policy     = "treatment_policy",
    seed       = 1234
  )

  # Check realistic values
  c(mean(df$follow_time),
    mean(df$event) * 100,
    mean(df$switched) * 100,
    mean(df$event) / mean(df$follow_time <= 85) * 100)

  write.csv(df, here::here("data/sim_hep_renal.csv"), row.names = FALSE)


  df_complex <- generate_hep_data(
    np_hazard  = TRUE,
    dep_censor = TRUE,
    complexity = TRUE,
    policy     = "treatment_policy",
    seed       = 1234
  )

  c(mean(df_complex$follow_time),
    mean(df_complex$event) * 100,
    mean(df_complex$switched) * 100,
    mean(df_complex$event) / mean(df_complex$follow_time <= 85) * 100)

  write.csv(df_complex, here::here("data/sim_hep_renal_complex.csv"), row.names = FALSE)
}
