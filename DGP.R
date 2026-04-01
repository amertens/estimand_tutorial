library(tidyverse)


generate_hep_data <- function(
    N               = 125000,
    ## treatment assignment
    p_trt_target    = 0.36,
    ## baseline hazard parameters
    h0              = 5e-5,
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
    lambda_sw0      = 2.5e-5,     # baseline switch hazard
    gamma_A         = 0.60,       # log-HR for treatment on switch
    gamma_ckd       = 0.40,       # log-HR per CKD on switch
    ## estimand policy
    policy          = c("treatment_policy", "no_switch"),
    ## misc
    censor_base     = 1/100,      # admin censoring rate (1/days), only used when dep_censor=TRUE
    treat_override  = c("simulate", "all_treated", "all_control"),
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
  ## 5. Event times (renal failure) ------------------------------------------
  if (!np_hazard) {
    ## proportional hazards
    rate <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    cohort$event_time <- rexp(nrow(cohort), rate = rate)
  } else {
    ## non-PH: change-point in HR at tau days
    rpexp_piece <- function(n, r1, r2, tau) {
      u  <- runif(n); p1 <- 1 - exp(-r1 * tau); t <- numeric(n)
      e  <- u <= p1
      t[e]  <- -log(1 - u[e]) / r1[e]
      t[!e] <- tau - log((1 - u[!e]) / (1 - p1[!e])) / r2[!e]
      t
    }
    r1 <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    r2 <- base_rate * ifelse(cohort$treatment == 1, HR_late,  1)
    cohort$event_time <- rpexp_piece(nrow(cohort), r1, r2, tau)
  }

  ##--------------------------------------------------------------------------
  ## 6. Administrative censoring ---------------------------------------------
  if (!dep_censor) {
    censor_admin <- rep(max_follow, nrow(cohort))
  } else {
    c_rate <- censor_base * exp(0.04 * lp_out + 0.03 * cohort$treatment)
    censor_admin <- rexp(nrow(cohort), rate = c_rate)
  }
  cohort$censor_admin <- pmin(censor_admin, max_follow)

  ##--------------------------------------------------------------------------
  ## 7. Treatment switching --------------------------------------------------
  # switch_on = FALSE -> no switching at all
  # switch_on = TRUE  -> informative switching modeled with its own hazard
  if (!switch_on) {
    cohort$switch_time   <- Inf
    cohort$switched      <- 0L
  } else {
    lambda_sw <- lambda_sw0 *
      exp(gamma_A * cohort$treatment + gamma_ckd * cohort$ckd)
    Sw_lat <- rexp(nrow(cohort), rate = lambda_sw)

    cohort$switch_time <- Sw_lat
    # A subject "switched" if their latent switch time falls before the
    # administrative censoring and before the event.
    cohort$switched <- as.integer(Sw_lat < cohort$censor_admin)
  }

  ##--------------------------------------------------------------------------
  ## 8. Observed follow-up & event indicator ---------------------------------
  # Under treatment_policy: switching does NOT censor. Subjects are
  # followed until the event or administrative censoring regardless of
  # whether they switched treatments.
  # Under no_switch: switch_on is already FALSE (set above), so
  # switch_time = Inf and has no effect.

  cohort$follow_time <- pmin(cohort$event_time, cohort$censor_admin)
  cohort$event       <- as.integer(cohort$event_time <= cohort$follow_time)

  # Refine switched indicator: only count switching before follow-up end
  cohort$switched <- as.integer(
    cohort$switch_time < cohort$follow_time & cohort$switched == 1L
  )

  ##--------------------------------------------------------------------------
  ## 9. Analysis dataset -----------------------------------------------------
  # Keep switch_time for downstream analyses (e.g., censor-at-switch Cox)
  ana <- cohort %>%
    dplyr::select(
      id, age, sex_male, race,
      ckd, diabetes, hypertension, cirrhosis, heart_failure,
      nsaid, acearb, statin,
      treatment, event_time, switch_time, follow_time, event, switched
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
