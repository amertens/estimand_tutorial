library(tidyverse)


generate_hcv_data <- function(
    N               = 125000,
    ## treatment assignment
    p_sof           = 0.36,
    ## baseline hazard parameters
    h0              = 5e-5,
    HR_early        = 1.25,
    HR_late         = 0.70,
    tau             = 45,          # change‑point for non‑PH (days)
    max_follow      = 720,         # admin cut‑off (days)
    risk_window     = 30,          # at‑risk window after switch (days)
    ## behaviour toggles
    np_hazard       = TRUE,       # TRUE → HR varies over time (non‑PH)
    dep_censor      = TRUE,       # TRUE → admin censoring depends on risk
    complexity      = TRUE,       # TRUE → nonlinear PS & outcome models
    ## informative switching
    switch_on       = TRUE,       # TRUE → censor by switch w/ hazard model
    lambda_sw0      = 2.5e-5,        # baseline switch hazard
    gamma_A         = 0.60,        # log‑HR for treatment on switch
    gamma_ckd       = 0.40,        # log‑HR per CKD on switch
    ## misc
    censor_base     = 1/100,       # admin censoring rate (1/days)
    treat_override  = c("simulate","all_treated","all_control"),
    add_missing     = FALSE,
    impute          = FALSE,
    seed            = NULL)
{
  if (!is.null(seed)) set.seed(seed)
  treat_override <- match.arg(treat_override)
  if (impute) requireNamespace("missForest")
  requireNamespace("tidyverse")


  ##--------------------------------------------------------------------------
  ## 1. Demography -----------------------------------------------------------
  raw <- tibble::tibble(
    id           = seq_len(N),
    age          = pmax(rnorm(N, 48, 13), 18),
    sex_male     = rbinom(N, 1, 0.58),
    race         = sample(c("white","black","hispanic","asian","other"), N, TRUE,
                          prob = c(.48,.14,.06,.02,.30)),
    region       = sample(c("NE","MW","S","W"), N, TRUE,
                          prob = c(.20,.18,.37,.25)),
    enroll_days  = rpois(N, 420))

  ## 2. Clinical history & concomitant meds ----------------------------------
  add_bin <- function(p) rbinom(N, 1, p)
  raw <- raw %>%
    mutate(
      ckd     = add_bin(.08),  prior_aki = add_bin(.05),
      heart_failure = add_bin(.07), sepsis = add_bin(.03),
      dehydration   = add_bin(.06), obstruction = add_bin(.04),
      cirrhosis      = add_bin(.18), portal_htn = add_bin(.04),
      esld   = add_bin(.02), hiv   = add_bin(.04),
      diabetes = add_bin(.20), hypertension = add_bin(.45),
      bmi       = rnorm(N, 28, 5), overweight_obese = add_bin(.20),
      smoking   = add_bin(.40), alcohol = add_bin(.18),
      substance_abuse = add_bin(.25), cancer = add_bin(.08),
      chemo      = add_bin(.01),
      nsaid      = add_bin(.25), acearb = add_bin(.30), diuretic = add_bin(.22),
      aminoglycoside = add_bin(.05), contrast    = add_bin(.08),
      statin     = add_bin(.15), aspirin = add_bin(.10),
      beta_blocker = add_bin(.14), ccb    = add_bin(.16), art = add_bin(.05),
      prior_sof    = add_bin(.05), prior_nonsof = add_bin(.05))

  ## 3. Baseline exclusions ---------------------------------------------------
  cohort <- raw %>%
    filter(enroll_days >= 365, age >= 18, prior_aki == 0,
           !(prior_sof == 1 | prior_nonsof == 1))

  ##--------------------------------------------------------------------------
  ## 4. Treatment assignment --------------------------------------------------
  if (treat_override == "simulate") {

    ## baseline PS linear predictor
    lp0 <- with(cohort,
                0.015*age + 0.30*cirrhosis + 0.25*ckd + 0.15*hiv + 0.10*diabetes -
                  0.10*cancer + rnorm(nrow(cohort), 0, 0.6))

    if (complexity) {
      lp0 <- with(cohort, lp0 +
                    0.02*(bmi^2)/100 - 0.3*sin(0.1*bmi) + 0.5*(age/50)^3 +
                    1.5*ckd*cancer + 0.8*hiv*log1p(age))
    }
    alpha0  <- qlogis(p_sof) - mean(lp0)
    p_trt   <- plogis(alpha0 + lp0) |> pmin(0.95) |> pmax(0.05)
    cohort$treatment <- rbinom(nrow(cohort), 1, p_trt)
  } else {
    cohort$treatment <- ifelse(treat_override == "all_treated", 1L, 0L)
  }

  ##--------------------------------------------------------------------------
  ## 5. Individual baseline hazard -------------------------------------------
  if (!complexity) {
    lp_out <- with(cohort,
                   -2.8 + 0.03*age + 0.7*ckd + 0.5*cirrhosis +
                     0.3*heart_failure + 0.25*nsaid + 0.20*contrast)
  } else {
    lp_out <- with(cohort,
                   -2.8 + 0.03*age + 0.0005*age^2 + 0.7*ckd + 0.5*cirrhosis +
                     0.02*(bmi^2)/100 - 0.3*sin(0.1*bmi) + 0.4*heart_failure*acearb +
                     0.6*nsaid*treatment + 0.3*contrast*log1p(age))
  }
  base_rate <- h0 * exp(lp_out)

  ##--------------------------------------------------------------------------
  ## 6. Event times (AKI) -----------------------------------------------------
  if (!np_hazard) {
    ## proportional hazards
    rate <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    cohort$event_time <- rexp(nrow(cohort), rate = rate)
  } else {
    ## non‑PH: change‑point in HR at tau days
    rpexp_piece <- function(n, r1, r2, tau) {
      u  <- runif(n); p1 <- 1 - exp(-r1*tau); t <- numeric(n)
      e  <- u <= p1
      t[e]  <- -log(1 - u[e]) / r1[e]
      t[!e] <- tau - log((1 - u[!e])/(1 - p1[!e])) / r2[!e]
      t
    }
    r1 <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    r2 <- base_rate * ifelse(cohort$treatment == 1, HR_late,  1)
    cohort$event_time <- rpexp_piece(nrow(cohort), r1, r2, tau)
  }

  ##--------------------------------------------------------------------------
  ## 7. Administrative censoring --------------------------------------------
  if (!dep_censor) {
    censor_admin <- rexp(nrow(cohort), rate = censor_base)
  } else {
    c_rate <- censor_base * exp(0.04*lp_out + 0.03*cohort$treatment)
    censor_admin <- rexp(nrow(cohort), rate = c_rate)
  }
  cohort$censor_admin <- pmin(censor_admin, max_follow)

  ##--------------------------------------------------------------------------
  ## 8. Treatment‑switch censoring -------------------------------------------
  if (!switch_on) {
    cohort$tx_days <- ifelse(cohort$treatment == 1,
                             rpois(nrow(cohort), 84), rpois(nrow(cohort), 70))
    cohort$switch  <- rbinom(nrow(cohort), 1, 0.03)
    cohort$censor_switch <- ifelse(cohort$switch == 1,
                                   cohort$tx_days + risk_window,
                                   max_follow)
  } else {
    lambda_sw <- lambda_sw0 *
      exp(gamma_A * cohort$treatment + gamma_ckd * cohort$ckd)
    Sw_lat <- rexp(nrow(cohort), rate = lambda_sw)

    cohort$tx_days <- Sw_lat                    # store raw switch time
    cohort$switch  <- as.integer(Sw_lat < max_follow)
    cohort$censor_switch <- pmin(Sw_lat + risk_window, max_follow)
  }

  ##--------------------------------------------------------------------------
  ## 9. Observed follow‑up & event indicator ---------------------------------
  cohort$follow_time <- pmin(cohort$event_time,
                             cohort$censor_admin,
                             cohort$censor_switch)
  cohort$event <- as.integer(cohort$event_time <= cohort$follow_time)

  ##--------------------------------------------------------------------------
  ## 10. Analysis dataset -----------------------------------------------------
  ana <- cohort %>%
    dplyr::select(-enroll_days, -prior_aki, -prior_sof, -prior_nonsof,
                  -tx_days, -event_time, -censor_admin, -censor_switch)

  ##--------------------------------------------------------------------------
  ## 11. Optional missingness / imputation -----------------------------------
  if (add_missing) {
    ana$region[sample(nrow(ana), 0.05*nrow(ana))] <- NA
    ana$ckd[sample(nrow(ana), 0.10*nrow(ana))]   <- NA
    if (impute) {
      imp_vars <- c("age","race","region","ckd","cirrhosis","hiv",
                    "diabetes","hypertension","bmi")
      imp_in   <- ana %>% dplyr::select(dplyr::all_of(imp_vars)) %>%
        dplyr::mutate(across(c(race, region), as.factor))
      ana[, imp_vars] <- missForest::missForest(as.data.frame(imp_in),
                                                verbose = FALSE)$ximp
    }
  }

  return(ana)
}



set.seed(1234)
df <- generate_hcv_data(    h0              = 3e-4,
                            np_hazard     = FALSE,
                            dep_censor    = FALSE,
                            complexity    = FALSE,
                            seed=1234)

#check realistic values
c(mean(df$follow_time) ,
  mean(df$event)*100,
  mean(df$switch)*100,
  mean(df$event) / mean(df$follow_time <= 85)*100)  # crude CI within follow-up


write.csv(df,here::here("data/sim_hcv_aki.csv"))


df_complex <- generate_hcv_data(    np_hazard     = TRUE,
                            dep_censor    = TRUE,
                            complexity    = TRUE,
                            seed=1234)

#check realistic values
c(mean(df_complex$follow_time) ,
  mean(df_complex$event)*100,
  mean(df_complex$switch)*100,
  mean(df_complex$event) / mean(df$follow_time <= 85)*100)  # crude CI within follow-up


write.csv(df_complex,here::here("data/sim_hcv_aki_complex.csv"))
