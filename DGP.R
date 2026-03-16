library(tidyverse)


generate_hcv_data <- function(
    N               = 125000,
    ## treatment assignment
    p_sof           = 0.36,
    ## baseline hazard parameters
    h0              = 5e-5,
    HR_early        = 1.25,
    HR_late         = 0.70,
    tau             = 45,          # change-point for non-PH (days)
    max_follow      = 720,         # admin cut-off (days)
    risk_window     = 30,          # at-risk window after switch (days)
    ## behaviour toggles
    np_hazard       = TRUE,       # TRUE -> HR varies over time (non-PH)
    dep_censor      = TRUE,       # TRUE -> censoring depends on risk; FALSE -> no censoring (only max_follow)
    complexity      = TRUE,       # TRUE -> nonlinear PS & outcome models
    ## informative switching
    switch_on       = TRUE,       # TRUE -> informative switching; FALSE -> no switching at all
    lambda_sw0      = 2.5e-5,     # baseline switch hazard
    gamma_A         = 0.60,       # log-HR for treatment on switch
    gamma_ckd       = 0.40,       # log-HR per CKD on switch
    ## misc
    censor_base     = 1/100,      # admin censoring rate (1/days), only used when dep_censor=TRUE
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
  ## 1. Demographics ---------------------------------------------------------
  raw <- tibble::tibble(
    id           = seq_len(N),
    age          = pmax(rnorm(N, 48, 13), 18),
    sex_male     = rbinom(N, 1, 0.58),
    race         = sample(c("white","black","hispanic","asian","other"), N, TRUE,
                          prob = c(.48,.14,.06,.02,.30)))

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
                0.015*age + 0.30*cirrhosis + 0.25*ckd + 0.10*diabetes +
                  0.05*heart_failure + rnorm(nrow(cohort), 0, 0.6))

    if (complexity) {
      lp0 <- with(cohort, lp0 +
                    0.5*(age/50)^3 + 1.5*ckd*cirrhosis)
    }
    alpha0  <- qlogis(p_sof) - mean(lp0)
    p_trt   <- plogis(alpha0 + lp0) |> pmin(0.95) |> pmax(0.05)
    cohort$treatment <- rbinom(nrow(cohort), 1, p_trt)
  } else {
    cohort$treatment <- ifelse(treat_override == "all_treated", 1L, 0L)
  }

  ##--------------------------------------------------------------------------
  ## 4. Individual baseline hazard -------------------------------------------
  if (!complexity) {
    lp_out <- with(cohort,
                   -2.8 + 0.03*age + 0.7*ckd + 0.5*cirrhosis +
                     0.3*heart_failure + 0.25*nsaid)
  } else {
    lp_out <- with(cohort,
                   -2.8 + 0.03*age + 0.0005*age^2 + 0.7*ckd + 0.5*cirrhosis +
                     0.4*heart_failure*acearb + 0.6*nsaid*cohort$treatment)
  }
  base_rate <- h0 * exp(lp_out)

  ##--------------------------------------------------------------------------
  ## 5. Event times (AKI) ----------------------------------------------------
  if (!np_hazard) {
    ## proportional hazards
    rate <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    cohort$event_time <- rexp(nrow(cohort), rate = rate)
  } else {
    ## non-PH: change-point in HR at tau days
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
  ## 6. Administrative censoring ---------------------------------------------
  ## dep_censor = FALSE  -> no random censoring, only admin max_follow
  ## dep_censor = TRUE   -> censoring depends on risk profile
  if (!dep_censor) {
    censor_admin <- rep(max_follow, nrow(cohort))
  } else {
    c_rate <- censor_base * exp(0.04*lp_out + 0.03*cohort$treatment)
    censor_admin <- rexp(nrow(cohort), rate = c_rate)
  }
  cohort$censor_admin <- pmin(censor_admin, max_follow)

  ##--------------------------------------------------------------------------
  ## 7. Treatment-switch censoring -------------------------------------------
  ## switch_on = FALSE -> no switching at all (switch=0, censor_switch=max_follow)
  ## switch_on = TRUE  -> informative switching modeled with its own hazard
  if (!switch_on) {
    cohort$tx_days       <- max_follow
    cohort$switch        <- 0L
    cohort$censor_switch <- max_follow
  } else {
    lambda_sw <- lambda_sw0 *
      exp(gamma_A * cohort$treatment + gamma_ckd * cohort$ckd)
    Sw_lat <- rexp(nrow(cohort), rate = lambda_sw)

    cohort$tx_days <- Sw_lat                    # store raw switch time
    cohort$switch  <- as.integer(Sw_lat < max_follow)
    cohort$censor_switch <- pmin(Sw_lat + risk_window, max_follow)
  }

  ##--------------------------------------------------------------------------
  ## 8. Observed follow-up & event indicator ---------------------------------

  cohort$follow_time <- pmin(cohort$event_time,
                             cohort$censor_admin,
                             cohort$censor_switch)
  cohort$event <- as.integer(cohort$event_time <= cohort$follow_time)

  ##--------------------------------------------------------------------------
  ## 9. Analysis dataset -----------------------------------------------------
  ana <- cohort %>%
    dplyr::select(-tx_days, -event_time, -censor_admin, -censor_switch)

  ##--------------------------------------------------------------------------
  ## 10. Optional missingness / imputation -----------------------------------
  if (add_missing) {
    ana$race[sample(nrow(ana), 0.05*nrow(ana))] <- NA
    ana$ckd[sample(nrow(ana), 0.10*nrow(ana))]  <- NA
    if (impute) {
      imp_vars <- c("age","race","ckd","cirrhosis","diabetes","hypertension")
      imp_in   <- ana %>% dplyr::select(dplyr::all_of(imp_vars)) %>%
        dplyr::mutate(across(c(race), as.factor))
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
