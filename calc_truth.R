# calc_truth.R
# Compute Monte Carlo ground truth under each DGP policy.
#
# Each truth function generates large counterfactual datasets (all-treated
# and all-control) under a specific DGP policy and computes 180-day marginal
# risks, risk differences, risk ratios, and marginal HRs.
#
# Important notes on interpretation:
# - treatment_policy: switching occurs and modifies the hazard, but does not
#   censor. Truth reflects outcomes under natural switching behaviour.
# - no_switch: switching is suppressed entirely. Truth reflects potential
#   outcomes under sustained treatment.
# - while_on_treatment: uses policy="no_switch" as a shortcut. Under uniform
#   counterfactual treatment (all-treated or all-control), nobody switches,
#   so the numerical truth is identical to no_switch. This does NOT exercise
#   the censoring-at-switch mechanism; the truth should not be interpreted as
#   a full validation of while-on-treatment estimators.
# - composite: event = AKI or switching. Under counterfactual uniform treatment,
#   switching still occurs (driven by treatment-dependent hazard).
# - principal_stratum: restricts to subjects who would not switch under either
#   treatment assignment. Requires generating both potential switching outcomes.

library(dplyr)
library(survival)
library(here)

source(here("DGP.R"))

# ── Helper: compute marginal HR from two counterfactual datasets ─────────────
#' Stack all-treated and all-control datasets, fit unadjusted Cox,
#' and extract the marginal HR. This is the "true" HR under each estimand
#' for comparing against Cox estimates from observed data.
#' @param df_a1 data.frame; all-treated counterfactual.
#' @param df_a0 data.frame; all-control counterfactual.
#' @param tau numeric; follow-up horizon.
#' @return numeric; marginal hazard ratio.
compute_true_hr <- function(df_a1, df_a0, tau) {
  df_a1$trt <- 1L
  df_a0$trt <- 0L
  combined <- bind_rows(
    df_a1 %>% mutate(time_use = pmin(follow_time, tau),
                     event_use = as.integer(event == 1 & follow_time <= tau)),
    df_a0 %>% mutate(time_use = pmin(follow_time, tau),
                     event_use = as.integer(event == 1 & follow_time <= tau))
  )
  fit <- coxph(Surv(time_use, event_use) ~ trt, data = combined)
  as.numeric(exp(coef(fit)["trt"]))
}

# ── Truth: Treatment-Policy Estimand ─────────────────────────────────────────
#' Compute true risks under the treatment-policy estimand.
#' Under treatment-policy, switching occurs naturally but is not censored.
#' We generate all-treated and all-control datasets with policy="treatment_policy",
#' no dependent censoring, and compute 180-day cumulative incidence.
#'
#' @param N integer; Monte Carlo sample size (large for precision).
#' @param tau integer; risk window in days.
#' @param seed integer; random seed.
#' @param ... additional DGP arguments (np_hazard, complexity, etc.).
#' @return list with true_risk_1, true_risk_0, true_rd, true_rr.
truth_treatment_policy <- function(N = 500000, tau = 180, seed = 9999, ...) {
  dots <- list(...)

  common <- list(
    N          = N,
    dep_censor = FALSE,
    policy     = "treatment_policy",
    seed       = seed
  )
  common <- modifyList(common, dots)

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  risk_1 <- mean(df_a1$event == 1 & df_a1$follow_time <= tau)
  risk_0 <- mean(df_a0$event == 1 & df_a0$follow_time <= tau)
  true_hr <- compute_true_hr(df_a1, df_a0, tau)

  list(
    estimand    = "treatment_policy",
    true_risk_1 = risk_1,
    true_risk_0 = risk_0,
    true_rd     = risk_1 - risk_0,
    true_rr     = risk_1 / risk_0,
    true_hr     = true_hr,
    N           = N,
    tau         = tau
  )
}


# ── Truth: Hypothetical No-Switch Estimand ───────────────────────────────────
#' Compute true risks under the hypothetical no-switch estimand.
#' Under no-switch, no switching occurs at all. This gives the counterfactual
#' outcome distribution in a world where nobody switches treatment.
#'
#' @param N integer; Monte Carlo sample size.
#' @param tau integer; risk window in days.
#' @param seed integer; random seed.
#' @param ... additional DGP arguments.
#' @return list with true_risk_1, true_risk_0, true_rd, true_rr.
truth_no_switch <- function(N = 500000, tau = 180, seed = 9999, ...) {
  dots <- list(...)

  common <- list(
    N          = N,
    dep_censor = FALSE,
    policy     = "no_switch",
    seed       = seed
  )
  common <- modifyList(common, dots)

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  risk_1 <- mean(df_a1$event == 1 & df_a1$follow_time <= tau)
  risk_0 <- mean(df_a0$event == 1 & df_a0$follow_time <= tau)
  true_hr <- compute_true_hr(df_a1, df_a0, tau)

  list(
    estimand    = "no_switch",
    true_risk_1 = risk_1,
    true_risk_0 = risk_0,
    true_rd     = risk_1 - risk_0,
    true_rr     = risk_1 / risk_0,
    true_hr     = true_hr,
    N           = N,
    tau         = tau
  )
}


# ── Truth: While-on-Treatment Estimand ───────────────────────────────────────
#' Compute true risks under the while-on-treatment estimand.
#' Under while-on-treatment, subjects are censored at their switch time.
#' Only events occurring before switching contribute. We generate large
#' counterfactual datasets with policy="while_on_treatment" and no dependent
#' censoring, then compute 180-day cumulative incidence among non-switchers.
#'
#' Note: for the truth calculation under treat_override, all subjects receive
#' the same treatment so switching patterns differ from the observed data.
#' We use policy="no_switch" with switch_on=FALSE to get the clean
#' potential outcome under sustained treatment (equivalent to while-on-treatment
#' when everyone stays on their assigned treatment).
#'
#' @param N integer; Monte Carlo sample size.
#' @param tau integer; risk window in days.
#' @param seed integer; random seed.
#' @param ... additional DGP arguments.
#' @return list with true_risk_1, true_risk_0, true_rd, true_rr.
truth_while_on_treatment <- function(N = 500000, tau = 180, seed = 9999, ...) {
  dots <- list(...)

  # While-on-treatment under all-treated or all-control means nobody switches
  # away (everyone is on the same drug). This is equivalent to the no-switch
  # truth for counterfactual populations, but the estimand interpretation
  # differs: WOT conditions on not switching rather than intervening to
  # prevent switching. The numerical truth is identical when treat_override
  # forces uniform treatment.
  common <- list(
    N          = N,
    dep_censor = FALSE,
    policy     = "no_switch",
    seed       = seed
  )
  common <- modifyList(common, dots)

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  risk_1 <- mean(df_a1$event == 1 & df_a1$follow_time <= tau)
  risk_0 <- mean(df_a0$event == 1 & df_a0$follow_time <= tau)
  true_hr <- compute_true_hr(df_a1, df_a0, tau)

  list(
    estimand    = "while_on_treatment",
    true_risk_1 = risk_1,
    true_risk_0 = risk_0,
    true_rd     = risk_1 - risk_0,
    true_rr     = risk_1 / risk_0,
    true_hr     = true_hr,
    N           = N,
    tau         = tau
  )
}


# ── Truth: Composite Estimand ────────────────────────────────────────────────
#' Compute true risks under the composite estimand.
#' Under composite, the outcome is AKI *or* switching (whichever comes first).
#' We generate all-treated and all-control datasets with policy="composite"
#' and compute 180-day cumulative incidence of the composite endpoint.
#'
#' @param N integer; Monte Carlo sample size.
#' @param tau integer; risk window in days.
#' @param seed integer; random seed.
#' @param ... additional DGP arguments.
#' @return list with true_risk_1, true_risk_0, true_rd, true_rr.
truth_composite <- function(N = 500000, tau = 180, seed = 9999, ...) {
  dots <- list(...)

  common <- list(
    N          = N,
    dep_censor = FALSE,
    policy     = "composite",
    seed       = seed
  )
  common <- modifyList(common, dots)

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  risk_1 <- mean(df_a1$event == 1 & df_a1$follow_time <= tau)
  risk_0 <- mean(df_a0$event == 1 & df_a0$follow_time <= tau)
  true_hr <- compute_true_hr(df_a1, df_a0, tau)

  list(
    estimand    = "composite",
    true_risk_1 = risk_1,
    true_risk_0 = risk_0,
    true_rd     = risk_1 - risk_0,
    true_rr     = risk_1 / risk_0,
    true_hr     = true_hr,
    N           = N,
    tau         = tau
  )
}


# ── Truth: Principal Stratum Estimand ────────────────────────────────────────
#' Compute true risks under the principal stratum estimand.
#' The principal stratum of "never-switchers" consists of subjects who
#' would not switch under either treatment assignment. Since we control
#' the DGP, we can identify this stratum by generating switching times
#' under both treatments for the same subjects and restricting to those
#' who would not switch under either.
#'
#' Implementation: generate a cohort with no treatment override, draw
#' switching times under both treatments using the same random seed for
#' baseline covariates, identify the never-switcher stratum, then compute
#' event risks for that subpopulation.
#'
#' @param N integer; Monte Carlo sample size.
#' @param tau integer; risk window in days.
#' @param seed integer; random seed.
#' @param ... additional DGP arguments.
#' @return list with true_risk_1, true_risk_0, true_rd, true_rr.
truth_principal_stratum <- function(N = 500000, tau = 180, seed = 9999, ...) {
  dots <- list(...)

  common <- list(
    N          = N,
    dep_censor = FALSE,
    policy     = "treatment_policy",
    seed       = seed
  )
  common <- modifyList(common, dots)

  # Generate under all-treated and all-control with the SAME seed
  # so baseline covariates are identical across the two worlds.
  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  # Never-switchers: would not switch under either treatment
  # would_switch is based on switch_time <= max_follow
  never_switcher <- df_a1$would_switch == 0 & df_a0$would_switch == 0

  risk_1 <- mean(df_a1$event[never_switcher] == 1 &
                   df_a1$follow_time[never_switcher] <= tau)
  risk_0 <- mean(df_a0$event[never_switcher] == 1 &
                   df_a0$follow_time[never_switcher] <= tau)
  true_hr <- compute_true_hr(df_a1[never_switcher, ], df_a0[never_switcher, ], tau)

  list(
    estimand      = "principal_stratum",
    true_risk_1   = risk_1,
    true_risk_0   = risk_0,
    true_rd       = risk_1 - risk_0,
    true_rr       = risk_1 / risk_0,
    true_hr       = true_hr,
    n_never_switch = sum(never_switcher),
    pct_never_switch = mean(never_switcher),
    N             = N,
    tau           = tau
  )
}


# ── Compute and save ─────────────────────────────────────────────────────────

if (sys.nframe() == 0) {
  # Only run when sourced as a script, not when sourced by other files
  set.seed(12345)

  message("Computing treatment-policy truth (N = 500,000)...")
  truth_tp <- truth_treatment_policy(
    N = 500000, tau = 180, seed = 9999,
    np_hazard = TRUE, complexity = TRUE,
    switch_on = TRUE  # switching happens but does not censor
  )

  message("Computing no-switch truth (N = 500,000)...")
  truth_ns <- truth_no_switch(
    N = 500000, tau = 180, seed = 9999,
    np_hazard = TRUE, complexity = TRUE
  )

  truth_all <- list(
    treatment_policy = truth_tp,
    no_switch        = truth_ns
  )

  # TODO(Joy): verify that these parameters match the hepatitis B
  #   renal failure incidence rates from the literature.

  message("\n=== Treatment-Policy Truth ===")
  message("  Risk (treated):  ", round(truth_tp$true_risk_1, 6))
  message("  Risk (control):  ", round(truth_tp$true_risk_0, 6))
  message("  Risk difference: ", round(truth_tp$true_rd, 6))
  message("  Risk ratio:      ", round(truth_tp$true_rr, 4))

  message("\n=== No-Switch Truth ===")
  message("  Risk (treated):  ", round(truth_ns$true_risk_1, 6))
  message("  Risk (control):  ", round(truth_ns$true_risk_0, 6))
  message("  Risk difference: ", round(truth_ns$true_rd, 6))
  message("  Risk ratio:      ", round(truth_ns$true_rr, 4))

  dir.create(here("results", "sim_results"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(truth_all, file = here("results", "sim_results", "ground_truth.rds"))
  message("\nSaved to results/sim_results/ground_truth.rds")
}
