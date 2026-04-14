# calc_truth.R
# Compute Monte Carlo ground truth under each estimand strategy.
#
# Each truth function generates large counterfactual datasets (all-treated
# and all-control) and computes 180-day marginal risks, RD, RR, and HR.
#
# The DGP always generates under treatment-policy (switching occurs and
# modifies the hazard). Estimand-specific outcomes are derived from the
# raw timing variables using derive_estimand().
#
# For the no-switch truth, switching is suppressed entirely (switch_on=FALSE).

library(dplyr)
library(survival)
library(here)

source(here("DGP.R"))

# ── Helper: compute marginal HR from two counterfactual datasets ─────────────
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

# ── Helper: compute risk at tau ──────────────────────────────────────────────
compute_risk <- function(df, tau) {
  mean(df$event == 1 & df$follow_time <= tau)
}

# ── Truth: Treatment-Policy ──────────────────────────────────────────────────
#' Generate all-treated and all-control with switching enabled (switch_on=TRUE).
#' Follow-up is not censored at switch (treatment-policy default).
truth_treatment_policy <- function(N = 500000, tau = 180, seed = 9999, ...) {
  common <- list(N = N, dep_censor = FALSE, switch_on = TRUE, seed = seed)
  common <- modifyList(common, list(...))

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  risk_1 <- compute_risk(df_a1, tau)
  risk_0 <- compute_risk(df_a0, tau)

  list(estimand = "treatment_policy",
       true_risk_1 = risk_1, true_risk_0 = risk_0,
       true_rd = risk_1 - risk_0, true_rr = risk_1 / risk_0,
       true_hr = compute_true_hr(df_a1, df_a0, tau),
       N = N, tau = tau)
}

# ── Truth: No-Switch ─────────────────────────────────────────────────────────
#' Generate all-treated and all-control with switch_on=FALSE.
#' No switching occurs; potential outcomes under sustained treatment.
truth_no_switch <- function(N = 500000, tau = 180, seed = 9999, ...) {
  common <- list(N = N, dep_censor = FALSE, switch_on = FALSE, seed = seed)
  common <- modifyList(common, list(...))

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  risk_1 <- compute_risk(df_a1, tau)
  risk_0 <- compute_risk(df_a0, tau)

  list(estimand = "no_switch",
       true_risk_1 = risk_1, true_risk_0 = risk_0,
       true_rd = risk_1 - risk_0, true_rr = risk_1 / risk_0,
       true_hr = compute_true_hr(df_a1, df_a0, tau),
       N = N, tau = tau)
}

# ── Truth: While-on-Treatment ────────────────────────────────────────────────
#' Compute true risks under the while-on-treatment (censor-at-switch) estimand.
#' Generate all-treated and all-control with switching enabled, then apply
#' derive_estimand("while_on_treatment") which censors follow-up at the
#' switch time. The truth is the 180-day cumulative incidence under this
#' censoring regime — NOT conditional on non-switching.
#'
#' Because dep_censor=FALSE and there is no administrative censoring,
#' the only censoring mechanism is switching. The KM estimator on this
#' data gives the correct WOT risk (no confounding under counterfactual
#' uniform treatment).
truth_while_on_treatment <- function(N = 500000, tau = 180, seed = 9999, ...) {
  common <- list(N = N, dep_censor = FALSE, switch_on = TRUE, seed = seed)
  common <- modifyList(common, list(...))

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  # Apply censor-at-switch: follow_time = min(event, switch, admin_censor)
  df_a1_wot <- derive_estimand(df_a1, "while_on_treatment")
  df_a0_wot <- derive_estimand(df_a0, "while_on_treatment")

  # KM-based risk at tau under switch-censoring
  # Under counterfactual uniform treatment with no dependent admin censoring,
  # the KM estimator is unbiased for the WOT risk.
  sf1 <- survfit(Surv(pmin(follow_time, tau),
                       as.integer(event == 1 & follow_time <= tau)) ~ 1,
                 data = df_a1_wot)
  sf0 <- survfit(Surv(pmin(follow_time, tau),
                       as.integer(event == 1 & follow_time <= tau)) ~ 1,
                 data = df_a0_wot)

  s1 <- summary(sf1, times = tau, extend = TRUE)
  s0 <- summary(sf0, times = tau, extend = TRUE)
  risk_1 <- 1 - s1$surv
  risk_0 <- 1 - s0$surv

  true_hr <- compute_true_hr(df_a1_wot, df_a0_wot, tau)

  list(estimand = "while_on_treatment",
       true_risk_1 = risk_1, true_risk_0 = risk_0,
       true_rd = risk_1 - risk_0, true_rr = risk_1 / risk_0,
       true_hr = true_hr,
       pct_switched_a1 = mean(df_a1_wot$switched),
       pct_switched_a0 = mean(df_a0_wot$switched),
       N = N, tau = tau)
}

# ── Truth: Composite ─────────────────────────────────────────────────────────
#' Generate all-treated and all-control with switching, then derive
#' composite outcome (renal failure or switching, whichever first).
truth_composite <- function(N = 500000, tau = 180, seed = 9999, ...) {
  common <- list(N = N, dep_censor = FALSE, switch_on = TRUE, seed = seed)
  common <- modifyList(common, list(...))

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  # Apply composite outcome definition
  df_a1 <- derive_estimand(df_a1, "composite")
  df_a0 <- derive_estimand(df_a0, "composite")

  risk_1 <- compute_risk(df_a1, tau)
  risk_0 <- compute_risk(df_a0, tau)

  list(estimand = "composite",
       true_risk_1 = risk_1, true_risk_0 = risk_0,
       true_rd = risk_1 - risk_0, true_rr = risk_1 / risk_0,
       true_hr = compute_true_hr(df_a1, df_a0, tau),
       N = N, tau = tau)
}

# ── Truth: Principal Stratum ─────────────────────────────────────────────────
#' Generate all-treated and all-control with paired switching draws
#' to identify the true never-switcher stratum.
truth_principal_stratum <- function(N = 500000, tau = 180, seed = 9999, ...) {
  common <- list(N = N, dep_censor = FALSE, switch_on = TRUE,
                 return_potential_switching = TRUE, seed = seed)
  common <- modifyList(common, list(...))

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  never_sw <- df_a1$never_switcher == 1 & df_a0$never_switcher == 1

  risk_1 <- mean(df_a1$event_time[never_sw] <= tau)
  risk_0 <- mean(df_a0$event_time[never_sw] <= tau)
  true_hr <- compute_true_hr(df_a1[never_sw, ], df_a0[never_sw, ], tau)

  list(estimand = "principal_stratum",
       true_risk_1 = risk_1, true_risk_0 = risk_0,
       true_rd = risk_1 - risk_0, true_rr = risk_1 / risk_0,
       true_hr = true_hr,
       n_never_switch = sum(never_sw),
       pct_never_switch = mean(never_sw),
       N = N, tau = tau)
}


# ── Compute and save ─────────────────────────────────────────────────────────
if (sys.nframe() == 0) {
  set.seed(12345)
  dgp_args <- list(np_hazard = TRUE, complexity = TRUE)

  message("Computing treatment-policy truth...")
  truth_tp <- do.call(truth_treatment_policy,
    c(list(N = 500000, tau = 180, seed = 9999), dgp_args))
  message("Computing no-switch truth...")
  truth_ns <- do.call(truth_no_switch,
    c(list(N = 500000, tau = 180, seed = 9999), dgp_args))
  message("Computing while-on-treatment truth...")
  truth_wot <- do.call(truth_while_on_treatment,
    c(list(N = 500000, tau = 180, seed = 9999), dgp_args))
  message("Computing composite truth...")
  truth_comp <- do.call(truth_composite,
    c(list(N = 500000, tau = 180, seed = 9999), dgp_args))
  message("Computing principal stratum truth...")
  truth_ps <- do.call(truth_principal_stratum,
    c(list(N = 500000, tau = 180, seed = 9999), dgp_args))

  truth_all <- list(
    treatment_policy   = truth_tp,
    no_switch          = truth_ns,
    while_on_treatment = truth_wot,
    composite          = truth_comp,
    principal_stratum  = truth_ps
  )

  for (nm in names(truth_all)) {
    t <- truth_all[[nm]]
    message(sprintf("  %-20s RD=%+.6f  RR=%.4f  HR=%.4f",
                    nm, t$true_rd, t$true_rr, t$true_hr))
  }

  dir.create(here("results", "sim_results"), recursive = TRUE,
             showWarnings = FALSE)
  saveRDS(truth_all, file = here("results", "sim_results",
                                  "ground_truth.rds"))
  message("Saved to results/sim_results/ground_truth.rds")
}
