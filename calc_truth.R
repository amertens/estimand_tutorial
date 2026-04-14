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
#' Compute true 180-day risk under baseline treatment with follow-up
#' truncated at switching.
#'
#' Under counterfactual uniform treatment (treat_override) with
#' dep_censor=FALSE, the only censoring mechanism is switching.
#' Within each arm, all subjects share the same treatment, so the
#' switching hazard depends only on CKD (not on treatment-by-covariate
#' interactions that differ across subjects' treatments). The truth
#' is computed directly as the Monte Carlo proportion of subjects
#' whose event occurs before both their switch time and tau.
truth_while_on_treatment <- function(N = 500000, tau = 180, seed = 9999, ...) {
  common <- list(N = N, dep_censor = FALSE, switch_on = TRUE, seed = seed)
  common <- modifyList(common, list(...))

  df_a1 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_treated")))
  df_a0 <- do.call(generate_hep_data,
                    modifyList(common, list(treat_override = "all_control")))

  # Direct Monte Carlo risk under censor-at-switch:
  # An event is observed if event_time <= min(switch_time, tau).
  # With N=500,000 this is effectively exact.
  risk_1 <- mean(df_a1$event_time <= pmin(df_a1$switch_time, tau))
  risk_0 <- mean(df_a0$event_time <= pmin(df_a0$switch_time, tau))

  # For the HR, apply derive_estimand to get the truncated follow-up
  df_a1_wot <- derive_estimand(df_a1, "while_on_treatment")
  df_a0_wot <- derive_estimand(df_a0, "while_on_treatment")
  true_hr <- compute_true_hr(df_a1_wot, df_a0_wot, tau)

  list(estimand = "while_on_treatment",
       true_risk_1 = risk_1, true_risk_0 = risk_0,
       true_rd = risk_1 - risk_0, true_rr = risk_1 / risk_0,
       true_hr = true_hr,
       pct_switched_a1 = mean(df_a1$switch_time <= tau),
       pct_switched_a0 = mean(df_a0$switch_time <= tau),
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
