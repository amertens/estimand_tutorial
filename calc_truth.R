# calc_truth.R
# Compute Monte Carlo ground truth for each estimand from the DGP.
# Generates large counterfactual datasets under each policy to obtain
# true 180-day risks, risk differences, and risk ratios.

library(dplyr)
library(here)

source(here("DGP.R"))

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

  list(
    estimand    = "treatment_policy",
    true_risk_1 = risk_1,
    true_risk_0 = risk_0,
    true_rd     = risk_1 - risk_0,
    true_rr     = risk_1 / risk_0,
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

  list(
    estimand    = "no_switch",
    true_risk_1 = risk_1,
    true_risk_0 = risk_0,
    true_rd     = risk_1 - risk_0,
    true_rr     = risk_1 / risk_0,
    N           = N,
    tau         = tau
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
