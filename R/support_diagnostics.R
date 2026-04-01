# support_diagnostics.R
# Functions for diagnosing data support / positivity for target interventions.

# TODO(Joy): review whether the support diagnostics chosen are the most
#   interpretable ones for the handbook audience.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ── Propensity score diagnostics ─────────────────────────────────────────────
#' Compute propensity scores and return diagnostics.
#' @param dat data.frame with treatment and baseline covariates.
#' @param covars character vector of covariate names for PS model.
#' @return data.frame with ps column added, plus diagnostic summary.
compute_ps_diagnostics <- function(dat,
                                   covars = c("age", "ckd", "cirrhosis",
                                              "diabetes", "hiv",
                                              "hypertension", "bmi")) {
  covars <- intersect(covars, names(dat))
  fml <- as.formula(paste("treatment ~", paste(covars, collapse = " + ")))
  ps_mod <- glm(fml, data = dat, family = binomial)
  dat$ps <- predict(ps_mod, type = "response")

  diag <- list(
    ps_summary    = summary(dat$ps),
    near_0        = mean(dat$ps < 0.05),
    near_1        = mean(dat$ps > 0.95),
    near_0_or_1   = mean(dat$ps < 0.05 | dat$ps > 0.95),
    ess_treated   = sum(dat$treatment == 1 & dat$ps > 0.05 & dat$ps < 0.95),
    ess_control   = sum(dat$treatment == 0 & dat$ps > 0.05 & dat$ps < 0.95)
  )

  list(data = dat, diagnostics = diag)
}


# ── Plot PS overlap ──────────────────────────────────────────────────────────
#' Plot propensity score distributions by treatment group.
#' @param dat data.frame with ps and treatment columns.
#' @param title character; plot title.
#' @return ggplot object.
plot_ps_overlap <- function(dat, title = "Propensity Score Distribution") {
  ggplot(dat, aes(x = ps, fill = factor(treatment,
                                         labels = c("Comparator", "Active")))) +
    geom_histogram(alpha = 0.5, bins = 60, position = "identity") +
    labs(x = "P(Treatment = Active | Baseline)", y = "Count",
         fill = "Treatment", title = title) +
    theme_minimal(base_size = 12) +
    scale_fill_manual(values = c("steelblue", "tomato"))
}


# ── Switching probability diagnostics ────────────────────────────────────────
#' Compute and summarize switching probabilities.
#' @param dat data.frame with switch, treatment, and covariates.
#' @param covars character vector for switch model.
#' @return list with model, predicted probabilities, and diagnostics.
compute_switch_diagnostics <- function(dat,
                                       covars = c("age", "ckd", "treatment",
                                                   "cirrhosis", "diabetes")) {
  covars <- intersect(covars, names(dat))
  # Use the 'switched' column from generate_hep_data()
  sw_col <- if ("switched" %in% names(dat)) "switched" else "switch"
  fml <- as.formula(paste(sw_col, "~", paste(covars, collapse = " + ")))
  sw_mod <- glm(fml, data = dat, family = binomial)
  dat$p_switch <- predict(sw_mod, type = "response")

  diag <- list(
    p_switch_summary = summary(dat$p_switch),
    # Proportion with very high predicted switching (near-deterministic)
    near_certain_switch = mean(dat$p_switch > 0.8),
    near_certain_stay   = mean(dat$p_switch < 0.02),
    # By CKD status
    mean_pswitch_ckd    = mean(dat$p_switch[dat$ckd == 1]),
    mean_pswitch_nockd  = mean(dat$p_switch[dat$ckd == 0])
  )

  list(data = dat, model = sw_mod, diagnostics = diag)
}


# ── Natural intervention compliance ─────────────────────────────────────────
#' Compute the proportion of subjects who naturally follow a target
#' intervention at each time point (proxy for support).
#' @param dat data.frame with treatment and switch columns.
#' @return tibble with natural compliance summaries.
compute_natural_compliance <- function(dat) {
  # For static "always treat" intervention: proportion who start on treatment
  # AND do not switch away before tau.
  sw_col <- if ("switched" %in% names(dat)) "switched" else "switch"
  treat_start <- mean(dat$treatment == 1)
  treat_no_switch <- mean(dat$treatment == 1 & dat[[sw_col]] == 0)
  ctrl_start <- mean(dat$treatment == 0)
  ctrl_no_switch <- mean(dat$treatment == 0 & dat[[sw_col]] == 0)

 tibble::tibble(
    intervention = c("Always active (no switch)",
                     "Always comparator (no switch)",
                     "Start active (any switching)",
                     "Start comparator (any switching)"),
    prop_natural = c(treat_no_switch, ctrl_no_switch,
                     treat_start, ctrl_start)
  )
}


# ── Combined support diagnostic table ────────────────────────────────────────
#' Build a summary table of support diagnostics for a given scenario.
#' @param dat data.frame from DGP.
#' @param scenario_label character label for the scenario.
#' @return tibble summarizing key diagnostics.
summarize_support <- function(dat, scenario_label = "Default") {
  ps_diag <- compute_ps_diagnostics(dat)
  sw_diag <- compute_switch_diagnostics(dat)
  nc      <- compute_natural_compliance(dat)

  tibble::tibble(
    scenario      = scenario_label,
    n             = nrow(dat),
    switch_rate   = mean(dat[[if ("switched" %in% names(dat)) "switched" else "switch"]]),
    ps_near_0     = ps_diag$diagnostics$near_0,
    ps_near_1     = ps_diag$diagnostics$near_1,
    ps_violation  = ps_diag$diagnostics$near_0_or_1,
    p_switch_ckd  = sw_diag$diagnostics$mean_pswitch_ckd,
    p_switch_nockd = sw_diag$diagnostics$mean_pswitch_nockd,
    natural_sof_nosw   = nc$prop_natural[1],
    natural_nonsof_nosw = nc$prop_natural[2]
  )
}


# ── Multi-scenario support plot ──────────────────────────────────────────────
#' Create a faceted PS overlap plot across multiple scenarios.
#' @param scenario_list named list of data.frames, each from DGP.
#' @return ggplot object.
plot_support_scenarios <- function(scenario_list) {
  combined <- bind_rows(lapply(names(scenario_list), function(nm) {
    d <- scenario_list[[nm]]
    covars <- intersect(c("age", "ckd", "cirrhosis", "diabetes", "hiv",
                          "hypertension", "bmi"), names(d))
    fml <- as.formula(paste("treatment ~", paste(covars, collapse = " + ")))
    ps_mod <- glm(fml, data = d, family = binomial)
    d$ps <- predict(ps_mod, type = "response")
    d$scenario <- nm
    d %>% select(id, treatment, ps, scenario)
  }))

  ggplot(combined, aes(x = ps, fill = factor(treatment,
                                              labels = c("Comparator", "Active")))) +
    geom_histogram(alpha = 0.5, bins = 50, position = "identity") +
    facet_wrap(~ scenario, scales = "free_y", ncol = 1) +
    labs(x = "Estimated P(Active | Baseline)", y = "Count",
         fill = "Treatment",
         title = "Data Support Diagnostic: Propensity Score Overlap") +
    theme_minimal(base_size = 11) +
    scale_fill_manual(values = c("steelblue", "tomato")) +
    theme(legend.position = "bottom")
}
