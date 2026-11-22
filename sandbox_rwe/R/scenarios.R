# scenarios.R -- the base case and the one-at-a-time deviations (spec section 2).
#   c1_ph_rare_switch   boundary: proportional hazards and a quarter of the base
#                       switching and gap rates
#   d2a_ph              proportional hazards: one HR for all time, chosen so the
#                       adjusted Cox HR matches the base target (1.41)
#   d2b_heavy_ltfu      loss to follow-up at the commercial-claims rate (21.5% per year)
#   d2c_renal_switch    switching strongly driven by measured eGFR (odds 3 per 10%
#                       decline, 15 below 60, on TDF); switching volume refit to the
#                       base targets. Built to make informative censoring visible.
#   d2d_common_outcome  1-year entecavir risk 10% (TARGETS row 23)
# Supplement only (SUPP_SCENARIOS):
#   s1_tdf_susceptibility  d2c plus 10% of patients susceptible to TDF nephrotoxicity
#                       (a 5x TDF dip), switching volume refit. Built to violate
#                       principal ignorability: the always-on stratum excludes
#                       susceptible patients for a post-baseline, prognostic reason
#                       baseline covariates miss (Andrew, 2026-10-12).

SCENARIOS <- c("base", "c1_ph_rare_switch", "d2a_ph", "d2b_heavy_ltfu",
               "d2c_renal_switch", "d2d_common_outcome")
SUPP_SCENARIOS <- c("s1_tdf_susceptibility")

params_fingerprint <- function(p) digest::digest(p[order(names(p))], algo = "md5")

scenario_params <- function(name, base = readRDS(rwe_path("sandbox_rwe/results/params_calibrated.rds")),
                            tg = rwe_targets()) {
  p <- base
  if (name == "base") return(p)
  if (name == "d2b_heavy_ltfu") { p$ltfu_rate <- -log(1 - 0.215) / 365; return(p) }
  if (name == "d2a_ph") return(fit_hr(p, tg$adj_hr, constant = TRUE))
  if (name == "d2c_renal_switch") {
    p$sw_dec[2] <- log(3); p$sw_low60[2] <- log(15)
    return(fit_switching(p, tg$switch_12m, tg$switch_36m))
  }
  if (name == "d2d_common_outcome") return(fit_h0(p, 0.10))
  # Supplementary: loss to follow-up 5 times as fast once measured eGFR falls below 60, so it
  # depends on the post-baseline kidney function that also drives the event (informative
  # given baseline covariates). Built to make informative censoring visible (review A12).
  if (name == "d2e_informative_ltfu") { p$ltfu_low60 <- log(5); return(p) }
  # Supplementary, stronger (Andrew 2026-11-22, after d2e changed no analysis's bias): loss to
  # follow-up 30 times as fast from the first visit with a measured eGFR decline of 10% or more
  # from baseline, a trigger reached in the first year by about 36% of Drug A and 8% of Drug B
  # starters (pilot, scratchpad pilot_ltfu2.R). Run on the base case's seeds (run/75), so the
  # change from the base case is paired.
  if (name == "d2f_informative_ltfu_strong") { p$ltfu_low60 <- log(30); p$ltfu_dec_thr <- 0.10; return(p) }
  if (name == "s1_tdf_susceptibility") {
    p <- scenario_params("d2c_renal_switch", base, tg); p$sus_prev <- 0.1; p$sus_mult <- 5
    return(fit_switching(p, tg$switch_12m, tg$switch_36m))
  }
  if (name == "c1_ph_rare_switch") {
    p$gap_r1 <- p$gap_r1 / 4; p$gap_r2 <- p$gap_r2 / 4
    p <- fit_switching(p, tg$switch_12m / 4, tg$switch_36m / 4)
    return(fit_hr(p, tg$adj_hr, constant = TRUE))
  }
  stop("unknown scenario: ", name)
}

get_scenario_params <- function(name) {
  base <- readRDS(rwe_path("sandbox_rwe/results/params_calibrated.rds"))
  f <- rwe_path(sprintf("sandbox_rwe/results/params_%s.rds", name))
  if (file.exists(f)) {
    x <- readRDS(f)
    if (identical(x$base_fingerprint, params_fingerprint(base))) return(x$params)
  }
  p <- scenario_params(name, base)
  saveRDS(list(base_fingerprint = params_fingerprint(base), params = p), f)
  p
}
