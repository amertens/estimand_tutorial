# params.R -- parameters of the continuous-time RWE safety cohort (sandbox).
#
# Clinical frame: adults with chronic hepatitis B starting tenofovir disoproxil
# fumarate (TDF, A = 1) or entecavir (ETV, A = 0) in routine care. Time is in days.
# Two-element vectors are indexed by drug: [1] = ETV, [2] = TDF.
# Targets and sources: sandbox_rwe/TARGETS.md. Parameters marked "fitted" are
# overwritten by smoke/00_calibrate.R; the values here are starting points.

# Absolute project paths: testthat (and parallel workers) may run with another
# working directory.
# Override with the RWE_ROOT environment variable when the repository is cloned elsewhere.
RWE_ROOT <- Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial")
rwe_path <- function(...) file.path(RWE_ROOT, ...)

rwe_params <- function(...) {
  p <- list(
    ## baseline covariates
    age_mean = 50, age_sd = 12,
    p_male = 0.66,                              # TARGETS row 3
    dm_b0 = qlogis(0.13), dm_age = 0.04,        # fitted (row 5)
    htn_b0 = qlogis(0.155), htn_age = 0.06,     # fitted (row 4)
    p_cirr = 0.23, p_nsaid = 0.073,             # rows 6-7
    egfr_b0 = 95, egfr_age = -0.8, egfr_dm = -4, egfr_htn = -4,
    egfr_sd = 14,                               # egfr_b0, egfr_sd fitted (rows 8-9)
    ## natural eGFR slope per year (row 12)
    slope_mean = -1.0, slope_dm = -1.0, slope_htn = -0.5, slope_sd = 1.5,
    ## drug choice (fitted, rows 1-2, 4-5, 8-9; cirrhosis and NSAID set, rows 6-7)
    trt_b0 = 0, trt_egfr = 0.01, trt_low60 = -1.5, trt_age = -0.02,
    trt_dm = -0.1, trt_htn = -0.3, trt_cirr = -0.12, trt_nsaid = -0.3,
    ## TDF effect on eGFR (rows 10, 11, 13)
    dip_early = 4.5, dip_tau = 30, dip_late = 0.6, recover_tau = 270,
    meas_sd = 5,                                # assumption
    ## clinic visits with an eGFR test (row 14; mean gap fitted by drug)
    visit_first_min = 60, visit_first_max = 120,
    visit_gap_mean = c(165, 120), visit_gap_shape = 4, J = 20,
    ## switching at a visit (rows 15-16: intercept and per-year trend fitted;
    ## odds ratios are assumptions, TARGETS.md)
    sw_b0 = c(-4.5, -3.0), sw_t = c(0, 0.3),
    sw_dec = log(c(1.2, 1.6)),                  # per 10% measured decline from baseline
    sw_low60 = log(c(1.5, 4)),                  # measured eGFR below 60
    sw_age60 = log(c(1.2, 2)),                  # age over 60
    ## gaps in supply of 90 days or more (rows 17-18): year-1 and later rates
    gap_r1 = -log(1 - c(0.124, 0.113)) / 365,
    gap_r2 = (-log(1 - c(0.215, 0.161)) + log(1 - c(0.124, 0.113))) / 730,
    gap_age = -0.02,
    ## renal event (rows 21-22: h0 and HR_late fitted)
    h0 = 8e-5, b_age = 0.02, b_dm = 0.3, b_htn = 0.2, b_cirr = 0.4, b_nsaid = 0.3,
    b_egfr = 0.035, egfr_ref = 90,
    HR_early = 2.0, HR_late = 1.15, t_change = 90,
    ## loss to follow-up (row 19) and staggered entry (row 20)
    ltfu_rate = -log(0.75) / (3 * 365), ltfu_age = -0.015, ltfu_A = 0,
    entry_window = 1460, max_follow = 1470, piece = 15
  )
  modifyList(p, list(...))
}

rwe_targets <- function() list(
  p_tdf = 0.40, age_gap = 2.2,                       # ETV minus TDF
  htn = c(etv = 0.177, tdf = 0.123), dm = c(etv = 0.139, tdf = 0.121),
  egfr_median = c(etv = 90.0, tdf = 91.7), egfr_lt60 = c(etv = 0.081, tdf = 0.018),
  test_6m = c(etv = 0.60, tdf = 0.80),
  switch_12m = c(etv = 0.037, tdf = 0.165), switch_24m = c(etv = 0.079, tdf = 0.354),
  switch_36m = c(etv = 0.120, tdf = 0.529),
  gap_12m = c(etv = 0.124, tdf = 0.113), gap_36m = c(etv = 0.215, tdf = 0.161),
  ltfu_3y = 0.25, etv_risk_1y = 0.043, adj_hr = 1.41)
