# tests/test_helpers.R
# Sanity tests for helpers.R and DGP.R. Run via:
#   Rscript tests/test_helpers.R
#
# These are not formal testthat tests; they are lightweight checks that
# exercise the most error-prone code paths and exit non-zero on failure.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(survival)
})
source(here("DGP.R"))
source(here("R", "helpers.R"))

ok <- function(cond, msg) {
  if (!isTRUE(cond)) {
    cat("FAIL: ", msg, "\n", sep = "")
    quit(status = 1)
  }
  cat("PASS: ", msg, "\n", sep = "")
}

# ── 1. DGP produces required columns ─────────────────────────────────────────
dat <- generate_hep_data(N = 1000, switch_on = TRUE, seed = 1,
                         lambda_sw0 = 1e-3, gamma_A = 0.8, gamma_ckd = 0.6)
required <- c("id", "treatment", "follow_time", "event",
              "switch_time", "switched", "event_time",
              "age", "ckd", "cirrhosis")
ok(all(required %in% names(dat)),
   "generate_hep_data returns required columns")
ok(nrow(dat) == 1000, "DGP returns N rows")
ok(all(dat$treatment %in% c(0L, 1L)), "treatment is 0/1")
ok(all(dat$follow_time > 0 & dat$follow_time <= 180 |
       dat$follow_time > 180),
   "follow_time is positive")

# ── 2. prepare_lmtp_data: default WOT (competing risk) shape ────────────────
prep_wot <- prepare_lmtp_data(dat, tau = 180, bin_width = 14,
                              competing_risk_at_switch = TRUE)
ok(prep_wot$n_bins == 13, "WOT prep produces 13 bins for tau=180, bw=14")
ok(is.null(prep_wot$A_cols),
   "WOT prep with default time_varying_trt = FALSE has A_cols = NULL")
ok(length(prep_wot$Y_cols) == 13 && length(prep_wot$C_cols) == 13,
   "Y_cols and C_cols have length n_bins")

# Switchers without a prior event: Y must be 0 in every bin
sw_no_event <- which(dat$switched == 1 &
                     dat$event_time > dat$switch_time &
                     dat$switch_time <= 180)
if (length(sw_no_event)) {
  Y_mat <- as.matrix(prep_wot$data[sw_no_event, prep_wot$Y_cols])
  ok(all(Y_mat == 0, na.rm = TRUE),
     "competing-risk WOT: switchers without prior event have Y=0 in all bins")
}

# ── 3. prepare_lmtp_data: time-varying treatment shape ──────────────────────
prep_ns <- prepare_lmtp_data(dat, tau = 180, bin_width = 14,
                             time_varying_trt = TRUE)
ok(!is.null(prep_ns$A_cols) && length(prep_ns$A_cols) == 13,
   "time-varying prep produces A_cols of length n_bins")

# ── 4. cox_gcomp_rd matches survfit-based g-comp on a small example ─────────
fit <- coxph(Surv(time_use, event_use) ~ treatment + age + ckd,
             data = dat %>% mutate(time_use = pmin(follow_time, 180),
                                   event_use = as.integer(event == 1 &
                                                          follow_time <= 180)))
gcomp_fast <- cox_gcomp_rd(fit, dat, tau = 180)

# Reference: survfit.coxph standardization
d0 <- dat; d0$treatment <- 0L
d1 <- dat; d1$treatment <- 1L
sf0 <- survfit(fit, newdata = d0)
sf1 <- survfit(fit, newdata = d1)
s0 <- as.numeric(summary(sf0, times = 180, extend = TRUE)$surv)
s1 <- as.numeric(summary(sf1, times = 180, extend = TRUE)$surv)
gcomp_ref <- mean(1 - s1) - mean(1 - s0)

ok(abs(gcomp_fast$risk_diff - gcomp_ref) < 1e-3,
   sprintf("cox_gcomp_rd matches survfit-based g-comp (fast=%.4f, ref=%.4f)",
           gcomp_fast$risk_diff, gcomp_ref))

# ── 5. M1 truth round-trip on a tiny cohort ─────────────────────────────────
small <- generate_hep_data(N = 100, switch_on = TRUE, seed = 2,
                           lambda_sw0 = 1e-3, gamma_A = 0.8, gamma_ckd = 0.6)
m1_truth <- mean(small$event_time <= pmin(small$switch_time, 180))
ok(is.finite(m1_truth) && m1_truth >= 0 && m1_truth <= 1,
   sprintf("M1 WOT truth is a valid probability (got %.4f)", m1_truth))

cat("\nAll tests passed.\n")
