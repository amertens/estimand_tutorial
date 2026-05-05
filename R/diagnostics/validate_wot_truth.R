# validate_wot_truth.R
# Cross-check the WOT truth computation against three independent estimators
# on a single huge counterfactual cohort. If all four methods agree, the
# +0.035 LMTP bias is not a truth-side artifact.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(survival)
})
source(here("DGP.R"))

tau <- 180
N <- 1e6  # huge to crush MC noise

dgp_args <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, complexity = FALSE, dep_censor = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60
)

cat(sprintf("Generating counterfactual cohorts (N = %d each)...\n", N))
df_a1 <- do.call(generate_hep_data,
                 c(list(N = N, switch_on = TRUE, seed = 9999,
                        treat_override = "all_treated"), dgp_args))
df_a0 <- do.call(generate_hep_data,
                 c(list(N = N, switch_on = TRUE, seed = 9999,
                        treat_override = "all_control"), dgp_args))

# ── Method 1: current calc_truth.R formula ──
risk1_m1 <- mean(df_a1$event_time <= pmin(df_a1$switch_time, tau))
risk0_m1 <- mean(df_a0$event_time <= pmin(df_a0$switch_time, tau))
rd_m1 <- risk1_m1 - risk0_m1

# ── Method 2: independent reformulation ──
# Construct WOT-derived (time, status) and compute risk = mean(time <= tau & status == 1)
wot_dat <- function(df) {
  t <- pmin(df$event_time, df$switch_time, tau)
  s <- as.integer(df$event_time <= df$switch_time & df$event_time <= tau)
  data.frame(t = t, s = s)
}
w1 <- wot_dat(df_a1)
w0 <- wot_dat(df_a0)
risk1_m2 <- mean(w1$t <= tau & w1$s == 1)
risk0_m2 <- mean(w0$t <= tau & w0$s == 1)
rd_m2 <- risk1_m2 - risk0_m2

# ── Method 3: Kaplan-Meier on WOT-derived data ──
# Treat switch as censoring; estimate 1 - S(tau).
km_risk_at_tau <- function(w, t = tau) {
  sf <- survfit(Surv(t, s) ~ 1, data = w)
  s_tau <- summary(sf, times = t, extend = TRUE)$surv
  as.numeric(1 - s_tau)
}
risk1_m3 <- km_risk_at_tau(w1)
risk0_m3 <- km_risk_at_tau(w0)
rd_m3 <- risk1_m3 - risk0_m3

# ── Method 4: KM on LMTP's discretization grid (bin_width = 14) ──
# Mimics how LMTP would aggregate within-bin events/switches.
discretize <- function(w, bin_width = 14, tau = 180) {
  edges <- seq(0, tau, by = bin_width)
  if (tail(edges, 1) < tau) edges <- c(edges, tau)
  # Snap each time to its bin's right edge, except ties at exact edges
  bin_idx <- findInterval(w$t, edges, rightmost.closed = TRUE)
  bin_idx <- pmax(1, pmin(length(edges) - 1, bin_idx))
  w$t_disc <- edges[bin_idx + 1]
  w
}
w1d <- discretize(w1)
w0d <- discretize(w0)
risk1_m4 <- as.numeric(1 - summary(survfit(Surv(t_disc, s) ~ 1, data = w1d),
                                   times = tau, extend = TRUE)$surv)
risk0_m4 <- as.numeric(1 - summary(survfit(Surv(t_disc, s) ~ 1, data = w0d),
                                   times = tau, extend = TRUE)$surv)
rd_m4 <- risk1_m4 - risk0_m4

# ── Sanity: switching rates and event rates ──
cat(sprintf("\nSanity checks:\n"))
cat(sprintf("  Switch rate by tau:  A=1: %.3f  A=0: %.3f\n",
            mean(df_a1$switch_time <= tau),
            mean(df_a0$switch_time <= tau)))
cat(sprintf("  Event rate by tau:   A=1: %.3f  A=0: %.3f\n",
            mean(df_a1$event_time <= tau),
            mean(df_a0$event_time <= tau)))

# ── Report ──
cat(sprintf("\n=== WOT TRUTH CROSS-CHECK (N = %d per arm) ===\n", N))
cat(sprintf("  %-50s  %s   %s   %s\n", "Method", "risk_1", "risk_0", "RD"))
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("  %-50s  %.4f   %.4f   %+.4f\n",
            "M1: current formula (calc_truth.R)",
            risk1_m1, risk0_m1, rd_m1))
cat(sprintf("  %-50s  %.4f   %.4f   %+.4f\n",
            "M2: reformulation (time + status)",
            risk1_m2, risk0_m2, rd_m2))
cat(sprintf("  %-50s  %.4f   %.4f   %+.4f\n",
            "M3: KM on continuous WOT data",
            risk1_m3, risk0_m3, rd_m3))
cat(sprintf("  %-50s  %.4f   %.4f   %+.4f\n",
            "M4: KM on bin_width=14 discretization",
            risk1_m4, risk0_m4, rd_m4))
cat(sprintf("\nCached truth: RD = %+.4f\n",
            readRDS(here("results", "sim_results", "ground_truth.rds"))$while_on_treatment$true_rd))
