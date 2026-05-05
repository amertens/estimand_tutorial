suppressPackageStartupMessages({ library(dplyr); library(here) })
source(here("R", "run_simulations.R"))
res <- readRDS(here("results", "sim_study_main.rds"))
truth_all <- readRDS(here("results", "sim_results", "ground_truth.rds"))

for (e in c("treatment_policy", "no_switch", "while_on_treatment",
            "composite", "principal_stratum")) {
  s <- summarize_simulation(
    res[res$estimand == e, ],
    truth_rd = truth_all[[e]]$true_rd,
    truth_rr = truth_all[[e]]$true_rr,
    truth_hr = truth_all[[e]]$true_hr
  )
  cat(sprintf("\n=== %s | RD=%+.4f  RR=%.3f  HR=%.3f ===\n",
              e, truth_all[[e]]$true_rd,
              truth_all[[e]]$true_rr, truth_all[[e]]$true_hr))
  cat(sprintf("  %-32s %3s %7s %7s %7s %7s %7s %7s\n",
              "Method", "n", "mn_HR", "HR_cov",
              "mn_RD", "RD_cov", "mn_RR", "RR_cov"))
  cat(paste(rep("-", 90), collapse = ""), "\n")
  for (i in seq_len(nrow(s))) {
    f <- function(x, fmt = "%+.3f") if (is.na(x)) "-" else sprintf(fmt, x)
    p <- function(x) if (is.na(x)) "-" else sprintf("%4.0f%%", 100 * x)
    cat(sprintf("  %-32s %3d %7s %7s %7s %7s %7s %7s\n",
                s$method[i], s$n_converged[i],
                f(s$mean_hr[i]), p(s$coverage_hr[i]),
                f(s$mean_rd[i]), p(s$coverage_rd[i]),
                f(s$mean_rr[i], "%.3f"), p(s$coverage_rr[i])))
  }
}
