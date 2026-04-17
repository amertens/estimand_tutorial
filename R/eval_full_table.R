# eval_full_table.R
# Print a full results table for the simulation study, with bias, RD coverage
# (bootstrap for Cox, analytic for LMTP), and HR coverage (Cox only).
#
# Prefers the production cache (sim_study_main.rds) and falls back to the
# 5-iteration test cache (sim_study_test.rds) if the production run has
# not yet completed.

library(dplyr)

sim_main <- "results/sim_study_main.rds"
sim_test <- "results/sim_study_test.rds"
sim_path <- if (file.exists(sim_main)) sim_main else sim_test
sim <- readRDS(sim_path)
truth <- readRDS("results/sim_results/ground_truth.rds")

n_iter <- length(unique(sim$iter))
cache_label <- if (identical(sim_path, sim_main)) "PRODUCTION" else "TEST"
cat(sprintf("\n=== FULL RESULTS TABLE [%s]: All Estimators x All Estimands (%d iterations) ===\n\n",
            cache_label, n_iter))
cat(sprintf("%-22s  %-32s  %4s  %7s  %8s  %8s  %8s  %8s\n",
    "Estimand", "Method", "n", "Truth", "Mean RD", "RD Bias", "RD Cov", "HR Cov"))
cat(paste(rep("-", 105), collapse = ""), "\n")

for (est in c("treatment_policy", "no_switch", "while_on_treatment", "composite", "principal_stratum")) {
  sub <- sim[sim$estimand == est, ]
  if (nrow(sub) == 0) next
  tr_rd <- truth[[est]]$true_rd
  tr_hr <- truth[[est]]$true_hr

  est_label <- c(treatment_policy = "Treatment-policy",
                 no_switch = "No-switch",
                 while_on_treatment = "While-on-treatment",
                 composite = "Composite",
                 principal_stratum = "Principal stratum")[[est]]

  first <- TRUE
  for (m in unique(sub$method)) {
    s <- sub[sub$method == m & sub$converged == TRUE, ]
    n <- nrow(s)
    if (n == 0) {
      cat(sprintf("%-22s  %-32s  %4d  %7s  %8s  %8s  %8s  %8s\n",
          ifelse(first, est_label, ""), m, 0, "", "FAILED", "", "", ""))
      first <- FALSE
      next
    }
    mean_rd <- mean(s$risk_diff, na.rm = TRUE)
    bias <- mean_rd - tr_rd

    # RD coverage (ci_low/ci_high are RD CIs for all methods now)
    if (!all(is.na(s$ci_low))) {
      rd_cov <- mean(s$ci_low <= tr_rd & s$ci_high >= tr_rd, na.rm = TRUE)
      rd_cov_str <- sprintf("%.0f%%", rd_cov * 100)
    } else {
      rd_cov_str <- "—"
    }

    # HR coverage (hr_ci_low/hr_ci_high if present)
    if ("hr_ci_low" %in% names(s) && !all(is.na(s$hr_ci_low))) {
      hr_cov <- mean(s$hr_ci_low <= tr_hr & s$hr_ci_high >= tr_hr, na.rm = TRUE)
      hr_cov_str <- sprintf("%.0f%%", hr_cov * 100)
    } else {
      hr_cov_str <- "—"
    }

    cat(sprintf("%-22s  %-32s  %4d  %+.4f  %+.4f   %+.4f   %6s    %6s\n",
        ifelse(first, est_label, ""), m, n, tr_rd, mean_rd, bias, rd_cov_str, hr_cov_str))
    first <- FALSE
  }
  cat("\n")
}
