suppressPackageStartupMessages({ library(dplyr); library(here) })
res <- readRDS(here("results", "sim_study_main.rds"))
truth_all <- readRDS(here("results", "sim_results", "ground_truth.rds"))
truth_rd <- sapply(truth_all, function(t) t$true_rd)

out <- res %>%
  group_by(estimand, method) %>%
  summarise(
    n       = sum(converged, na.rm = TRUE),
    mean_rd = mean(risk_diff, na.rm = TRUE),
    bias_rd = mean_rd - truth_rd[first(estimand)],
    cov_rd  = mean(ci_low <= truth_rd[first(estimand)] &
                   ci_high >= truth_rd[first(estimand)], na.rm = TRUE),
    .groups = "drop"
  )

for (e in c("treatment_policy", "no_switch", "while_on_treatment",
            "composite", "principal_stratum")) {
  cat(sprintf("\n=== %s (truth RD = %+.4f) ===\n", e, truth_rd[e]))
  sub <- out[out$estimand == e, ]
  cat(sprintf("  %-32s %4s %9s %9s %9s\n",
              "Method", "n", "mean_RD", "bias", "cov_RD"))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-32s %4d %+9.4f %+9.4f %8.0f%%\n",
                sub$method[i], sub$n[i],
                sub$mean_rd[i], sub$bias_rd[i],
                100 * sub$cov_rd[i]))
  }
}
