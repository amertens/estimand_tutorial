# diagnose_no_switch.R
#
# Per-iteration diagnostics for the no-switch estimand:
#   - predicted switching probabilities at each biweekly bin
#     (pooled-logistic working model fit on each simulated dataset)
#   - effective sample size (ESS) of the no-switching IP weights
#
# Outputs:
#   results/sim_diagnostics/no_switch_probs.rds (long: iter x bin x trt x p)
#   results/sim_diagnostics/no_switch_ess.rds   (iter x bin x ESS)
#   results/sim_diagnostics/fig_switch_probs.png
#   results/sim_diagnostics/fig_ess.png
#
# Usage: Rscript R/diagnose_no_switch.R

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2)
})
source(here("DGP.R"))

tau       <- 180
bin_width <- 14
n_bins    <- ceiling(tau / bin_width)
N_ITER    <- 20
N_SAMPLE  <- 5000

dgp_args <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, dep_censor = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60,
  return_potential_switching = TRUE
)

out_dir <- here("results", "sim_diagnostics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------ helpers
build_long <- function(dat) {
  # One row per subject-bin while still on assigned treatment.
  # NA-safe: a non-switcher (switched == 0) is treated as "still on" forever.
  bins <- seq_len(n_bins)
  rep_cov <- function(x) rep(x, each = n_bins)
  n <- nrow(dat)
  long <- tibble(
    id        = rep_cov(seq_len(n)),
    bin       = rep(bins, n),
    treatment = rep_cov(dat$treatment),
    age       = rep_cov(dat$age),
    ckd       = rep_cov(dat$ckd),
    cirrhosis = rep_cov(dat$cirrhosis),
    diabetes  = rep_cov(dat$diabetes),
    heart_failure = rep_cov(dat$heart_failure),
    switched  = rep_cov(dat$switched),
    switch_time = rep_cov(dat$switch_time)
  )
  long <- long %>%
    mutate(
      bin_end = bin * bin_width,
      bin_start = (bin - 1) * bin_width,
      still_on = (switched == 0L) |
                 (!is.na(switch_time) & switch_time > bin_start),
      switch_now = ifelse(
        is.na(switch_time), 0L,
        as.integer(switched == 1L &
                   switch_time > bin_start &
                   switch_time <= bin_end)
      )
    ) %>%
    filter(still_on)
  long
}

ess_of <- function(w) {
  s1 <- sum(w, na.rm = TRUE)
  s2 <- sum(w^2, na.rm = TRUE)
  if (s2 == 0 || !is.finite(s2)) return(NA_real_)
  (s1 ^ 2) / s2
}

run_one <- function(i) {
  set.seed(2000 + i)
  dat <- do.call(generate_hep_data,
                 c(list(N = N_SAMPLE, seed = 2000 + i), dgp_args))
  long <- build_long(dat)
  if (nrow(long) == 0) return(NULL)
  fit <- tryCatch(
    glm(switch_now ~ factor(bin) + treatment + age + ckd + cirrhosis +
          diabetes + heart_failure,
        data = long, family = binomial()),
    error = function(e) NULL,
    warning = function(w) {
      # Refit without the warning-trigger by just suppressing it
      suppressWarnings(
        glm(switch_now ~ factor(bin) + treatment + age + ckd + cirrhosis +
              diabetes + heart_failure,
            data = long, family = binomial())
      )
    }
  )
  if (is.null(fit)) return(NULL)
  long$p_switch <- predict(fit, type = "response")

  probs <- long %>% group_by(bin, treatment) %>%
    summarise(
      iter   = i,
      n      = n(),
      p_med  = median(p_switch, na.rm = TRUE),
      p_p95  = quantile(p_switch, 0.95, na.rm = TRUE),
      p_max  = max(p_switch, na.rm = TRUE),
      p_lt01 = mean(p_switch < 0.01, na.rm = TRUE),
      p_gt99 = mean(p_switch > 0.99, na.rm = TRUE),
      .groups = "drop"
    )

  ess <- long %>% arrange(id, bin) %>%
    group_by(id) %>%
    mutate(stay_prob = pmax(1 - p_switch, 1e-6),
           w         = 1 / cumprod(stay_prob)) %>%
    ungroup() %>%
    group_by(bin) %>%
    summarise(
      iter      = i,
      n_at_risk = n(),
      ess       = ess_of(w),
      ess_frac  = ess_of(w) / n(),
      .groups   = "drop"
    )
  list(probs = probs, ess = ess)
}

# ------------------------------------------------------------------ run
message(sprintf("Running %d diagnostic iterations...", N_ITER))
res_list <- vector("list", N_ITER)
for (i in seq_len(N_ITER)) {
  res_list[[i]] <- tryCatch(run_one(i),
                            error = function(e) {
                              message("  iter ", i, " errored: ",
                                      conditionMessage(e))
                              NULL
                            })
  if ((i %% 10) == 0) message("  done ", i, "/", N_ITER)
}
res_list <- res_list[!vapply(res_list, is.null, logical(1))]
message("usable iterations: ", length(res_list))
if (length(res_list) == 0) stop("no usable iterations")

probs_df <- bind_rows(lapply(res_list, `[[`, "probs"))
ess_df   <- bind_rows(lapply(res_list, `[[`, "ess"))

saveRDS(probs_df, file.path(out_dir, "no_switch_probs.rds"))
saveRDS(ess_df,   file.path(out_dir, "no_switch_ess.rds"))

# ------------------------------------------------------------------ figures
p_probs <- ggplot(probs_df,
                  aes(x = factor(bin), y = p_p95, fill = factor(treatment))) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "tomato"),
                    labels = c("Comparator (ETV)", "Active (TDF)"),
                    name = "Treatment arm") +
  labs(x = "Biweekly bin", y = "95th percentile of P(switch | history)",
       title = "Predicted switching probabilities under the no-switch intervention",
       subtitle = sprintf("Across %d diagnostic iterations (N = %d each)",
                          length(res_list), N_SAMPLE)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig_switch_probs.png"), p_probs,
       width = 8, height = 4.5, dpi = 200, bg = "white")

p_ess <- ggplot(ess_df, aes(x = factor(bin), y = ess_frac)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "tomato") +
  labs(x = "Biweekly bin", y = "ESS / N at risk",
       title = "Effective sample size of no-switching IP weights",
       subtitle = "Dashed line: 50% efficiency threshold. Below this the no-switch counterfactual is poorly supported.") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "fig_ess.png"), p_ess,
       width = 8, height = 4.5, dpi = 200, bg = "white")

message("Wrote diagnostics to ", out_dir)
