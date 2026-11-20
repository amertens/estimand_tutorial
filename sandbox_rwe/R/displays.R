# displays.R -- main and supplementary displays from summaries and truths.
# Colors follow the paper: SDR blue #0072B2, Cox vermillion #D55E00.
STRAT_LABEL <- c(treatment_policy = "Treatment policy", sustained_use = "Sustained use",
                 while_on_treatment = "While on treatment", composite = "Composite",
                 principal_stratum = "Principal stratum")
SCEN_LABEL <- c(base = "Base case", c1_ph_rare_switch = "Proportional hazards, rare switching",
                d2a_ph = "Proportional hazards", d2b_heavy_ltfu = "Heavier loss to follow-up",
                d2c_renal_switch = "Renal-driven switching", d2d_common_outcome = "Common outcome",
                s1_tdf_susceptibility = "Susceptibility to the nephrotoxicity of Drug A (supplement)")

make_table1 <- function(perf, truth) {
  b <- perf[perf$scenario == "base", ]
  st <- truth$strategies
  out <- do.call(rbind, lapply(seq_len(nrow(PAIRS)), function(i) {
    s <- PAIRS$strategy[i]; m <- PAIRS$method[i]
    cx1 <- b[b$family == "cox" & b$method == m & b$window %in% "1y" & b$strategy == s, ]
    cxa <- b[b$family == "cox" & b$method == m & b$window %in% "all" & b$strategy == s, ]
    sm <- if (s == "principal_stratum") "SDR, true stratum" else "SDR"
    sd <- b[b$family == "sdr" & b$strategy == s & b$method == sm, ]
    if (!nrow(cx1) || !nrow(sd)) return(NULL)
    tr <- st[st$strategy == s, ]
    data.frame(strategy = STRAT_LABEL[[s]], cox_analysis = m, risk_tdf = tr$risk_tdf, risk_etv = tr$risk_etv,
               true_rd = tr$rd, true_rr = tr$rr, true_hr = tr$hr_marginal,
               hr_1y = cx1$mean_est, cover_1y = cx1$cover, hr_all = cxa$mean_est, cover_all = cxa$cover,
               sdr_rd_bias = sd$bias, sdr_rd_mcse = sd$mcse_bias, sdr_rd_cover = sd$cover,
               sdr_logrr_bias = sd$bias_logrr, sdr_rr_cover = sd$cover_rr)
  }))
  rownames(out) <- NULL
  out
}

make_figure1 <- function(perf) {
  suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))   # "built under" notes
  perf <- perf[perf$scenario %in% SCENARIOS, ]                # main scenarios; the supplement has its own table
  keep <- perf$strategy %in% c("treatment_policy", "sustained_use", "while_on_treatment") &
    ((perf$family == "cox" & perf$window %in% "1y" & perf$method %in% c("itt", "astreated")) |
       (perf$family == "sdr" & perf$method == "SDR"))
  f <- perf[keep & !(perf$strategy == "sustained_use" & perf$method == "itt") &
              !(perf$strategy == "treatment_policy" & perf$method == "astreated"), ]
  # Cox: mean(log HR) - log(true RR). SDR: mean(log RR) - log(true RR) (summarise_runs'
  # bias_logrr; pre-pilot review I7), so both sit on the log risk-ratio scale.
  f$log_bias <- f$bias_logrr
  f$mc_lo <- f$log_bias - 1.96 * f$mcse_logrr; f$mc_hi <- f$log_bias + 1.96 * f$mcse_logrr
  f$cover <- f$cover_rr
  f$who <- ifelse(f$family == "cox", "Cox analysis paired with the strategy", "SDR")
  sc <- intersect(SCENARIOS, f$scenario)
  f$scenario <- factor(SCEN_LABEL[f$scenario], levels = rev(SCEN_LABEL[sc]))
  # the text layer is grouped by estimator so each coverage label dodges with its point
  ggplot(f, aes(x = log_bias, y = scenario, colour = who)) +
    geom_vline(xintercept = 0, colour = "grey55") +
    geom_linerange(aes(xmin = mc_lo, xmax = mc_hi), position = position_dodge(width = 0.5),
                   linewidth = 0.6, show.legend = FALSE) +              # 95% Monte Carlo intervals
    geom_point(position = position_dodge(width = 0.5), size = 2.4) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * cover), group = who), position = position_dodge(width = 0.5),
              hjust = -0.3, size = 2.9, colour = "grey25", show.legend = FALSE) +
    facet_wrap(~ factor(STRAT_LABEL[strategy], levels = STRAT_LABEL[1:3]), nrow = 1) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    scale_colour_manual(values = c("Cox analysis paired with the strategy" = "#D55E00", SDR = "#0072B2"), name = NULL) +
    scale_y_discrete(labels = function(x) sapply(x, function(s) paste(strwrap(s, 20), collapse = "
"))) +
    labs(x = "Bias, log risk ratio", y = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank())
}

# Principal-stratum approaches in each scenario that ran them (the base case and the
# renal-driven-switching deviation); truths is a list of truths named by scenario.
make_table2 <- function(perf, truths) {
  sc <- intersect(c(SCENARIOS, SUPP_SCENARIOS), names(truths))
  do.call(rbind, lapply(sc, function(s) {
    b <- perf[perf$scenario == s & perf$strategy == "principal_stratum" & perf$family != "cox", ]
    if (!nrow(b)) return(NULL)
    col <- function(nm) if (nm %in% names(b)) b[[nm]] else rep(NA_real_, nrow(b))
    data.frame(scenario = SCEN_LABEL[[s]], approach = b$method, true_rd = b$true_rd, mean_rd = b$mean_est,
               bias = b$bias, mcse = b$mcse_bias, cover = b$cover, logrr_bias = col("bias_logrr"),
               rr_cover = col("cover_rr"), est_share = col("mean_share"),
               true_share = truths[[s]]$strategies$stratum_share[1])
  }))
}

# Each target recomputed on a fresh draw, defined as in calibration.R: switching
# and gaps are net of each other (cumulative probability on the drug's schedule,
# no censoring), the ETV 1-year risk is weighted to ETV starters, and the hazard
# ratio is the adjusted intention-to-treat Cox HR over all follow-up.
make_calibration_table <- function(p, tg, n = 2e5, seed = 11) {
  sim <- sim_rwe(n, p, seed = seed); d <- sim$data
  arm <- c("ETV", "TDF")
  by_arm <- function(x) as.numeric(tapply(x, d$A, mean))
  row <- function(q, target, simulated, a = arm) data.frame(quantity = q, arm = a, target = as.numeric(target),
                                                             simulated = as.numeric(simulated))
  tested <- sapply(0:1, function(a) { v <- sim$vis[d$A == a, , drop = FALSE]
    mean(sapply(0:3, function(w) rowSums(v >= w * 182.5 & v < (w + 1) * 182.5) > 0)) })
  set.seed(seed + 1); b <- sim_baseline(n, p); vd <- sim_visit_draws(n, p)
  eps <- matrix(rnorm(n * p$J, 0, p$meas_sd), n); u <- matrix(runif(n * p$J), n)
  ps <- plogis(trt_lp(b, p)); ug <- runif(n)
  sw <- sapply(0:1, function(a) { w <- if (a == 1) ps else 1 - ps
    t <- leave_time(b, a, visits_for(vd, a, p), eps, u, rep(Inf, n), p)$time
    c(.wmean(t <= 365, w), .wmean(t <= 1095, w)) })
  gp <- sapply(0:1, function(a) { w <- if (a == 1) ps else 1 - ps
    t <- gap_time(ug, b$age, a, p); c(.wmean(t <= 365, w), .wmean(t <= 1095, w)) })
  etv <- sim_rwe(n, p, arm = "etv", censoring = FALSE, seed = seed + 2)$data
  dd <- with_stages(d); dd <- cbind(dd[, c("id", "A", COX_X)], time = dd$follow, status = dd$event)
  rbind(
    row("share starting TDF", tg$p_tdf, mean(d$A), "all"),
    row("age, ETV minus TDF (years)", tg$age_gap, -diff(as.numeric(tapply(d$age, d$A, mean))), "all"),
    row("hypertension", tg$htn, by_arm(d$htn)),
    row("diabetes", tg$dm, by_arm(d$dm)),
    row("median baseline eGFR", tg$egfr_median, as.numeric(tapply(d$egfr0, d$A, median))),
    row("baseline eGFR < 60", tg$egfr_lt60, by_arm(d$egfr0 < 60)),
    row("eGFR test in a 6-month window", tg$test_6m, tested),
    row("switched by 12 months", tg$switch_12m, sw[1, ]),
    row("switched by 36 months", tg$switch_36m, sw[2, ]),
    row("gap in supply by 12 months", tg$gap_12m, gp[1, ]),
    row("gap in supply by 36 months", tg$gap_36m, gp[2, ]),
    row("lost to follow-up by 3 years", tg$ltfu_3y, mean(d$ltfu <= 1095), "all"),
    row("renal event by 1 year, ETV starters", tg$etv_risk_1y, .wmean(etv$T_nat <= 365, 1 - etv$ps), "ETV"),
    row("adjusted Cox HR, all follow-up", tg$adj_hr, suppressWarnings(cox_hr(dd))[["hr"]], "all"))
}

write_numbers <- function(perf, truths, path) {
  con <- file(path, "w"); on.exit(close(con))
  writeLines("# Numbers for the manuscript (generated by run/50_displays.R; do not edit by hand)\n", con)
  for (s in names(truths)) {
    writeLines(sprintf("## Truths: %s\n", s), con)
    writeLines(knitr::kable(truths[[s]]$strategies, digits = 4), con)
    writeLines(knitr::kable(truths[[s]]$reduced, digits = 4), con)
  }
  writeLines("\n## Performance\n", con)
  writeLines(knitr::kable(perf, digits = 4), con)
}
