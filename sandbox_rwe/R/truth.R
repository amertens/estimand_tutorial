# truth.R -- Monte Carlo truths for the five strategies at the horizon (sandbox).
#
# Counterfactual cohorts under TDF and under entecavir share every random draw
# (same seed per chunk), with no loss to follow-up or administrative censoring.
#   treatment policy   P(T <= h) under the natural course of leaving
#   sustained use      P(T <= h) with leaving switched off
#   while on treatment P(T <= min(leave, h)): the event before leaving the drug
#   composite          P(min(T, leave) <= h)
#   principal stratum  P(T <= h) among patients who would stay on either drug
#                      through h (leave1 > h and leave0 > h)
# The marginal HR to the horizon for each strategy is the coefficient of a
# treatment-only Cox model on the stacked counterfactual cohorts (first chunk).

rwe_truth <- function(p = rwe_params(), N = 1e6, chunk = 2.5e5, h = 365, seed = 1) {
  suppressPackageStartupMessages(library(survival))
  n_chunks <- ceiling(N / chunk)
  acc <- NULL
  hr_rows <- NULL
  for (cidx in seq_len(n_chunks)) {
    s <- seed * 1000 + cidx
    d1 <- sim_rwe(chunk, p, arm = "tdf", censoring = FALSE, sustained = TRUE, seed = s)$data
    d0 <- sim_rwe(chunk, p, arm = "etv", censoring = FALSE, sustained = TRUE, seed = s)$data
    stopifnot(identical(d1$leave1, d0$leave1), identical(d1$age, d0$age))
    stratum <- d1$leave1 > h & d1$leave0 > h
    risks <- function(d) c(
      treatment_policy   = mean(d$T_nat <= h),
      sustained_use      = mean(d$T_sus <= h),
      while_on_treatment = mean(d$T_nat <= pmin(d$leave, h)),
      composite          = mean(pmin(d$T_nat, d$leave) <= h),
      principal_stratum  = mean(d$T_nat[stratum] <= h))
    acc <- rbind(acc, data.frame(chunk = cidx, arm = c(1, 0),
                                 rbind(risks(d1), risks(d0)),
                                 stratum_share = mean(stratum)))
    if (cidx == 1) {
      st <- rbind(d1, d0)
      cut <- function(t, e) list(time = pmin(t, h), status = as.integer(e & t <= h))
      mk <- list(
        treatment_policy   = cut(st$T_nat, TRUE),
        sustained_use      = cut(st$T_sus, TRUE),
        while_on_treatment = cut(pmin(st$T_nat, st$leave), st$T_nat <= st$leave),
        composite          = cut(pmin(st$T_nat, st$leave), TRUE))
      hr <- sapply(mk, function(z)
        unname(exp(coef(coxph(Surv(z$time, z$status) ~ st$A, ties = "breslow")))))
      ss <- rep(stratum, 2)
      z <- cut(st$T_nat[ss], TRUE)
      hr["principal_stratum"] <- exp(coef(coxph(Surv(z$time, z$status) ~ st$A[ss],
                                                ties = "breslow")))
      hr_rows <- hr
    }
  }
  strategies <- c("treatment_policy", "sustained_use", "while_on_treatment",
                  "composite", "principal_stratum")
  r1 <- colMeans(acc[acc$arm == 1, strategies]); r0 <- colMeans(acc[acc$arm == 0, strategies])
  data.frame(strategy = strategies, risk_tdf = r1, risk_etv = r0,
             rd = r1 - r0, rr = r1 / r0, hr_marginal = unname(hr_rows[strategies]),
             stratum_share = mean(acc$stratum_share), N_per_arm = n_chunks * chunk,
             row.names = NULL)
}

# Truths under stochastic reductions of leaving the drug: each would-be switch or
# gap happens with probability delta (delta = 1 is treatment policy, delta = 0 is
# sustained use). Same seeds as rwe_truth(), so the endpoints can be checked.
rwe_truth_reduced <- function(p, deltas, N = 5e5, chunk = 2.5e5, h = 365, seed = 1) {
  do.call(rbind, lapply(deltas, function(dl) {
    acc <- t(sapply(seq_len(ceiling(N / chunk)), function(cidx) {
      s <- seed * 1000 + cidx
      d1 <- sim_rwe(chunk, p, arm = "tdf", censoring = FALSE, seed = s, delta = dl)$data
      d0 <- sim_rwe(chunk, p, arm = "etv", censoring = FALSE, seed = s, delta = dl)$data
      c(mean(d1$T_nat <= h), mean(d0$T_nat <= h), mean(d1$leave <= h), mean(d0$leave <= h))
    }))
    m <- colMeans(acc)
    data.frame(delta = dl, risk_tdf = m[1], risk_etv = m[2], rd = m[1] - m[2], rr = m[1] / m[2],
               left_1y_tdf = m[3], left_1y_etv = m[4])
  }))
}
