# principal_score.R -- real-world estimate of the principal stratum: patients who
# would stay on either drug through day h. Assumes monotonicity (anyone who would
# stay on TDF would stay on entecavir), so the stratum is {S*(1) = 1}, and
# principal ignorability given the baseline covariates X.
#
# Adherence through h is truncated by the event: after an event we never see
# whether the patient would have left later, so "did not leave before the event"
# is not membership (13% of such TDF patients would have left before h in the
# test cohort). Both pieces are therefore survival models fitted on observed data:
#   e1(X)  principal score P(S*(1) = 1 | X): time to leaving the drug among TDF
#          starters, censored at the event and at loss to follow-up (Cox)
#   m_a(X) risk by h under staying on drug a: time to the event among starters of
#          a, censored at leaving and at loss to follow-up (Cox)
# Stratum risk on drug a: sum(e1 * m_a) / sum(e1). Share of the cohort: mean(e1).

.cox_surv_at <- function(time, status, X, newX, h) {
  fit <- survival::coxph(survival::Surv(time, status) ~ ., data = data.frame(time, status, X),
                         ties = "breslow")
  H0 <- survival::basehaz(fit, centered = FALSE)
  H0h <- if (any(H0$time <= h)) max(H0$hazard[H0$time <= h]) else 0
  lp <- as.vector(as.matrix(newX[, names(coef(fit))]) %*% coef(fit))
  exp(-H0h * exp(lp))
}

.psw_once <- function(d, h) {
  X <- d[, COX_X]
  t_off <- ifelse(d$left_obs == 1, d$leave, d$follow)       # leaving, or end of follow-up
  i1 <- d$A == 1
  e1 <- suppressWarnings(.cox_surv_at(t_off[i1], d$left_obs[i1], X[i1, ], X, h))
  ev_on <- as.integer(d$event == 1 & d$left_obs == 0)       # event while still on the drug
  m <- sapply(c(1, 0), function(a) {
    i <- d$A == a
    1 - suppressWarnings(.cox_surv_at(t_off[i], ev_on[i], X[i, ], X, h))
  })
  r1 <- sum(e1 * m[, 1]) / sum(e1); r0 <- sum(e1 * m[, 2]) / sum(e1)
  c(risk1 = r1, risk0 = r0, rd = r1 - r0, rr = r1 / r0, share = mean(e1))
}

ps_weighting <- function(sim, h = 365, B = 0) {
  d <- with_stages(sim$data)
  est <- .psw_once(d, h)
  ci <- c(rd_lo = NA_real_, rd_hi = NA_real_, rr_lo = NA_real_, rr_hi = NA_real_)
  if (B > 0) {
    bs <- replicate(B, .psw_once(d[sample.int(nrow(d), replace = TRUE), ], h)[c("rd", "rr")])
    q <- function(x, pr) unname(quantile(x, pr, na.rm = TRUE))
    ci <- c(rd_lo = q(bs["rd", ], 0.025), rd_hi = q(bs["rd", ], 0.975),
            rr_lo = q(bs["rr", ], 0.025), rr_hi = q(bs["rr", ], 0.975))
  }
  c(est, ci)
}
