# analysis.R -- Cox analyses and SDR (lmtp) for the RWE sandbox.
#
# Cox analyses (baseline covariates COX_X), each fit to all follow-up and to the
# 1-year window:
#   itt          intention to treat: all follow-up regardless of leaving the drug
#   astreated    as treated: follow-up ends when the patient leaves the drug
#                (as an estimate of the cause-specific HR, this is the
#                while-on-treatment analysis)
#   astreated_w  as treated, weighted by the inverse probability of not having
#                switched, from a visit-level model with the measured eGFR
#   composite    renal event or leaving the drug, whichever first
#   adherers     intention to treat among patients not seen to leave the drug in
#                their first year of follow-up
# SDR: lmtp_sdr on K monthly intervals, baseline treatment, time-varying measured
# eGFR (last value carried forward), one coding per strategy.

suppressPackageStartupMessages({ library(survival); library(lmtp); library(dplyr) })

# Adjustment set shared by Cox and SDR. Baseline eGFR enters linearly and as CKD
# stage indicators (60-89, 45-59, < 45; reference >= 90), the usual coding in
# renal safety analyses.
COX_X <- c("age", "male", "dm", "htn", "cirr", "nsaid", "egfr0",
           "egfr_60_89", "egfr_45_59", "egfr_lt45")
cox_rhs <- paste("A +", paste(COX_X, collapse = " + "))
with_stages <- function(d) {
  d$egfr_60_89 <- as.integer(d$egfr0 >= 60 & d$egfr0 < 90)
  d$egfr_45_59 <- as.integer(d$egfr0 >= 45 & d$egfr0 < 60)
  d$egfr_lt45  <- as.integer(d$egfr0 < 45)
  d
}

window <- function(time, status, h) list(time = pmin(time, h), status = as.integer(status & time <= h))

cox_hr <- function(dd, weighted = FALSE) {
  f <- if ("start" %in% names(dd)) {
    as.formula(paste("Surv(start, stop, status) ~", cox_rhs, if (weighted) "+ cluster(id)"))
  } else as.formula(paste("Surv(time, status) ~", cox_rhs))
  fit <- if (weighted) coxph(f, data = dd, weights = w, ties = "breslow")
         else coxph(f, data = dd, ties = "breslow")
  b <- unname(coef(fit)["A"]); se <- sqrt(vcov(fit)["A", "A"])
  if (!all(is.finite(c(b, se))) || se > 10) stop("treatment effect not estimable")
  c(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))
}

# Counting-process data for the weighted as-treated analysis
# Switching model at visits: decline from baseline (per 10%), measured eGFR below
# 60, age over 60 and calendar time on the drug, each by drug (an analyst can fit
# this from the visit records). oracle_p = DGP parameters gives the true switching
# probabilities instead (semi-oracle weights).
ipcw_long <- function(sim, h = Inf, trunc = c(0.01, 0.99), oracle_p = NULL) {
  d <- with_stages(sim$data); vis <- sim$vis; meas <- sim$meas
  end <- pmin(d$follow, d$leave)
  status <- as.integer(d$event == 1 & d$T_nat <= end)
  J <- ncol(vis)
  vv <- data.frame(id = rep(d$id, J), j = rep(seq_len(J), each = nrow(d)),
                   vis = as.vector(vis), meas = as.vector(meas),
                   A = rep(d$A, J), age = rep(d$age, J),
                   end = rep(end, J), follow = rep(d$follow, J),
                   leave = rep(d$leave, J), type = rep(d$leave_type, J),
                   egfr0 = rep(d$egfr0, J))
  # visits at which the patient was on the drug and under observation
  vv <- vv[vv$vis <= vv$end & vv$vis < vv$follow, ]
  vv$sw <- as.integer(vv$type == "switch" & abs(vv$vis - vv$leave) < 1e-8)
  vv$dec10 <- 10 * pmax(0, (vv$egfr0 - vv$meas) / vv$egfr0)
  vv$lt60 <- as.integer(vv$meas < 60); vv$old <- as.integer(vv$age > 60); vv$yr <- vv$vis / 365
  p_den <- if (is.null(oracle_p)) {
    predict(glm(sw ~ A * (dec10 + lt60 + old + yr), family = binomial, data = vv), type = "response")
  } else {
    k <- vv$A + 1
    with(oracle_p, plogis(sw_b0[k] + sw_t[k] * vv$yr + sw_dec[k] * vv$dec10 +
                            sw_low60[k] * vv$lt60 + sw_age60[k] * vv$old))
  }
  num <- glm(sw ~ A, family = binomial, data = vv)
  vv$ratio <- (1 - predict(num, type = "response")) / (1 - p_den)
  vv <- vv[vv$sw == 0 & vv$vis < vv$end, ]           # visits that start a new interval
  vv <- vv[order(vv$id, vv$vis), ]
  vv$cw <- ave(vv$ratio, vv$id, FUN = cumprod)
  starts <- rbind(data.frame(id = d$id, start = 0, w = 1),
                  data.frame(id = vv$id, start = vv$vis, w = vv$cw))
  starts <- starts[order(starts$id, starts$start), ]
  starts$stop <- ave(starts$start, starts$id, FUN = function(s) c(s[-1], NA))
  idx <- match(starts$id, d$id)
  last <- is.na(starts$stop)
  starts$stop[last] <- end[idx[last]]
  starts$status <- as.integer(last & status[idx] == 1)
  if (is.finite(h)) {
    starts <- starts[starts$start < h, ]
    starts$status <- as.integer(starts$status == 1 & starts$stop <= h)
    starts$stop <- pmin(starts$stop, h)
  }
  starts <- starts[starts$stop > starts$start, ]
  q <- quantile(starts$w, trunc); starts$w <- pmin(pmax(starts$w, q[1]), q[2])
  starts[, c("id", "start", "stop", "status", "w")]
}

cox_suite <- function(sim, h = 365) {
  d <- with_stages(sim$data)
  as_end <- pmin(d$follow, d$leave)
  as_status <- as.integer(d$event == 1 & d$T_nat <= as_end)
  comp_status <- as.integer(pmin(d$T_nat, d$leave) <= d$follow)
  adh <- d$leave >= pmin(d$follow, h)
  specs <- list(
    itt       = list(d, d$follow, d$event),
    astreated = list(d, as_end, as_status),
    composite = list(d, as_end, comp_status),
    adherers  = list(d[adh, ], d$follow[adh], d$event[adh]))
  # One degenerate analysis (e.g. no events in an arm of the adherer subset) gives
  # a failed row; the others still run.
  safe_hr <- function(dd, ...) tryCatch(c(cox_hr(dd, ...), ok = 1, error = NA),
                                        error = function(e) c(hr = NA, lo = NA, hi = NA, ok = 0,
                                                              error = conditionMessage(e)))
  row_of <- function(nm, win, r) data.frame(
    analysis = nm, window = win, hr = as.numeric(r[["hr"]]), lo = as.numeric(r[["lo"]]),
    hi = as.numeric(r[["hi"]]), ok = r[["ok"]] == "1", error = r[["error"]], stringsAsFactors = FALSE)
  out <- list()
  for (nm in names(specs)) {
    s <- specs[[nm]]
    for (win in c("all", "1y")) {
      tt <- if (win == "all") list(time = s[[2]], status = s[[3]]) else window(s[[2]], s[[3]], h)
      dd <- cbind(s[[1]][, c("id", "A", COX_X)], time = tt$time, status = tt$status)
      out[[length(out) + 1]] <- row_of(nm, win, safe_hr(dd))
    }
  }
  for (win in c("all", "1y")) {
    r <- tryCatch({
      long <- ipcw_long(sim, h = if (win == "all") Inf else h)
      long <- cbind(long[, c("id", "start", "stop", "status", "w")],
                    d[match(long$id, d$id), c("A", COX_X)])
      safe_hr(long, weighted = TRUE)
    }, error = function(e) c(hr = NA, lo = NA, hi = NA, ok = 0, error = conditionMessage(e)))
    out[[length(out) + 1]] <- row_of("astreated_w", win, r)
  }
  do.call(rbind, out)
}

# ---------------------------------------------------------------- SDR (lmtp)
# Smoke-test learner library beyond GLM: MARS (earth defaults, degree 2) and
# boosted trees (depth 3, 200 trees, learning rate 0.1, no minimum leaf weight).
# Untuned; the paper's library is set with the supervised-learning rules later.
SL.xgb_s <- function(...) SuperLearner::SL.xgboost(..., ntrees = 200, max_depth = 3,
                                                   shrinkage = 0.1, minobspernode = 0)
LIB_FLEX <- c("SL.glm", "SL.earth", "SL.xgb_s")

# Logistic regression with treatment-by-covariate interactions, for lmtp's
# density-ratio classifier: the ratio for a fixed treatment is 1/P(A | history)
# in the treated and 0 otherwise, which a main-effects model cannot represent.
SL.glm.Aint <- function(Y, X, newX, family, obsWeights, ...) {
  rhs <- if ("A" %in% names(X)) {
    paste("A * (", paste(setdiff(names(X), "A"), collapse = " + "), ")")
  } else "."
  fit <- suppressWarnings(glm(as.formula(paste("Y ~", rhs)), data = cbind(Y = Y, X),
                              family = family, weights = obsWeights))
  pred <- predict(fit, newdata = newX, type = "response")
  out <- list(object = fit); class(out) <- "SL.glm"
  list(pred = pred, fit = out)
}

prep_lmtp <- function(sim, strategy, K = 12, h = 365, rows = NULL, tv_features = FALSE,
                      L_at = NULL, switch_terms = strategy == "reduced_leaving") {
  # The history L_k is the measured eGFR carried forward to a cut point in interval
  # k. Where leaving the drug is censored or intervened on, the cut matches how
  # leaving is coded, so the switch model at interval k sees the measurement from
  # the visit where the switch was decided: both codings place leaving at the
  # nearest interval edge, so leaving coded at interval k happened before the
  # middle of interval k ("mid"). With the interval start, these models missed that
  # measurement for about half of the switches. Elsewhere leaving is not modelled
  # and the history is the interval start.
  if (is.null(L_at)) L_at <- switch(strategy, sustained_use = "mid", reduced_leaving = "mid", "start")
  L_at <- match.arg(L_at, c("start", "mid", "end"))
  d <- with_stages(sim$data); vis <- sim$vis; meas <- sim$meas
  if (is.null(rows)) rows <- rep(TRUE, nrow(d))
  edges <- h * (0:K) / K
  # Follow-up, and so visits, stop at the event; for sustained use and reduced
  # leaving the history also stops at leaving (the switch visit included), so no
  # measurement taken after leaving enters a switch model.
  stop_at <- if (strategy %in% c("sustained_use", "reduced_leaving")) d$leave else Inf
  L <- sapply(2:K, function(k) {
    cut <- switch(L_at, start = edges[k], mid = (edges[k] + edges[k + 1]) / 2, end = edges[k + 1])
    n_seen <- rowSums(vis <= pmin(cut, stop_at) & vis < d$follow)
    ifelse(n_seen == 0, d$egfr0, meas[cbind(seq_len(nrow(d)), pmax(n_seen, 1))])
  })
  colnames(L) <- paste0("L", 2:K)
  # Optional switching-model features for follow-up eGFR: below 60, and the
  # fractional decline from baseline (the drivers of switching).
  L60 <- (L < 60) * 1L; colnames(L60) <- paste0("L60_", 2:K)
  Ldec <- (d$egfr0 - L) / d$egfr0; Ldec[Ldec < 0] <- 0; colnames(Ldec) <- paste0("Ldec", 2:K)
  T <- d$T_nat; fol <- d$follow; lv <- d$leave
  ev_time <- rep(Inf, nrow(d)); cens_time <- rep(Inf, nrow(d)); comp_time <- rep(Inf, nrow(d))
  if (strategy %in% c("treatment_policy", "principal_stratum", "reduced_leaving")) {
    ev <- d$event == 1
    ev_time[ev] <- T[ev]; cens_time[!ev] <- fol[!ev]
  } else if (strategy == "sustained_use") {
    end <- pmin(fol, lv); ev <- d$event == 1 & T <= end
    ev_time[ev] <- T[ev]; cens_time[!ev] <- end[!ev]
  } else if (strategy == "while_on_treatment") {
    ev <- d$event == 1 & T <= lv
    left <- lv < fol & !ev
    ev_time[ev] <- T[ev]; comp_time[left] <- lv[left]
    cens_time[!ev & !left] <- fol[!ev & !left]
  } else if (strategy == "composite") {
    ct <- pmin(T, lv); ev <- ct <= fol
    ev_time[ev] <- ct[ev]; cens_time[!ev] <- fol[!ev]
  }
  Y <- sapply(1:K, function(k) as.integer(ev_time <= edges[k + 1]))
  # Censoring at the nearest interval edge (actuarial convention): a patient
  # censored in the second half of interval k counts as followed through k.
  # Dropping every patient censored mid-interval from that interval inflated each
  # monthly hazard by about half its censoring probability (pre-pilot review, C2).
  last_obs <- ifelse(is.finite(cens_time), round(cens_time / (h / K)), Inf)
  C <- sapply(1:K, function(k) as.integer(k <= last_obs))
  colnames(Y) <- paste0("Y", 1:K); colnames(C) <- paste0("C", 1:K)
  wide <- if (tv_features) data.frame(d[, c("id", "A", COX_X)], L, L60, Ldec, C, Y)
          else data.frame(d[, c("id", "A", COX_X)], L, C, Y)
  D_cols <- NULL
  if (strategy == "while_on_treatment") {
    D <- sapply(1:K, function(k) as.integer(comp_time <= edges[k + 1]))
    colnames(D) <- D_cols <- paste0("D", 1:K)
    wide <- data.frame(wide, D)
  }
  # Sustained use censors at leaving the drug, which depends on the measured eGFR.
  # With a single baseline treatment column lmtp's censoring models see only the
  # baseline covariates, so the treatment is given as copies A1..AK of A: the
  # censoring model at interval k then conditions on the eGFR history, and for a
  # patient who stays on the assigned drug the treatment part of the ratio is 1.
  trt <- "A"
  if (strategy == "sustained_use") {
    for (k in 1:K) wide[[paste0("A", k)]] <- wide$A
    trt <- paste0("A", 1:K)
  }
  # Reduced leaving (stochastic intervention): outcome and censoring as for
  # treatment policy, plus treatment nodes S_k = 1 once the patient has left the
  # initial drug (switch or gap), with leaving placed at the nearest interval edge
  # like censoring (S_k = 1 from interval round(leave / width) + 1). Counting
  # leaving anywhere in interval k as S_k = 1 dropped leavers from the whole
  # interval and raised the no-departure (sustained-use) risk difference by about
  # 0.0005 (final review, I1). The shift is reduced_leave_shift().
  if (strategy == "reduced_leaving") {
    last_on <- round(lv / (h / K))
    for (k in 1:K) wide[[paste0("S", k)]] <- as.integer(last_on < k & lv < fol)
    trt <- c(list(c("A", "S1")), as.list(paste0("S", 2:K)))
  }
  tv <- c(list(NULL), if (tv_features || switch_terms)
    lapply(2:K, function(k) c(paste0("L", k), paste0("L60_", k), paste0("Ldec", k))) else as.list(colnames(L)))
  W <- COX_X
  # The switching model gets the covariates of the Cox weight model (ipcw_long):
  # measured eGFR below 60, fractional decline and age over 60, each by drug (time
  # through lmtp's interval-specific fits). The drug products enter the time-varying
  # nodes, after the treatment node, so the propensity model never sees them.
  # Without them the pooled main-effects model under-fitted switching off TDF and
  # the sustained-use TDF risk was biased downward (pilot iteration 12).
  if (switch_terms) {
    if (!tv_features) wide <- data.frame(wide, L60, Ldec)
    wide$age60 <- as.integer(wide$age > 60); W <- c(W, "age60")
    for (k in 2:K) {
      wide[[paste0("AL60_", k)]] <- wide$A * wide[[paste0("L60_", k)]]
      wide[[paste0("ALdec", k)]] <- wide$A * wide[[paste0("Ldec", k)]]
      tv[[k]] <- c(tv[[k]], paste0(c("AL60_", "ALdec"), k))
    }
    wide$Aage60 <- wide$A * wide$age60; tv[[2]] <- c(tv[[2]], "Aage60")
  }
  wide <- wide[rows, ]
  wide <- as.data.frame(lmtp::event_locf(wide, outcomes = colnames(Y)))
  list(data = wide, trt = trt, Y = colnames(Y), C = colnames(C), D = D_cols, tv = tv, W = W)
}

# Density ratios are not trimmed (lmtp's default trims at the 99.9th percentile,
# which here clips the legitimate large weights of high-risk TDF patients with low
# eGFR and biased the risk difference downward; overlap is good, propensity
# 0.06-0.56).
# Stochastic reduction of leaving the drug, for the "reduced_leaving" coding: the
# initial drug is fixed at a, a patient who has already left stays off, and a new
# departure in an interval is kept with probability delta (lmtp's ipsi() logic,
# without letting patients restart the drug).
# With a list-style treatment, lmtp passes all of an interval's columns at once
# (c("A", "S1") at interval 1) and expects a data frame of shifted columns.
reduced_leave_shift <- function(a, delta) {
  force(a); force(delta)
  one_col <- function(data, col) {
    if (col == "A") return(rep(as.integer(a), nrow(data)))
    k <- as.integer(sub("^S", "", col))
    prev <- if (k > 1) data[[paste0("S", k - 1)]] else rep(0L, nrow(data))
    keep <- stats::runif(nrow(data)) < delta
    as.integer(ifelse(prev == 1, 1L, ifelse(keep, data[[col]], 0L)))
  }
  function(data, trt) {
    as.data.frame(setNames(lapply(trt, function(col) one_col(data, col)), trt))
  }
}

# Largest cumulative weight over every patient and interval (explosion check); a
# patient can reach a high weight and then be censored or have the event.
max_cum_weight <- function(dr) max(apply(dr, 1, function(r) max(cumprod(r))))

# Two Super Learner folds: with a single learner the inner cross-validated risk
# is never used, and lmtp's default of 10 folds took about 70% of the time.
# The flexible-library sensitivity passes its own control with 10 folds.
run_sdr <- function(pr, learners = "SL.glm", folds = 2,
                    learners_trt = learners,
                    control = lmtp_control(.trim = 0.9999, .learners_outcome_folds = 2,
                                           .learners_trt_folds = 2),
                    shift1 = static_binary_on, shift0 = static_binary_off, mtp = TRUE) {
  fit <- function(shift) lmtp_sdr(
    pr$data, trt = pr$trt, outcome = pr$Y, cens = pr$C, compete = pr$D,
    baseline = pr$W, time_vary = pr$tv, outcome_type = "survival", shift = shift, mtp = mtp,
    learners_trt = learners_trt, learners_outcome = learners, folds = folds,
    control = control)
  f1 <- suppressWarnings(fit(shift1)); f0 <- suppressWarnings(fit(shift0))
  maxw <- function(fx) max_cum_weight(fx$density_ratios)
  r1 <- 1 - f1$estimate@x; r0 <- 1 - f0$estimate@x
  e1 <- -f1$estimate@eif; e0 <- -f0$estimate@eif        # EIFs of the risks
  n <- length(e1)
  se_rd <- sd(e1 - e0) / sqrt(n)
  se_lrr <- sd(e1 / r1 - e0 / r0) / sqrt(n)
  c(risk1 = r1, risk0 = r0, rd = r1 - r0, rd_lo = r1 - r0 - 1.96 * se_rd,
    rd_hi = r1 - r0 + 1.96 * se_rd, rr = r1 / r0,
    rr_lo = exp(log(r1 / r0) - 1.96 * se_lrr), rr_hi = exp(log(r1 / r0) + 1.96 * se_lrr),
    max_w = max(maxw(f1), maxw(f0)))
}

sdr_suite <- function(sim, K = 12, h = 365, learners = "SL.glm", folds = 2) {
  d <- with_stages(sim$data)
  out <- list()
  for (s in c("treatment_policy", "sustained_use", "while_on_treatment", "composite")) {
    out[[s]] <- data.frame(strategy = s, method = "SDR",
                           t(run_sdr(prep_lmtp(sim, s, K, h), learners, folds)))
  }
  oracle <- d$leave1 > h & d$leave0 > h
  observed_adh <- d$leave >= pmin(d$follow, h)
  out$ps_oracle <- data.frame(strategy = "principal_stratum", method = "SDR, true stratum",
    t(run_sdr(prep_lmtp(sim, "principal_stratum", K, h, rows = oracle), learners, folds)))
  out$ps_adh <- data.frame(strategy = "principal_stratum", method = "SDR, observed adherers",
    t(run_sdr(prep_lmtp(sim, "principal_stratum", K, h, rows = observed_adh), learners, folds)))
  do.call(rbind, out)
}
