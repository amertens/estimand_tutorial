# dgp.R -- continuous-time RWE safety cohort (sandbox).
#
# One call simulates N patients. Each has baseline covariates, an initial drug
# (TDF, A = 1, or entecavir, A = 0), a latent eGFR that declines faster on TDF
# (early dip, small persistent excess, slow recovery after leaving TDF), clinic
# visits at irregular times with eGFR measured with error, and three clocks:
#   - leaving the initial drug: a switch decided at a visit from the measured eGFR
#     and age, or a gap in supply in continuous time, whichever comes first;
#   - the renal event: piecewise-constant hazard updated every `piece` days from
#     the latent eGFR and current drug, with an early-toxicity window on TDF;
#   - loss to follow-up and administrative censoring (staggered entry).
# Decisions depend only on measured history and baseline covariates, so the
# hypothetical strategies are identified from the observed data.
#
# Common random numbers: everything except the treatment draw uses the same
# random stream whatever `arm` is, so counterfactual cohorts simulated with the
# same seed are coupled (needed for the principal stratum).

dip_on <- function(t, p) p$dip_early * (1 - exp(-t / p$dip_tau)) + p$dip_late * t / 365

# TDF-induced eGFR loss at time t for a patient who leaves TDF at `leave`
tdf_dip <- function(t, leave, p) {
  ifelse(t < leave, dip_on(t, p),
         dip_on(pmin(leave, t), p) * exp(-(t - leave) / p$recover_tau))
}

sim_baseline <- function(N, p) {
  age   <- pmin(pmax(rnorm(N, p$age_mean, p$age_sd), 20), 85)
  male  <- rbinom(N, 1, p$p_male)
  dm    <- rbinom(N, 1, plogis(p$dm_b0 + p$dm_age * (age - p$age_mean)))
  htn   <- rbinom(N, 1, plogis(p$htn_b0 + p$htn_age * (age - p$age_mean)))
  cirr  <- rbinom(N, 1, p$p_cirr)
  nsaid <- rbinom(N, 1, p$p_nsaid)
  egfr0 <- pmin(pmax(p$egfr_b0 + p$egfr_age * (age - p$age_mean) + p$egfr_dm * dm +
                       p$egfr_htn * htn + rnorm(N, 0, p$egfr_sd), 15), 140)
  slope <- rnorm(N, p$slope_mean + p$slope_dm * dm + p$slope_htn * htn, p$slope_sd)
  # Optional susceptibility to TDF nephrotoxicity (supplementary scenario): a share
  # sus_prev has a TDF dip sus_mult times larger. Absent or 0: no draw, multiplier 1.
  dip_mult <- if (!is.null(p$sus_prev) && p$sus_prev > 0)
    ifelse(stats::runif(N) < p$sus_prev, p$sus_mult, 1) else rep(1, N)
  data.frame(age, male, dm, htn, cirr, nsaid, egfr0, slope, dip_mult)
}

trt_lp <- function(b, p) {
  with(b, p$trt_b0 + p$trt_egfr * (egfr0 - 90) + p$trt_low60 * (egfr0 < 60) +
         p$trt_age * (age - p$age_mean) + p$trt_dm * dm + p$trt_htn * htn +
         p$trt_cirr * cirr + p$trt_nsaid * nsaid)
}

# Unit draws for the visit schedule: first visit time and mean-1 gamma gaps.
# The schedule under drug a is first + visit_gap_mean[a] * cumsum(gaps), so the
# two drugs' schedules are coupled.
sim_visit_draws <- function(N, p) {
  list(first = runif(N, p$visit_first_min, p$visit_first_max),
       G = matrix(rgamma(N * (p$J - 1), p$visit_gap_shape, rate = p$visit_gap_shape), N))
}
visits_for <- function(vd, a, p) {
  v <- matrix(0, length(vd$first), p$J)
  v[, 1] <- vd$first
  m <- p$visit_gap_mean[a + 1]
  for (j in 2:p$J) v[, j] <- v[, j - 1] + m * vd$G[, j - 1]
  v
}

# Time of a gap in supply under drug a: piecewise-constant hazard (year 1, later)
gap_time <- function(u, age, a, p, delta = 1) {
  m <- exp(p$gap_age * (age - p$age_mean))
  r1 <- p$gap_r1[a + 1] * m * delta; r2 <- p$gap_r2[a + 1] * m * delta
  E <- -log(u)
  ifelse(E <= r1 * 365, E / r1, 365 + (E - r1 * 365) / r2)
}

# Time of leaving the initial drug under arm a (0 or 1), from the coupled draws
leave_time <- function(b, a, vis, eps, u_sw, g_time, p, delta = 1) {
  k <- a + 1
  meas_on <- b$egfr0 + b$slope * vis / 365 - a * b$dip_mult * dip_on(vis, p) + eps
  dec10 <- 10 * pmax(0, (b$egfr0 - meas_on) / b$egfr0)
  lp <- p$sw_b0[k] + p$sw_t[k] * vis / 365 + p$sw_dec[k] * dec10 +
    p$sw_low60[k] * (meas_on < 60) + p$sw_age60[k] * (b$age > 60)
  sw <- u_sw < delta * plogis(lp)
  any_sw <- rowSums(sw) > 0
  first <- max.col(sw, ties.method = "first")
  t_sw <- ifelse(any_sw, vis[cbind(seq_len(nrow(vis)), first)], Inf)
  data.frame(time = pmin(t_sw, g_time),
             type = ifelse(is.infinite(pmin(t_sw, g_time)), "none",
                           ifelse(t_sw < g_time, "switch", "gap")))
}

# Renal event time by inverse cumulative hazard over `piece`-day intervals
event_time <- function(b, A, leave, E, p) {
  N <- nrow(b); Tev <- rep(Inf, N); H <- numeric(N); open <- rep(TRUE, N)
  lp <- with(b, p$b_age * (age - p$age_mean) + p$b_dm * dm + p$b_htn * htn +
               p$b_cirr * cirr + p$b_nsaid * nsaid)
  starts <- seq(0, p$max_follow - p$piece, by = p$piece)
  for (s in starts) {
    e <- s + p$piece
    egfr_s <- b$egfr0 + b$slope * s / 365 - A * b$dip_mult * tdf_dip(s, leave, p)
    base <- p$h0 * exp(lp + p$b_egfr * (p$egfr_ref - egfr_s))
    # optional link to the decline from baseline (log-HR b_dec per 10% decline);
    # absent or 0 in the base case
    if (!is.null(p$b_dec) && p$b_dec != 0)
      base <- base * exp(p$b_dec * 10 * pmax(0, (b$egfr0 - egfr_s) / b$egfr0))
    hr <- ifelse(A == 1, if (s < p$t_change) p$HR_early else p$HR_late, 1)
    r_on <- base * hr; r_off <- base
    a1 <- r_on * pmin(pmax(leave - s, 0), p$piece)        # on-drug part
    a2 <- r_off * (p$piece - pmin(pmax(leave - s, 0), p$piece))
    h <- a1 + a2
    hit <- open & (H + h >= E)
    if (any(hit)) {
      need <- (E - H)[hit]
      in_on <- need <= a1[hit]
      t_on  <- s + need / r_on[hit]
      t_off <- pmax(s, pmin(leave[hit], e)) + (need - a1[hit]) / r_off[hit]
      Tev[hit] <- ifelse(in_on, t_on, t_off)
      open[hit] <- FALSE
    }
    H <- H + h
    if (!any(open)) break
  }
  Tev
}

# Simulate a cohort. arm: "observed" (drug chosen by the propensity model), "tdf"
# or "etv" (counterfactual). delta < 1 thins leaving the drug (each would-be switch
# or gap happens with probability delta), for the stochastic-intervention truths.
# Returns list(data, vis, meas).
sim_rwe <- function(N, p = rwe_params(), arm = c("observed", "tdf", "etv"),
                    censoring = TRUE, sustained = FALSE, seed = NULL, delta = 1) {
  if (!is.null(seed)) set.seed(seed)
  arm <- match.arg(arm)
  b <- sim_baseline(N, p)
  vd <- sim_visit_draws(N, p)
  eps <- matrix(rnorm(N * p$J, 0, p$meas_sd), N)
  u_sw <- matrix(runif(N * p$J), N)
  u_gap <- runif(N)
  E_nat <- rexp(N); E_sus <- E_nat
  ltfu <- rexp(N, 1)   # unit draw, scaled below (rate depends on arm)
  admin <- runif(N, 0, p$entry_window)
  u_trt <- runif(N)    # drawn for every arm so the streams stay coupled
  ps <- plogis(trt_lp(b, p))
  A <- switch(arm, observed = as.integer(u_trt < ps),
              tdf = rep(1L, N), etv = rep(0L, N))

  vis1 <- visits_for(vd, 1, p); vis0 <- visits_for(vd, 0, p)
  lv1 <- leave_time(b, 1, vis1, eps, u_sw, gap_time(u_gap, b$age, 1, p, delta), p, delta)
  lv0 <- leave_time(b, 0, vis0, eps, u_sw, gap_time(u_gap, b$age, 0, p, delta), p, delta)
  vis <- vis0; vis[A == 1, ] <- vis1[A == 1, ]
  leave <- ifelse(A == 1, lv1$time, lv0$time)
  leave_type <- ifelse(A == 1, lv1$type, lv0$type)

  T_nat <- event_time(b, A, leave, E_nat, p)
  T_sus <- if (sustained) event_time(b, A, rep(Inf, N), E_sus, p) else NA_real_

  if (censoring) {
    ltfu <- ltfu / (p$ltfu_rate * exp(p$ltfu_age * (b$age - p$age_mean) + p$ltfu_A * A))
  } else {
    ltfu <- rep(Inf, N); admin <- rep(Inf, N)
  }
  follow <- pmin(T_nat, ltfu, admin, p$max_follow)
  d <- data.frame(id = seq_len(N), b, A = A, ps = ps,
                  T_nat = T_nat, T_sus = T_sus,
                  leave = leave, leave_type = leave_type,
                  leave1 = lv1$time, leave0 = lv0$time,
                  ltfu = ltfu, admin = admin, follow = follow,
                  event = as.integer(T_nat <= follow),
                  left_obs = as.integer(leave < follow))
  # measured eGFR at visits along the patient's actual drug path
  meas <- b$egfr0 + b$slope * vis / 365 - A * b$dip_mult * tdf_dip(vis, leave, p) + eps
  list(data = d, vis = vis, meas = meas)
}
