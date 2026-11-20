# calibration.R -- fit the cohort to its targets (TARGETS.md). Each function
# returns an updated parameter list; calibrate_all() runs them in dependency order.

.wmean <- function(x, w) sum(w * x) / sum(w)
.wmed <- function(x, w) { o <- order(x); x[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]] }
.base_at <- function(p, n, seed) { set.seed(seed); sim_baseline(n, p) }

fit_baseline_dist <- function(p, tg, n = 1e5, seed = 7) {
  w <- tg$p_tdf
  overall <- function(v) w * v[["tdf"]] + (1 - w) * v[["etv"]]
  p$dm_b0 <- uniroot(function(x) mean(.base_at(modifyList(p, list(dm_b0 = x)), n, seed)$dm) -
                       overall(tg$dm), c(-6, 1))$root
  p$htn_b0 <- uniroot(function(x) mean(.base_at(modifyList(p, list(htn_b0 = x)), n, seed)$htn) -
                        overall(tg$htn), c(-6, 1))$root
  med_target <- overall(tg$egfr_median)
  fa <- function(th) {
    b <- .base_at(modifyList(p, list(egfr_b0 = th[1], egfr_sd = th[2])), n, seed)
    ((median(b$egfr0) - med_target) / 0.3)^2 + ((mean(b$egfr0 < 60) - overall(tg$egfr_lt60)) / 0.003)^2
  }
  o <- optim(c(p$egfr_b0, p$egfr_sd), fa)
  p$egfr_b0 <- o$par[1]; p$egfr_sd <- o$par[2]
  p
}

fit_drug_choice <- function(p, tg, n = 1e5, seed = 7) {
  b <- .base_at(p, n, seed)
  nm <- c("trt_b0", "trt_age", "trt_dm", "trt_htn", "trt_egfr", "trt_low60")
  stats_of <- function(pp) {
    ps <- plogis(trt_lp(b, pp)); w1 <- ps; w0 <- 1 - ps
    c(mean(ps), .wmean(b$age, w0) - .wmean(b$age, w1), .wmean(b$htn, w1), .wmean(b$htn, w0),
      .wmean(b$dm, w1), .wmean(b$dm, w0), .wmed(b$egfr0, w1), .wmed(b$egfr0, w0),
      .wmean(b$egfr0 < 60, w1), .wmean(b$egfr0 < 60, w0))
  }
  tgt <- c(tg$p_tdf, tg$age_gap, tg$htn[["tdf"]], tg$htn[["etv"]], tg$dm[["tdf"]], tg$dm[["etv"]],
           tg$egfr_median[["tdf"]], tg$egfr_median[["etv"]], tg$egfr_lt60[["tdf"]], tg$egfr_lt60[["etv"]])
  sc <- c(0.005, 0.1, 0.005, 0.005, 0.005, 0.005, 0.3, 0.3, 0.002, 0.002)
  f <- function(th) sum(((stats_of(modifyList(p, as.list(setNames(th, nm)))) - tgt) / sc)^2)
  o <- optim(unlist(p[nm]), f, control = list(maxit = 3000))
  p[nm] <- as.list(o$par)
  p
}

fit_visits <- function(p, tg, n = 1e5, seed = 7) {
  set.seed(seed); vd <- sim_visit_draws(n, p)
  share <- function(pp, a) {
    v <- visits_for(vd, a, pp)
    mean(sapply(0:3, function(w) rowSums(v >= w * 182.5 & v < (w + 1) * 182.5) > 0))
  }
  for (a in 0:1) {
    p$visit_gap_mean[a + 1] <- uniroot(function(m) {
      pp <- p; pp$visit_gap_mean[a + 1] <- m; share(pp, a) - tg$test_6m[a + 1]
    }, c(40, 800))$root
  }
  p
}

fit_switching <- function(p, t12, t36, n = 1e5, seed = 7) {
  set.seed(seed); b <- sim_baseline(n, p); vd <- sim_visit_draws(n, p)
  eps <- matrix(rnorm(n * p$J, 0, p$meas_sd), n); u <- matrix(runif(n * p$J), n)
  ps <- plogis(trt_lp(b, p))
  for (a in 0:1) {
    w <- if (a == 1) ps else 1 - ps
    v <- visits_for(vd, a, p)
    cum <- function(pp, t) .wmean(leave_time(b, a, v, eps, u, rep(Inf, n), pp)$time <= t, w)
    f <- function(th) {
      pp <- p; pp$sw_b0[a + 1] <- th[1]; pp$sw_t[a + 1] <- th[2]
      ((cum(pp, 365) - t12[[a + 1]]) / 0.005)^2 + ((cum(pp, 1095) - t36[[a + 1]]) / 0.005)^2
    }
    o <- optim(c(p$sw_b0[a + 1], p$sw_t[a + 1]), f)
    p$sw_b0[a + 1] <- o$par[1]; p$sw_t[a + 1] <- o$par[2]
  }
  p
}

fit_h0 <- function(p, target, n = 1e5, seed = 7) {
  f <- function(lh0) {
    pp <- p; pp$h0 <- exp(lh0)
    d <- sim_rwe(n, pp, arm = "etv", censoring = FALSE, seed = seed)$data
    .wmean(d$T_nat <= 365, 1 - d$ps) - target
  }
  p$h0 <- exp(uniroot(f, log(c(1e-7, 1e-2)), tol = 1e-5)$root)
  p
}

# Late HR (early fixed) or, with constant = TRUE, one HR for all time, chosen so
# the adjusted intention-to-treat Cox HR over all follow-up equals `target`.
fit_hr <- function(p, target, constant = FALSE, n = 2e5, seed = 7) {
  f <- function(h) {
    pp <- p
    if (constant) { pp$HR_early <- h; pp$HR_late <- h; pp$np_hazard <- FALSE } else pp$HR_late <- h
    d <- with_stages(sim_rwe(n, pp, seed = seed)$data)
    dd <- cbind(d[, c("id", "A", COX_X)], time = d$follow, status = d$event)
    cox_hr(dd)[["hr"]] - target
  }
  h <- uniroot(f, c(0.5, 3), tol = 0.005)$root
  if (constant) { p$HR_early <- h; p$HR_late <- h } else p$HR_late <- h
  p
}

calibrate_all <- function(p, tg, n = 1e5, seed = 7) {
  p <- fit_baseline_dist(p, tg, n, seed)
  p <- fit_drug_choice(p, tg, n, seed)
  p <- fit_visits(p, tg, n, seed)
  p <- fit_switching(p, tg$switch_12m, tg$switch_36m, n, seed)
  p <- fit_h0(p, tg$etv_risk_1y, n, seed)
  fit_hr(p, tg$adj_hr, constant = FALSE, n = 2 * n, seed = seed)
}
