# 75_job5.R -- desktop job 5 (Andrew 2026-11-22, after the sci-review of the job 1-4 results):
#   points   Cox-family risk estimators not in run/65, in the 6 main scenarios x 200 cohorts:
#              while on treatment  fg_std   Fine-Gray model (renal event, leaving the drug competing,
#                                           first year), standardized 1-year cumulative incidence
#              sustained use       wat_std  weighted as-treated Cox model (the ipcw_long weights of the
#                                           main analysis), standardized 1-year risks, fit to the first
#                                           year and to all follow-up
#            Point estimates only (no analytic interval exists for either).
#   boot     percentile intervals for fg_std and wat_std (first year) from B = 100 bootstrap samples
#            of patients, in the 200 base-case cohorts (weights and models refit in every sample)
#   ltfu     scenario d2f_informative_ltfu_strong (loss to follow-up 30 times as fast after a measured
#            eGFR decline of 10% from baseline) on the BASE CASE's seeds, so its difference from the
#            base case is paired; the one-feature estimator set (Cox analyses, SDR for treatment
#            policy and sustained use); computes the scenario's truth first (about an hour)
# Fine-Gray times are grouped into 7-day bins before the finegray() expansion (10 times fewer
# rows; in a test cohort the subdistribution log HR agreed to 4 decimals with daily times).
# One checkpoint per cohort under sandbox_rwe/results/runs/job5/<part>/.
#
#   Rscript --vanilla sandbox_rwe/run/75_job5.R pilot          # 1 cohort of each part, timed
#   Rscript --vanilla sandbox_rwe/run/75_job5.R points|boot|ltfu
setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
for (f in list.files("sandbox_rwe/R", pattern = "[.]R$", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(survival))
H <- 365; TSTAR <- 364.9; B <- 100
ROOT <- "sandbox_rwe/results/runs/job5"
S_LTFU <- "d2f_informative_ltfu_strong"
cores <- as.integer(Sys.getenv("RWE_CORES", "4"))

# ---- estimators: each returns c(risk1, risk0)
fg_std <- function(sim) {
  d <- with_stages(sim$data)
  end <- pmin(d$follow, d$leave, H)
  st <- ifelse(d$event == 1 & d$T_nat <= end, 1L, ifelse(d$leave <= pmin(d$follow, H), 2L, 0L))
  dc <- cbind(d[, c("A", COX_X)], time = pmin(ceiling(end / 7) * 7, H), status = st)
  dc$ev <- factor(st, 0:2, c("censor", "renal", "left"))
  fg <- finegray(as.formula(paste("Surv(time, ev) ~", paste(c("A", COX_X), collapse = " + "))), data = dc, etype = "renal")
  fit <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", cox_rhs)), data = fg, weights = fgwt, ties = "breslow")
  bh <- basehaz(fit, centered = FALSE)
  L0 <- approx(bh$time, bh$hazard, TSTAR, method = "constant", rule = 2)$y
  X <- model.matrix(as.formula(paste("~", cox_rhs)), dc)[, -1]
  r <- function(a) { X[, "A"] <- a; mean(1 - exp(-L0 * exp(drop(X %*% coef(fit))))) }
  c(r(1), r(0))
}
wat_std <- function(sim, h) {
  d <- with_stages(sim$data)
  long <- ipcw_long(sim, h = h)
  long <- cbind(long, d[match(long$id, d$id), c("A", COX_X)])
  fit <- coxph(as.formula(paste("Surv(start, stop, status) ~", cox_rhs)), data = long, weights = w, ties = "breslow")
  bh <- basehaz(fit, centered = FALSE)
  L0 <- approx(bh$time, bh$hazard, TSTAR, method = "constant", rule = 2)$y
  X <- model.matrix(as.formula(paste("~", cox_rhs)), d)[, -1]
  r <- function(a) { X[, "A"] <- a; mean(1 - exp(-L0 * exp(drop(X %*% coef(fit))))) }
  c(r(1), r(0))
}
ESTS <- list(fg_std_1y = function(s) fg_std(s), wat_std_1y = function(s) wat_std(s, H),
             wat_std_all = function(s) wat_std(s, Inf))

# bootstrap sample of patients: rows of data, vis and meas together, ids renumbered
boot_sim <- function(sim) {
  i <- sample.int(nrow(sim$data), replace = TRUE)
  out <- sim; out$data <- sim$data[i, ]; out$data$id <- seq_along(i)
  out$vis <- sim$vis[i, , drop = FALSE]; out$meas <- sim$meas[i, , drop = FALSE]
  out
}

one_points <- function(s, r) {
  p <- get_scenario_params(s); sim <- sim_rwe(20000, p, seed = seed_for(s, r))
  do.call(rbind, lapply(names(ESTS), function(nm) {
    t0 <- Sys.time()
    x <- tryCatch(ESTS[[nm]](sim), error = function(e) structure(c(NA, NA), msg = conditionMessage(e)))
    data.frame(scenario = s, rep = r, method = nm, risk1 = x[1], risk0 = x[2],
               ok = all(is.finite(x)), error = if (is.null(attr(x, "msg"))) NA else attr(x, "msg"),
               secs = as.numeric(Sys.time() - t0, units = "secs"), fingerprint = params_fingerprint(p))
  }))
}
one_boot <- function(s, r) {
  p <- get_scenario_params(s); sim <- sim_rwe(20000, p, seed = seed_for(s, r))
  set.seed(estimator_seed(s, r, "job5_boot"))
  t0 <- Sys.time()
  bs <- t(sapply(seq_len(B), function(b) {
    sb <- boot_sim(sim)
    c(fg = tryCatch(fg_std(sb), error = function(e) c(NA, NA)),
      wat = tryCatch(wat_std(sb, H), error = function(e) c(NA, NA)))
  }))
  q <- function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE)
  data.frame(scenario = s, rep = r, method = c("fg_std_1y", "wat_std_1y"),
             rd_lo = c(q(bs[, 1] - bs[, 2])[1], q(bs[, 3] - bs[, 4])[1]),
             rd_hi = c(q(bs[, 1] - bs[, 2])[2], q(bs[, 3] - bs[, 4])[2]),
             rr_lo = c(q(bs[, 1] / bs[, 2])[1], q(bs[, 3] / bs[, 4])[1]),
             rr_hi = c(q(bs[, 1] / bs[, 2])[2], q(bs[, 3] / bs[, 4])[2]),
             n_boot_ok = c(sum(is.finite(bs[, 1] + bs[, 2])), sum(is.finite(bs[, 3] + bs[, 4]))),
             secs = as.numeric(Sys.time() - t0, units = "secs"), fingerprint = params_fingerprint(p))
}
# the one-feature estimator set of run_replicate, on the base case's seeds
one_ltfu <- function(s, r) {
  p <- get_scenario_params(s)
  sim <- sim_rwe(20000, p, seed = seed_for("base", r))
  ests <- estimator_set("deviation", sim)
  rows <- do.call(rbind, lapply(names(ests), function(nm) {
    set.seed(estimator_seed("base", r, nm))
    safe_rows(ests[[nm]]$fun, ests[[nm]]$meta)
  }))
  rows$scenario <- s; rows$rep <- r; rows$set <- "deviation_baseseeds"
  rows$fingerprint <- params_fingerprint(p)
  rows[, ROW_COLS]
}

cached <- function(part, s, r) {
  f <- file.path(ROOT, part, s, sprintf("rep_%03d.rds", r))
  if (file.exists(f)) { x <- tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(x)) return(x) }
  x <- switch(part, points = one_points(s, r), boot = one_boot(s, r), ltfu = one_ltfu(s, r))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(f, ".tmp"); saveRDS(x, tmp); if (file.exists(f)) file.remove(f); file.rename(tmp, f)
  x
}
run_part <- function(part, jobs) {
  todo <- jobs[!file.exists(file.path(ROOT, part, jobs$s, sprintf("rep_%03d.rds", jobs$r))), ]
  cat(format(Sys.time(), "%m-%d %H:%M"), "-", part, ":", nrow(todo), "cohorts to run\n")
  for (ch in split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / (2 * cores)))) {
    cl <- parallel::makeCluster(cores)
    tryCatch({
      parallel::clusterEvalQ(cl, {
        setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
        for (f in list.files("sandbox_rwe/R", pattern = "[.]R$", full.names = TRUE)) source(f)
        suppressPackageStartupMessages(library(survival))
      })
      parallel::clusterExport(cl, c("H", "TSTAR", "B", "ROOT", "fg_std", "wat_std", "ESTS", "boot_sim",
                                    "one_points", "one_boot", "one_ltfu", "cached"))
      parallel::parLapply(cl, ch, function(i, todo, part) tryCatch(cached(part, todo$s[i], todo$r[i]),
                                                                   error = function(e) NULL), todo = todo, part = part)
    }, finally = parallel::stopCluster(cl))
    cat(format(Sys.time(), "%m-%d %H:%M"), "- done", max(ch), "of", nrow(todo), "\n")
  }
  rows <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
    f <- file.path(ROOT, part, jobs$s[i], sprintf("rep_%03d.rds", jobs$r[i])); if (file.exists(f)) readRDS(f)
  }))
  rows
}

mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode)) stop("mode: pilot, points, boot or ltfu")
if (mode == "pilot") {
  print(one_points("base", 1), digits = 4)
  B <- 3; t0 <- Sys.time(); print(one_boot("base", 1)[, c("method", "rd_lo", "rd_hi", "n_boot_ok", "secs")], digits = 4)
  cat(sprintf("bootstrap: %.0f s for 3 samples -> about %.0f min per cohort at B = 100\n",
              as.numeric(Sys.time() - t0, units = "secs"), as.numeric(Sys.time() - t0, units = "mins") / 3 * 100))
  t0 <- Sys.time(); x <- one_ltfu(S_LTFU, 1)
  print(x[, c("family", "strategy", "method", "window", "est", "lo", "hi", "ok")], digits = 3)
  cat(sprintf("ltfu cohort: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
} else if (mode == "points") {
  rows <- run_part("points", expand.grid(r = 1:200, s = SCENARIOS, stringsAsFactors = FALSE))
  saveRDS(rows, "sandbox_rwe/results/job5_points_rows.rds")
  cat("rows:", nrow(rows), " failed:", sum(!rows$ok %in% TRUE), "\n")
} else if (mode == "boot") {
  rows <- run_part("boot", data.frame(s = "base", r = 1:200))
  saveRDS(rows, "sandbox_rwe/results/job5_boot_rows.rds")
  cat("rows:", nrow(rows), " fewer than 95 good samples:", sum(rows$n_boot_ok < 95), "\n")
} else if (mode == "ltfu") {
  p <- get_scenario_params(S_LTFU)
  tf <- sprintf("sandbox_rwe/results/truth_%s.rds", S_LTFU)
  if (!file.exists(tf) || !identical(readRDS(tf)$fingerprint, params_fingerprint(p))) saveRDS(scenario_truth(S_LTFU, p), tf)
  tr <- get_truth(S_LTFU)
  print(cbind(tr$strategies[, c("strategy", "rd", "rr")], base_rd = get_truth("base")$strategies$rd), digits = 4)
  rows <- run_part("ltfu", data.frame(s = S_LTFU, r = 1:200))
  perf <- summarise_runs(rows, tr)
  saveRDS(list(rows = rows, perf = perf, code_version = code_version()), "sandbox_rwe/results/informative_ltfu_strong_rows.rds")
  cat("rows:", nrow(rows), " failed:", sum(!rows$ok %in% TRUE), "\n")
}
