# 65_comparators.R -- Cox-family analyses that target the strategies' own risk-scale or
# competing-risk quantities, so the paired Cox HRs are not the only Cox-based comparison
# (independent review 2026-11-20: "add Fine-Gray and standardized flexible-parametric risk
# comparators so the failures can't be dismissed as straw men").
#
# In the same 200 cohorts per main scenario as the main run (seed_for(scenario, rep)):
#   treatment policy   cox_std     Cox model (baseline covariates COX_X), standardized 1-year
#                                  risks under each drug (g-formula; riskRegression::ate,
#                                  influence-function intervals); fit to the first year and to
#                                  all follow-up
#                      rp_std      Royston-Parmar model (rstpm2::stpm2, 4 df baseline, 2 df
#                                  time-varying drug effect), all follow-up, standardized 1-year
#                                  risks (delta-method intervals)
#   composite          cox_std     as for treatment policy, composite event, first year
#   while on treatment csc_std     cause-specific Cox models for the renal event and for leaving
#                                  the drug, standardized 1-year cumulative incidence of the
#                                  renal event before leaving (riskRegression::ate)
#                      finegray    Fine-Gray subdistribution HR for the renal event, leaving as
#                                  the competing event, first year (survival::finegray + coxph)
# Risk-scale rows carry the RD in est/lo/hi and the RR in rr/rr_lo/rr_hi; the Fine-Gray row
# carries the subdistribution HR in est/lo/hi. Its true marginal subdistribution HR is fit on
# the first truth chunk (seed 1001, 250,000 per drug), as for the other true HRs.
# Not in R/run.R, so the main run's code version and checkpoints are untouched.
#
#   Rscript --vanilla sandbox_rwe/run/65_comparators.R pilot     # rep 1 of base, timings
#   Rscript --vanilla sandbox_rwe/run/65_comparators.R truth     # true subdistribution HRs
#   Rscript --vanilla sandbox_rwe/run/65_comparators.R           # 6 scenarios x 200 reps
# Launch sandbox_rwe/run/no_ecoqos.ps1 -Watch alongside the full run.
setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
for (f in list.files("sandbox_rwe/R", pattern = "\\.R$", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(survival); library(riskRegression); library(prodlim); library(rstpm2) })
H <- 365; TSTAR <- 364.9           # risks at 1 year; just inside the 1-year window
ROOT <- "sandbox_rwe/results/runs/comparators"

.row <- function(strategy, method, window, est, lo, hi, rr = NA, rr_lo = NA, rr_hi = NA)
  data.frame(family = if (method == "finegray") "cox" else "cox_risk", strategy = strategy, method = method,
             window = window, est = est, lo = lo, hi = hi, rr = rr, rr_lo = rr_lo, rr_hi = rr_hi,
             ok = all(is.finite(c(est, lo, hi))), error = NA_character_, stringsAsFactors = FALSE)
.fail <- function(strategy, method, window, e) {
  r <- .row(strategy, method, window, NA, NA, NA); r$ok <- FALSE; r$error <- conditionMessage(e); r
}

# Standardized 1-year risks from a Cox model (riskRegression::ate). The RR interval is
# formed on the log scale.
std_cox <- function(d, time, status, strategy, window) {
  dd <- cbind(d[, c("A", COX_X)], time = time, status = status); dd$A <- factor(dd$A, levels = c(0, 1))
  fit <- coxph(as.formula(paste("Surv(time, status) ~", cox_rhs)), data = dd, ties = "breslow", x = TRUE, y = TRUE)
  a <- ate(fit, treatment = "A", times = TSTAR, data = dd, se = TRUE, verbose = FALSE)
  rd <- a$diffRisk; rr <- a$ratioRisk
  lrr_se <- rr$se / rr$estimate                                  # delta method, log scale
  .row(strategy, "cox_std", window, rd$estimate, rd$lower, rd$upper, rr$estimate,
       exp(log(rr$estimate) - 1.96 * lrr_se), exp(log(rr$estimate) + 1.96 * lrr_se))
}

std_rp <- function(d) {
  dr <- cbind(d[, c("A", COX_X)], time = d$follow, status = d$event)
  rp <- stpm2(as.formula(paste("Surv(time, status) ~", cox_rhs)), data = dr, df = 4, tvc = list(A = 2))
  risk <- function(obj, newdata, a) 1 - predict(obj, newdata = transform(newdata, time = H, A = a), type = "meansurv")
  f_rd <- function(obj, newdata) risk(obj, newdata, 1) - risk(obj, newdata, 0)
  f_lrr <- function(obj, newdata) log(risk(obj, newdata, 1) / risk(obj, newdata, 0))
  x <- predictnl(rp, f_rd, newdata = dr); y <- predictnl(rp, f_lrr, newdata = dr)
  .row("treatment_policy", "rp_std", "all", x$fit, x$fit - 1.96 * x$se.fit, x$fit + 1.96 * x$se.fit,
       exp(y$fit), exp(y$fit - 1.96 * y$se.fit), exp(y$fit + 1.96 * y$se.fit))
}

# Competing-risk data on the first year: 1 = renal event on the drug, 2 = left the drug first
cr_data <- function(d) {
  end <- pmin(d$follow, d$leave, H)
  st <- ifelse(d$event == 1 & d$T_nat <= end, 1L, ifelse(d$leave <= pmin(d$follow, H), 2L, 0L))
  dc <- cbind(d[, c("A", COX_X)], time = end, status = st); dc$A <- factor(dc$A, levels = c(0, 1)); dc
}
std_csc <- function(d) {
  dc <- cr_data(d)
  csc <- CSC(as.formula(paste("Hist(time, status) ~", cox_rhs)), data = dc)
  a <- ate(csc, treatment = "A", times = TSTAR, cause = 1, data = dc, se = TRUE, verbose = FALSE, product.limit = FALSE)
  rd <- a$diffRisk; rr <- a$ratioRisk; lrr_se <- rr$se / rr$estimate
  .row("while_on_treatment", "csc_std", "1y", rd$estimate, rd$lower, rd$upper, rr$estimate,
       exp(log(rr$estimate) - 1.96 * lrr_se), exp(log(rr$estimate) + 1.96 * lrr_se))
}
fine_gray <- function(d) {
  dc <- cr_data(d); dc$ev <- factor(dc$status, levels = 0:2, labels = c("censor", "renal", "left"))
  fg <- finegray(as.formula(paste("Surv(time, ev) ~", paste(c("A", COX_X), collapse = " + "))), data = dc, etype = "renal")
  fit <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", cox_rhs)), data = fg, weights = fgwt)
  b <- unname(coef(fit)["A1"]); se <- sqrt(vcov(fit)["A1", "A1"])
  .row("while_on_treatment", "finegray", "1y", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
}

comp_rep <- function(scenario, rep) {
  p <- get_scenario_params(scenario)
  d <- with_stages(sim_rwe(20000, p, seed = seed_for(scenario, rep))$data)
  as_end <- pmin(d$follow, d$leave)
  comp_status <- as.integer(pmin(d$T_nat, d$leave) <= d$follow)
  w_tp <- window(d$follow, d$event, H); w_comp <- window(as_end, comp_status, H)
  jobs <- list(
    list("treatment_policy", "cox_std", "1y", function() std_cox(d, w_tp$time, w_tp$status, "treatment_policy", "1y")),
    list("treatment_policy", "cox_std", "all", function() std_cox(d, d$follow, d$event, "treatment_policy", "all")),
    list("treatment_policy", "rp_std", "all", function() std_rp(d)),
    list("composite", "cox_std", "1y", function() std_cox(d, w_comp$time, w_comp$status, "composite", "1y")),
    list("while_on_treatment", "csc_std", "1y", function() std_csc(d)),
    list("while_on_treatment", "finegray", "1y", function() fine_gray(d)))
  rows <- do.call(rbind, lapply(jobs, function(j) {
    t0 <- Sys.time()
    r <- tryCatch(suppressWarnings(j[[4]]()), error = function(e) .fail(j[[1]], j[[2]], j[[3]], e))
    r$secs <- as.numeric(Sys.time() - t0, units = "secs"); r
  }))
  rows$scenario <- scenario; rows$rep <- rep; rows$fingerprint <- params_fingerprint(p)
  rows
}

cached <- function(scenario, rep) {
  f <- file.path(ROOT, scenario, sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) { x <- tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(x)) return(x) }
  rows <- comp_rep(scenario, rep)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(f, ".tmp"); saveRDS(rows, tmp); if (file.exists(f)) file.remove(f); file.rename(tmp, f)
  rows
}

# True marginal subdistribution HR of the renal event (leaving competing), first year: on
# the uncensored counterfactual cohorts a patient who leaves first stays in the risk set.
sub_truth <- function(scenario) {
  p <- get_scenario_params(scenario); s <- 1 * 1000 + 1
  d1 <- sim_rwe(2.5e5, p, arm = "tdf", censoring = FALSE, seed = s)$data
  d0 <- sim_rwe(2.5e5, p, arm = "etv", censoring = FALSE, seed = s)$data
  st <- rbind(d1, d0)
  ev <- st$T_nat <= pmin(st$leave, H)
  dd <- data.frame(A = st$A, time = ifelse(ev, st$T_nat, H), status = as.integer(ev))
  data.frame(scenario = scenario,
             true_shr = exp(coef(coxph(Surv(time, status) ~ A, data = dd, ties = "breslow"))[["A"]]))
}

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "full"
if (mode == "pilot") {
  x <- comp_rep("base", 1)
  print(x[, c("strategy", "method", "window", "est", "lo", "hi", "rr", "rr_lo", "rr_hi", "ok", "error", "secs")], digits = 4)
  print(get_truth("base")$strategies[, c("strategy", "rd", "rr", "hr_marginal")], digits = 4)
} else if (mode == "truth") {
  tr <- do.call(rbind, lapply(SCENARIOS, sub_truth))
  saveRDS(tr, "sandbox_rwe/results/truth_subdistribution_hr.rds"); print(tr, digits = 4)
} else {
  cores <- as.integer(Sys.getenv("RWE_CORES", "4"))   # workers; each peaks at a few GB (memory, not CPU, is the limit)
  todo <- expand.grid(rep = 1:200, scenario = SCENARIOS, stringsAsFactors = FALSE)
  todo <- todo[!file.exists(file.path(ROOT, todo$scenario, sprintf("rep_%03d.rds", todo$rep))), ]
  cat(format(Sys.time(), "%m-%d %H:%M"), "-", nrow(todo), "replicates to run\n")
  for (ch in split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / (2 * cores)))) {
    cl <- parallel::makeCluster(cores)
    tryCatch({
      parallel::clusterEvalQ(cl, {
        setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
        for (f in list.files("sandbox_rwe/R", pattern = "\\.R$", full.names = TRUE)) source(f)
        suppressPackageStartupMessages({ library(survival); library(riskRegression); library(prodlim); library(rstpm2) })
      })
      parallel::clusterExport(cl, c("comp_rep", "cached", "ROOT", "H", "TSTAR", ".row", ".fail", "std_cox",
                                    "std_rp", "cr_data", "std_csc", "fine_gray"))
      parallel::parLapply(cl, ch, function(i, todo) tryCatch(cached(todo$scenario[i], todo$rep[i]),
                                                             error = function(e) NULL), todo = todo)
    }, finally = parallel::stopCluster(cl))
    cat(format(Sys.time(), "%m-%d %H:%M"), "- done", max(ch), "of", nrow(todo), "\n")
  }
  rows <- do.call(rbind, lapply(seq_len(nrow(expand.grid(rep = 1:200, scenario = SCENARIOS))), function(i) {
    g <- expand.grid(rep = 1:200, scenario = SCENARIOS, stringsAsFactors = FALSE)
    f <- file.path(ROOT, g$scenario[i], sprintf("rep_%03d.rds", g$rep[i])); if (file.exists(f)) readRDS(f)
  }))
  saveRDS(rows, "sandbox_rwe/results/comparators_rows.rds")
  cat("rows:", nrow(rows), " failed:", sum(!rows$ok %in% TRUE), "\n")
}
