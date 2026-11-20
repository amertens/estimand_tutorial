# run.R -- one replicate of one scenario: simulate, run an estimator set, return
# rows in one long schema; checkpoint each replicate with the parameter
# fingerprint and the code version.

ROW_COLS <- c("scenario", "rep", "set", "family", "strategy", "method", "window", "delta",
              "est", "lo", "hi", "rr", "rr_lo", "rr_hi", "share", "max_w", "ok", "error",
              "secs", "fingerprint")

seed_for <- function(scenario, rep) {
  u <- utf8ToInt(scenario)
  as.integer((sum(u * seq_along(u)) %% 100000L) * 1000L + rep)
}

# Each estimator gets its own seed (cross-fitting folds, stochastic shifts), so it
# gives the same result alone as within its set.
estimator_seed <- function(scenario, rep, name) {
  u <- utf8ToInt(name)
  as.integer((seed_for(scenario, rep) * 31 + sum(u * seq_along(u))) %% .Machine$integer.max)
}

# Version of the code that produces rows: a checkpoint from different code is
# recomputed. Scenario parameters and truths carry their own fingerprints, so
# scenarios, calibration, truths and summaries can change without invalidating
# replicates.
code_version_files <- function()
  rwe_path("sandbox_rwe/R", c("analysis.R", "dgp.R", "principal_score.R", "run.R"))
code_version <- function() digest::digest(lapply(code_version_files(), readLines), algo = "md5")

.blank_row <- function(meta) {
  r <- as.list(setNames(rep(NA, length(ROW_COLS)), ROW_COLS))
  for (k in names(meta)) r[[k]] <- meta[[k]]
  as.data.frame(r, stringsAsFactors = FALSE)
}

safe_rows <- function(fun, meta) {
  t0 <- Sys.time()
  out <- tryCatch(fun(), error = function(e) e)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(out, "error")) {
    r <- .blank_row(meta); r$ok <- FALSE; r$error <- conditionMessage(out); r$secs <- secs
    return(r)
  }
  rows <- do.call(rbind, lapply(seq_len(nrow(out)), function(i) {
    r <- .blank_row(meta)
    for (k in intersect(names(out), ROW_COLS)) r[[k]] <- out[[k]][i]
    r
  }))
  if (!"ok" %in% names(out)) rows$ok <- TRUE        # estimators may report per-row status
  bad <- rows$ok %in% TRUE & !is.finite(suppressWarnings(as.numeric(rows$est)))
  rows$ok[bad] <- FALSE; rows$error[bad] <- "estimate not finite"
  rows$secs <- secs
  rows
}

# An estimator: the function and the metadata that identifies its rows, even when it fails.
.e <- function(fun, family, strategy = NA, method = NA, delta = NA)
  list(fun = fun, meta = list(family = family, strategy = strategy, method = method, delta = delta))

# An unstable fit can return a risk outside [0, 1]; that row is a failure, not a result.
.sdr_row <- function(x, strategy, method, delta = NA, share = NA) {
  in01 <- all(is.finite(c(x[["risk1"]], x[["risk0"]]))) && x[["risk1"]] >= 0 && x[["risk1"]] <= 1 &&
    x[["risk0"]] >= 0 && x[["risk0"]] <= 1
  data.frame(strategy = strategy, method = method, delta = delta, est = x[["rd"]], lo = x[["rd_lo"]],
             hi = x[["rd_hi"]], rr = x[["rr"]], rr_lo = x[["rr_lo"]], rr_hi = x[["rr_hi"]],
             share = share, max_w = x[["max_w"]], ok = in01,
             error = if (in01) NA_character_ else "risk outside [0, 1]")
}

.cox_fun <- function(sim) function() {
  cx <- suppressWarnings(cox_suite(sim))
  data.frame(method = cx$analysis, window = cx$window, est = cx$hr, lo = cx$lo, hi = cx$hi,
             ok = cx$ok, error = cx$error)
}
.sdr_fun <- function(sim, strategy, method = "SDR", rows = NULL, ...) {
  force(strategy); force(method)
  share <- if (is.null(rows)) NA else mean(rows)
  function() .sdr_row(run_sdr(prep_lmtp(sim, strategy, rows = rows), ...), strategy, method,
                      share = share)
}
# Reduced leaving at delta, or, with delta = 0 (no departures), sustained use. Sustained
# use is an intervention on leaving rather than censoring at leaving: censoring
# needed copies of the treatment node, whose classifier missed the eGFR-driven
# switching and biased the TDF risk downward (pilot iterations 8 and 11).
.reduced_fun <- function(sim, delta, strategy = "reduced_leaving", method = "SDR", ...) {
  force(delta); force(strategy); force(method)
  function() .sdr_row(run_sdr(prep_lmtp(sim, "reduced_leaving"), shift1 = reduced_leave_shift(1, delta),
                              shift0 = reduced_leave_shift(0, delta), ...),
                      strategy, method, if (strategy == "reduced_leaving") delta else NA)
}
.su_e <- function(sim, method = "SDR", ...)
  .e(.reduced_fun(sim, 0, "sustained_use", method, ...), "sdr", "sustained_use", method)
.sdr_e <- function(sim, strategy, method = "SDR", rows = NULL, ...)
  .e(.sdr_fun(sim, strategy, method, rows, ...), "sdr", strategy, method)

CTL_UNTRIMMED <- lmtp_control(.trim = 1, .learners_outcome_folds = 2, .learners_trt_folds = 2)
CTL_FLEX <- lmtp_control(.trim = 0.9999)                      # 10 Super Learner folds

# The three principal-stratum estimators: SDR in the true stratum (oracle), SDR in
# observed adherers, principal-score weighting. Used by the base set and by the
# "principal" set (renal-driven switching, where principal ignorability is strained).
.ps_ests <- function(sim) {
  d <- sim$data
  oracle <- d$leave1 > 365 & d$leave0 > 365
  adh <- d$leave >= pmin(d$follow, 365)
  list(
    ps_oracle = .sdr_e(sim, "principal_stratum", "SDR, true stratum", rows = oracle),
    ps_adh = .sdr_e(sim, "principal_stratum", "SDR, observed adherers", rows = adh),
    psw = .e(function() {
      x <- ps_weighting(sim, B = 100)
      data.frame(est = x[["rd"]], lo = x[["rd_lo"]], hi = x[["rd_hi"]], rr = x[["rr"]],
                 rr_lo = x[["rr_lo"]], rr_hi = x[["rr_hi"]], share = x[["share"]])
    }, "psw", "principal_stratum", "Principal-score weighting"))
}

estimator_set <- function(set, sim) {
  d <- sim$data
  probe_cox <- .e(function() {
    dd <- with_stages(d); dd <- cbind(dd[, c("id", "A", COX_X)], time = dd$follow, status = dd$event)
    h <- suppressWarnings(cox_hr(dd))
    data.frame(method = "itt", window = "all", est = h[["hr"]], lo = h[["lo"]], hi = h[["hi"]])
  }, "cox")
  switch(set,
    probe = list(cox = probe_cox),
    probe_fail = list(cox = probe_cox,
                      tp = .e(function() stop("forced failure"), "sdr", "treatment_policy", "SDR")),
    probe_rng = list(a = .e(function() data.frame(est = stats::runif(1)), "sdr", "x", "a"),
                     b = .e(function() data.frame(est = stats::runif(1)), "sdr", "x", "b")),
    deviation = list(
      cox = .e(.cox_fun(sim), "cox"),
      tp = .sdr_e(sim, "treatment_policy"),
      su = .su_e(sim)),
    base = c(list(
      cox = .e(.cox_fun(sim), "cox"),
      tp = .sdr_e(sim, "treatment_policy"),
      su = .su_e(sim),
      su_untrimmed = .su_e(sim, "SDR, untrimmed", control = CTL_UNTRIMMED),
      wot = .sdr_e(sim, "while_on_treatment"),
      comp = .sdr_e(sim, "composite")),
      .ps_ests(sim)[c("ps_oracle", "ps_adh")],
      list(reduced = .e(.reduced_fun(sim, 0.5), "sdr", "reduced_leaving", "SDR", 0.5)),
      .ps_ests(sim)["psw"]),
    principal = .ps_ests(sim),
    # SDR for while on treatment and the composite in the one-feature scenarios (run/72);
    # same estimator names, so the same per-estimator seeds as the base set
    wotcomp = list(
      wot = .sdr_e(sim, "while_on_treatment"),
      comp = .sdr_e(sim, "composite")),
    grid = setNames(lapply(c(0.75, 0.25, 0.1), function(dl)
      .e(.reduced_fun(sim, dl), "sdr", "reduced_leaving", "SDR", dl)),
      paste0("reduced_", c(0.75, 0.25, 0.1))),
    flex = list(
      tp = .sdr_e(sim, "treatment_policy", "SDR, flexible library", learners = LIB_FLEX, control = CTL_FLEX),
      su = .su_e(sim, "SDR, flexible library", learners = LIB_FLEX, control = CTL_FLEX)),
    stop("unknown estimator set: ", set))
}

run_replicate <- function(scenario, rep, p, set, only = NULL) {
  n <- if (startsWith(set, "probe")) 2000 else 20000
  sim <- sim_rwe(n, p, seed = seed_for(scenario, rep))
  ests <- estimator_set(set, sim)
  if (!is.null(only)) ests <- ests[only]
  rows <- do.call(rbind, lapply(names(ests), function(nm) {
    set.seed(estimator_seed(scenario, rep, nm))
    safe_rows(ests[[nm]]$fun, ests[[nm]]$meta)
  }))
  rows$scenario <- scenario; rows$rep <- rep; rows$set <- set
  rows$fingerprint <- params_fingerprint(p)
  rows[, ROW_COLS]
}

# cv: the code version the run started with (run_scenario passes it, so every
# checkpoint of one run carries one version)
run_replicate_cached <- function(scenario, rep, set, root = rwe_path("sandbox_rwe/results/runs"),
                                 cv = code_version()) {
  p <- get_scenario_params(scenario)
  fp <- params_fingerprint(p)
  f <- file.path(root, scenario, set, sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)       # unreadable = missing
    if (!is.null(x) && identical(x$fingerprint, fp) && identical(x$code_version, cv)) return(x$rows)
  }
  rows <- run_replicate(scenario, rep, p, set)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(f, ".tmp")                                     # write, then move into place
  saveRDS(list(fingerprint = fp, code_version = cv, rows = rows), tmp)
  if (file.exists(f)) file.remove(f)
  file.rename(tmp, f)
  rows
}

.failed_replicate <- function(scenario, rep, set, e) {
  x <- .blank_row(list(scenario = scenario, rep = rep, set = set))
  x$ok <- FALSE; x$error <- conditionMessage(e)
  x[, ROW_COLS]
}

# Replicates run in chunks with a fresh cluster per chunk (lmtp workers grow in
# memory); an error in one replicate gives a failed row and the run continues.
# The code version is fixed when the run starts: a change to the row-producing code
# before a chunk stops the run, and a worker whose sourced code differs records a
# failed replicate (final review C2).
run_scenario <- function(scenario, reps, set, cores = 6, root = rwe_path("sandbox_rwe/results/runs"),
                         chunk = 2 * cores, cv_fun = code_version) {
  get_scenario_params(scenario)                    # build the parameter cache once
  cv <- cv_fun()
  one <- function(r, s, st, rt, cv) tryCatch({
    if (!identical(code_version(), cv)) stop("code changed since the run started (", cv, ")")
    run_replicate_cached(s, r, st, rt, cv)
  }, error = function(e) .failed_replicate(s, r, st, e))
  if (cores <= 1) return(do.call(rbind, lapply(reps, one, s = scenario, st = set, rt = root, cv = cv)))
  out <- list()
  for (ch in split(reps, ceiling(seq_along(reps) / chunk))) {
    if (!identical(cv_fun(), cv))
      stop("row-producing code changed during the run (started as ", cv, "); restart it")
    cl <- parallel::makeCluster(cores)
    res <- tryCatch({
      parallel::clusterEvalQ(cl, {
        setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
        for (f in list.files("sandbox_rwe/R", pattern = "\\.R$", full.names = TRUE)) source(f)
      })
      parallel::parLapply(cl, ch, one, s = scenario, st = set, rt = root, cv = cv)
    }, finally = parallel::stopCluster(cl))
    out <- c(out, res)
  }
  do.call(rbind, out)
}
