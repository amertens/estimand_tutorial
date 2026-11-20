# 70_psw_boot1000.R -- principal-score weighting with 1,000 bootstrap samples instead of 100
# (review 2026-11-19, item A15): its coverage of 91.8-93.5% may partly reflect bootstrap error.
# Same cohorts and the same estimator seed ("psw") as the main runs, so the point estimates
# reproduce full_rows.rds and supp_s1_rows.rds and only the percentile intervals change.
# Scenarios: base (200 cohorts), renal-driven switching (200), susceptibility to the
# nephrotoxicity of Drug A (400), as in the main runs. One checkpoint per cohort.
#
#   Rscript --vanilla sandbox_rwe/run/70_psw_boot1000.R pilot   # 1 cohort, timed, checked
#   Rscript --vanilla sandbox_rwe/run/70_psw_boot1000.R         # 800 cohorts, RWE_CORES workers
setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
for (f in list.files("sandbox_rwe/R", pattern = "[.]R$", full.names = TRUE)) source(f)
B <- 1000
ROOT <- "sandbox_rwe/results/runs/psw1000"
cores <- as.integer(Sys.getenv("RWE_CORES", "4"))

psw_rep <- function(s, r) {
  p <- get_scenario_params(s)
  sim <- sim_rwe(20000, p, seed = seed_for(s, r))
  set.seed(estimator_seed(s, r, "psw"))
  t0 <- Sys.time()
  x <- tryCatch(ps_weighting(sim, B = B), error = function(e) NULL)
  if (is.null(x)) return(data.frame(scenario = s, rep = r, ok = FALSE))
  data.frame(scenario = s, rep = r, t(x), ok = TRUE, secs = as.numeric(Sys.time() - t0, units = "secs"),
             fingerprint = params_fingerprint(p))
}
cached <- function(s, r) {
  f <- file.path(ROOT, s, sprintf("rep_%03d.rds", r))
  if (file.exists(f)) { x <- tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(x)) return(x) }
  x <- psw_rep(s, r)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(f, ".tmp"); saveRDS(x, tmp); if (file.exists(f)) file.remove(f); file.rename(tmp, f)
  x
}
jobs <- rbind(data.frame(s = "base", r = 1:200), data.frame(s = "d2c_renal_switch", r = 1:200),
              data.frame(s = "s1_tdf_susceptibility", r = 1:400))

if (length(commandArgs(trailingOnly = TRUE)) && commandArgs(trailingOnly = TRUE)[1] == "pilot") {
  x <- psw_rep("base", 1)
  print(x, digits = 4)
  main <- readRDS("sandbox_rwe/results/full_rows.rds")
  m <- main[main$scenario == "base" & main$rep == 1 & main$method == "Principal-score weighting", ]
  if (nrow(m)) cat("point estimate reproduces the main run:", isTRUE(all.equal(x$rd, as.numeric(m$est[1]))), "\n")
} else {
  todo <- jobs[!file.exists(file.path(ROOT, jobs$s, sprintf("rep_%03d.rds", jobs$r))), ]
  cat(format(Sys.time(), "%m-%d %H:%M"), "-", nrow(todo), "cohorts to run\n")
  for (ch in split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / (2 * cores)))) {
    cl <- parallel::makeCluster(cores)
    tryCatch({
      parallel::clusterEvalQ(cl, {
        setwd(Sys.getenv("RWE_ROOT", "C:/Users/andre/OneDrive/Documents/estimand_tutorial"))
        for (f in list.files("sandbox_rwe/R", pattern = "[.]R$", full.names = TRUE)) source(f)
      })
      parallel::clusterExport(cl, c("psw_rep", "cached", "ROOT", "B"))
      parallel::parLapply(cl, ch, function(i, todo) tryCatch(cached(todo$s[i], todo$r[i]), error = function(e) NULL),
                          todo = todo)
    }, finally = parallel::stopCluster(cl))
    cat(format(Sys.time(), "%m-%d %H:%M"), "- done", max(ch), "of", nrow(todo), "\n")
  }
  rows <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
    f <- file.path(ROOT, jobs$s[i], sprintf("rep_%03d.rds", jobs$r[i])); if (file.exists(f)) readRDS(f)
  }))
  saveRDS(rows, "sandbox_rwe/results/psw_boot1000_rows.rds")
  cat("rows:", nrow(rows), " failed:", sum(!rows$ok %in% TRUE), "\n")
}
