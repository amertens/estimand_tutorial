# rerun_sim_extra140.R
# Run iterations 61..200 for each of the five estimands (140 new iters per
# estimand) using the same DGP and estimator pipeline as the 60-iter run.
#
# Chunked caching: each estimand's iterations are run in chunks of CHUNK_SIZE,
# and each chunk is saved to disk immediately. If the script crashes or is
# interrupted, restarting it picks up at the next un-cached chunk. Per-chunk
# files live under results/sim_partial_extra/<estimand>/iter<NNN>-<MMM>.rds.
#
# When all chunks complete, the script concatenates the original 60-iter
# results (results/sim_study_main_60iter.rds) with the new 140-iter results
# to give a 200-iter combined file at results/sim_study_main.rds.
# The original 60-iter cache is preserved at results/sim_study_main_60iter.rds.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(survival); library(broom)
  library(parallel); library(lmtp); library(SuperLearner); library(arm)
})
source(here("DGP.R"))
source(here("R", "helpers.R"))
source(here("R", "run_simulations.R"))

tau <- 180
BIN_WIDTH <- 14
ITER_START <- 61
ITER_END <- 200
N_CORES <- 4
CHUNK_SIZE <- 10  # save after every CHUNK_SIZE iterations

dgp_args <- list(
  h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, dep_censor = TRUE, complexity = FALSE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60,
  return_potential_switching = TRUE
)

ESTIMANDS <- c("treatment_policy", "no_switch", "while_on_treatment",
               "composite", "principal_stratum")

partial_dir <- here("results", "sim_partial_extra")
dir.create(partial_dir, showWarnings = FALSE, recursive = TRUE)

# Build the chunk index: list of (start, end) pairs covering ITER_START..ITER_END
chunks <- list()
i <- ITER_START
while (i <= ITER_END) {
  j <- min(i + CHUNK_SIZE - 1, ITER_END)
  chunks[[length(chunks) + 1]] <- c(i, j)
  i <- j + 1
}

chunk_filename <- function(est_dir, chunk) {
  file.path(est_dir, sprintf("iter%03d-%03d.rds", chunk[1], chunk[2]))
}

run_chunk <- function(est, iters, n_cores, dgp_args, tau, bin_width) {
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(dplyr); library(survival); library(broom)
      library(here); library(lmtp); library(SuperLearner); library(arm)
    })
    source(here("DGP.R"))
    source(here("R", "helpers.R"))
  })
  worker_args <- list(
    est = est, sample_size = 10000, tau = tau,
    bin_width = bin_width,
    lmtp_learners = c("SL.mean", "SL.glm", "SL.bayesglm"),
    run_lmtp = TRUE, run_cox_td = FALSE,
    dgp_args = dgp_args
  )
  clusterExport(cl, c("run_one_iter", "worker_args"), envir = environment())

  iter_results <- clusterApply(cl, iters, function(i) {
    run_one_iter(
      i, estimand = worker_args$est,
      sample_size = worker_args$sample_size,
      tau = worker_args$tau,
      bin_width = worker_args$bin_width,
      lmtp_learners = worker_args$lmtp_learners,
      run_lmtp = worker_args$run_lmtp,
      run_cox_td = worker_args$run_cox_td,
      dgp_args = worker_args$dgp_args
    )
  })
  dplyr::bind_rows(iter_results)
}

cat(sprintf("Running iterations %d..%d in chunks of %d on %d cores\n",
            ITER_START, ITER_END, CHUNK_SIZE, N_CORES))
cat(sprintf("Total chunks per estimand: %d\n", length(chunks)))

for (est in ESTIMANDS) {
  est_dir <- file.path(partial_dir, est)
  dir.create(est_dir, showWarnings = FALSE, recursive = TRUE)

  est_t0 <- Sys.time()
  est_completed_chunks <- 0
  cat(sprintf("\n=== %s ===\n", est))

  for (chunk in chunks) {
    cf <- chunk_filename(est_dir, chunk)
    if (file.exists(cf)) {
      cat(sprintf("  [skip] iter %d..%d cached\n", chunk[1], chunk[2]))
      est_completed_chunks <- est_completed_chunks + 1
      next
    }
    iters <- seq(chunk[1], chunk[2])
    t0 <- Sys.time()
    res <- run_chunk(est, iters, N_CORES, dgp_args, tau, BIN_WIDTH)
    saveRDS(res, cf)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    est_completed_chunks <- est_completed_chunks + 1
    cat(sprintf("  iter %d..%d: %d rows in %.1f min (%d/%d chunks done)\n",
                chunk[1], chunk[2], nrow(res), elapsed,
                est_completed_chunks, length(chunks)))
  }
  est_elapsed <- round(as.numeric(difftime(Sys.time(), est_t0, units = "mins")), 1)
  cat(sprintf("  %s complete in %.1f min\n", est, est_elapsed))
}

# ── Concatenate original 60-iter + new 140-iter into 200-iter combined ──
cat("\nConcatenating 60-iter + 140-iter into 200-iter combined file...\n")
old_60 <- readRDS(here("results", "sim_study_main_60iter.rds"))

new_140 <- bind_rows(lapply(ESTIMANDS, function(e) {
  est_dir <- file.path(partial_dir, e)
  cfs <- list.files(est_dir, pattern = "^iter\\d{3}-\\d{3}\\.rds$",
                    full.names = TRUE)
  bind_rows(lapply(cfs, readRDS))
}))

combined <- bind_rows(old_60, new_140)
saveRDS(combined, here("results", "sim_study_main.rds"))
cat(sprintf("Wrote %d rows (60+140 = 200 iter, %d estimands) to results/sim_study_main.rds\n",
            nrow(combined), length(ESTIMANDS)))
