# Prompt for Claude Code on the desktop

Start Claude Code in this folder (a clone of the `desktop-comparators` branch) and say:
"Follow DESKTOP_PROMPT.md." Background and details are in `sandbox_rwe/run/DESKTOP_RUN.md`.

---

Run the comparator simulation in this repository and report back. Do not edit any code or input
file; the only new files are the ones the run writes.

1. **Setup.** Confirm R 4.4.x is installed (use `Rscript --vanilla`; a user profile can crash R).
   Install any of these packages that are missing: lmtp, survival, riskRegression, prodlim, rstpm2,
   dplyr, digest, SuperLearner. Set the environment variable `RWE_ROOT` to this folder's absolute path
   with forward slashes (for example `C:/Users/<you>/estimand_comparators`), in the shell that runs R,
   so that every R process and worker inherits it.

2. **Pilot and memory check.** Run
   `Rscript --vanilla sandbox_rwe/run/65_comparators.R pilot`
   and watch the peak memory of the Rscript process (Task Manager, or `Get-Process Rscript` in
   PowerShell). The pilot fits six analyses to one cohort of 20,000 patients and takes 2-3 minutes.
   Expected: six rows, all `ok = TRUE`, and treatment-policy RDs near 0.016-0.020 against a true RD of
   0.019. Stop and report if any row fails.

3. **Choose the number of workers.** On a laptop, 8 workers ran out of memory and crashed the
   machine. Set `RWE_CORES` to the smaller of (free RAM in GB) / (1.5 x the pilot's peak GB) and
   (physical cores - 2), and at most 12.

4. **Full run.** In the background, run
   `Rscript --vanilla sandbox_rwe/run/65_comparators.R > sandbox_rwe/results/comparators.log 2>&1`
   It runs 6 scenarios x 200 cohorts = 1,200 cohorts, saves one checkpoint per cohort under
   `sandbox_rwe/results/runs/comparators/`, and skips cohorts already done, so the same command
   restarts it after any interruption. Check memory after the first chunk; if free RAM drops below
   about 4 GB, stop it, lower `RWE_CORES` and restart. At the end it writes
   `sandbox_rwe/results/comparators_rows.rds` and prints the number of rows and failures.

5. **Report and push.** Report the run time, `RWE_CORES`, the rows and failures line, and the failed
   rows by method with their error messages. Then commit only
   `sandbox_rwe/results/comparators_rows.rds` and `sandbox_rwe/results/comparators.log` to this
   branch and push it (`git push origin desktop-comparators`).
