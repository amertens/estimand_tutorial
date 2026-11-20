# Prompt for Claude Code on the desktop

Clone the `desktop-comparators` branch, start Claude Code in that folder, and say
"Follow DESKTOP_PROMPT.md" (the same text sits at the branch root as DESKTOP_PROMPT.md):

```
git clone --branch desktop-comparators --single-branch https://github.com/amertens/estimand_tutorial.git estimand_comparators
cd estimand_comparators
claude
```

---

Run four simulation jobs in this repository, one after another, and report back after each.
Do not edit any code or input file; the only new files are the ones the runs write.

**Setup (once).** Confirm R 4.4.x is installed and always use `Rscript --vanilla` (a user profile
can crash R). Install any of these packages that are missing: lmtp, survival, riskRegression,
prodlim, rstpm2, dplyr, digest, SuperLearner. Set the environment variable `RWE_ROOT` to this
folder's absolute path with forward slashes (for example `C:/Users/<you>/estimand_comparators`) in
the shell that runs R, so every R process and worker inherits it.

**Memory rule for every job.** On a laptop, 8 workers ran out of memory and crashed the machine.
Before each full run, run its pilot and watch the peak memory of the Rscript process (Task Manager,
or `Get-Process Rscript` in PowerShell). Set `RWE_CORES` to the smaller of (free RAM in GB) / (1.5 x
the pilot's peak GB) and (physical cores - 2), and at most 12. Check memory after the first chunk
of the full run; if free RAM drops below about 4 GB, stop it, lower `RWE_CORES` and restart with the
same command (every job saves one checkpoint per cohort and skips cohorts already done).

**Job 1: Cox-family comparators** (`sandbox_rwe/run/65_comparators.R`).
- Pilot: `Rscript --vanilla sandbox_rwe/run/65_comparators.R pilot` (2-3 minutes). Expect six rows,
  all `ok = TRUE`, and treatment-policy RDs near 0.016-0.020 against a true RD of 0.019.
- Full run, in the background:
  `Rscript --vanilla sandbox_rwe/run/65_comparators.R > sandbox_rwe/results/comparators.log 2>&1`
  (6 scenarios x 200 cohorts; the cause-specific model takes about a minute per cohort).
- Output: `sandbox_rwe/results/comparators_rows.rds`.

**Job 2: informative loss to follow-up** (`sandbox_rwe/run/69_informative_ltfu.R`).
- Pilot: `Rscript --vanilla sandbox_rwe/run/69_informative_ltfu.R pilot`. It first computes the
  scenario's true values (about an hour, one process); the printed `rd` column should match
  `base_rd` closely, because loss to follow-up does not change the true values. Then it runs one
  cohort and prints its rows; expect all `ok = TRUE`.
- Full run, in the background:
  `Rscript --vanilla sandbox_rwe/run/69_informative_ltfu.R > sandbox_rwe/results/informative_ltfu.log 2>&1`
  (200 cohorts with SDR; several hours).
- Output: `sandbox_rwe/results/informative_ltfu_rows.rds`, plus the new
  `sandbox_rwe/results/params_d2e_informative_ltfu.rds` and `truth_d2e_informative_ltfu.rds`.

**Job 3: principal-score weighting with 1,000 bootstrap samples** (`sandbox_rwe/run/70_psw_boot1000.R`).
- Pilot: `Rscript --vanilla sandbox_rwe/run/70_psw_boot1000.R pilot` (about 11 minutes). Expect
  `ok = TRUE` and "point estimate reproduces the main run: TRUE".
- Full run, in the background:
  `Rscript --vanilla sandbox_rwe/run/70_psw_boot1000.R > sandbox_rwe/results/psw_boot1000.log 2>&1`
  (800 cohorts, about 11 minutes each per worker, so roughly 15 hours with 10 workers; Cox models only, so memory per worker is small).
- Output: `sandbox_rwe/results/psw_boot1000_rows.rds`.

**Job 4: SDR for while on treatment and the composite in the five one-feature scenarios** (`sandbox_rwe/run/72_wot_comp_scenarios.R`).
- Pilot: `Rscript --vanilla sandbox_rwe/run/72_wot_comp_scenarios.R pilot` (one cohort of the common-outcome
  scenario). Expect two rows, both `ok = TRUE`, with RDs near the printed true values.
- Full run, in the background:
  `Rscript --vanilla sandbox_rwe/run/72_wot_comp_scenarios.R > sandbox_rwe/results/wotcomp.log 2>&1`
  (5 scenarios x 200 cohorts, two SDR fits per cohort, about 1 minute per cohort per worker; SDR uses more memory than Cox, so apply the memory rule).
- Output: `sandbox_rwe/results/wotcomp_rows.rds`.

**After each job: report and push.** Report the run time, `RWE_CORES`, the final "rows / failed"
line, and any failed rows by method with their error messages. Then commit only that job's output
files listed above and its log, and push: `git push origin desktop-comparators`. Do not commit the
checkpoint folders under `sandbox_rwe/results/runs/` (they are gitignored).
