# Prompt for Claude Code on the desktop: job 5

In the existing clone of the `desktop-comparators` branch (from jobs 1-4), pull the new code,
start Claude Code in that folder, and say "Follow DESKTOP_PROMPT.md":

```
cd estimand_comparators
git pull origin desktop-comparators
claude
```

(If there is no clone yet: `git clone --branch desktop-comparators --single-branch
https://github.com/amertens/estimand_tutorial.git estimand_comparators`.)

---

Run job 5 in this repository: three parts, one after another, and report back after each.
Do not edit any code or input file; the only new files are the ones the runs write.

**Setup.** As for jobs 1-4: R 4.4.x, always `Rscript --vanilla`, packages lmtp, survival,
riskRegression, prodlim, dplyr, digest, SuperLearner, nnls. Set `RWE_ROOT` to this folder's
absolute path with forward slashes in the shell that runs R. Launch
`sandbox_rwe/run/no_ecoqos.ps1 -Watch` alongside each full run if Windows throttles background R.

**Memory rule.** As before: set `RWE_CORES` to the smaller of (free RAM in GB) / (1.5 x the pilot's
peak GB) and (physical cores - 2), at most 12; check memory after the first chunk and lower
`RWE_CORES` and restart (checkpoints are kept per cohort) if free RAM drops below about 4 GB.

**Pilot (once).** `Rscript --vanilla sandbox_rwe/run/75_job5.R pilot` (about 10 minutes). Expect three
`points` rows with `ok = TRUE` (risks near 0.054/0.038 and 0.063/0.041), two bootstrap rows with
`n_boot_ok = 3`, and twelve rows for the informative loss-to-follow-up cohort, all `ok = TRUE`.

**Part A: points** (6 scenarios x 200 cohorts, about 25 s per cohort per worker):
`Rscript --vanilla sandbox_rwe/run/75_job5.R points > sandbox_rwe/results/job5_points.log 2>&1`
Output: `sandbox_rwe/results/job5_points_rows.rds`.

**Part B: ltfu** (computes the new scenario's true values first, about an hour in one process, then
200 cohorts with SDR, about 7 minutes per cohort per worker; SDR needs more memory, apply the rule):
`Rscript --vanilla sandbox_rwe/run/75_job5.R ltfu > sandbox_rwe/results/job5_ltfu.log 2>&1`
Output: `sandbox_rwe/results/informative_ltfu_strong_rows.rds`, `truth_d2f_informative_ltfu_strong.rds`,
`params_d2f_informative_ltfu_strong.rds`. The printed `rd` column should match `base_rd` closely.

**Part C: boot** (200 base-case cohorts, 100 bootstrap samples each, about 19 minutes per cohort per
worker, so about 5-6 hours with 12 workers; Cox models only, small memory):
`Rscript --vanilla sandbox_rwe/run/75_job5.R boot > sandbox_rwe/results/job5_boot.log 2>&1`
Output: `sandbox_rwe/results/job5_boot_rows.rds`.

**After each part: report and push.** Report the run time, `RWE_CORES`, the final "rows / failed"
line, and any failed rows with their error messages. Then commit only that part's output files
listed above and its log, and push: `git push origin desktop-comparators`. Do not commit the
checkpoint folders under `sandbox_rwe/results/runs/` (they are gitignored).
