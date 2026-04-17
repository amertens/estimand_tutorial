# Estimand–Estimator Mismatch in Time-to-Event Safety Analyses

Simulation study and manuscript comparing Cox-based analyses (naive, censor-at-switch, IPCW) with LMTP SDR across five ICH E9(R1) estimand strategies, using a hepatitis B renal safety scenario with informative treatment switching.

## Repository layout

### Top-level documents
- [`draft-estimands-results-outline.qmd`](draft-estimands-results-outline.qmd) — main manuscript (Quarto → Word). Renders from cached simulation results.
- [`presentation-estimands.qmd`](presentation-estimands.qmd) — RevealJS slide deck. Reads the same caches as the manuscript.
- [`tutorial-overview.qmd`](tutorial-overview.qmd) — plain-language companion document.
- [`references.bib`](references.bib) / [`vancouver-superscript.csl`](vancouver-superscript.csl) — citation database and style.

### Core simulation code
- [`DGP.R`](DGP.R) — `generate_hep_data()` data-generating process and `derive_estimand()` outcome transformation.
- [`calc_truth.R`](calc_truth.R) — Monte Carlo ground-truth functions for each of the five estimands (N = 500,000 per arm).
- [`R/helpers.R`](R/helpers.R) — Cox wrappers (`fit_cox_naive`, `fit_cox_censor_switch`, `fit_cox_ipcw`, `fit_cox_td`), LMTP data prep (`prepare_lmtp_data`), LMTP estimation wrapper (`run_lmtp_analysis`), and bootstrap KM-RD confidence-interval utility (`boot_km_rd_ci`).
- [`R/run_simulations.R`](R/run_simulations.R) — `run_one_iter()` applies all estimators to one dataset; `run_simulation_study()` orchestrates repeated iterations across estimands via parallel `clusterApply`.
- [`R/support_diagnostics.R`](R/support_diagnostics.R) — positivity diagnostic plots used in the manuscript.

### Entry-point scripts
- [`R/run_cached_analyses.R`](R/run_cached_analyses.R) — **the main entry point**. Generates everything the manuscript and presentation consume: ground truth, single-dataset LMTP results, support scenarios, and the 60-iteration repeated simulation study. Skips cached steps.
- [`R/run_sim_test.R`](R/run_sim_test.R) — 5-iteration smoke test (separate cache, `sim_study_test.rds`) for verifying the pipeline works without committing ~5 hours of compute.
- [`R/run_wot_rich.R`](R/run_wot_rich.R) — WOT sensitivity analysis with richer SuperLearner libraries (`SL.glmnet`, `SL.ranger`); tests whether a stronger censoring model improves the WOT LMTP coverage beyond the 48% achieved with the minimal library.
- [`R/eval_full_table.R`](R/eval_full_table.R) — prints a bias / RD-coverage / HR-coverage table across all estimator × estimand cells. Prefers `sim_study_main.rds`; falls back to `sim_study_test.rds`.

### Outputs
- `results/sim_results/ground_truth.rds` — true risk difference, risk ratio, and marginal hazard ratio for each estimand.
- `results/lmtp_main.rds` — single-dataset LMTP treatment-policy fit (used in `tbl-lmtp` in the manuscript).
- `results/support_estimation.rds` — Cox + LMTP estimates under three support scenarios (good / strained / poor).
- `results/sim_study_main.rds` — 60-iteration simulation results (4 estimators × 5 estimands). Contains `hr`, `hr_ci_low`, `hr_ci_high`, `risk_diff`, `ci_low`, `ci_high`, `converged` per row.
- `results/sim_study_test.rds` — 5-iteration test cache (created by `run_sim_test.R`).
- `results/wot_rich_sl.rds` — WOT rich-learner comparison results (created by `run_wot_rich.R`).

### Archive
- [`archive/`](archive/) — historical variants (old presentation, alternate CSL, salvaged utility functions).

## Workflow

### 1. Compute ground truth and full simulation
```r
# ~5 hours on 8 cores for the full 60-iteration study.
# Uses cached files where available; only recomputes missing steps.
Rscript R/run_cached_analyses.R
```
This produces the four `.rds` files the manuscript and presentation depend on. To force recomputation, delete the relevant cache files first.

### 2. (Optional) Quick pipeline check
```r
# 5 iterations; ~25 min. Writes to sim_study_test.rds, so it does not
# overwrite production results.
Rscript R/run_sim_test.R
```

### 3. (Optional) WOT sensitivity to learner library
```r
# 15 iterations × 4 SL configs, parallel on 4 cores. ~1-2 hours.
Rscript R/run_wot_rich.R
```

### 4. Inspect results
```r
Rscript R/eval_full_table.R
```
Prints bias, RD coverage (bootstrap for Cox, EIF for LMTP), and HR coverage (Cox only) per estimator × estimand.

### 5. Render manuscript and presentation
```bash
quarto render draft-estimands-results-outline.qmd    # → .docx
quarto render presentation-estimands.qmd             # → .html (RevealJS)
```
Both documents load from the cached `.rds` files and will fail cleanly with a `stopifnot` message if a cache is missing.

## Data-generating process

The simulation uses `generate_hep_data()` with the following pedagogical parameters (see `R/run_cached_analyses.R`):

| Parameter | Value | Meaning |
|---|---|---|
| `h0` | 5e-3 | Baseline hazard (~16% 180-day event rate) |
| `HR_early` | 2.0 | Treatment HR before day 90 (early harm) |
| `HR_late` | 0.90 | Treatment HR after day 90 (modest attenuation) |
| `lambda_sw0` | 1e-3 | Baseline switching rate (~11% switch by day 180) |
| `gamma_A` | 0.80 | log-HR for treatment on switching |
| `gamma_ckd` | 0.60 | log-HR for CKD on switching |
| `complexity` | `FALSE` | Linear confounding only (simpler for teaching) |

Under these parameters, the five estimands produce distinct true risk differences:
- Treatment-policy: +7.1 pp
- No-switch: +9.1 pp
- While-on-treatment: +6.4 pp
- Composite: +18.3 pp
- Principal stratum: +9.0 pp (~19% of subjects are never-switchers)

## Performance metrics

Each estimator reports a point estimate plus a 95% CI. We score performance on two scales:

- **Risk-difference scale.** Every estimator produces an RD point estimate. We report bias (mean − truth) and RD coverage (proportion of 95% CIs containing the true RD). LMTP CIs come from the efficient influence function; Cox methods use bootstrap (200 replicates) on the KM risk difference.
- **Hazard-ratio scale.** Cox methods additionally produce a model-based HR CI. We report HR coverage against the true *marginal* HR (computed via Cox regression on the counterfactual potential outcomes, N = 500,000 per arm). LMTP does not target an HR.

Reporting both metrics exposes two failure modes: (1) the Cox HR reports a parameter on a scale that does not correspond to the estimand, and (2) even on its own HR scale, Cox fails to recover the correct marginal HR when confounding or informative switching is present.

## Key findings (60-iteration production run)

| Estimand | LMTP RD bias | LMTP RD cov. | Cox RD cov. (best) | Cox HR cov. (best) |
|---|---|---|---|---|
| Treatment-policy | +0.006 | 98% | 48% (naive) | 10% (naive) |
| No-switch | +0.040 | 86% | 52% | 8% |
| While-on-treatment | +0.036 | 48% | 8% | 7% |
| Composite | +0.006 | 95% | 78% (IPCW) | 38% |
| Principal stratum | +0.006 | 95% | 97% (naive, coincidental) | 60% (naive) |

LMTP achieves nominal coverage for four of five estimands. While-on-treatment remains the weak spot (48% coverage with the minimal learner library); `run_wot_rich.R` tests whether a richer library closes the gap.
