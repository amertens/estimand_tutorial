# Estimand–Estimator Mismatch in Time-to-Event Safety Analyses

Simulation study and manuscript comparing Cox-based analyses (naive, censor-at-switch, IPCW) with LMTP SDR across five ICH E9(R1) estimand strategies, using a hepatitis B renal safety scenario with informative treatment switching.

## Repository layout

### Top-level documents
- [`draft-estimands-manuscript.qmd`](draft-estimands-manuscript.qmd) — main manuscript (Quarto → Word). Renders from cached simulation results.
- [`presentation-estimands.qmd`](presentation-estimands.qmd) — RevealJS slide deck. Reads the same caches as the manuscript.
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
# 20 iterations × 4 SL configs, parallel on 4 cores. ~1-2 hours.
Rscript R/run_wot_rich.R
```

### 4. Inspect results
```r
Rscript R/eval_full_table.R
```
Prints bias, RD coverage (bootstrap for Cox, EIF for LMTP), and HR coverage (Cox only) per estimator × estimand.

### 5. Render manuscript and presentation
```bash
quarto render draft-estimands-manuscript.qmd         # → .docx
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

## Current implementation: all five estimands

### Overview

One dataset is generated under the treatment-policy regime (`generate_hep_data()` with `switch_on = TRUE`). Switching occurs and modifies the post-switch hazard, but does not censor follow-up. All raw timing variables (`event_time`, `switch_time`, `censor_admin`) are retained. Each estimand analysis derives its own outcome/censoring structure from this common dataset. The no-switch estimand additionally requires a separate DGP call with `switch_on = FALSE` for the Cox analyses; the no-switch LMTP uses the switching dataset and intervenes to prevent switching.

LMTP analyses use baseline covariates `c("age", "sex_male", "ckd", "cirrhosis", "nsaid", "diabetes", "hypertension", "heart_failure")` by default (see `prepare_lmtp_data()`). Cox analyses use `c("age", "ckd", "cirrhosis", "diabetes", "heart_failure")`.

---

### 1. Treatment-policy

**Scientific question:** What is the 180-day risk comparing populations initiated on active vs comparator, regardless of subsequent switching?

**Target parameter:** ψ^TP = E[Y^(a=1)(180)] − E[Y^(a=0)(180)], where switching occurs naturally post-baseline.

**Truth** (`calc_truth.R: truth_treatment_policy`):
- Generate 500K all-treated and 500K all-control with `switch_on = TRUE`, `dep_censor = FALSE`
- Switching occurs and modifies the post-switch hazard (switchers move to the other regimen's hazard, including the comparator's HR = 1.0)
- `follow_time = min(event_time, censor_admin)` — switching does NOT censor
- Compute `mean(event == 1 & follow_time <= 180)` in each arm

**LMTP specification** (`run_simulations.R`):
- Data: `dat_raw` (treatment-policy, switching present, not censored)
- `time_varying_trt = FALSE` — single baseline `"treatment"` column
- `censor_at_switch = FALSE` — C columns encode only administrative censoring
- Shift: `static_binary_on` sets A = 1; `static_binary_off` sets A = 0
- **Mechanism:** LMTP assigns treatment at baseline. Post-baseline switching is part of the natural course. The censoring model handles administrative censoring only; switching affects outcomes through the Y columns (events may occur after switching, under the post-switch hazard).

**Cox analyses** (applied to `dat_raw`):
- **Cox naive:** baseline treatment only, follow-up not censored at switch. This is the correct Cox model *for this estimand*, though its conditional HR does not correspond to the marginal RD.
- **Cox censor-at-switch:** targets a different estimand (WOT/no-switch); reported for comparison.
- **Cox IPCW:** censors at switch and reweights non-switchers by the inverse probability of remaining unswitched.

---

### 2. Hypothetical no-switch

**Scientific question:** What would the 180-day risk be if no patient could switch?

**Target parameter:** ψ^HYP = E[Y^(a=1, s̄=0)(180)] − E[Y^(a=0, s̄=0)(180)], where s̄ = 0 means switching is prevented at all times.

**Truth** (`calc_truth.R: truth_no_switch`):
- Generate 500K all-treated and 500K all-control with `switch_on = FALSE`
- No switching occurs; event times are drawn under sustained assigned treatment only
- `follow_time = min(event_time, censor_admin)` with `censor_admin = max_follow` (no dependent censoring)
- Compute `mean(event == 1 & follow_time <= 180)` in each arm

**LMTP specification** (`run_simulations.R`):
- Data: `dat_raw` (the **switching** dataset with `switch_on = TRUE`). LMTP needs to observe the switching process in order to model it and then intervene to prevent it.
- `time_varying_trt = TRUE` — creates A_1, A_2, …, A_13 columns (one per biweekly bin) that reflect the observed switching trajectory (A_j = baseline treatment before switch, 1 − baseline after)
- `censor_at_switch = FALSE` — no additional censoring
- Shift: `static_binary_on` sets ALL A_j = 1; `static_binary_off` sets ALL A_j = 0
- **Mechanism:** The static binary intervention overrides the observed switching trajectory, forcing treatment to stay constant from baseline through the full 180 days — directly implementing the counterfactual "no switching." Identification requires sequential exchangeability (no unmeasured confounders of the switching decision at each time point).

**Cox analyses** (applied to `dat_noswitch` generated with `switch_on = FALSE`):
- Cox naive, censor-at-switch, and IPCW all produce identical results because there is no switching to distinguish them; IPCW falls back to naive Cox internally when `sum(did_switch) == 0`.

---

### 3. While-on-treatment (WOT)

**Scientific question:** What is the 180-day risk under baseline treatment with follow-up truncated at the time of switching?

**Target parameter:** 180-day cumulative incidence under baseline treatment assignment with switching treated as a censoring event.

**Truth** (`calc_truth.R: truth_while_on_treatment`):
- Generate 500K all-treated and 500K all-control with `switch_on = TRUE`, `dep_censor = FALSE`
- Compute **direct Monte Carlo**: `mean(event_time <= pmin(switch_time, tau))` in each arm
- Under counterfactual uniform treatment with no dependent admin censoring, this is the exact WOT risk — no KM machinery needed
- **This is not the conditional risk among non-switchers.** It is the risk under follow-up that is truncated at switching, summed over the whole arm.

**LMTP specification** (`run_simulations.R`):
- Data: `dat_raw` (treatment-policy data, switching present)
- `time_varying_trt = FALSE` — single baseline `"treatment"` column
- `censor_at_switch = TRUE` — **this is the key difference from treatment-policy**
  - Follow-up is redefined as `min(event_time, switch_time, follow_time, tau)`
  - Events after switching are excluded from Y
  - C columns encode BOTH administrative censoring AND switch-induced censoring
- Shift: `static_binary_on/off` at baseline only
- **Mechanism:** Treatment is assigned at baseline. Switching enters through the censoring nodes (C columns), not the treatment nodes. The LMTP censoring model adjusts for the informative component of switch-censoring (CKD drives both switching and the outcome). This is a proper censoring-based WOT analysis, not a treatment-node intervention to prevent switching.

**Cox analyses** (applied to `dat_cox = derive_estimand(dat_raw, "while_on_treatment")`):
- After `derive_estimand`, `follow_time` is already min(event, switch, admin), so Cox naive, censor-at-switch, and IPCW all produce identical results on this data. IPCW falls back to naive Cox because `switched` in the derived data is 0 for everyone (switching has been absorbed into the new follow_time).

---

### 4. Composite

**Scientific question:** What is the 180-day risk of renal failure OR treatment switching (whichever comes first)?

**Target parameter:** ψ^COMP = E[Y_comp^(a=1)(180)] − E[Y_comp^(a=0)(180)], where Y_comp^a(t) = I(min(T_event^a, T_switch^a) ≤ t).

**Truth** (`calc_truth.R: truth_composite`):
- Generate 500K all-treated and 500K all-control with `switch_on = TRUE`, `dep_censor = FALSE`
- Apply `derive_estimand("composite")`: `follow_time = min(event_time, switch_time, censor_admin)`, `event = I(min(event_time, switch_time) <= censor_admin)`
- Compute `mean(event == 1 & follow_time <= 180)` in each arm

**LMTP specification** (`run_simulations.R`):
- Data: `dat_cox = derive_estimand(dat_raw, "composite")` — Y columns encode the composite outcome
- `time_varying_trt = FALSE` — single baseline `"treatment"` column
- `censor_at_switch = FALSE` — switching is in the Y columns (outcome), not the C columns
- Shift: `static_binary_on/off` at baseline
- **Mechanism:** Switching is part of the outcome, not a treatment change or a censoring event. The LMTP censoring model adjusts only for administrative censoring. No time-varying treatment is needed.

**Cox analyses** (applied to `dat_cox` from `derive_estimand("composite")`):
- All three Cox approaches produce identical results because switching has been absorbed into the composite event indicator. IPCW falls back to naive internally.

---

### 5. Principal stratum

**Scientific question:** What is the 180-day risk among patients who would not switch under either treatment assignment?

**Target parameter:** ψ^PS = E[Y^(a=1)(180) | S^(a=1) > 180, S^(a=0) > 180] − E[Y^(a=0)(180) | S^(a=1) > 180, S^(a=0) > 180].

**Truth** (`calc_truth.R: truth_principal_stratum`):
- Generate 500K all-treated and 500K all-control with `return_potential_switching = TRUE` and the **same seed** for baseline covariates
- The DGP draws potential switching times under both treatments using paired random number streams (save/restore `.Random.seed`)
- `never_switcher = I(switch_time_a1 > max_follow AND switch_time_a0 > max_follow)`
- Restrict to subjects where BOTH `df_a1$never_switcher == 1` AND `df_a0$never_switcher == 1`
- Compute `mean(event_time[never_sw] <= 180)` in each arm

**LMTP specification — oracle (primary)** (`run_simulations.R`):
- Data: `dat_raw[dat_raw$never_switcher == 1, ]` — restricted to true never-switchers
- `time_varying_trt = FALSE`, `censor_at_switch = FALSE`
- Shift: `static_binary_on/off` at baseline
- Label: `"LMTP SDR (oracle PS)"`
- **Mechanism:** Since never-switchers by definition do not switch under either treatment, their treatment is constant throughout. Baseline-only treatment is appropriate. This is an oracle analysis using simulation-only information (the `never_switcher` indicator requires knowing both potential switching outcomes).

**LMTP specification — observed non-switcher approximation (supplementary):**
- Data: `dat_raw[dat_raw$switched == 0, ]` — restricted to subjects who did not switch under their observed treatment
- Same LMTP specification as oracle
- Label: `"LMTP SDR (obs. non-sw, approx)"`
- **Limitation:** Observed non-switchers are NOT the same as the principal stratum. A subject who did not switch on active might have switched on comparator (and vice versa). This proxy systematically selects lower-risk patients and does not recover the oracle truth.

**Cox analyses** (applied to full `dat_raw`, same as treatment-policy):
- Cox naive, censor-at-switch, and IPCW all run on the full dataset (not restricted to the stratum). Compared against the PS truth. Note that Cox naive achieves high coverage here not because it targets the principal stratum but because, under moderate switching selection (`gamma_A` = 0.80), the never-switcher subpopulation is compositionally similar to the overall cohort.

---

### Summary table

| Estimand | LMTP data | `time_varying_trt` | `censor_at_switch` | LMTP shift | Switching handled via |
|---|---|---|---|---|---|
| Treatment-policy | `dat_raw` (switching) | FALSE | FALSE | Baseline on/off | Natural course (in Y) |
| No-switch | `dat_raw` (switching) | **TRUE** | FALSE | Hold all A_j constant | Treatment intervention (prevented) |
| While-on-treatment | `dat_raw` (switching) | FALSE | **TRUE** | Baseline on/off | Censoring event (in C) |
| Composite | `dat_cox` (composite-derived) | FALSE | FALSE | Baseline on/off | Part of outcome (in Y) |
| Principal stratum | `dat_raw` subset (oracle) | FALSE | FALSE | Baseline on/off | Restricted to never-switchers |

---

## Key findings (60-iteration production run)

| Estimand | LMTP RD bias | LMTP RD cov. | Cox RD cov. (best) | Cox HR cov. (best) |
|---|---|---|---|---|
| Treatment-policy | +0.006 | 98% | 48% (naive) | 10% (naive) |
| No-switch | +0.040 | 86% | 52% | 8% |
| While-on-treatment | +0.036 | 48% | 8% | 7% |
| Composite | +0.006 | 95% | 78% (IPCW) | 38% |
| Principal stratum | +0.006 | 95% | 97% (naive, coincidental) | 60% (naive) |

LMTP achieves nominal coverage for four of five estimands. While-on-treatment remains the weak spot (48% coverage with the minimal learner library); `run_wot_rich.R` tests whether a richer library closes the gap.
