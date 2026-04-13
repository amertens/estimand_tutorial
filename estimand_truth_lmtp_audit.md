# Estimand, Truth, and LMTP Implementation Audit

**Date:** April 2026
**Auditor:** Claude (automated methodological review)
**Scope:** DGP.R, calc_truth.R, R/helpers.R, R/run_simulations.R

---

## Executive Summary

1. **CRITICAL: treatment-policy LMTP targets the wrong estimand.** `prepare_lmtp_data()` always creates time-varying A_j columns reflecting observed switching. `static_binary_on/off` then overrides ALL A_j to a constant. This holds treatment fixed across all time points — which is a **no-switch** intervention, not treatment-policy. Treatment-policy should assign baseline treatment and let switching occur naturally.

2. **treatment-policy and no-switch LMTP are structurally identical.** Both apply the same static binary intervention on time-varying A_j columns. The only difference is the input data's censoring distribution. Both target "hold treatment constant at every time point" = hypothetical no-switch.

3. **Truth calculations are correct for all 5 estimands.** Each truth function uses the appropriate DGP configuration, outcome derivation, and conditioning.

4. **while-on-treatment LMTP uses the wrong data.** It receives `dat_raw` (treatment-policy outcomes without censoring at switch) but should use `dat_cox` (WOT-derived outcomes censored at switch). The Y/C columns do not reflect the WOT estimand.

5. **composite LMTP intervenes on time-varying treatment when switching IS the outcome.** Holding treatment constant (preventing switching) contradicts the composite definition. Should use baseline-only treatment.

6. **principal-stratum LMTP is a reasonable oracle analysis.** For true never-switchers, A_j are already constant, so the static intervention is consistent.

7. **`intervention_no_switch()` is defined but never called.** It would be the correct shift function for no-switch/WOT with time-varying A_j.

8. **The LMTP method label "LMTP SDR (no-switch)" is used for both no_switch and WOT.** WOT should have a distinct label.

---

## Estimand-by-Estimand Audit

| Estimand | Formal Intervention | Truth | LMTP | Main Issue | Fix |
|---|---|---|---|---|---|
| Treatment-policy | Assign baseline A; let switching occur naturally | **Correct** | **WRONG**: static_binary on all A_j = no-switch | LMTP holds trt constant = no-switch | Use baseline-only trt (`"treatment"`) |
| No-switch | Hold treatment constant at assigned value | **Correct** | **Correct** (A_j constant since switch_on=FALSE) | Works by coincidence | Document equivalence |
| While-on-treatment | Follow only while on assigned trt (censor at switch) | **Correct** (conditional on non-switching) | **WRONG**: uses dat_raw (TP outcomes) not dat_cox (WOT outcomes) | Y/C columns don't reflect switch censoring | Use dat_cox for LMTP prep |
| Composite | Outcome = min(event, switch); assign baseline trt | **Correct** | **Partially wrong**: intervenes on A_j when switching IS the outcome | Should use baseline-only trt | Use baseline-only trt |
| Principal stratum | Restrict to latent never-switchers under both trt | **Correct** | **Reasonable**: A_j constant for never-switchers | Oracle analysis; acceptable | Document limitations |

---

## Detailed Findings

### treatment_policy

**Truth (calc_truth.R):** Generates all-treated/all-control with `switch_on=TRUE`. Switching occurs and modifies post-switch hazard. `follow_time = pmin(event_time, censor_admin)`. Correct: marginal risk under "assign treatment, let nature take its course."

**LMTP (run_simulations.R lines 212-227):** Uses `dat_raw` → `prepare_lmtp_data()` → time-varying A_j columns reflecting observed switching → `static_binary_on/off` overrides all A_j to constant. **This is a no-switch intervention, not treatment-policy.**

**Fix:** For treatment-policy, use baseline-only treatment. Add `time_varying_trt = FALSE` parameter to `prepare_lmtp_data()`. When FALSE, skip A_j creation and use `trt = "treatment"` in the LMTP call.

### no_switch

**Truth:** `switch_on=FALSE`. No switching generated. Correct.

**LMTP:** Data has `switch_time = Inf`, so `make_A()` produces constant A_j = baseline treatment. `static_binary_on/off` overrides to the same constant. **Works correctly** because the data already has no switching.

### while_on_treatment

**Truth:** Conditions on `switch_time > tau`. Computes risk among non-switchers. Standard WOT definition.

**LMTP:** Uses `dat_raw` (treatment-policy outcomes). Y and C columns do NOT reflect censoring at switch. LMTP sees full treatment-policy follow-up, not WOT.

**Fix:** Use `dat_cox` (from `derive_estimand("while_on_treatment")`) for LMTP data preparation. This censors Y/C at switch. Label as "LMTP SDR (WOT)".

### composite

**Truth:** Uses `derive_estimand("composite")`. `event = I(min(event_time, switch_time) <= censor_admin)`. Correct.

**LMTP:** Uses `dat_cox` with composite outcome. Y columns correctly encode the composite event. But A_j columns model the switching trajectory, and `static_binary_on/off` intervenes to hold treatment constant — preventing the switching that IS part of the composite outcome.

**Fix:** Use baseline-only treatment for composite. The estimand is "what is the 180-day composite risk under baseline assignment A=1 vs A=0?" — no need to model time-varying treatment.

### principal_stratum

**Truth:** Paired RNG draws via `.Random.seed` save/restore. Correct identification of the latent never-switcher stratum.

**LMTP:** Two variants: observed non-switchers (approximation) and true never-switchers (oracle). For both subsets, A_j are approximately constant (non-switchers don't switch), so static intervention is consistent. Reasonable.

---

## Cross-Cutting Issues

1. **Structural identity of TP and NS LMTP.** The code dispatches differently (TP goes to `else` branch, NS goes to `no_switch/WOT` branch) but both end up calling `static_binary_on/off` on A_j columns. The TP branch uses `dat_raw` with switching; the NS branch uses `dat_raw` without switching. But LMTP overrides all A_j regardless, so the target intervention is identical.

2. **`intervention_no_switch()` is orphaned.** Defined in helpers.R but never invoked. It returns `data[["treatment"]]` — preserving baseline assignment at each time point — which is exactly the no-switch shift function.

3. **Label collision.** Both no_switch and WOT use the label `"LMTP SDR (no-switch)"` in `run_simulations.R` line 175. WOT should be labeled distinctly.

4. **Comment accuracy.** The comment block in `run_simulations.R` lines 132-143 claims treatment_policy uses "time-fixed treatment" but the code path uses time-varying A_j.

---

## Recommended Code Changes

### Priority 1: Fix treatment-policy LMTP

Add `time_varying_trt` parameter to `prepare_lmtp_data()`:

```r
prepare_lmtp_data <- function(dat, tau = 180, bin_width = 1,
                              trt_var = "treatment",
                              switch_var = "switch_time",
                              time_varying_trt = TRUE, ...) {
  ...
  if (time_varying_trt) {
    # create A_j columns via make_A() [existing code]
    A_cols <- paste0("A", seq_len(n_bins))
  } else {
    A_cols <- NULL  # use baseline trt only
  }
  ...
  list(..., A_cols = A_cols, ...)
}
```

In `run_lmtp_analysis()`, when `A_cols` is NULL, use `trt = "treatment"`.

In `run_simulations.R`:
- treatment_policy: `prepare_lmtp_data(dat_raw, time_varying_trt = FALSE)`
- no_switch: `prepare_lmtp_data(dat_raw, time_varying_trt = TRUE)` (or FALSE — equivalent since A_j are constant)
- WOT: `prepare_lmtp_data(dat_cox, time_varying_trt = TRUE)` — use WOT-derived data
- composite: `prepare_lmtp_data(dat_cox, time_varying_trt = FALSE)` — baseline only
- PS: `prepare_lmtp_data(subset, time_varying_trt = FALSE)` — never-switchers don't switch

### Priority 2: Fix WOT LMTP data source

Change `dat_raw` to `dat_cox` for WOT in `run_simulations.R` line 173.

### Priority 3: Add distinct WOT label

Change `"LMTP SDR (no-switch)"` to `"LMTP SDR (WOT)"` for the WOT branch.

### Priority 4: Update comments

Fix the comment block in `run_simulations.R` to accurately describe the per-estimand dispatch.

---

## Recommended Wording Changes

1. **Presentation slide on estimand-aligned estimation:** Should note that treatment-policy LMTP uses baseline-only treatment (not time-varying A_j), because intervening on time-varying treatment to hold it constant is a no-switch (hypothetical) intervention.

2. **Manuscript Limitations section:** Update to reflect that treatment-policy and no-switch are now correctly distinguished via different LMTP specifications (baseline-only vs time-varying).

3. **Code comments in run_simulations.R:** Replace the current comment block (lines 132-143) with an accurate description of the per-estimand dispatch logic.
