# WOT and Principal Stratum Revision Note

**Date:** April 2026

## What was wrong before

### While-on-treatment (WOT)
- **Truth:** Computed as the conditional risk among natural non-switchers (`switch_time > tau`). This is a subgroup-conditional quantity, not a censor-at-switch estimand.
- **LMTP:** Used time-varying A_j columns with `static_binary_on/off` to hold treatment constant. This is a no-switch intervention, not a WOT censoring-based analysis. Switching was handled through treatment nodes rather than censoring nodes.
- **Result:** WOT LMTP was identical to the no-switch LMTP and did not target the WOT estimand.

### Principal stratum
- **Primary estimator:** Observed non-switcher approximation (`switched == 0`) was listed first, with the oracle analysis (`never_switcher == 1`) listed second.
- **Label:** The oracle analysis was called "LMTP SDR (true PS)" — should be primary.
- **Problem:** Simulation results compared primarily to the observed-non-switcher proxy rather than the oracle.

## What was changed

### WOT — redesigned as a censoring-based estimand

**calc_truth.R:**
- `truth_while_on_treatment()` now applies `derive_estimand("while_on_treatment")` to the counterfactual datasets, which censors follow-up at switch time.
- The truth is computed as the KM estimate of 180-day cumulative incidence under this censoring regime — not as the conditional risk among non-switchers.
- Under counterfactual uniform treatment with no dependent admin censoring, the KM estimator is unbiased for the WOT risk.

**R/helpers.R:**
- Added `censor_at_switch` parameter to `prepare_lmtp_data()`.
- When `censor_at_switch = TRUE`: follow-up is redefined as `min(event_time, switch_time, follow_time, tau)`. Events after switching are excluded. Censoring (C columns) encodes both administrative censoring AND switch-induced censoring.
- When `censor_at_switch = FALSE` (default): existing behavior unchanged.

**R/run_simulations.R:**
- WOT dispatch now uses `time_varying_trt = FALSE` (baseline treatment only) with `censor_at_switch = TRUE`.
- Switching enters through C columns (censoring nodes), not A columns (treatment nodes).
- The LMTP censoring model adjusts for informative switch-censoring.

### Principal stratum — oracle estimator is now primary

**R/run_simulations.R:**
- The oracle analysis (`never_switcher == 1`) is now the primary result, labeled "LMTP SDR (oracle PS)".
- The observed non-switcher approximation is moved to a supplementary slot, relabeled "LMTP SDR (obs. non-sw, approx)".

**QMD:**
- PS section updated to state that the primary simulation analysis uses the oracle indicator.
- Observed non-switcher approximation is described as supplementary.

## Files changed

| File | Change |
|---|---|
| `calc_truth.R` | Rewrote `truth_while_on_treatment()` — KM under censor-at-switch, not conditional on non-switching |
| `R/helpers.R` | Added `censor_at_switch` param to `prepare_lmtp_data()`. Updated header comments. |
| `R/run_simulations.R` | WOT: baseline trt + censor_at_switch. PS: oracle primary, obs approx supplementary. Updated dispatch comments. |
| `draft-estimands-results-outline.qmd` | Updated WOT and PS descriptions in Estimand Framework and LMTP specification sections. |

## Remaining limitations

1. **WOT truth uses KM under switch-censoring.** Under counterfactual uniform treatment (all-treated or all-control) with `dep_censor = FALSE`, the KM estimator is unbiased because switching is the only censoring mechanism and there are no confounders. In observed data with mixed treatment and dependent admin censoring, the LMTP censoring model must correctly specify both switching and admin censoring.

2. **WOT LMTP censoring model specification.** The LMTP censoring model uses the same Super Learner library (`SL.mean`, `SL.glm`, `SL.bayesglm`) for both admin censoring and switch-censoring. A richer library and more CV folds would improve performance.

3. **Principal stratum oracle is simulation-only.** The `never_switcher` indicator requires paired potential switching draws under both treatments — not available in practice. The oracle analysis validates the estimation framework but is not applicable to real data.

4. **Observed non-switcher approximation.** Restricting to `switched == 0` selects a different subpopulation than the true never-switcher stratum. In practice, principal stratum identification requires structural assumptions (monotonicity) or sensitivity analysis.
