# Language Precision and Code Audit Note

**Date:** April 2026

## Key Findings

1. **WOT truth used unnecessary KM machinery.** Under counterfactual uniform treatment with `dep_censor=FALSE`, the only censoring is switching. The KM estimator adjusts for this censoring, but a direct Monte Carlo calculation — `mean(event_time <= min(switch_time, tau))` — gives the same result more transparently. Replaced KM with direct computation.

2. **WOT implementation is correctly censoring-based.** `time_varying_trt=FALSE` with `censor_at_switch=TRUE` encodes switching in C columns. Treatment is baseline-only. This is a proper censor-at-switch WOT analysis, not a treatment-node intervention.

3. **Stale WOT language persisted in the manuscript.** Line 190 described WOT as "conditions on not switching" — a remnant of the old conditional-on-non-switching implementation. Corrected to "follow-up is truncated at switching."

4. **Principal stratum oracle vs approximation are clearly separated in code.** Oracle (`never_switcher==1`) is `results[["lmtp"]]` (primary); observed non-switcher (`switched==0`) is `results[["lmtp_approx"]]` (supplementary). Labels are distinct.

5. **Principal stratum distinction was not explicitly stated in the manuscript.** Added two paragraphs explaining that the principal stratum is defined by joint potential switching outcomes (latent), while observed non-switchers are defined by realised switching under one treatment (observed post-treatment event). These are not the same population.

## WOT-Specific Findings

### Truth computation
- **Before:** `survfit()` + `summary(sf, times=tau)` on `derive_estimand("while_on_treatment")` data
- **After:** `mean(event_time <= pmin(switch_time, tau))` — direct Monte Carlo proportion
- **Why:** Under counterfactual uniform treatment with no dependent admin censoring, the switching hazard is the same for all subjects within each arm. The direct calculation is exact for the truncation-based WOT risk and avoids unnecessary KM estimation.

### Estimator
- **Code:** `prepare_lmtp_data(dat_raw, time_varying_trt=FALSE, censor_at_switch=TRUE)`
- **Status:** Correct. Switching enters through C columns. Treatment is baseline-only.

### Truth-estimator alignment
- **Truth:** 180-day risk = P(event before both switch and 180 days)
- **Estimator:** LMTP with baseline treatment, switch-censoring in C columns
- **Aligned:** Yes, both target the 180-day risk under follow-up truncated at switching.

## Principal-Stratum-Specific Findings

### Oracle target
- `never_switcher = I(switch_time_a1 > max_follow AND switch_time_a0 > max_follow)`
- Constructed via paired RNG draws: save `.Random.seed`, draw under A=1, restore, draw under A=0.
- Truth restricts to subjects where BOTH `df_a1$never_switcher == 1` AND `df_a0$never_switcher == 1`.
- **Correct.**

### Main estimator
- `dat_raw[dat_raw$never_switcher == 1, ]` with `time_varying_trt=FALSE`
- Label: `"LMTP SDR (oracle PS)"`
- Key name: `results[["lmtp"]]` (primary slot)
- **Correct.**

### Approximation
- `dat_raw[dat_raw$switched == 0, ]` with `time_varying_trt=FALSE`
- Label: `"LMTP SDR (obs. non-sw, approx)"`
- Key name: `results[["lmtp_approx"]]` (supplementary slot)
- **Correctly separated.**

## Files Changed

| File | Change |
|---|---|
| `calc_truth.R` | Replaced KM-based WOT truth with direct Monte Carlo: `mean(event_time <= pmin(switch_time, tau))` |
| `draft-estimands-results-outline.qmd` | Fixed "conditions on not switching" → "follow-up truncated at switching". Updated WOT truth description. Added PS vs observed non-switcher distinction. Updated supplementary results caveat. |
| `presentation-estimands.qmd` | Updated WOT speaker note: "direct Monte Carlo, not KM" |

## Exact Language Changes

1. **WOT truth description (manuscript):**
   - Before: "Compute the Kaplan-Meier estimate... Under counterfactual uniform treatment with no dependent administrative censoring, the KM estimator is unbiased."
   - After: "The 180-day risk is computed directly as the proportion of subjects whose event occurs before both their switch time and 180 days... this Monte Carlo average is the exact WOT risk."

2. **WOT interpretation (manuscript):**
   - Before: "The while-on-treatment truth differs from the no-switch truth because it conditions on not switching before 180 days."
   - After: "The while-on-treatment truth differs from the no-switch truth because follow-up is truncated at switching. Subjects who switch before 180 days contribute only their pre-switch person-time."

3. **Principal stratum (manuscript, new):**
   - Added: "The principal stratum is defined by joint potential switching outcomes... latent... cannot be identified from observed switching behaviour... observed non-switchers... should be interpreted as a descriptive subgroup analysis or proxy, not as identification of the principal stratum estimand."

## Remaining Limitations

1. The direct WOT Monte Carlo truth (`mean(event_time <= pmin(switch_time, tau))`) is correct under counterfactual uniform treatment but does not account for dependent administrative censoring. In truth calculations, `dep_censor=FALSE`, so this is not an issue.

2. The WOT LMTP estimator (`censor_at_switch=TRUE`) combines switching and administrative censoring in the same C columns. A richer specification might model these separately.

3. The principal stratum oracle requires `return_potential_switching=TRUE` and paired RNG draws — simulation-only. In practice, identification requires monotonicity or sensitivity analysis.
