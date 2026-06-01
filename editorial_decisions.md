# Editorial decisions and revision changelog

## Base and working file

- **Base:** `draft-estimands-manuscript.qmd` (the May-4 source rendered to the
  `5-4` docx that Joy commented on). Confirmed by content: it holds all 15
  `[CO-AUTHOR NOTE]` blocks and the exact phrasings flagged for rewrite, and it
  lacks the PH-satisfied and time-varying-confounder sections that M1/M2 ask to
  port from v5.
- **Donor:** `draft-estimands-manuscript-v5.qmd` (already contains the
  PH-satisfied, TV-confounder, variance-diagnostic, IPCW, and Practical-
  Implications material).
- **Working revision:** `draft-estimands-manuscript-5-4-rev.qmd` on git branch
  `revise-estimands-5-4`. The canonical base and v5 are left untouched.

## Revision status

### Phase A - global style and mechanical sweep (COMPLETE, committed 17230ac)

- All 15 `[CO-AUTHOR NOTE]` blocks removed; open decisions captured below.
- Estimator naming unified to SDR (124 sites); abstract defines the family
  (longitudinal TMLE and the closely related SDR estimator; SDR used for all
  results). [M6, C1, C17]
- `active arm` / `active-arm` -> `TDF` / `TDF arm` (16 sites). [C4]
- ICH attribute list now reads "...time horizon, and population-level summary
  measure" in the Introduction, the estimand-definition paragraph, and Reporting
  Recommendations item 1. [inline edit]
- Phrasing: removed "Headline", "bake-off", "escape valve", "box to tick";
  rewrote the three named antithesis sentences; "aligns with" -> "follows" (the
  one remaining "aligns" is the legitimate estimator-estimand sense); reduced
  "deliberately"/"intentionally" to zero. Em-dashes and contractions were
  already absent in the base.
- Cox results reframed as an estimand-mismatch diagnostic with explicit
  two-target framing. [M1, partial]
- No-switch diagnosis sharpened: variance diagnostic (empirical SD 0.062 vs mean
  EIF SE 0.060), bias-driven, nuisance under-fitting, richer-library confirmation
  flagged as compute-blocked. [M3]
- IPCW limitation rewritten as an implemented same-data result (pooled
  time-varying IPCW; ESS ~ 1320; Cox IPCW misses the marginal HR by +0.28;
  50-iteration MC SE noted). [M4]
- Two explicit contributions added to the Introduction. [M7]
- Abstract trimmed from 612 to ~300 words; method-level granularity removed; the
  C0 worked-example sentence retained.

### Phase B - content ports from v5 (PENDING)

These are written and Joy-reviewed in v5 but not yet grafted onto the 5-4 base:

- **M1:** promote the proportional-hazards-satisfied decomposition (`tbl-ph-cox`,
  `tbl-ph-lmtp`) from supplement to MAIN; state the +0.07 to +0.38
  non-collapsibility overshoot under PH; add the "DGP chosen so Cox loses"
  rebuttal.
- **M2:** port the time-varying-confounder DGP section (`tbl-tv-lmtp`,
  `tbl-tv-cox`; `R/dgp_tv_confounder.R`); add the time-structure clarification
  paragraph; cite the feedback-DGP no-switch result (-0.032 -> -0.010).
- **M4/C21/C23:** add the same-data IPCW vs SDR Supplementary section body.
- **M8/C11:** Practical Implications early section; per-estimand Cox-implementation
  paragraphs; "When Cox regression is the right tool" subsection.
- **C6:** move the DGP parameter table to Supplementary Methods.
- Display equations for every estimand definition and counterfactual contrast.
- Compress the two motivating-context subsections (C7/C9/C10).
- Sensitivity placement: Markov result stated; dep_censor kept-with-explanation
  or moved here; stochastic arm cited from Discussion; MSCM labelled a pilot.

### Phase C - render and full audit (PENDING)

- Update `_render_v2.R` xref map for the 5-4-rev table/figure set; render to docx;
  run the full verification audit (0 "Table Table", 0 unresolved xref tokens,
  every estimand definition as a display equation, each coverage caption carries
  the Monte Carlo SE, etc.).

## Deferred decisions for co-authors

1. **Monte Carlo iterations:** 200 retained; captions to report MC SE (~1.5 pp at
   95% coverage for 200 iterations). Rerun at 500 if compute permits. [M5]
2. **Target journal:** undecided (JCE / Biometrics / Pharmaceutical Statistics);
   methods-journal prescriptive tone retained. [C30]
3. **dep_censor sensitivity:** counter-intuitive (non-informative censoring makes
   no-switch coverage worse, 94% -> 74%). Keep only with a mechanistic
   explanation a reviewer will accept; otherwise cut. Co-authors to decide. [C27]
4. **Stochastic half-switching arm:** kept in supplement, cited briefly from the
   Discussion; can demote to a future-work footnote. [C29]
5. **MSCM pilot (10 iters):** describe as an independent confirmatory check
   without citing 10-iteration coverage, or scale to >=200.

## Joy comment crosswalk (C0-C34)

Resolved in 5-4-rev (Phase A): C0 (worked-example note), C1/C17 (SDR naming),
C4 (TDF naming), C18-C22 (Cox-naive naming), inline ICH attribute edit. Reframed
in prose: M1 Cox-diagnostic, M3 no-switch, M4 IPCW. Pending port from v5
(Phase B): C3/C5/C16/C28 (baseline-CKD wording + TV clarification), C6 (DGP table
to supplement), C7/C9/C10 (compress motivating context), C8 (RR in
tables/figures), C11/C13 (per-estimand Cox paragraphs), C12/C14 (informative
loss-to-follow-up terminology), C15 (principal-stratum latent wording), C19
(informative process), C24/C26 (figure-caption boundary), C25 (formatting), C29
(stochastic), C30 (journal tone).
