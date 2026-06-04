# Response to Joy's comments — round 4 (2026-05-20)
# File: draft-estimands-manuscript 5-4.docx
# Revised manuscript: draft-estimands-manuscript-v5.qmd (all edits in this round applied to v5)

This document covers two categories of change in this revision round:

1. **The time-varying-confounder cluster** (C3, C5, C16, C28, C31 and associated
   framing questions): a combined substantive response describing what was done, what
   the 200-iteration simulation shows, and where the open placement decision stands.

2. **All other comments C0 through C34 and the 2 inline edits**: one row per
   comment with the resolution location and status.

---

## 1. Combined response to the time-varying-confounder cluster (C3, C5, C16, C28, C31)

**Joy's concern, stated cumulatively across the four rounds:**
The switch hazard appeared to be time-constant and driven only by baseline CKD.
If L-TMLE is a longitudinal estimator, should the DGP contain a time-varying
confounder? Does baseline CKD capture what is really a time-varying marker?

**What the manuscript now says (§"Data-generating process", clarification paragraph):**

The DGP is explicit that it contains no time-varying covariate L(t). What is
time-varying is the switching trajectory itself, which depends on baseline CKD
and the baseline arm — not on an evolving marker. The switch hazard equation
makes this precise: it depends only on the baseline CKD indicator and the
baseline treatment assignment. The manuscript states this as a design choice with
a rationale: holding confounding at baseline isolates the two questions the paper
is about — the choice of summary measure (hazard ratio versus risk difference)
and the choice of intercurrent-event strategy — from the distinct and
well-studied problem of treatment-confounder feedback. Baseline adjustment for
CKD suffices for the confounding structure encoded here.

**Why a longitudinal estimator is needed without a time-varying confounder
(§"L-TMLE", "Why L-TMLE for all five estimands?" paragraph):**

The no-switch (hypothetical) estimand requires a sequential intervention on the
switching process at every biweekly bin. Identification requires a product of
bin-level outcome and switching nuisances — the longitudinal modified treatment
policy structure — regardless of whether confounding is baseline or time-varying.
The other four estimands admit simpler estimators in this DGP; L-TMLE is used
uniformly across all five strategies so that cross-estimand comparisons reflect
estimand structure rather than estimator switching. The time-varying censoring
process also requires the longitudinal machinery. The paper does not use the
phrase "time-varying confounding" when describing the primary DGP.

**What the 200-iteration feedback DGP supplement shows
(§"Supplementary Methods: Time-varying-confounder DGP sensitivity"):**

A supplementary sensitivity was run at two calibrations (200 iterations each,
N = 3,000):

- *Pure exogenous L(t)*: no treatment-confounder feedback; L(t) evolves through
  its baseline-CKD anchor, mean reversion, and noise. L-TMLE baseline-only:
  no-switch bias +0.003 (MC SE 0.007), coverage 92%. L-TMLE +L(t): bias +0.004
  (MC SE 0.007), coverage 91%. The two specifications are indistinguishable.

- *A(t−1) → L(t) feedback* (the canonical Robins setting): past treatment drives
  the marker up; the full feedback loop is active. L-TMLE baseline-only: no-switch
  bias −0.001 (MC SE 0.008), coverage 90%. L-TMLE +L(t): bias +0.012 (MC SE
  0.008), coverage 91%. Again indistinguishable within Monte Carlo noise
  (MC SE on coverage ≈ 0.02). Convergence rates were high for both (193/200 and
  190/200 respectively).

**What this means for Joy's concern:**

The result is a boundary of the paper's claims, not the positive payoff of
time-varying adjustment. Baseline adjustment recovered the no-switch estimand
even under the canonical feedback DGP, because L(t) in this calibration is
largely predictable from baseline CKD and treatment history. The part of L(t)
that would require explicit time-varying adjustment — an independent innovation
not captured by past treatment history — was not strong enough in this
calibration to break baseline adjustment. The paper therefore does not claim to
demonstrate L-TMLE's handling of irreducible treatment-confounder feedback. The
longitudinal estimator earns its place through the longitudinal no-switch
intervention and the censoring process, not through a time-varying confounder
that the primary DGP does not contain.

**Open decision for Joy and co-authors:**

Should the feedback DGP remain in the supplement or move to the headline
comparison? The recommendation is to keep it in the supplement. The finding is
a robustness check (baseline adjustment suffices even under feedback), not a
load-bearing demonstration of time-varying-confounding handling. Promoting it to
the headline would invite the criticism that the paper claims to address
treatment-confounder feedback when the calibration was not strong enough to
require it. If the group wants a positive payoff demonstration of time-varying
adjustment (adding L(t) materially improves performance over baseline-only), a
stronger calibration with larger independent innovation variance in L(t) would
be needed; that would be a different paper. Either way, the 200-iteration runs
are now in place for both calibrations, so the supplement result is
well-powered regardless of where it lives.

---

## 2. Comment-by-comment table (C0 to C34 plus 2 inline edits)

### Resolved by substantive change in this round

| # | Joy's comment | What was done | Location in v5 |
|---|---|---|---|
| C3, C5, C16, C28 | "evolving CKD?", "switch hazard seems constant?", "needs to be time-varying?" | DGP section now states plainly: no L(t) in the primary DGP; switch hazard depends only on *baseline* CKD indicator (equation given); design choice stated with rationale. "A clarification on the time structure of confounding" paragraph added. | §"Data-generating process", clarification paragraph after the switch-hazard equation |
| C31 | "I suspect this might be flagged in review — if possible we could do it; if not, maintain as a limitation and remove this section." | Production-quality pooled IPCW vs L-TMLE comparison implemented (50 iterations, same observed data, stabilised weights, ESS ≈ 1320). Cox IPCW still misses the marginal HR by +0.28. The Limitations paragraph no longer calls the input asymmetry the largest unresolved limitation. | §"Supplementary Methods: Direct IPCW vs L-TMLE comparison" + Limitations |

### Resolved in prior rounds (unchanged in v5)

| # | Joy's comment | Status |
|---|---|---|
| C0 | "are we only looking at treatment switching?" | Abstract states treatment switching is the worked example; the framework applies to other intercurrent events. |
| C1, C15 | TMLE vs SDR / LTMLE naming | "L-TMLE" in prose, "LMTP SDR" in table labels; convention stated once in the L-TMLE paragraph. |
| C2 | "what does this mean?" (competing-risks Monte Carlo truth) | Unpacked in Abstract Results. |
| C4 (new, this round) | "since both arms are active, would not it be more clear to say TDF?" | "active arm"/"active-arm" replaced by "TDF"/"TDF-treated" throughout wherever the phrase refers to the TDF group. |
| C6 (new, this round) | "Can move this to supplementary?" (DGP parameter table) | DGP parameter table moved to §"Supplementary Methods: DGP parameter specification". Main text retains qualitative DGP description and pointer to `DGP.R`. |
| C7, C9, C10 | Compress motivating context; oncology paragraph | Two motivating-context subsections merged into one shorter passage; Joy's compressed oncology paragraph integrated. |
| C8 | "risk ratios also" | Truth table has RR column; simulation tables annotate RD, RR, HR truths in captions; fig-sim-supp shows RD (top) and RR (bottom) with truth dashes on both. |
| C11, C13 | Juxtapose L-TMLE/Cox per estimand; "not observed did not switch?" | Each of the five estimand subsections has a *Cox implementation:* paragraph alongside *L-TMLE implementation:*. Principal-stratum row of tbl-alignment states: "this is the latent stratum, not the observed-did-not-switch subgroup." |
| C12, C14, C17 | "administrative censoring" terminology | Renamed to "informative loss to follow-up" wherever the DGP makes it informative (all five L-TMLE-implementation bullets and the DGP description). |
| C18, C19, C20 | Baseline-treatment Cox vs Cox naive naming | Global pass: "Cox naive (ignoring switching)" at first mention, "Cox naive" downstream. |
| C21, C23 | "why a separate no-switch dataset?"; "would not Cox implementation censor at switch?" | No-switch results paragraph explains: censor-at-switch is the natural Cox route to a no-switch quantity on the observed data; the separately-generated switching-off dataset is used only to compute the marginal-HR truth; the same-data IPCW vs L-TMLE comparison is in the supplement. |
| C22, C24 | Figure/table formatting; figure-caption boundary | No "Table Table"/"Figure Figure" artifacts; support-diagnostic caption reads "P(switch \| history) bounded away from 1 (equivalently, 1 − P(switch \| history) bounded above 0)." |
| C26 | Markov sensitivity | Markov vs full-history results integrated (200 iterations; essentially identical coverage → rules out long-history overfitting). |

### Resolved by minor edit or confirmed as no-change needed

| # | Joy's comment | Status |
|---|---|---|
| C2 | Competing-risks truth explanation | Plain-language Abstract Results sentence added. |
| C5 | "Separate dataset for no-switch" | Clarified: truth-only Monte Carlo evaluation on a separately-generated switching-off cohort; the observed-data estimator uses one dataset. |
| C25 | "Beautify" (table formatting) | Cross-reference artifacts resolved in render script. Further cosmetic formatting deferred to final submission pass. |
| C30 | Target journal | CO-AUTHOR NOTE removed. Methods-journal prescriptive tone retained as appropriate for JCE / Biometrics / Pharmaceutical Statistics / Statistics in Medicine. |

### Decisions still open for Joy and co-authors

| # | Question | Current default | Notes |
|---|---|---|---|
| TV-DGP placement | Keep feedback DGP in supplement or move to headline? | **Supplement** (robustness check) | See combined TV response above. The 200-iter result is in place either way. |
| Iteration count | 200 (current) vs 500 for main simulation tables? | 200 | Joy prefers 500 if feasible. Depends on available compute time before submission. |
| Target journal | JCE, Biometrics, Pharmaceutical Statistics, Statistics in Medicine? | Undecided | Tone already calibrated for a methods journal. |
| C27 (dep_censor sensitivity) | Keep the counter-intuitive non-informative-censoring result or cut? | Keep, with mechanistic explanation | Joy defers to co-authors. A mechanistic explanation (censoring selectively removes high-risk TDF person-time; removing that selection makes the no-switch estimand harder) is needed if it stays. |
| C29 (stochastic-arm placement) | Supplement (current) vs demote to pure future work? | Supplement | Joy defers to co-authors. Current placement: cited briefly from Discussion §"Beyond the five strategies." |

---

## Inline edits incorporated

- **Reporting Recommendations item 1** (Conclusions): Joy's inline addition of
  "population-level summary measure" as the sixth ICH E9(R1) attribute is
  incorporated. The list now reads: "population, treatment, outcome, intercurrent
  events, strategy per intercurrent event, time horizon, and population-level
  summary measure."

- **Sixth attribute also added** in the Introduction's estimand-definition
  paragraph and in the estimand-framework opening for consistency.
