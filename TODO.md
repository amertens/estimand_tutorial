# Estimand manuscript — outstanding items
# Working file: draft-estimands-manuscript-v6.qmd
# Branch: revise-estimands-5-4
# Last updated: 2026-06-05

---

## BLOCKING — must resolve before submission

### B1. L-TMLE / SDR naming decision and global sweep
**What:** v6 has 174 prose instances of "L-TMLE" but all table labels say "LMTP SDR".
The paper includes a paragraph acknowledging the mismatch, which reads as an admission of
inconsistency to a reviewer.
**Decision needed:** Either (a) keep L-TMLE in prose and switch the running estimator to
`lmtp_tmle` so the name is accurate, or (b) sweep prose to "SDR / sequentially doubly
robust estimator" and define the family once. Option (b) matches what the code actually
runs. Option (a) requires a re-run.
**Who decides:** Co-authors. Once decided, one deliberate global pass (~1 hour).

### B2. Display equations for every estimand definition
**What:** The risk-difference, risk-ratio, and counterfactual contrast definitions and the
five estimand formalisations are all inline math. In the rendered docx they collapse into
unreadable run-on text.
**Fix:** Wrap each in `$$...$$` display blocks.
**Effort:** ~30 min.

### B3. Render v6 to docx and run the full cross-reference audit
**What:** v6 has never been rendered. The `_render_v2.R` xref map needs updating for v6's
table/figure labels. Until the docx exists, "Table Table"/"Figure Figure" artifacts and
broken cross-references cannot be confirmed clean.
**Steps:** (1) Audit v6 for all `@tbl-` and `@fig-` tokens; (2) update the xref map in
`_render_v2.R`; (3) run `Rscript --vanilla _render_v2.R draft-estimands-manuscript-v6.qmd`;
(4) open the docx and check every cross-reference resolved.
**Effort:** ~1 hour.

### B4. Figure 2 caption direction error
**What:** The caption says censoring at switch removes high-risk CKD patients, "leaving
more healthier patients in the risk set" — but the KM risk difference *rises* from 0.087
to 0.112 under censor-at-switch, which is the opposite direction. The stated mechanism
contradicts the numbers.
**Fix:** Re-derive the selection mechanism (censoring selectively removes high-risk TDF
person-time → the remaining TDF risk set is depleted of fast-events → the observed TDF
hazard drops → the HR-scale gap widens) and rewrite the caption so direction matches the
numbers.
**Effort:** ~15 min once the mechanism is confirmed.

---

## SHOULD FIX — would be flagged in peer review

### S1. PH-satisfied decomposition: move from supplement to main paper
**What:** The proportional-hazards-satisfied sensitivity (which shows Cox misses the
marginal HR by +0.07–0.38 from non-collapsibility alone, even when PH holds) is in the
supplement. Reviewers who object that "the DGP was tuned so Cox fails" will not read the
supplement. The rebuttal needs to be in the main paper.
**Fix:** Add ~2 paragraphs in the main Results (§"The hazard ratio does not answer the
risk question") citing the supplement for the full table. State the +0.07–0.38 numbers.
**Effort:** ~30 min.

### S2. dep_censor sensitivity: add mechanistic explanation or cut
**What:** Non-informative censoring makes no-switch coverage *worse* (94% → 74%). This
counter-intuitive result sits in the paper without a mechanism. A reviewer will demand an
explanation or flag it as a potential bug.
**Explanation to add:** Under informative censoring, the L-TMLE censoring model is doing
useful identification work that partially offsets the longitudinal nuisance pressure on the
no-switch estimand; removing informative censoring removes that scaffolding and exposes the
underlying difficulty of the sustained-treatment counterfactual.
**Effort:** ~20 min. Alternatively: move to `editorial_decisions.md` and cut from the
manuscript if the group does not want to defend it.

### S3. MSCM pilot: soften or scale up
**What:** The 10-iteration MSCM pilot (mean HR 1.74, 10% coverage) is cited alongside
50- and 200-iteration results. At 10 iterations the MC SE on coverage is ~10 pp; the
coverage figure is essentially noise.
**Fix:** Either (a) scale to ≥50 iterations, or (b) describe the MSCM result qualitatively
only ("consistent with the IPCW Cox finding") without citing the coverage number.
**Effort:** (a) ~2h compute; (b) ~10 min rewrite.

### S4. No-switch bias explanation: lean harder on the richer-learner negative result
**What:** The paper says no-switch confirmation with richer learners is "outstanding."
But the 25-paired-iteration richer-library sensitivity (lean vs lasso vs ridge vs XGBoost)
already shows *no improvement*, which is fairly strong evidence against the learner-
flexibility hypothesis. The writing treats this as inconclusive when it could be
conclusive.
**Fix:** Reframe: the richer-learner sensitivity rules out linear-learner limitation as
the mechanism; the residual bias is structural, not learner-dependent. Remove the
"confirmation outstanding" hedge.
**Effort:** ~20 min.

---

## NICE TO HAVE

### N1. Scale main simulation to 500 iterations
**What:** Joy prefers 500; current count is 200. At 200, MC SE on a 95% coverage estimate
is ~1.5 pp — borderline for Biometrics / Statistics in Medicine.
**Note:** Each rep is ~18s, so 300 additional reps ≈ 90 min on the current machine.
The run is checkpointed and resumable. Set `N_ITER=500` in the run script env.
**Decision:** Co-authors; depends on timeline.

### N2. Abstract structure: consider a short sensitivity-analysis table
**What:** v6 abstract lists "seven sensitivity analyses" in a long parenthetical. A
two-column mini-table (sensitivity name | headline result) in the Methods section would
be cleaner and make the paper look more comprehensive.
**Effort:** ~15 min.

### N3. Conclusions synthesis subsection: tighten
**What:** The "How the simulation supports each recommendation" subsection in Conclusions
is listy and repeats Discussion material. A reviewer reading linearly feels déjà vu.
**Fix:** Cut to 3–4 sentences that synthesize rather than enumerate.
**Effort:** ~20 min.

---

## HOUSEKEEPING

### H1. Remove the "Possible extensions" internal section before submission
The section (§"Possible extensions [INTERNAL: for co-authors only; remove before
submission]") is ~50 lines in v6. Easy to forget. Remove it as the final pre-submission
step.

### H2. Author affiliations, contributions, and acknowledgements
Placeholders only. Needed before any journal submission. Up to co-authors.

### H3. YODA proposal: fill in PI/institutional/funding fields
`yoda_proposal.md` is drafted but has placeholder text. Complete separately before
submitting to YODA.

### H4. Target journal decision
Affects word limits, supplement policy, and whether Practical Implications section reads
as appropriate for the outlet. Options: JCE, Pharmaceutical Statistics, Statistics in
Medicine, Biometrics, Pharmacoepidemiology and Drug Safety.

### H5. Archive the old working files
`draft-estimands-manuscript.qmd` (5-4 base), `draft-estimands-manuscript-5-4-rev.qmd`,
and v2–v5 drafts are in the repo root alongside v6. Before submission, either archive to
a subfolder or leave as git history and add to .gitignore for the submission bundle.

### H6. References audit
Check for duplicate citations, stale preprint DOIs that now have published versions, and
any `[@key]` tokens that do not resolve in `references.bib`.

---

## OPEN CO-AUTHOR DECISIONS (cannot be resolved unilaterally)

| # | Question | Current default |
|---|---|---|
| D1 | L-TMLE vs SDR naming (B1 above) | Undecided |
| D2 | 200 vs 500 main-simulation iterations (N1 above) | 200 |
| D3 | Target journal | Undecided |
| D4 | dep_censor sensitivity: keep with explanation or cut? | Keep (Joy defers) |
| D5 | TV feedback DGP: supplement vs headline? | Supplement (recommended) |
| D6 | Stochastic arm: supplement vs pure future-work footnote? | Supplement (Joy defers) |

---

## COMPLETED (reference)

- Phase A style sweep on 5-4-rev: CO-AUTHOR NOTEs, em-dashes, antithesis rewrites,
  abstract trim, M3/M4 rewrites, M7 contributions, ICH attribute list
- TV-confounder 200-iteration run (both calibrations, fully checkpointed)
- v5 TV framing updated with real 200-iter numbers + Option-1 boundary statement
- v6 created: v5 + Phase A style rules applied (0 em-dashes, 0 deliberately, etc.)
- joy3-comments-response.md rewritten with combined TV answer and C0-C34 crosswalk
- All changes pushed to branch `revise-estimands-5-4` on GitHub
