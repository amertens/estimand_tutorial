# Editorial decisions for co-author review

## Decisions needed before submission

1. Monte Carlo iterations
Current manuscript uses 200 production iterations. Decide whether to rerun at 500–1000 iterations for greater stability of coverage estimates.

2. Markov vs full-history sensitivity
Currently specified but not run. Decide whether to implement, move to future work, or delete.

3. Informative censoring sensitivity
Currently specified but not run. Decide whether to implement, move to future work, or delete.

4. Production IPCW comparison
Current manuscript acknowledges that the fair no-switch comparison would require time-varying IPCW and L-TMLE on the same observed switching data. Decide whether to implement before submission or retain as a limitation/future-work item.

5. No-switch diagnostic interpretation
Confirm whether the diagnostic evidence supports the claim that longitudinal nuisance-estimation burden, rather than severe weight instability, is the main driver of no-switch undercoverage in the primary DGP.

6. Stochastic intervention arm
Currently discussed as future work. Decide whether to add a small stochastic-switching simulation arm or keep the paper focused on the five ICH E9(R1) strategies.

7. Targeted-learning tone
Confirm that the paper reads as an estimand-alignment tutorial rather than a one-sided L-TMLE advocacy piece.

## Suggested default path

Proceed with the current focused manuscript unless co-authors strongly prefer additional simulations. Highest-value optional addition is the Markov vs full-history sensitivity for no-switch. The production IPCW comparison is valuable but may become a separate methodological paper.
