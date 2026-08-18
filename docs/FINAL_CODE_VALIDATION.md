# Final code-validation record

This update consolidates the code changes validated during the processed reviewer-data rerun.

## Validated analytical blocks

- DEP primary analysis and sensitivity/robustness workflows
- WGCNA network construction, module characterization, trait associations, stability, structural preservation, and reporting
- ML primary nested CV, p-tau217, APOE, LOCO, fixed-panel, matched sensitivity, clinical regression, LOCO meta-analysis, and publication-output auditing

## Locked input quantities checked by `RUN_PRECHECK`

- 653 total participants: 327 CN and 326 AD
- 639 primary complete cases: 313 CN and 326 AD
- 579 APOE analytic participants
- 343 AT(N) analytic participants
- 382 selected matched participants: 191 CN and 191 AD
- 10,751 human SOMAmer features
- 9,638 gene-level DEP/WGCNA features
- 9,638 historical ML representative features

## ML concordance anchors

The processed-data rerun recovered the submitted seven-protein panel:

`SPC25`, `CPLX2`, `TCP11L1`, `ACHE`, `ODC1`, `SPON1`, `RTN4RL1`.

Country-held-out AUCs and the random-effects LOCO synthesis were also reproduced from the processed-data interface. The exact submitted matched cohort was restored from its pseudonymized inclusion flag rather than regenerated through a tie-sensitive rematching step.

## Code/static checks

The final reviewer code package passed project prechecks for personal absolute paths, active ADAT reads, runtime package installation, transient patch files, data dimensions, and locked manuscript counts. A fresh R syntax parse checked 43 R files with zero failures. Python source files in the validated package were syntax-compiled successfully.

This record concerns the processed-data code interface and does not publish or embed the governed reviewer dataset.
