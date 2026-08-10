# Legacy script status

The `scripts/legacy_exploratory/` directory preserves all historical R scripts recovered from the project. These files are included for auditability, not as a second authoritative pipeline.

Do not use the following outputs to support the revised manuscript:

- `REVIEW_R6_nomogram_mixedmodel_side.R`: the nodal nomogram was developed in 89 patients without independent validation and is removed.
- `11_A7_ECM_TSR_LN_risk.R` and related TSR scripts: transcriptomic nodal analyses are retained only as negative sensitivity analyses; unreproducible pathological TSR data are excluded.
- `REVIEW_R13_senescence_expansion.R`: superseded by `scripts/review_corrected/REVIEW_R13_senescence_expansion_tumor566.R`; the original target-probe extraction did not include the full CellAge gene set.
- cell-level mixed-model P values calculated with `df = n_cells - 2`: not used for inference.
- FAP differential expression when FAP itself defines the positive group: classifier leakage; FAP is excluded from the displayed differential-expression result.

Historical scripts were path-scrubbed to remove workstation-specific information. Some still expect files or packages used in earlier development and are not guaranteed to run without adaptation.
