# Verified evidence ledger (R 4.6.1)

Generated: 2026-07-15 18:38:34 CST

All values below were read from the completed R 4.6.1 result tables. Claim levels distinguish observed, validated/reproduced, sensitivity, contextual, exploratory and negative findings.

## TCGA bulk RNA: Original versus ligand-excluded FAP-CAF score

- Estimate: rho=0.987; n=434
- Uncertainty: 95% CI 0.984 to 0.989; P=0.00e+00
- Multiplicity/unit: Not part of a confirmatory test family; score robustness check
- Interpretation: Removing four candidate ligands does not materially alter the stromal/ECM score.
- Claim level: `validated`
- Source: `L0_TCGA/score_correlations.csv`

## TCGA bulk RNA: Ligand-excluded FAP-CAF score versus CLDN core across all profiles

- Estimate: rho=-0.054; n=434
- Uncertainty: 95% CI -0.147 to 0.040; P=0.262
- Multiplicity/unit: Unadjusted prespecified correlation; primary-tumour sensitivity reported separately
- Interpretation: No positive bulk association was detected; the interval excludes a moderate positive correlation.
- Claim level: `observed`
- Source: `L0_TCGA/score_correlations.csv`

## TCGA bulk RNA: Primary-tumour-only ligand-excluded FAP-CAF versus CLDN core

- Estimate: rho=-0.152; n=380
- Uncertainty: P=0.003; BH-FDR=0.005
- Multiplicity/unit: BH correction across three prespecified sensitivity correlations
- Interpretation: A weak inverse, not positive, association was detected in primary tumours.
- Claim level: `observed`
- Source: `L0_TCGA_primary_sensitivity/score_correlations.csv`

## TCGA bulk RNA: Tumour versus adjacent-normal expression

- Estimate: FAP: median difference 4.20; CLDN1: median difference 5.00; CLDN2: median difference 6.55; CLDN4: median difference -0.30
- Uncertainty: FAP FDR=0.00e+00; CLDN1 FDR=0.00e+00; CLDN2 FDR=0.00e+00; CLDN4 FDR=1.52e-04
- Multiplicity/unit: BH correction across FAP, CLDN1, CLDN2 and CLDN4
- Interpretation: FAP, CLDN1 and CLDN2 were higher in tumour; CLDN4 was modestly lower.
- Claim level: `observed`
- Source: `L0_TCGA/tumor_normal_gene_comparisons.csv`

## TCGA CMS: FAP-CAF score across single-sample predictor CMS groups

- Estimate: n=314 classified tumours; CMS4 median=1.126
- Uncertainty: Kruskal-Wallis P=1.51e-28
- Multiplicity/unit: Overall test followed by Bonferroni-adjusted pairwise tests
- Interpretation: The FAP-CAF programme is enriched in CMS4; this is contextual rather than novel subtype discovery.
- Claim level: `validated`
- Source: `L0_TCGA/cms_kruskal_wallis.csv; cms_score_summary.csv`

## TCGA MSI: MANTIS-defined MSI context

- Estimate: n=359; ligand-excluded FAP-CAF FDR=0.105; CLDN core FDR=1.44e-05
- Uncertainty: MSI-H=61; MSS=298
- Multiplicity/unit: BH correction across three score comparisons
- Interpretation: The stromal score was not associated with MSI after correction, whereas the CLDN core was.
- Claim level: `observed`
- Source: `L0_TCGA/msi_score_comparisons.csv`

## CPTAC proteomics: FAP protein versus ECM proteins

- Estimate: COL1A1 rho=0.720 (n=97); COL1A2 rho=0.729 (n=97); FN1 rho=0.814 (n=97)
- Uncertainty: COL1A1 FDR=1.98e-16; COL1A2 FDR=8.08e-17; FN1 FDR=2.37e-23
- Multiplicity/unit: BH correction across six prespecified protein correlations
- Interpretation: FAP protein covaries strongly with collagen/fibronectin proteins.
- Claim level: `validated`
- Source: `CPTAC_protein/cptac_prespecified_correlations.csv`

## CPTAC proteomics: FAP/ECM protein measures versus available CLDN proteins

- Estimate: FAP vs CLDN4: rho=-0.090 (n=71); ECM_protein_score vs CLDN4: rho=-0.179 (n=71); ECM_protein_score vs CLDN1: rho=-0.088 (n=18)
- Uncertainty: FDR=0.549; FDR=0.202; FDR=0.729
- Multiplicity/unit: BH correction across six prespecified protein correlations; pairwise complete cases
- Interpretation: No positive association was detected, but CLDN1 coverage was sparse and CLDN2 was unavailable.
- Claim level: `observed`
- Source: `CPTAC_protein/cptac_prespecified_correlations.csv`

## GSE166555 single-cell pseudobulk: Paired tumour-normal compartment-specific expression

- Estimate: CAF_like FAP: median difference 2.264 (n=8); Epithelial CLDN1: median difference 2.648 (n=11); Epithelial CLDN2: median difference 3.259 (n=11); Epithelial CLDN4: median difference -0.433 (n=11)
- Uncertainty: FAP FDR=0.014; CLDN1 FDR=0.005; CLDN2 FDR=0.005; CLDN4 FDR=0.005
- Multiplicity/unit: BH correction across four paired tests
- Interpretation: CAF-like FAP, epithelial CLDN1 and CLDN2 increased; epithelial CLDN4 decreased.
- Claim level: `validated`
- Source: `L1_TISCH_GSE166555/paired_tumor_normal_tests.csv`

## GSE166555 single-cell pseudobulk: Primary cross-compartment patient-level association

- Estimate: rho=-0.382; n=10
- Uncertainty: bootstrap 95% CI -0.962 to 0.447; P=0.276; BH-FDR=0.750
- Multiplicity/unit: BH correction across 16 predictor-outcome combinations at the 20-cell threshold
- Interpretation: No positive patient-level coupling was detected; the small cohort leaves a wide interval.
- Claim level: `observed`
- Source: `L1_TISCH_GSE166555/patient_level_associations.csv`

## GSE166555 single-cell pseudobulk: Lower CAF-cell-threshold sensitivity association

- Estimate: rho=-0.391; n=11
- Uncertainty: bootstrap 95% CI -0.944 to 0.355; P=0.235; BH-FDR=0.716
- Multiplicity/unit: BH correction across 16 predictor-outcome combinations at the 5-cell threshold
- Interpretation: The direction remained inverse and non-significant after including one additional patient.
- Claim level: `sensitivity`
- Source: `L1_TISCH_GSE166555/patient_level_associations.csv`

## Qi spatial cohort: FAP spatial autocorrelation and FAP-high/epithelial-high co-occurrence

- Estimate: Moran I range 0.461-0.653; co-occurrence 3.34%, 0.67% and 0.24%
- Uncertainty: Each Moran permutation P=0.001 with 999 permutations
- Multiplicity/unit: Per-patient spatial tests; no pooled spot-level inference
- Interpretation: FAP-rich regions were organised, but upper-quartile overlap with epithelial-high spots was limited.
- Claim level: `validated`
- Source: `L2_spatial/spatial_morans_i.csv; spatial_fap_epithelial_cooccurrence.csv`

## Valdeolivas raw spatial: raw_bivariate_moran

- Estimate: median=-0.0464; n=7; positive=3; negative=4
- Uncertainty: two-sided exact Wilcoxon P=0.688; BH-FDR=0.688
- Multiplicity/unit: BH correction across three patient-level spatial metrics
- Interpretation: Raw scores did not show a consistent direction.
- Claim level: `observed`
- Source: `L2_Valdeolivas_spatial/spatial_patient_summary.csv`

## Valdeolivas raw spatial: adjusted_bivariate_moran

- Estimate: median=-0.0256; n=7; positive=2; negative=5
- Uncertainty: two-sided exact Wilcoxon P=0.078; BH-FDR=0.117
- Multiplicity/unit: BH correction across three patient-level spatial metrics
- Interpretation: No positive spatial coupling was detected; negative estimates are suggestive, not confirmatory, of compartmental separation.
- Claim level: `observed`
- Source: `L2_Valdeolivas_spatial/spatial_patient_summary.csv`

## Valdeolivas raw spatial: high_high_enrichment_minus_one

- Estimate: median=-0.1178; n=6; positive=1; negative=5
- Uncertainty: two-sided exact Wilcoxon P=0.062; BH-FDR=0.117
- Multiplicity/unit: BH correction across three patient-level spatial metrics
- Interpretation: No positive spatial coupling was detected; negative estimates are suggestive, not confirmatory, of compartmental separation.
- Claim level: `observed`
- Source: `L2_Valdeolivas_spatial/spatial_patient_summary.csv`

## GSE132465 CellChat: Full pooled FAP-high myofibroblast-to-epithelial analysis

- Estimate: 80 significant pairs; 29 collagen; 3 FN1
- Uncertainty: 100 bootstrap iterations; CellChat threshold P<0.05
- Multiplicity/unit: Exploratory pooled cell-level inference; no patient-level replication
- Interpretation: Co-expression-compatible communication was dominated by ECM-receptor pairs.
- Claim level: `exploratory`
- Source: `work/cellchat_reanalysis/results/cellchat_reanalysis_summary.csv`

## GSE132465 CellChat: SMC20-excluded sensitivity analysis

- Estimate: 93 significant pairs; 35 COLLAGEN-pathway; 4 FN1-pathway
- Uncertainty: 162 sender cells and 16,945 epithelial cells; 100 bootstrap iterations
- Multiplicity/unit: Sensitivity exclusion because SMC20 contributed 122/284 candidate senders
- Interpretation: All six prespecified collagen/fibronectin-SDC4/CD44 pairs remained among the strongest signals; no TGFB-family, INHBA or CCL2 ligand was detected in this direction.
- Claim level: `sensitivity`
- Source: `work/cellchat_reanalysis/results/cellchat_FAPhigh_to_epithelial_significant_excluding_SMC20.csv`

## Valdeolivas spatial CellChat: Prespecified pair reproducibility across paired sections

- Estimate: FN1 - CD44: both 5/6; any 6/6; COL1A1 - CD44: both 4/6; any 6/6; COL1A2 - CD44: both 4/6; any 6/6; COL1A1 - SDC4: both 6/6; any 6/6; COL1A2 - SDC4: both 6/6; any 6/6; FN1 - SDC4: both 6/6; any 6/6
- Uncertainty: 12 sections from six patients; 100 bootstrap iterations per section
- Multiplicity/unit: BH correction across all six prespecified pairs within each section
- Interpretation: SDC4 pairs reproduced in both sections for 6/6 patients; CD44 pairs were less consistent.
- Claim level: `computationally_reproduced`
- Source: `L3_spatial_CellChat/spatial_cellchat_pair_prevalence.csv`

## NicheNet: Receiver differential-expression gate

- Estimate: 12,391 genes tested; 0 positive genes at FDR<0.10
- Uncertainty: eligible patients=13 (6 high, 7 low)
- Multiplicity/unit: edgeR BH-FDR; ligand ranking used top 200 unadjusted positive genes only
- Interpretation: Ligand ranking is hypothesis-generating and cannot support a direct mechanism.
- Claim level: `negative`
- Source: `L3_NicheNet/epithelial_DE_summary.csv`

## PROGENy/DoRothEA: Patient-level pathway and EMT-related TF contrasts

- Estimate: PROGENy minimum group FDR=0.824; selected DoRothEA minimum group FDR=0.709
- Uncertainty: n=13 patients
- Multiplicity/unit: BH correction performed separately for group contrasts and sender-burden correlations
- Interpretation: No tested pathway or selected EMT-related TF survived correction.
- Claim level: `negative`
- Source: `L3_pathway_TF/progeny_high_vs_low_and_burden.csv; dorothea_selected_emt_tfs.csv`

## TCGA RPPA: TGF-beta RNA score versus SMAD3, SMAD4 and E-cadherin protein

- Estimate: SMAD3_RPPA rho=0.222; SMAD4_RPPA rho=0.224; CDH1_RPPA rho=-0.346
- Uncertainty: FDR=6.81e-05; FDR=6.81e-05; FDR=7.91e-10
- Multiplicity/unit: BH correction across three proteins; n=316
- Interpretation: The expected TGF-beta/SMAD context was internally consistent, but direction from FAP-high CAFs was not tested.
- Claim level: `contextual`
- Source: `RPPA/rppa_tgf_beta_correlations_from_raw_api.csv`

## UALCAN CPTAC interface: Primary tumour versus normal display summaries

- Estimate: FAP: Higher in primary tumour (median difference 1.375); COL1A1: Lower in primary tumour (median difference -1.531); COL1A2: Lower in primary tumour (median difference -1.385); FN1: Higher in primary tumour (median difference 0.109); CLDN1: Higher in primary tumour (median difference 0.889); CLDN4: Higher in primary tumour (median difference 0.444)
- Uncertainty: P=6.23e-26; P=3.51e-25; P=7.25e-25; P=0.914; P=2.83e-04; P=0.015
- Multiplicity/unit: Website-reported unpaired t tests; same CPTAC source as the primary protein analysis
- Interpretation: A reproducibility check of the public interface, not an independent cohort.
- Claim level: `contextual`
- Source: `UALCAN_CPTAC/ualcan_cptac_primary_vs_normal.csv`

