# FAP–CAF–CRC manuscript: frozen hypothesis, outline and evidence ledger

Date frozen: 2026-08-14  
Purpose: governing blueprint for the rebuilt AJCR manuscript. Earlier manuscripts (v5.5, v5.7, v6.0 and candidate v6.1) are source materials only; none is an authoritative narrative or numerical source.

## 1. Scientific question and claim hierarchy

### Central question

How are senescence-associated transcriptional features distributed between tumour fibroblast-lineage and epithelial compartments in colorectal cancer, and are those features specifically attributable to FAP-detected fibroblasts or coupled to the FAP-associated matrix program?

### Primary supported claim

Tumour fibroblast-lineage cells show higher SenMayo transcriptional scores than patient-paired epithelial cells in GSE132465, and the direction is reproduced in independent GSE166555 pseudobulks. This is evidence for a fibroblast-enriched senescence-associated transcriptional pattern, not a diagnosis of cellular senescence.

### Secondary supported claim

Within tumour fibroblasts, FAP detection is not independently associated with SenMayo or SASP scores after sequencing-depth and subtype adjustment. FAP instead tracks a matrix-rich fibroblast state.

### Composition-aware claim

The strong bulk FAP13–matrix4 association is highly sensitive to fibroblast-content adjustment. SenMayo–FAP13 covariation is more stable, whereas SenMayo–matrix4 coupling is not supported. The FAP–matrix and senescence-associated transcriptional axes are therefore partially separable.

### Explicit evidence ceiling

Transcriptomic signatures do not establish durable arrest, senescent-cell identity, causal SASP signalling, immune dysfunction, drug resistance or a stromal–epithelial “division of labour”. These remain hypotheses requiring orthogonal tissue and functional validation.

## 2. Manuscript architecture

1. **Title** — lead with the reproducible stromal transcriptional pattern and separation from FAP–matrix biology; avoid “FAP+ senescent CAFs” and “compartmentalized senescence”.
2. **Abstract** — problem; tumour-only paired and external evidence; FAP-specific null result; composition-aware matrix result; bounded conclusion.
3. **Introduction** — CRC stromal heterogeneity; why tissue senescence requires convergent markers; recent CRC senescent-CAF studies; unresolved distinction among lineage, FAP state and senescence-associated transcription; prespecified questions.
4. **Methods** — datasets and cell universes; frozen signatures; patient-level paired analysis; expression-matched nulls; external pseudobulk; depth/subtype-adjusted FAP models; composition-sensitive bulk analysis; protein and ligand–immune boundary analyses.
5. **Results** — begin with tumour-only compartment signal; test specificity and external reproducibility; then show why it is not FAP-specific; establish FAP–matrix axis; quantify bulk composition sensitivity; report negative protein/receptor/immune boundaries.
6. **Discussion** — positive finding first; distinction from FAP-specific senescence; implications of two separable stromal axes; relation to 2026 CRC CAF-senescence studies; limitations in one consolidated section; minimum decisive validation experiment.
7. **Conclusion** — one positive, one boundary and one next-step sentence; no rhetorical overreach.

## 3. Frozen evidence ledger

| ID | Claim | Frozen data universe / model | Frozen result | Allowed wording | Prohibited wording | Main placement |
|---|---|---|---|---|---|---|
| C1 | Fibroblast compartment has higher SenMayo transcription | GSE132465 tumour only; 1,501 fibroblast-lineage and 17,469 epithelial cells; 23 patient-paired summaries; frozen overlap-removed SenMayo119; exact paired Wilcoxon; patient bootstrap CI | 23/23 higher; median standardized difference 0.466; 95% CI 0.410–0.553; P=2.38×10^-7 | fibroblast-enriched senescence-associated transcription | compartmentalized senescence; stromal senescence established | Abstract, Results 3.2, Discussion 4.1 |
| C2 | Cell-count sensitivity | GSE132465 tumour only; patients with ≥20 cells in both compartments; n=15 | median difference 0.466; 95% CI 0.410–0.551; P=6.10×10^-5 | direction retained under cell-count threshold | proof independent of sampling | Results 3.2 |
| C3 | SASP and proliferation directions | Same tumour-only paired universe | SASP25 higher in 22/23, mean difference 0.236, P=4.77×10^-7; MKI67 lower in 21/23, mean difference -0.104, P=1.81×10^-4 | concordant supporting transcriptional pattern | stable arrest; functional SASP secretion | Abstract (brief), Results 3.2 |
| C4 | Matched-null specificity | GSE132465; 5,000 sets matched for mean, SD and detection | SenMayo observed 0.499; null median 0.219; 97.5th percentile 0.340; empirical P=2.00×10^-4. SASP observed 0.739; null median 0.405; 97.5th percentile 0.759; P=0.0336 | SenMayo magnitude exceeded matched transcriptional background | null eliminates lineage biology; SASP is highly specific | Results 3.3 |
| C5 | External compartment direction | GSE166555 raw counts; 10 eligible tumour patients; patient-compartment pseudobulks; frozen SenMayo119 without re-selection | SenMayo 9/10 higher; median difference 0.1485; 95% CI 0.0438–0.3364; P=0.01953125. SASP 10/10 higher, difference 0.4754, P=0.001953. MKI67 9/10 lower, difference -2.033, P=0.003906 | independent reproduction/evaluation of direction | validation of cellular senescence; interchangeable effect size | Abstract, Results 3.4, Discussion 4.1 |
| C6 | FAP-specific SenMayo/SASP association is null | GSE132465 tumour fibroblasts only, n=1,501; log UMI depth + fibroblast subtype + MKI67 + patient clustering | SenMayo beta 0.0238, 95% CI -0.0051–0.0527, P=0.106; SASP beta 0.0323, 95% CI -0.0429–0.1074, P=0.400 | no independent association after adjustment | FAP+ CAFs carry a senescence program | Abstract, Results 3.5, Discussion 4.2 |
| C7 | Individual markers do not rescue FAP-specific claim | Tumour fibroblast marker-detection logistic mixed models adjusted for depth/subtype; patient clustering | CDKN2A OR 0.89 P=0.529; CDKN2B 0.98 P=0.911; CDKN1A 0.96 P=0.788; LMNB1 0.86 P=0.601 (singular); MKI67 0.84 P=0.614 | adjusted marker models were null; LMNB1 model singular | marker concordance establishes senescence | Results 3.5, Limitations |
| C8 | FAP tracks matrix within tumour fibroblasts | GSE132465 tumour fibroblasts; patient summaries and full-library mixed model | FAP–matrix4 rho=0.649 P=8.03×10^-4; matrix4~FAP slope 0.253, 95% CI 0.213–0.294, P=9.17×10^-33 | FAP-associated matrix state | FAP causes ECM remodelling | Results 3.6, Discussion 4.3 |
| C9 | Receptor co-induction is unsupported | Same full-library mixed-model universe | receptor2~FAP slope -0.021, 95% CI -0.062–0.019, P=0.302 | no SDC4/CD44 co-induction in this analysis | receptor function excluded | Results 3.6 |
| C10 | SenMayo and matrix axes are separable in fibroblasts | GSE132465 patient fibroblast summaries | FAP–SenMayo rho=0.331 P=0.123; SenMayo–matrix4 rho=-0.048 P=0.826 | patient-level axes were not correlated | biological independence proven | Results 3.6, Discussion 4.3 |
| C11 | Bulk FAP–matrix association is composition sensitive | TCGA-COAD/READ; rank partial correlations; V5.7 full-library pipeline | marginal 0.930; fib5 0.650; MCP-counter 0.113; EPIC 0.288; joint 0.150 | association attenuated and proxy dependent | survived fibroblast adjustment (without qualification) | Results 3.7, Discussion 4.4 |
| C12 | Bulk SenMayo–FAP covariation is more stable than SenMayo–matrix | TCGA; overlap-removed signatures; four fibroblast specifications | SenMayo–FAP13 0.378–0.575; SenMayo–matrix4 -0.098 to 0.133, with zero included under three of four specifications | more stable covariation; no independent senescence–matrix axis | FAP is a senescence biomarker; causal coupling | Results 3.7 |
| C13 | Protein evidence is restricted to FAP–ECM | CPTAC colon proteomics; FAP-excluded ECM score | supported FAP–ECM covariation; SASP–ECM not supported | protein-level support for FAP–ECM only | proteomic validation of senescence | Results 3.8 |
| C14 | No frozen SASP–immune pair | Prespecified cognate SASP ligand–immune receptor screen; 14 comparisons | minimum adjusted q=0.598; no pair survived | exploratory screen was null | CellChat validation; immune mechanism established | Results 3.8, Limitations |
| C15 | External separation of axes | GSE166555 eligible fibroblast pseudobulks, n=10 | FAP13–matrix4 rho=0.8788 P=0.0008139; FAP13–SenMayo rho=-0.0667 P=0.8548; SenMayo–matrix rho=-0.4061 P=0.2443 | external data reproduced FAP–matrix coupling and lack of FAP–SenMayo coupling | external causal validation | Results 3.4/3.6 |

## 4. Version-trace rules

| Source version | Permitted use | Excluded by default |
|---|---|---|
| v5.5 | historical wording and identification of analyses to audit | 0.112/+0.051 receptor-era mixed-model values; fib5-only interpretation |
| v5.7 | authoritative full-library matrix/receptor estimates and composition-aware bulk sensitivity | conclusions that ignore the new tumour-only senescence analysis |
| v6.0 | hypotheses, candidate section architecture and list of analyses | all-sample 3,462/18,539 B1/B2 coefficients; 120-gene count; “FAP+ senescent CAF” framing |
| candidate v6.1 | isolated wording improvements and submission metadata after verification | title, abstract, Methods 2.13, Results 3.12 and Cover Letter claims derived from mixed tumour/normal cells |
| rebuilt manuscript | current narrative source | any statement not traceable to C1–C15 or a verified background citation |

## 5. Voice and convergence rules

1. Lead each Results subsection with the supported observation; place the evidential boundary once at its end.
2. Use “senescence-associated transcription” for score-based findings and reserve “senescent cells” for cited experimental literature or proposed validation.
3. Avoid repeated defensive strings. Prefer neutral contrasts such as “The adjusted estimate was compatible with no FAP-specific enrichment” and “This separates the two transcriptional axes.”
4. Keep audit history out of the manuscript body. Explain the final cell universe and model directly; document version corrections only in the audit report.
5. Use one term per construct: tumour fibroblast-lineage cells; epithelial cells; FAP-detected/FAP-undetected; SenMayo119; SASP25; matrix4; receptor2.
6. Abstract, Results, tables, legends and Cover Letter must draw numbers from this ledger, not from prose copied from an earlier version.
7. References are added only when they support a specific sentence in the frozen outline; citation count is not an objective.

## 6. Final consistency gates

- Every numerical claim maps to C1–C15 and one frozen output file.
- No 3,462/18,539 all-sample cell universe is described as tumour-only CRC analysis.
- No beta=0.083, P=6.96×10^-37, paired difference 0.376 or adjusted beta=0.361 remains.
- No matrix slope 0.119 or receptor slope -0.005 remains in the final full-library model description.
- SenMayo source arithmetic is 125 minus four overlaps = 121 candidates; two absent genes = 119 represented genes.
- Title, Abstract, Discussion, Conclusion, figures and Cover Letter express the same claim hierarchy.
- Citations are verified against primary bibliographic records and appear in first-use order.
- The final DOCX is rendered page by page after regeneration.
