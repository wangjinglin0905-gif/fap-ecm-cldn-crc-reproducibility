#!/usr/bin/env python3
"""Build source-linked numeric and claim locks for the PLOS ONE v7 rerun.

The script reads only frozen rerun tables and writes compact ledgers used for
manuscript reconstruction. It does not recalculate statistics.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RERUN = PROJECT / "results"
LEDGERS = PROJECT / "ledgers"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def choose(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    hits = [row for row in rows if all(row.get(key) == value for key, value in criteria.items())]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one row for {criteria}; found {len(hits)}")
    return hits[0]


def rel(path: Path) -> str:
    return path.relative_to(PROJECT).as_posix()


def fmt(value: str | float | int | None, digits: int = 6) -> str:
    if value in (None, "", "NA"):
        return ""
    number = float(value)
    return f"{number:.{digits}g}"


single_dir = RERUN / "single_cell"
bulk_dir = RERUN / "bulk_composition"
common_dir = RERUN / "common_core"
null_dir = RERUN / "spatial_null"

aligned_path = single_dir / "scrna_paired_compartment_summary_aligned.csv"
cluster_path = single_dir / "gse132465_cluster_inference.csv"
contrast_path = single_dir / "dependent_correlation_contrasts.csv"
bulk_path = bulk_dir / "target_purged_partial_correlations.csv"
overlap_path = bulk_dir / "composition_proxy_feature_overlap_audit.csv"
spatial_null_path = null_dir / "spatial_matched_null_corrected.csv"
g132_core_path = common_dir / "gse132_common_core_summary.json"
g166_core_path = common_dir / "gse166_common_core_summary.csv"
spatial_core_path = common_dir / "spatial_common_core_summary.csv"
coverage_path = common_dir / "SenMayo_common_core_coverage.csv"
normalization_path = common_dir / "gse132_public_normalization_calibration.csv"

aligned = read_csv(aligned_path)
cluster = read_csv(cluster_path)
contrasts = read_csv(contrast_path)
bulk = read_csv(bulk_path)
overlaps = read_csv(overlap_path)
spatial_null = read_csv(spatial_null_path)
g132_core = json.loads(g132_core_path.read_text(encoding="utf-8"))
g166_core = read_csv(g166_core_path)[0]
spatial_core = read_csv(spatial_core_path)
coverage = read_csv(coverage_path)
normalization = read_csv(normalization_path)

numeric_rows: list[dict[str, str]] = []


def add_numeric(
    lock_id: str,
    domain: str,
    dataset: str,
    estimand: str,
    n_unit: str,
    estimate: str,
    ci_low: str = "",
    ci_high: str = "",
    p_value: str = "",
    direction_count: str = "",
    resampling: str = "",
    source: Path | None = None,
    status: str = "locked",
    interpretation_ceiling: str = "",
) -> None:
    numeric_rows.append(
        {
            "lock_id": lock_id,
            "domain": domain,
            "dataset": dataset,
            "estimand": estimand,
            "n_unit": n_unit,
            "estimate": estimate,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "p_value": p_value,
            "direction_count": direction_count,
            "resampling": resampling,
            "source_artifact": rel(source) if source else "",
            "status": status,
            "interpretation_ceiling": interpretation_ceiling,
        }
    )


for index, cohort in enumerate(("GSE132465", "GSE166555"), start=1):
    row = choose(aligned, cohort=cohort, endpoint="SenMayo_score")
    add_numeric(
        f"L{index:03d}",
        "single-cell primary",
        cohort,
        "available-gene SenMayo fibroblast-minus-epithelial paired median",
        f"{row['n_patients']} patients",
        fmt(row["median_paired_difference"]),
        p_value=fmt(row["p_two_sided"]),
        direction_count=f"{row['positive_direction']}/{row['n_patients']} positive",
        source=aligned_path,
        interpretation_ceiling="compartment-level senescence-associated transcription; not a senescence diagnosis",
    )

add_numeric(
    "L003",
    "common-core sensitivity",
    "GSE132465",
    "111-gene common-core fibroblast-minus-epithelial paired mean",
    f"{g132_core['patients']} patients",
    fmt(g132_core["paired_mean_difference"]),
    fmt(g132_core["mean_difference_ci_low"]),
    fmt(g132_core["mean_difference_ci_high"]),
    fmt(g132_core["wilcoxon_signed_rank_p"]),
    f"{g132_core['positive_differences']}/{g132_core['patients']} positive",
    "10,000 patient bootstraps for mean CI",
    g132_core_path,
    interpretation_ceiling="sensitivity using the same 111 genes across cohorts",
)
add_numeric(
    "L004",
    "common-core sensitivity",
    "GSE166555",
    "111-gene common-core fibroblast-minus-epithelial paired mean",
    f"{g166_core['patients']} patients",
    fmt(g166_core["paired_mean_difference"]),
    fmt(g166_core["mean_difference_ci_low"]),
    fmt(g166_core["mean_difference_ci_high"]),
    fmt(g166_core["wilcoxon_signed_rank_p"]),
    f"{g166_core['positive_differences']}/{g166_core['patients']} positive",
    "10,000 patient bootstraps for mean CI",
    g166_core_path,
    interpretation_ceiling="sensitivity using the same 111 genes across cohorts",
)

cluster_models = (
    ("L005", "GSE132_FAPdet_SenMayo_depth_subtype", "FAP detection coefficient for SenMayo", "no supported additional enrichment; not equivalence"),
    ("L006", "GSE132_FAPdet_SASP_depth_subtype", "FAP detection coefficient for SASP", "no supported additional enrichment; not equivalence"),
    ("L007", "GSE132_matrix4_on_FAP_depth_subtype", "FAP-expression coefficient for matrix4", "adjusted within-fibroblast covariation; not causality"),
)
for lock_id, model_id, label, ceiling in cluster_models:
    row = choose(
        cluster,
        model_id=model_id,
        estimator="patient fixed effects with patient-cluster percentile bootstrap",
    )
    add_numeric(
        lock_id,
        "cluster-aware single-cell model",
        "GSE132465",
        label,
        f"{row['clusters']} patient clusters",
        fmt(row["estimate"]),
        fmt(row["ci_low"]),
        fmt(row["ci_high"]),
        fmt(row["p_value"]),
        resampling=f"{row['bootstrap_valid']} valid patient-cluster bootstraps",
        source=cluster_path,
        interpretation_ceiling=ceiling,
    )

contrast_specs = (
    ("L008", "GSE132465", "rho_FAP_matrix", "FAP-matrix Spearman rho"),
    ("L009", "GSE132465", "FAP_matrix_minus_FAP_SenMayo", "rho(FAP,matrix) minus rho(FAP,SenMayo)"),
    ("L010", "GSE132465", "FAP_matrix_minus_SenMayo_matrix", "rho(FAP,matrix) minus rho(SenMayo,matrix)"),
    ("L011", "GSE166555", "rho_FAP_matrix", "FAP-matrix Spearman rho"),
    ("L012", "GSE166555", "FAP_matrix_minus_FAP_SenMayo", "rho(FAP,matrix) minus rho(FAP,SenMayo)"),
    ("L013", "GSE166555", "FAP_matrix_minus_SenMayo_matrix", "rho(FAP,matrix) minus rho(SenMayo,matrix)"),
)
for lock_id, cohort, estimand, label in contrast_specs:
    row = choose(contrasts, cohort=cohort, estimand=estimand)
    ceiling = (
        "positive FAP-matrix association"
        if estimand == "rho_FAP_matrix"
        else "direct contrast; CI crossing zero does not establish a difference"
        if estimand == "FAP_matrix_minus_FAP_SenMayo"
        else "FAP-matrix exceeds SenMayo-matrix on this patient-level scale"
    )
    add_numeric(
        lock_id,
        "dependent-correlation bootstrap",
        cohort,
        label,
        f"{row['n_patients']} patients",
        fmt(row["estimate"]),
        fmt(row["ci_low"]),
        fmt(row["ci_high"]),
        fmt(row["bootstrap_sign_p"]),
        resampling=f"{row['bootstrap_valid']} paired patient bootstraps",
        source=contrast_path,
        interpretation_ceiling=ceiling,
    )

bulk_specs = (
    ("L014", "none", "unadjusted FAP13-matrix4 rho"),
    ("L015", "fib5 transcript score", "fib5-adjusted partial rho"),
    ("L016", "MCP+EPIC recomputed full", "joint full-proxy partial rho"),
    ("L017", "MCP+EPIC pair-purged", "joint pair-purged partial rho"),
    ("L018", "MCP+EPIC global-disjoint", "joint globally disjoint partial rho"),
)
for lock_id, adjustment, label in bulk_specs:
    row = choose(bulk, pair="FAP13_matrix4", adjustment=adjustment)
    add_numeric(
        lock_id,
        "bulk composition sensitivity",
        "TCGA COAD/READ",
        label,
        f"{row['n']} primary tumours",
        fmt(row["rho"]),
        fmt(row["ci_low"]),
        fmt(row["ci_high"]),
        fmt(row["p_value"]),
        resampling="5,000 bootstrap CIs for adjusted estimates",
        source=bulk_path,
        interpretation_ceiling="positive residual association after target purging; magnitude remains proxy-sensitive",
    )

row = choose(bulk, pair="SenMayo_matrix4", adjustment="MCP+EPIC global-disjoint")
add_numeric(
    "L019",
    "bulk composition sensitivity",
    "TCGA COAD/READ",
    "SenMayo-matrix4 globally disjoint partial rho",
    f"{row['n']} primary tumours",
    fmt(row["rho"]),
    fmt(row["ci_low"]),
    fmt(row["ci_high"]),
    fmt(row["p_value"]),
    resampling="5,000 bootstrap CIs",
    source=bulk_path,
    interpretation_ceiling="estimate compatible with no residual monotonic association; not equivalence",
)

for offset, cohort in enumerate(("GSE280315", "Valdeolivas_Visium", "GSE334323"), start=20):
    row = choose(spatial_core, cohort=cohort)
    add_numeric(
        f"L{offset:03d}",
        "spatial common-core sensitivity",
        cohort,
        "111-gene stromal-minus-tumour patient effect",
        f"{row['patients']} patients/samples",
        fmt(row["mean_patient_effect"]),
        p_value=fmt(row["wilcoxon_signed_rank_p"]),
        direction_count=f"{row['positive_patients']}/{row['patients']} positive",
        source=spatial_core_path,
        interpretation_ceiling="cohort-specific heterogeneous stress test; no cross-platform pooling",
    )

for offset, cohort in enumerate(("GSE280315", "Valdeolivas_Visium"), start=23):
    row = choose(spatial_null, cohort=cohort)
    add_numeric(
        f"L{offset:03d}",
        "spatial matched-null correction",
        cohort,
        "observed minus empirical-null median",
        "5,000 matched sets",
        fmt(row["observed_minus_null_center"]),
        p_value=fmt(row["empirical_centered_absolute_p"]),
        resampling=f"{row['draws']} matched random sets",
        source=spatial_null_path,
        interpretation_ceiling="superiority to the tested centered matched-null distribution; supplementary only",
    )

norm = choose(normalization, candidate="log1p_CPTT_1e4")
add_numeric(
    "L025",
    "public-data reconstruction",
    "GSE132465",
    "correlation of reconstructed and archived SenMayo scores",
    f"{norm['n_cells']} cells",
    fmt(norm["pearson_vs_archived_SenMayo_zmean"]),
    source=normalization_path,
    interpretation_ceiling="exact score reconstruction under the selected normalization",
)

overlap_mcp = choose(overlaps, proxy="MCP-counter Fibroblasts", target_set="FAP13_matrix4")
overlap_epic = choose(overlaps, proxy="EPIC TRef signature genes", target_set="FAP13_matrix4")
add_numeric(
    "L026",
    "composition incorporation-bias audit",
    "TCGA COAD/READ",
    "proxy features overlapping FAP13 or matrix4",
    "MCP-counter and EPIC definitions",
    f"MCP {overlap_mcp['overlap_count']}; EPIC {overlap_epic['overlap_count']}",
    source=overlap_path,
    interpretation_ceiling=(
        f"MCP overlap: {overlap_mcp['overlap_genes']}; EPIC overlap: {overlap_epic['overlap_genes']}"
    ),
)

numeric_out = LEDGERS / "numeric_lock_v01_2026-09-04.csv"
with numeric_out.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(numeric_rows[0]))
    writer.writeheader()
    writer.writerows(numeric_rows)

claim_rows = [
    {
        "claim_id": "CL01",
        "status": "locked",
        "locked_wording": "Senescence-associated transcription was enriched in fibroblast-lineage relative to epithelial compartments in both scRNA-seq cohorts.",
        "prohibited_extension": "Do not call the score bona fide cellular senescence or stable growth arrest.",
        "evidence_anchor": "L001-L004",
    },
    {
        "claim_id": "CL02",
        "status": "locked null boundary",
        "locked_wording": "Patient-cluster analyses did not support additional SenMayo or SASP enrichment among FAP-detected tumour fibroblasts.",
        "prohibited_extension": "Do not claim equivalence, absence, or a uniquely senescent FAP-positive subgroup.",
        "evidence_anchor": "L005-L006",
    },
    {
        "claim_id": "CL03",
        "status": "locked",
        "locked_wording": "FAP expression covaried with matrix4 within tumour fibroblasts after depth, subtype, patient and cluster-aware uncertainty were considered.",
        "prohibited_extension": "Do not claim causality or a discrete FAP-matrix cell state.",
        "evidence_anchor": "L007-L008;L011",
    },
    {
        "claim_id": "CL04",
        "status": "locked qualified",
        "locked_wording": "FAP-matrix correlation exceeded SenMayo-matrix correlation in both scRNA cohorts, but its difference from FAP-SenMayo correlation remained imprecise.",
        "prohibited_extension": "Do not claim two statistically separable or independent biological states.",
        "evidence_anchor": "L009-L010;L012-L013",
    },
    {
        "claim_id": "CL05",
        "status": "locked qualified",
        "locked_wording": "TCGA FAP13-matrix4 covariation remained positive after target-purged composition adjustment, while its magnitude varied materially across proxy definitions.",
        "prohibited_extension": "Do not call the association composition-driven or composition-independent.",
        "evidence_anchor": "L014-L019;L026",
    },
    {
        "claim_id": "CL06",
        "status": "locked heterogeneous",
        "locked_wording": "The 111-gene spatial sensitivity showed positive effects in GSE280315 and Valdeolivas but mixed/negative effects in GSE334323.",
        "prohibited_extension": "Do not claim universal spatial validation, replication, or a pooled cross-platform effect.",
        "evidence_anchor": "L020-L024",
    },
    {
        "claim_id": "CL07",
        "status": "repository only",
        "locked_wording": "The epithelial-only organoid analysis is retained as a post hoc boundary control and is not part of the article's evidentiary chain.",
        "prohibited_extension": "Do not use epithelial organoids to validate CAF senescence or compartment localization.",
        "evidence_anchor": "Nature2026_organoid_boundary/analysis_summary.json",
    },
    {
        "claim_id": "CL08",
        "status": "discussion/future experiment only",
        "locked_wording": "Lipid transfer/peroxidation, trogocytosis and mitoxyperilysis may motivate event-resolving co-culture experiments but are not findings of this transcriptomic study.",
        "prohibited_extension": "Do not add inferred mechanism signatures or name these processes in Results without direct event-level measurements.",
        "evidence_anchor": "mechanism routing report",
    },
]

claim_out = LEDGERS / "claim_lock_v01_2026-09-04.csv"
with claim_out.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(claim_rows[0]))
    writer.writeheader()
    writer.writerows(claim_rows)

print(numeric_out)
print(claim_out)
