from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def split_genes(value: str | float) -> set[str]:
    if pd.isna(value) or not str(value).strip():
        return set()
    return {item.strip() for item in str(value).split(";") if item.strip()}


def relative_to_project(path: Path, project: Path) -> str:
    try:
        return path.resolve().relative_to(project.resolve()).as_posix()
    except ValueError:
        return path.name


def build_s1(project: Path, qa_root: Path, route_root: Path, out_dir: Path) -> Path:
    analysis = project / "results"
    key_files = {
        "TCGA-COAD/READ": analysis / "bulk_composition" / "tcga_primary_target_scores.csv",
        "GSE39582": analysis / "bulk_composition" / "gse39582_tumour_sample_scores.csv",
        "GSE132465": analysis / "single_cell" / "scrna_patient_compartment_scores_aligned.csv",
        "GSE166555": analysis / "single_cell" / "scrna_patient_compartment_scores_aligned.csv",
        "CPTAC-COAD": qa_root / "repro_cptac" / "cptac_fap_ecm_and_receptor_frozen.csv",
        "GSE280315": route_root / "GSE280315" / "patient_H1_compartment.csv",
        "Valdeolivas": route_root / "Valdeolivas_Visium" / "patient_H1_compartment.csv",
        "GSE334323": route_root / "GSE334323" / "patient_H1_compartment.csv",
    }
    rows = [
        {
            "dataset": "TCGA-COAD/READ",
            "repository_record": "TCGA-COAD and TCGA-READ",
            "modality": "bulk RNA sequencing",
            "analysis_role": "composition-proxy sensitivity",
            "source_population": "primary colon and rectal adenocarcinomas",
            "analysis_population": "380 barcode-reviewed primary tumours; one per patient",
            "independent_unit": "patient/tumour",
            "eligibility_or_exclusion": "primary-tumour barcode review and common-gene availability",
            "compartment_or_label_source": "not applicable",
        },
        {
            "dataset": "GSE39582",
            "repository_record": "GEO GSE39582",
            "modality": "microarray",
            "analysis_role": "external bulk marginal replication",
            "source_population": "585 profiled samples",
            "analysis_population": "566 tumours",
            "independent_unit": "tumour sample",
            "eligibility_or_exclusion": "19 non-tumour mucosa samples excluded by frozen sample map",
            "compartment_or_label_source": "discovery/validation sample map",
        },
        {
            "dataset": "GSE132465",
            "repository_record": "GEO GSE132465",
            "modality": "single-cell RNA sequencing",
            "analysis_role": "primary paired compartment localization and FAP-specific models",
            "source_population": "63,689 annotated cells from 23 patients",
            "analysis_population": "1,501 tumour fibroblast-lineage and 17,469 tumour epithelial cells",
            "independent_unit": "patient",
            "eligibility_or_exclusion": "source tumour labels; paired patient compartments",
            "compartment_or_label_source": "source Myofibroblasts/Stromal 1-3 and epithelial annotations",
        },
        {
            "dataset": "GSE166555",
            "repository_record": "GEO GSE166555",
            "modality": "single-cell RNA sequencing",
            "analysis_role": "external paired compartment evaluation",
            "source_population": "68,702 annotated cells from 12 tumour patients",
            "analysis_population": "10 eligible paired patients",
            "independent_unit": "patient",
            "eligibility_or_exclusion": "at least 20 CAF-like cells and 100 epithelial cells per patient",
            "compartment_or_label_source": "source CAF-like and epithelial annotations",
        },
        {
            "dataset": "CPTAC-COAD",
            "repository_record": "CPTAC colon cancer via cBioPortal",
            "modality": "mass-spectrometry proteomics",
            "analysis_role": "orthogonal FAP-matrix protein triangulation",
            "source_population": "colon tumours with available protein profiles",
            "analysis_population": "97 tumours with all required proteins",
            "independent_unit": "tumour sample",
            "eligibility_or_exclusion": "complete FAP, COL1A1, COL1A2, FN1, SDC4 and CD44 measurements",
            "compartment_or_label_source": "not applicable",
        },
        {
            "dataset": "GSE280315",
            "repository_record": "GEO GSE280315",
            "modality": "Visium HD spatial transcriptomics",
            "analysis_role": "spatial generalizability stress test",
            "source_population": "three colorectal tumour sections",
            "analysis_population": "three patients/sections meeting compartment-count criteria",
            "independent_unit": "patient/section",
            "eligibility_or_exclusion": "at least 30 assigned observations per compartment",
            "compartment_or_label_source": "author deconvolution singlet labels",
        },
        {
            "dataset": "Valdeolivas",
            "repository_record": "Zenodo 10.5281/zenodo.7760264",
            "modality": "Visium spatial transcriptomics",
            "analysis_role": "spatial generalizability stress test",
            "source_population": "replicate colorectal tumour sections",
            "analysis_population": "12 eligible sections from six patients",
            "independent_unit": "patient",
            "eligibility_or_exclusion": "replicate sections averaged within patient; at least 30 spots per compartment",
            "compartment_or_label_source": "pathologist Fibroblastic_stroma and Tumor labels",
        },
        {
            "dataset": "GSE334323",
            "repository_record": "GEO GSE334323; associated bioRxiv preprint",
            "modality": "FFPE Visium CytAssist spatial transcriptomics",
            "analysis_role": "independent marker-defined spatial sensitivity",
            "source_population": "three colorectal tumour sections",
            "analysis_population": "three tumour samples after host-feature denominator correction",
            "independent_unit": "tumour sample",
            "eligibility_or_exclusion": "at least 30 assigned spots per compartment; human ENSG features retained",
            "compartment_or_label_source": "frozen non-overlapping lineage-marker margin",
        },
    ]
    for row in rows:
        path = key_files[row["dataset"]]
        if not path.exists():
            raise FileNotFoundError(path)
        row["key_result_artifact"] = relative_to_project(path, project)
        row["key_artifact_sha256"] = sha256(path)
    output = out_dir / "S1_Table_dataset_accessions_eligibility_and_checksums.csv"
    pd.DataFrame(rows).to_csv(output, index=False)
    return output


def build_s2(project: Path, qa_root: Path, out_dir: Path) -> Path:
    availability = pd.read_csv(qa_root / "repro_gse132465" / "gene_availability.csv")
    availability_map = {row["set"]: split_genes(row["genes"]) for _, row in availability.iterrows()}
    coverage = pd.read_csv(
        project / "results" / "common_core" / "SenMayo_common_core_coverage.csv"
    )
    common_core = set(
        pd.read_csv(
            project / "results" / "common_core" / "SenMayo_common_core_gene_ledger.csv"
        )["gene"]
    )
    senmayo_source = availability_map["SenMayo_source"]
    senmayo_nonoverlap = availability_map["SenMayo_nonoverlap_source"]
    overlap = availability_map["SenMayo_FAP13_overlap"]
    sasp25 = availability_map["SASP25_source"]
    fap13 = {
        "FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9",
        "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A",
    }
    matrix4 = {"COL1A1", "COL1A2", "COL3A1", "FN1"}
    receptor2 = {"SDC4", "CD44"}
    represented: dict[str, set[str]] = {}
    for _, row in coverage.iterrows():
        represented[row["cohort"]] = senmayo_nonoverlap - split_genes(row["missing_from_source"])
    universe = sorted(senmayo_source | sasp25 | fap13 | matrix4 | receptor2)
    records = []
    for gene in universe:
        record = {
            "gene": gene,
            "SenMayo125_source": gene in senmayo_source,
            "excluded_from_SenMayo_for_FAP13_overlap": gene in overlap,
            "SenMayo121_nonoverlap_candidate": gene in senmayo_nonoverlap,
            "SASP25": gene in sasp25,
            "FAP13": gene in fap13,
            "matrix4": gene in matrix4,
            "receptor2": gene in receptor2,
            "common_core111": gene in common_core,
        }
        for cohort in coverage["cohort"]:
            record[f"SenMayo_represented_{cohort}"] = gene in represented[cohort]
        records.append(record)
    output = out_dir / "S2_Table_signature_membership_overlap_and_representation.csv"
    pd.DataFrame(records).to_csv(output, index=False)
    coverage.to_csv(out_dir / "S2_Data_SenMayo_cohort_coverage_summary.csv", index=False)
    return output


def build_s3(project: Path, out_dir: Path) -> Path:
    lock = pd.read_csv(project / "ledgers" / "scrna_endpoint_estimand_lock_v01_2026-09-04.csv")
    rows = []
    for _, row in lock.iterrows():
        rows.append(
            {
                "cohort": row["cohort"],
                "definition": "source-available",
                "endpoint": row["endpoint"],
                "score_definition": row["score_definition"],
                "n_patients": row["n_patients"],
                "primary_estimand": row["primary_summary"],
                "estimate": row["median_difference"],
                "ci_low": row["median_ci_low"],
                "ci_high": row["median_ci_high"],
                "fibroblast_higher": row["fibroblast_higher"],
                "epithelial_higher": row["epithelial_higher"],
                "p_value": row["exact_wilcoxon_p"],
                "bootstrap_resamples": row["bootstrap_resamples"],
            }
        )
    common_dir = project / "results" / "common_core"
    with (common_dir / "gse132_common_core_summary.json").open(encoding="utf-8") as handle:
        g132 = json.load(handle)
    g166 = pd.read_csv(common_dir / "gse166_common_core_summary.csv").iloc[0].to_dict()
    for summary in (g132, g166):
        rows.append(
            {
                "cohort": summary["cohort"],
                "definition": "111-gene common core",
                "endpoint": "SenMayo",
                "score_definition": summary["normalization"],
                "n_patients": summary["patients"],
                "primary_estimand": "mean fibroblast-minus-epithelial patient difference",
                "estimate": summary["paired_mean_difference"],
                "ci_low": summary["mean_difference_ci_low"],
                "ci_high": summary["mean_difference_ci_high"],
                "fibroblast_higher": summary["positive_differences"],
                "epithelial_higher": summary["negative_differences"],
                "p_value": summary["wilcoxon_signed_rank_p"],
                "bootstrap_resamples": 10000,
            }
        )
    output = out_dir / "S3_Table_single_cell_compartment_statistics.csv"
    pd.DataFrame(rows).to_csv(output, index=False)
    single_dir = project / "results" / "single_cell"
    pd.read_csv(single_dir / "scrna_patient_compartment_scores_aligned.csv").to_csv(
        out_dir / "S3_Data_patient_compartment_scores.csv", index=False
    )
    paired = pd.concat(
        [
            pd.read_csv(common_dir / "gse132_common_core_paired_values.csv"),
            pd.read_csv(common_dir / "gse166_common_core_paired_values.csv"),
        ],
        ignore_index=True,
    )
    paired.to_csv(out_dir / "S3_Data_common_core_patient_effects.csv", index=False)
    return output


def build_s4(project: Path, out_dir: Path) -> Path:
    single_dir = project / "results" / "single_cell"
    cluster = pd.read_csv(single_dir / "gse132465_cluster_inference.csv")
    cluster_rows = []
    endpoint_map = {
        "GSE132_FAPdet_SenMayo_depth_subtype": ("SenMayo", "binary FAP-detected"),
        "GSE132_FAPdet_SASP_depth_subtype": ("SASP25", "binary FAP-detected"),
        "GSE132_matrix4_on_FAP_depth_subtype": ("matrix4", "continuous FAP expression"),
    }
    for _, row in cluster.iterrows():
        endpoint, predictor = endpoint_map[row["model_id"]]
        cluster_rows.append(
            {
                "analysis_family": "signature or matrix score model",
                "model_id_or_gene": row["model_id"],
                "estimator": row["estimator"],
                "term": row["term"],
                "endpoint": endpoint,
                "predictor": predictor,
                "estimate": row["estimate"],
                "std_error": row["std_error"],
                "df": row["df"],
                "statistic": row["statistic"],
                "p_value": row["p_value"],
                "ci_low": row["ci_low"],
                "ci_high": row["ci_high"],
                "patient_clusters": row["clusters"],
                "bootstrap_requested": row["bootstrap_requested"],
                "bootstrap_valid": row["bootstrap_valid"],
                "diagnostic": "converged",
                "source_note": "depth, subtype, MKI67 and patient-adjusted model",
            }
        )
    marker_path = (
        project
        / "figures"
        / "supporting"
        / "source_data"
        / "S3_adjusted_marker_detection_ORs.csv"
    )
    marker = pd.read_csv(marker_path)
    for _, row in marker.iterrows():
        cluster_rows.append(
            {
                "analysis_family": "single-gene detection model",
                "model_id_or_gene": row["gene"],
                "estimator": "patient random-intercept logistic model",
                "term": "FAP-detected",
                "endpoint": f"{row['gene']} transcript detection",
                "predictor": "binary FAP-detected",
                "estimate": row["estimate"],
                "std_error": pd.NA,
                "df": pd.NA,
                "statistic": pd.NA,
                "p_value": row["p_value"],
                "ci_low": row["ci_low"],
                "ci_high": row["ci_high"],
                "patient_clusters": 23,
                "bootstrap_requested": pd.NA,
                "bootstrap_valid": pd.NA,
                "diagnostic": row["diagnostic"],
                "source_note": row["source"],
            }
        )
    output = out_dir / "S4_Table_FAP_specific_models_and_diagnostics.csv"
    pd.DataFrame(cluster_rows).to_csv(output, index=False)
    return output


def build_s5(project: Path, out_dir: Path) -> Path:
    source = (
        project
        / "results"
        / "single_cell"
        / "dependent_correlation_contrasts.csv"
    )
    data = pd.read_csv(source)
    output = out_dir / "S5_Table_dependent_patient_level_correlations.csv"
    data.to_csv(output, index=False)
    return output


def build_s6(project: Path, out_dir: Path) -> Path:
    bulk_dir = project / "results" / "bulk_composition"
    partial = pd.read_csv(bulk_dir / "target_purged_partial_correlations.csv")
    overlap = pd.read_csv(bulk_dir / "composition_proxy_feature_overlap_audit.csv")
    marginal = pd.read_csv(bulk_dir / "bulk_marginal_correlations_bootstrap.csv")
    rows = []
    for _, row in marginal.iterrows():
        rows.append(
            {
                "record_type": "marginal correlation",
                "cohort_or_proxy": row["cohort"],
                "target_pair": "FAP13_matrix4",
                "adjustment_or_target_set": "none",
                "proxy_definition": "not applicable",
                "n": row["n"],
                "estimate": row["rho"],
                "ci_low": row["ci_low"],
                "ci_high": row["ci_high"],
                "p_value": row["p_value"],
                "bootstrap_replicates": row["bootstrap_replicates"],
                "overlap_count": pd.NA,
                "overlap_genes": pd.NA,
                "retained_after_purge": pd.NA,
            }
        )
    for _, row in partial.iterrows():
        rows.append(
            {
                "record_type": "TCGA partial correlation",
                "cohort_or_proxy": "TCGA-COAD/READ",
                "target_pair": row["pair"],
                "adjustment_or_target_set": row["adjustment"],
                "proxy_definition": row["proxy_definition"],
                "n": row["n"],
                "estimate": row["rho"],
                "ci_low": row["ci_low"],
                "ci_high": row["ci_high"],
                "p_value": row["p_value"],
                "bootstrap_replicates": 5000 if pd.notna(row["ci_low"]) else pd.NA,
                "overlap_count": pd.NA,
                "overlap_genes": pd.NA,
                "retained_after_purge": pd.NA,
            }
        )
    for _, row in overlap.iterrows():
        rows.append(
            {
                "record_type": "proxy feature overlap",
                "cohort_or_proxy": row["proxy"],
                "target_pair": pd.NA,
                "adjustment_or_target_set": row["target_set"],
                "proxy_definition": pd.NA,
                "n": pd.NA,
                "estimate": pd.NA,
                "ci_low": pd.NA,
                "ci_high": pd.NA,
                "p_value": pd.NA,
                "bootstrap_replicates": pd.NA,
                "overlap_count": row["overlap_count"],
                "overlap_genes": row["overlap_genes"],
                "retained_after_purge": row["retained_after_purge"],
            }
        )
    output = out_dir / "S6_Table_bulk_composition_overlap_and_sensitivity.csv"
    pd.DataFrame(rows).to_csv(output, index=False)
    return output


def build_s7(project: Path, route_root: Path, out_dir: Path) -> Path:
    common_dir = project / "results" / "common_core"
    null_dir = project / "results" / "spatial_null"
    h1 = pd.read_csv(route_root / "cross_cohort" / "all_cohorts_H1_patient_effects.csv")
    h2 = pd.read_csv(route_root / "cross_cohort" / "all_cohorts_H2_patient_effects.csv")
    mki = pd.read_csv(route_root / "cross_cohort" / "all_cohorts_MKI67_patient_effects.csv")
    core = pd.read_csv(common_dir / "spatial_common_core_patient_effects.csv")
    null = pd.read_csv(null_dir / "spatial_matched_null_corrected.csv")
    definition = pd.read_csv(route_root / "cross_cohort" / "GSE280315_H1_definition_sensitivity.csv")
    pre = pd.read_csv(route_root / "GSE334323_pre_v01e_all_gene_expression" / "patient_H1_compartment.csv")
    corrected = pd.read_csv(route_root / "GSE334323" / "patient_H1_compartment.csv")

    rows: list[dict] = []
    for _, row in h1.iterrows():
        rows.append(
            {
                "analysis": "H1 source-available SenMayo",
                "cohort": row["cohort"],
                "patient": row["patient"],
                "definition": row.get("definition", pd.NA),
                "effect_name": "stromal-minus-tumour patient median",
                "effect_value": row["median_difference"],
                "n_stromal": row["n_stromal"],
                "n_tumour": row["n_tumour"],
                "n_sections": row.get("n_sections", pd.NA),
            }
        )
    for _, row in core.iterrows():
        rows.append(
            {
                "analysis": "H1 111-gene common core",
                "cohort": row["cohort"],
                "patient": row["patient"],
                "definition": "111-gene common core",
                "effect_name": "stromal-minus-tumour patient mean",
                "effect_value": row["stromal_minus_tumour_mean_score"],
                "n_stromal": pd.NA,
                "n_tumour": pd.NA,
                "n_sections": row["n_sections"],
            }
        )
    for _, row in h2.iterrows():
        rows.append(
            {
                "analysis": "H2 FAP-matrix comparative correlation",
                "cohort": row["cohort"],
                "patient": row["patient"],
                "definition": row.get("definition", pd.NA),
                "effect_name": "Fisher-z FAP-matrix contrast",
                "effect_value": row["fisher_z_axis_contrast"],
                "n_stromal": row["n_stromal"],
                "n_tumour": pd.NA,
                "n_sections": row.get("n_sections", pd.NA),
                "rho_fap_matrix4": row["rho_fap_matrix4"],
                "rho_fap_senmayo": row["rho_fap_senmayo"],
                "rho_senmayo_matrix4": row["rho_senmayo_matrix4"],
            }
        )
    for _, row in mki.iterrows():
        rows.append(
            {
                "analysis": "MKI67 polarity control",
                "cohort": row["cohort"],
                "patient": row["patient"],
                "definition": pd.NA,
                "effect_name": "epithelial-minus-stromal mean",
                "effect_value": row["epithelial_minus_stromal_mean"],
                "n_stromal": row["n_stromal"],
                "n_tumour": row["n_tumour"],
                "n_sections": row.get("n_sections", pd.NA),
            }
        )
    for _, row in null.iterrows():
        rows.append(
            {
                "analysis": "source-available spatial matched-null",
                "cohort": row["cohort"],
                "patient": "cohort-level statistic",
                "definition": "5,000 expression/detection-matched gene sets",
                "effect_name": "observed median patient mean difference",
                "effect_value": row["observed_median_patient_mean_difference"],
                "null_center": row["null_center_median"],
                "null_q025": row["null_q025"],
                "null_q975": row["null_q975"],
                "observed_minus_null_center": row["observed_minus_null_center"],
                "empirical_p": row["empirical_centered_absolute_p"],
                "monte_carlo_draws": row["draws"],
            }
        )
    for _, row in definition.iterrows():
        rows.append(
            {
                "analysis": "GSE280315 compartment-assignment sensitivity",
                "cohort": "GSE280315",
                "patient": row["patient"],
                "definition": str(row["definition"]).replace("\n", " ").replace("\r", " "),
                "effect_name": "stromal-minus-tumour patient median",
                "effect_value": row["median_difference"],
            }
        )
    for label, data in (("pre-correction", pre), ("human-feature corrected", corrected)):
        for _, row in data.iterrows():
            rows.append(
                {
                    "analysis": "GSE334323 denominator provenance",
                    "cohort": "GSE334323",
                    "patient": row["patient"],
                    "definition": label,
                    "effect_name": "stromal-minus-tumour patient median",
                    "effect_value": row["median_difference"],
                    "n_stromal": row["n_stromal"],
                    "n_tumour": row["n_tumour"],
                    "technical_note": "323 microbial, viral or control probes removed only in corrected denominator"
                    if label == "human-feature corrected"
                    else "original all-feature library-size denominator retained for provenance",
                }
            )

    primary_tests = []
    for cohort in ("GSE280315", "Valdeolivas_Visium"):
        table = pd.read_csv(route_root / cohort / "primary_tests.csv")
        table.insert(0, "cohort", cohort)
        primary_tests.append(table)
    tests = pd.concat(primary_tests, ignore_index=True)
    test_map = {
        (row["cohort"], "H1" if str(row["hypothesis"]).startswith("H1") else "H2"): row
        for _, row in tests.iterrows()
    }
    for row in rows:
        key = None
        if row["analysis"] == "H1 source-available SenMayo":
            key = (row["cohort"], "H1")
        elif row["analysis"] == "H2 FAP-matrix comparative correlation":
            key = (row["cohort"], "H2")
        if key in test_map:
            test = test_map[key]
            row["cohort_raw_p"] = test["raw_p"]
            row["cohort_holm_p"] = test["holm_p"]

    output = out_dir / "S7_Table_spatial_results_and_provenance.csv"
    pd.DataFrame(rows).to_csv(output, index=False)
    tests.to_csv(out_dir / "S7_Data_spatial_primary_tests.csv", index=False)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("qa64_root", type=Path)
    parser.add_argument("route_root", type=Path)
    args = parser.parse_args()

    project = args.project_root.resolve()
    qa_root = args.qa64_root.resolve()
    route_root = args.route_root.resolve()
    out_dir = project / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = [
        build_s1(project, qa_root, route_root, out_dir),
        build_s2(project, qa_root, out_dir),
        build_s3(project, out_dir),
        build_s4(project, out_dir),
        build_s5(project, out_dir),
        build_s6(project, out_dir),
        build_s7(project, route_root, out_dir),
    ]
    manifest_rows = []
    for path in sorted(out_dir.glob("*.csv")):
        frame = pd.read_csv(path)
        manifest_rows.append(
            {
                "file": path.name,
                "rows": len(frame),
                "columns": len(frame.columns),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "role": "primary supplementary table" if path in outputs else "companion source-data table",
            }
        )
    pd.DataFrame(manifest_rows).to_csv(out_dir / "Supplementary_table_manifest.csv", index=False)
    print(f"Wrote {len(manifest_rows)} supplementary table files to {out_dir}")


if __name__ == "__main__":
    main()
