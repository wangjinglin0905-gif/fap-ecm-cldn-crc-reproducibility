from __future__ import annotations

import gzip
import json
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


INPUT_DIR = Path("work/reproducibility/inputs/GSE166555")
OUTPUT_DIR = Path("work/reproducibility/results/L1_TISCH_GSE166555")
METADATA_PATH = INPUT_DIR / "GSE166555_meta_data.tsv.gz"
COUNTS_DIR = INPUT_DIR / "counts"

FAP_CAF_GENES = [
    "FAP",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "FN1",
    "POSTN",
    "THY1",
    "PDPN",
    "TAGLN",
    "ACTA2",
    "MMP2",
    "MMP9",
    "CXCL12",
    "TGFB1",
    "INHBA",
    "WNT2",
    "WNT5A",
]
DE_LIGAND_GENES = [
    gene for gene in FAP_CAF_GENES if gene not in {"TGFB1", "INHBA", "WNT2", "WNT5A"}
]
CLDN_CORE_GENES = ["CLDN1", "CLDN2", "CLDN4"]
LOCALIZATION_GENES = ["CLDN7", "EPCAM", "KRT19", "PTPRC", "PECAM1", "VWF"]
TARGET_GENES = sorted(set(FAP_CAF_GENES + CLDN_CORE_GENES + LOCALIZATION_GENES))


def zmean(frame: pd.DataFrame, genes: list[str]) -> pd.Series:
    available = [gene for gene in genes if gene in frame.columns]
    standardized = frame[available].apply(
        lambda column: (column - column.mean()) / column.std(ddof=1)
        if column.std(ddof=1) > 0
        else np.nan
    )
    return standardized.mean(axis=1, skipna=True)


def load_metadata() -> pd.DataFrame:
    columns = [
        "cell",
        "sample_id",
        "sample_origin",
        "main_cell_type",
        "sct_cell_type",
        "cell_type_str_custom",
        "cell_type_epi_custom",
        "nCount_RNA",
        "nFeature_RNA",
    ]
    metadata = pd.read_csv(
        METADATA_PATH,
        sep="\t",
        usecols=columns,
        low_memory=False,
    )
    metadata["patient"] = metadata["sample_id"].str.extract(r"(p\d{3})")
    metadata["caf_like"] = metadata["cell_type_str_custom"].fillna("").str.contains(
        r"FBs|CAFs|MyoFBs", regex=True
    )
    metadata["analysis_group"] = "Other"
    metadata.loc[metadata["main_cell_type"].eq("Immune"), "analysis_group"] = "Immune"
    metadata.loc[
        metadata["main_cell_type"].eq("Stromal"), "analysis_group"
    ] = "Other_stromal"
    metadata.loc[
        metadata["cell_type_str_custom"].isin(["Endothelial", "Pericytes"]),
        "analysis_group",
    ] = "Vascular"
    metadata.loc[metadata["caf_like"], "analysis_group"] = "CAF_like"
    metadata.loc[
        metadata["main_cell_type"].eq("Epithelial"), "analysis_group"
    ] = "Epithelial"
    metadata = metadata.set_index("cell", drop=False)
    return metadata


def read_selected_gene_counts(path: Path) -> tuple[list[str], dict[str, np.ndarray]]:
    selected = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        for line in handle:
            gene, separator, values = line.partition("\t")
            if separator and gene in TARGET_GENES:
                counts = np.fromstring(values, dtype=np.int64, sep="\t")
                if counts.size != len(cells):
                    raise RuntimeError(
                        f"Count length mismatch in {path.name} for {gene}: "
                        f"{counts.size} != {len(cells)}"
                    )
                selected[gene] = counts
    return cells, selected


def aggregate_counts(metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    count_records = []
    group_records = []
    for path in sorted(COUNTS_DIR.glob("*.tsv.gz")):
        sample_match = re.search(r"_(p\d{3}[nt]\d?)\.tsv\.gz$", path.name)
        if not sample_match:
            raise RuntimeError(f"Sample identifier could not be parsed from {path.name}")
        sample_id = sample_match.group(1)
        cells, selected_counts = read_selected_gene_counts(path)
        missing_cells = [cell for cell in cells if cell not in metadata.index]
        if missing_cells:
            raise RuntimeError(
                f"{path.name} contains {len(missing_cells)} cells absent from metadata"
            )
        sample_metadata = metadata.loc[cells].copy()
        if not sample_metadata["sample_id"].eq(sample_id).all():
            raise RuntimeError(f"Sample metadata mismatch for {path.name}")

        for analysis_group, group_metadata in sample_metadata.groupby("analysis_group"):
            positions = sample_metadata.index.get_indexer(group_metadata.index)
            total_umi = float(group_metadata["nCount_RNA"].sum())
            group_records.append(
                {
                    "sample_id": sample_id,
                    "patient": group_metadata["patient"].iloc[0],
                    "sample_origin": group_metadata["sample_origin"].iloc[0],
                    "analysis_group": analysis_group,
                    "cell_count": len(group_metadata),
                    "total_umi": total_umi,
                }
            )
            for gene in TARGET_GENES:
                counts = selected_counts.get(gene, np.zeros(len(cells), dtype=np.int64))[
                    positions
                ]
                count_records.append(
                    {
                        "sample_id": sample_id,
                        "patient": group_metadata["patient"].iloc[0],
                        "sample_origin": group_metadata["sample_origin"].iloc[0],
                        "analysis_group": analysis_group,
                        "gene": gene,
                        "count_sum": int(counts.sum()),
                        "positive_cells": int((counts > 0).sum()),
                        "cell_count": len(group_metadata),
                        "total_umi": total_umi,
                    }
                )
        print(f"processed {path.name}: {len(cells)} cells", flush=True)
    return pd.DataFrame(count_records), pd.DataFrame(group_records)


def patient_pseudobulk(counts: pd.DataFrame) -> pd.DataFrame:
    aggregated = (
        counts.groupby(
            ["patient", "sample_origin", "analysis_group", "gene"], as_index=False
        )
        .agg(
            count_sum=("count_sum", "sum"),
            positive_cells=("positive_cells", "sum"),
            cell_count=("cell_count", "sum"),
            total_umi=("total_umi", "sum"),
        )
    )
    aggregated["log_cpm"] = np.log1p(
        aggregated["count_sum"] / aggregated["total_umi"].replace(0, np.nan) * 1e6
    )
    aggregated["positive_fraction"] = (
        aggregated["positive_cells"] / aggregated["cell_count"].replace(0, np.nan)
    )
    return aggregated


def build_patient_scores(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    tumor = pseudobulk[pseudobulk["sample_origin"].eq("Tumor")]
    caf = tumor[tumor["analysis_group"].eq("CAF_like")]
    epithelial = tumor[tumor["analysis_group"].eq("Epithelial")]

    caf_expression = caf.pivot(index="patient", columns="gene", values="log_cpm")
    caf_positive = caf.pivot(
        index="patient", columns="gene", values="positive_fraction"
    )
    epithelial_expression = epithelial.pivot(
        index="patient", columns="gene", values="log_cpm"
    )
    caf_cells = caf.groupby("patient")["cell_count"].first()
    epithelial_cells = epithelial.groupby("patient")["cell_count"].first()

    common = sorted(set(caf_expression.index) & set(epithelial_expression.index))
    scores = pd.DataFrame(index=common)
    scores.index.name = "patient"
    scores["caf_cells"] = caf_cells.reindex(common)
    scores["epithelial_cells"] = epithelial_cells.reindex(common)
    scores["FAP_log_cpm"] = caf_expression.reindex(common)["FAP"]
    scores["FAP_positive_fraction"] = caf_positive.reindex(common)["FAP"]
    scores["FAP_CAF"] = zmean(caf_expression.reindex(common), FAP_CAF_GENES)
    scores["FAP_CAF_de_ligand"] = zmean(
        caf_expression.reindex(common), DE_LIGAND_GENES
    )
    for gene in CLDN_CORE_GENES:
        scores[f"epithelial_{gene}_log_cpm"] = epithelial_expression.reindex(common)[
            gene
        ]
    scores["epithelial_CLDN_core"] = zmean(
        epithelial_expression.reindex(common), CLDN_CORE_GENES
    )
    return scores.reset_index()


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    metadata = load_metadata()
    counts, group_summary = aggregate_counts(metadata)
    pseudobulk = patient_pseudobulk(counts)
    scores = build_patient_scores(pseudobulk)

    metadata_summary = (
        metadata.reset_index(drop=True)
        .groupby(
            ["patient", "sample_id", "sample_origin", "analysis_group"],
            as_index=False,
        )
        .agg(cells=("cell", "size"))
    )
    counts.to_csv(OUTPUT_DIR / "selected_gene_counts_by_sample_group.csv", index=False)
    group_summary.to_csv(OUTPUT_DIR / "sample_group_qc.csv", index=False)
    metadata_summary.to_csv(OUTPUT_DIR / "metadata_cell_counts.csv", index=False)
    pseudobulk.to_csv(OUTPUT_DIR / "patient_group_pseudobulk.csv", index=False)
    scores.to_csv(OUTPUT_DIR / "patient_cross_compartment_scores.csv", index=False)

    availability = (
        counts.groupby("gene", as_index=False)["count_sum"].sum().rename(
            columns={"count_sum": "total_counts"}
        )
    )
    availability["detected"] = availability["total_counts"].gt(0)
    availability.to_csv(OUTPUT_DIR / "target_gene_availability.csv", index=False)

    qc = {
        "dataset": "GSE166555",
        "metadata_cells": int(len(metadata)),
        "patients": int(metadata["patient"].nunique()),
        "samples": int(metadata["sample_id"].nunique()),
        "tumor_patients": int(
            metadata.loc[metadata["sample_origin"].eq("Tumor"), "patient"].nunique()
        ),
        "primary_eligible_patients": int(
            scores["caf_cells"].ge(20).mul(scores["epithelial_cells"].ge(100)).sum()
        ),
        "sensitivity_eligible_patients": int(
            scores["caf_cells"].ge(5).mul(scores["epithelial_cells"].ge(100)).sum()
        ),
        "caf_definition": "cell_type_str_custom matches FBs, CAFs, or MyoFBs",
        "primary_minimum_cells": {"CAF_like": 20, "Epithelial": 100},
        "normalization": "patient-group pseudobulk log1p(count_sum/total_UMI*1e6)",
    }
    (OUTPUT_DIR / "gse166555_qc.json").write_text(
        json.dumps(qc, indent=2), encoding="utf-8"
    )
    (OUTPUT_DIR / "sessionInfo.txt").write_text(
        "\n".join(
            [
                f"Python: {sys.version}",
                f"Platform: {platform.platform()}",
                f"numpy: {np.__version__}",
                f"pandas: {pd.__version__}",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
