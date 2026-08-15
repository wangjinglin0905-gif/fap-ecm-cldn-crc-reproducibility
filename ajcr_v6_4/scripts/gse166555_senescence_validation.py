from __future__ import annotations

import gzip
import json
import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


INPUT_DIR = Path()
COUNTS_DIR = Path()
METADATA_PATH = Path()
OUTPUT_DIR = Path()
GENE_AVAILABILITY_PATH = Path()

FAP13 = [
    "FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9",
    "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A",
]
MATRIX4 = ["COL1A1", "COL1A2", "COL3A1", "FN1"]
SASP25 = [
    "IL6", "CXCL8", "IL1A", "IL1B", "CCL2", "CCL5", "CXCL1", "CXCL2",
    "CXCL3", "CXCL10", "MMP1", "MMP3", "MMP9", "MMP10", "MMP13",
    "SERPINE1", "PLAU", "TIMP2", "VEGFA", "GDF15", "IGFBP3", "TNF",
    "CSF2", "HGF", "FAS",
]
MARKERS = ["CDKN2A", "CDKN2B", "CDKN1A", "LMNB1", "MKI67"]


def load_frozen_senmayo(path: Path) -> list[str]:
    availability = pd.read_csv(path)
    value = availability.loc[
        availability["set"].eq("SenMayo_nonoverlap_represented"), "genes"
    ].iloc[0]
    genes = value.split(";")
    if len(genes) != 119:
        raise RuntimeError(f"Expected frozen 119-gene SenMayo set, observed {len(genes)}")
    return genes


SENMAYO119: list[str] = []
TARGET_GENES: list[str] = []


def load_metadata() -> pd.DataFrame:
    columns = [
        "cell", "sample_id", "sample_origin", "main_cell_type",
        "cell_type_str_custom", "nCount_RNA", "nFeature_RNA",
    ]
    frame = pd.read_csv(
        METADATA_PATH, sep="\t", usecols=columns, low_memory=False
    )
    frame["patient"] = frame["sample_id"].str.extract(r"(p\d{3})")
    frame["caf_like"] = frame["cell_type_str_custom"].fillna("").str.contains(
        r"FBs|CAFs|MyoFBs", regex=True
    )
    frame["compartment"] = "Other"
    frame.loc[frame["caf_like"], "compartment"] = "Fibroblast"
    frame.loc[frame["main_cell_type"].eq("Epithelial"), "compartment"] = "Epithelial"
    return frame.set_index("cell", drop=False)


def read_target_counts(path: Path) -> tuple[list[str], dict[str, np.ndarray]]:
    selected: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        for line in handle:
            gene, separator, values = line.partition("\t")
            if separator and gene in TARGET_GENES:
                vector = np.fromstring(values, dtype=np.int64, sep="\t")
                if vector.size != len(cells):
                    raise RuntimeError(
                        f"{path.name} {gene}: {vector.size} values for {len(cells)} cells"
                    )
                selected[gene] = vector
    return cells, selected


def aggregate(metadata: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    paths = sorted(COUNTS_DIR.glob("*.tsv.gz"))
    if not paths:
        raise RuntimeError(f"No count matrices found under {COUNTS_DIR}")
    for path in paths:
        sample_match = re.search(r"_(p\d{3}[nt]\d?)\.tsv\.gz$", path.name)
        if not sample_match:
            raise RuntimeError(f"Could not parse sample identifier from {path.name}")
        sample_id = sample_match.group(1)
        cells, selected = read_target_counts(path)
        missing = [cell for cell in cells if cell not in metadata.index]
        if missing:
            raise RuntimeError(f"{path.name}: {len(missing)} cells absent from metadata")
        sample_meta = metadata.loc[cells]
        if not sample_meta["sample_id"].eq(sample_id).all():
            raise RuntimeError(f"Sample metadata mismatch for {path.name}")
        tumor_meta = sample_meta.loc[sample_meta["sample_origin"].eq("Tumor")]
        for compartment in ("Fibroblast", "Epithelial"):
            group = tumor_meta.loc[tumor_meta["compartment"].eq(compartment)]
            if group.empty:
                continue
            positions = sample_meta.index.get_indexer(group.index)
            total_umi = float(group["nCount_RNA"].sum())
            for gene in TARGET_GENES:
                vector = selected.get(gene)
                if vector is None:
                    count_sum = 0
                    positive_cells = 0
                else:
                    values = vector[positions]
                    count_sum = int(values.sum())
                    positive_cells = int((values > 0).sum())
                records.append(
                    {
                        "patient": group["patient"].iloc[0],
                        "sample_id": sample_id,
                        "compartment": compartment,
                        "gene": gene,
                        "count_sum": count_sum,
                        "positive_cells": positive_cells,
                        "cell_count": len(group),
                        "total_umi": total_umi,
                    }
                )
        print(f"processed {path.name}: {len(cells)} cells", flush=True)
    result = pd.DataFrame(records)
    return (
        result.groupby(["patient", "compartment", "gene"], as_index=False)
        .agg(
            count_sum=("count_sum", "sum"),
            positive_cells=("positive_cells", "sum"),
            cell_count=("cell_count", "sum"),
            total_umi=("total_umi", "sum"),
        )
    )


def zmean(frame: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    available = [gene for gene in genes if gene in frame.columns and frame[gene].std(ddof=1) > 0]
    standardized = frame[available].apply(
        lambda column: (column - column.mean()) / column.std(ddof=1)
    )
    return standardized.mean(axis=1), available


def build_scores(pseudobulk: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    cell_counts = (
        pseudobulk.groupby(["patient", "compartment"], as_index=False)
        .agg(cell_count=("cell_count", "first"))
        .pivot(index="patient", columns="compartment", values="cell_count")
    )
    eligible = cell_counts.index[
        cell_counts["Fibroblast"].ge(20) & cell_counts["Epithelial"].ge(100)
    ]
    data = pseudobulk.loc[pseudobulk["patient"].isin(eligible)].copy()
    data["log_cpm"] = np.log1p(data["count_sum"] / data["total_umi"] * 1e6)
    wide = data.pivot(index=["patient", "compartment"], columns="gene", values="log_cpm")
    result = pd.DataFrame(index=wide.index)
    result["SenMayo119"], sen_available = zmean(wide, SENMAYO119)
    result["SASP25"], sasp_available = zmean(wide, SASP25)
    result["FAP13"], fap_available = zmean(wide, FAP13)
    result["matrix4"], matrix_available = zmean(wide, MATRIX4)
    for marker in MARKERS:
        result[marker] = wide[marker] if marker in wide.columns else np.nan
    result = result.reset_index()
    result = result.merge(
        cell_counts.stack().rename("cell_count").reset_index(),
        on=["patient", "compartment"], how="left"
    )
    ledger = {
        "metadata_cells": int(pd.read_csv(METADATA_PATH, sep="\t", usecols=["cell"]).shape[0]),
        "tumor_patients": int(pseudobulk["patient"].nunique()),
        "eligible_patients": int(len(eligible)),
        "thresholds": {"Fibroblast": 20, "Epithelial": 100},
        "SenMayo119_available": sen_available,
        "SASP25_available": sasp_available,
        "FAP13_available": fap_available,
        "matrix4_available": matrix_available,
        "normalization": "patient-compartment pseudobulk log1p(CPM), followed by gene-wise z scoring across eligible patient-compartments",
        "caf_definition": "cell_type_str_custom regex FBs|CAFs|MyoFBs; tumor samples only",
    }
    return result, ledger


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recompute the frozen GSE166555 tumour-compartment evaluation."
    )
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--gene-availability", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    global INPUT_DIR, COUNTS_DIR, METADATA_PATH, OUTPUT_DIR
    global GENE_AVAILABILITY_PATH, SENMAYO119, TARGET_GENES
    INPUT_DIR = args.input_dir.resolve(strict=True)
    COUNTS_DIR = INPUT_DIR / "counts"
    METADATA_PATH = INPUT_DIR / "GSE166555_meta_data.tsv.gz"
    GENE_AVAILABILITY_PATH = args.gene_availability.resolve(strict=True)
    OUTPUT_DIR = args.output_dir
    if not COUNTS_DIR.is_dir() or not METADATA_PATH.is_file():
        raise FileNotFoundError(
            "Input directory must contain counts/ and GSE166555_meta_data.tsv.gz"
        )
    SENMAYO119 = load_frozen_senmayo(GENE_AVAILABILITY_PATH)
    TARGET_GENES = sorted(set(SENMAYO119 + SASP25 + MARKERS + FAP13 + MATRIX4))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    metadata = load_metadata()
    pseudobulk = aggregate(metadata)
    scores, ledger = build_scores(pseudobulk)
    pseudobulk.to_csv(OUTPUT_DIR / "target_gene_pseudobulk.csv", index=False)
    scores.to_csv(OUTPUT_DIR / "patient_compartment_scores.csv", index=False)
    (OUTPUT_DIR / "analysis_ledger.json").write_text(
        json.dumps(ledger, indent=2), encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in ledger.items() if not k.endswith("_available")}, indent=2))


if __name__ == "__main__":
    main()
