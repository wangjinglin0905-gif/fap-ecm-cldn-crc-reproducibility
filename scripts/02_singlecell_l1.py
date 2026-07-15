from __future__ import annotations

import csv
import gzip
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


TASK_ROOT = Path.cwd() / "work" / "reproducibility"
OUTPUT_DIR = TASK_ROOT / "results" / "L1_single_cell"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MATRIX_PATH = Path(
    os.environ.get(
        "FAP_GSE132465_MATRIX",
        TASK_ROOT / "inputs" / "GSE132465_raw_UMI_matrix.txt.gz",
    )
)
ANNOTATION_PATH = Path(
    os.environ.get(
        "FAP_GSE132465_METADATA",
        TASK_ROOT / "inputs" / "GSE132465_cell_annotation.txt.gz",
    )
)

MARKERS = [
    "FAP",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "FN1",
    "POSTN",
    "ACTA2",
    "TAGLN",
    "CXCL12",
    "CCL2",
    "CD74",
    "HLA-DRA",
    "CLDN1",
    "CLDN2",
    "CLDN3",
    "CLDN4",
    "CLDN7",
    "CLDN18",
    "OCLN",
    "TJP1",
    "CDH1",
    "EPCAM",
    "VIM",
    "ZEB1",
    "SNAI1",
    "SNAI2",
    "TWIST1",
]


def load_annotations(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row["Index"] for row in rows], {row["Index"]: row for row in rows}


def parse_values(raw_line: bytes, expected: int) -> tuple[str, np.ndarray]:
    gene_raw, values_raw = raw_line.rstrip(b"\r\n").split(b"\t", 1)
    values = np.fromstring(values_raw, dtype=np.int32, sep="\t")
    if values.size != expected:
        raise ValueError(f"{gene_raw!r}: expected {expected} values, found {values.size}")
    return gene_raw.decode("utf-8"), values


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize_expression(
    indices: np.ndarray,
    label_fields: dict[str, str],
    normalized: dict[str, np.ndarray],
    raw_counts: dict[str, np.ndarray],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for gene in MARKERS:
        values = normalized[gene][indices]
        counts = raw_counts[gene][indices]
        output.append(
            {
                **label_fields,
                "gene": gene,
                "n_cells": int(indices.size),
                "mean_log_normalized": float(np.mean(values)),
                "median_log_normalized": float(np.median(values)),
                "positive_cells": int(np.count_nonzero(counts > 0)),
                "positive_percent": float(100 * np.mean(counts > 0)),
            }
        )
    return output


def main() -> None:
    annotation_ids, annotations = load_annotations(ANNOTATION_PATH)
    with gzip.open(MATRIX_PATH, "rb") as handle:
        header = handle.readline().decode("utf-8").rstrip("\r\n").split("\t")
    cell_ids = header[1:]
    if set(cell_ids) != set(annotation_ids):
        raise ValueError("Expression and annotation cell identifiers do not match")

    library_sizes = np.zeros(len(cell_ids), dtype=np.int64)
    raw_counts = {gene: np.zeros(len(cell_ids), dtype=np.int32) for gene in MARKERS}
    genes_seen: set[str] = set()

    with gzip.open(MATRIX_PATH, "rb") as handle:
        handle.readline()
        for line_number, raw_line in enumerate(handle, start=2):
            gene, values = parse_values(raw_line, len(cell_ids))
            library_sizes += values
            if gene in raw_counts:
                raw_counts[gene] = values.copy()
                genes_seen.add(gene)
            if line_number % 5000 == 0:
                print(f"Processed {line_number - 1} genes", flush=True)

    missing_markers = sorted(set(MARKERS) - genes_seen)
    normalized = {
        gene: np.log1p(
            np.divide(
                counts * 10000.0,
                library_sizes,
                out=np.zeros_like(counts, dtype=np.float64),
                where=library_sizes > 0,
            )
        )
        for gene, counts in raw_counts.items()
    }

    group_indices: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    patient_counts = Counter()
    for index, cell_id in enumerate(cell_ids):
        row = annotations[cell_id]
        group_indices[(row["Class"], row["Cell_type"], row["Cell_subtype"])].append(index)
        patient_counts[(row["Patient"], row["Class"], row["Cell_type"], row["Cell_subtype"])] += 1

    annotation_rows = [
        {
            "class": key[0],
            "cell_type": key[1],
            "cell_subtype": key[2],
            "n_cells": len(indices),
        }
        for key, indices in sorted(group_indices.items())
    ]
    write_rows(OUTPUT_DIR / "gse132465_annotation_counts.csv", annotation_rows)

    patient_rows = [
        {
            "patient": key[0],
            "class": key[1],
            "cell_type": key[2],
            "cell_subtype": key[3],
            "n_cells": count,
        }
        for key, count in sorted(patient_counts.items())
    ]
    write_rows(OUTPUT_DIR / "gse132465_patient_cell_counts.csv", patient_rows)

    expression_rows: list[dict[str, object]] = []
    for key, indices_list in sorted(group_indices.items()):
        indices = np.asarray(indices_list, dtype=np.int32)
        expression_rows.extend(
            summarize_expression(
                indices,
                {"class": key[0], "cell_type": key[1], "cell_subtype": key[2]},
                normalized,
                raw_counts,
            )
        )
    write_rows(OUTPUT_DIR / "gse132465_marker_expression_by_subtype.csv", expression_rows)

    epithelial_all = np.asarray(
        [i for i, cell_id in enumerate(cell_ids) if annotations[cell_id]["Cell_type"] == "Epithelial cells"],
        dtype=np.int32,
    )
    epithelial_tumor = np.asarray(
        [
            i
            for i, cell_id in enumerate(cell_ids)
            if annotations[cell_id]["Cell_type"] == "Epithelial cells" and annotations[cell_id]["Class"] == "Tumor"
        ],
        dtype=np.int32,
    )
    myofibroblast_all = np.asarray(
        [i for i, cell_id in enumerate(cell_ids) if annotations[cell_id]["Cell_subtype"] == "Myofibroblasts"],
        dtype=np.int32,
    )
    myofibroblast_tumor = np.asarray(
        [
            i
            for i, cell_id in enumerate(cell_ids)
            if annotations[cell_id]["Cell_subtype"] == "Myofibroblasts" and annotations[cell_id]["Class"] == "Tumor"
        ],
        dtype=np.int32,
    )

    focused_rows: list[dict[str, object]] = []
    for label, indices in [
        ("epithelial_all", epithelial_all),
        ("epithelial_tumor", epithelial_tumor),
        ("myofibroblast_all", myofibroblast_all),
        ("myofibroblast_tumor", myofibroblast_tumor),
    ]:
        focused_rows.extend(
            summarize_expression(indices, {"population": label}, normalized, raw_counts)
        )
    write_rows(OUTPUT_DIR / "gse132465_focused_marker_summaries.csv", focused_rows)

    fap_threshold = float(np.quantile(normalized["FAP"][myofibroblast_tumor], 0.75, method="linear"))
    strict_sender = myofibroblast_tumor[
        (normalized["FAP"][myofibroblast_tumor] >= fap_threshold)
        & (raw_counts["COL1A1"][myofibroblast_tumor] > 0)
    ]
    sender_patient_counts = Counter(annotations[cell_ids[index]]["Patient"] for index in strict_sender)
    write_rows(
        OUTPUT_DIR / "strict_fap_high_sender_by_patient.csv",
        [{"patient": patient, "n_cells": count} for patient, count in sorted(sender_patient_counts.items())],
    )

    tumor_patients = sorted({annotations[cell_ids[index]]["Patient"] for index in epithelial_tumor})
    patient_to_code = {patient: code for code, patient in enumerate(tumor_patients)}
    epithelial_patient_codes = np.asarray(
        [patient_to_code[annotations[cell_ids[index]]["Patient"]] for index in epithelial_tumor], dtype=np.int16
    )
    myofibroblast_patient_codes = np.asarray(
        [patient_to_code[annotations[cell_ids[index]]["Patient"]] for index in myofibroblast_tumor], dtype=np.int16
    )
    myofibroblast_patient_counts = Counter(
        annotations[cell_ids[index]]["Patient"] for index in myofibroblast_tumor
    )
    burden_rows = []
    for patient in tumor_patients:
        total_myofibroblasts = myofibroblast_patient_counts.get(patient, 0)
        sender_count = sender_patient_counts.get(patient, 0)
        burden_rows.append(
            {
                "patient": patient,
                "tumor_myofibroblasts": total_myofibroblasts,
                "strict_fap_high_senders": sender_count,
                "strict_sender_fraction": sender_count / total_myofibroblasts if total_myofibroblasts else 0.0,
                "tumor_epithelial_cells": int(np.sum(epithelial_patient_codes == patient_to_code[patient])),
            }
        )
    write_rows(OUTPUT_DIR / "patient_sender_burden.csv", burden_rows)

    pseudobulk_path = OUTPUT_DIR / "epithelial_pseudobulk_counts.tsv.gz"
    myofibroblast_pseudobulk_path = OUTPUT_DIR / "myofibroblast_pseudobulk_counts.tsv.gz"
    sender_path = OUTPUT_DIR / "strict_sender_gene_expression.tsv.gz"
    with gzip.open(pseudobulk_path, "wt", encoding="utf-8", newline="") as pseudobulk_handle, gzip.open(
        myofibroblast_pseudobulk_path, "wt", encoding="utf-8", newline=""
    ) as myofibroblast_handle, gzip.open(sender_path, "wt", encoding="utf-8", newline="") as sender_handle:
        pseudobulk_writer = csv.writer(pseudobulk_handle, delimiter="\t", lineterminator="\n")
        myofibroblast_writer = csv.writer(myofibroblast_handle, delimiter="\t", lineterminator="\n")
        sender_writer = csv.writer(sender_handle, delimiter="\t", lineterminator="\n")
        pseudobulk_writer.writerow(["gene", *tumor_patients])
        myofibroblast_writer.writerow(["gene", *tumor_patients])
        sender_writer.writerow(["gene", "total_count", "positive_cells", "mean_count_per_cell", "positive_percent"])
        with gzip.open(MATRIX_PATH, "rb") as matrix_handle:
            matrix_handle.readline()
            for line_number, raw_line in enumerate(matrix_handle, start=2):
                gene, values = parse_values(raw_line, len(cell_ids))
                epithelial_sums = np.bincount(
                    epithelial_patient_codes,
                    weights=values[epithelial_tumor],
                    minlength=len(tumor_patients),
                ).astype(np.int64)
                myofibroblast_sums = np.bincount(
                    myofibroblast_patient_codes,
                    weights=values[myofibroblast_tumor],
                    minlength=len(tumor_patients),
                ).astype(np.int64)
                sender_values = values[strict_sender]
                pseudobulk_writer.writerow([gene, *epithelial_sums.tolist()])
                myofibroblast_writer.writerow([gene, *myofibroblast_sums.tolist()])
                sender_writer.writerow(
                    [
                        gene,
                        int(sender_values.sum()),
                        int(np.count_nonzero(sender_values)),
                        float(sender_values.mean()),
                        float(100 * np.mean(sender_values > 0)),
                    ]
                )
                if line_number % 5000 == 0:
                    print(f"Pseudobulk pass: {line_number - 1} genes", flush=True)

    qc = {
        "dataset": "GSE132465",
        "total_cells": len(cell_ids),
        "epithelial_all": int(epithelial_all.size),
        "epithelial_tumor": int(epithelial_tumor.size),
        "myofibroblast_all": int(myofibroblast_all.size),
        "myofibroblast_tumor": int(myofibroblast_tumor.size),
        "strict_fap_high_sender": int(strict_sender.size),
        "strict_fap_high_threshold": fap_threshold,
        "pseudobulk_patients": len(tumor_patients),
        "missing_markers": missing_markers,
    }
    (OUTPUT_DIR / "gse132465_l1_qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(json.dumps(qc, indent=2), flush=True)


if __name__ == "__main__":
    main()
