#!/usr/bin/env python3
"""Rescore GSE132465 from the public raw UMI matrix using an archived cell ledger.

Only requested gene rows are parsed numerically. The full gzip stream is read so
that the source remains the official GEO matrix, while exported nCount_RNA values
provide the exact per-cell library sizes used in the archived Seurat object.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def read_gene_set(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def read_archived_represented(path: Path) -> list[str]:
    if path.suffix.lower() != ".csv":
        return read_gene_set(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    hits = [row for row in rows if row.get("set") == "SenMayo_nonoverlap_represented"]
    if len(hits) != 1:
        raise RuntimeError("Archived represented-gene row was not found uniquely")
    return [gene for gene in hits[0]["genes"].split(";") if gene]


def read_archived_cells(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = reader.fieldnames or []
    required = {
        "cell", "Patient", "Class", "Sample", "Cell_type", "Cell_subtype",
        "nCount_RNA", "compartment", "SenMayo_zmean",
    }
    if not required.issubset(fields):
        raise RuntimeError(f"Archived cell ledger is missing columns: {sorted(required - set(fields))}")
    if len({row["cell"] for row in rows}) != len(rows):
        raise RuntimeError("Archived cell identifiers are not unique")
    return rows, fields


def verify_public_annotation(path: Path, archived: list[dict[str, str]]) -> dict[str, object]:
    wanted = {row["cell"]: row for row in archived}
    observed = 0
    mismatches: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cell = row["Index"]
            if cell not in wanted:
                continue
            observed += 1
            old = wanted[cell]
            for field in ("Patient", "Class", "Sample", "Cell_type", "Cell_subtype"):
                if row[field] != old[field]:
                    mismatches.append(f"{cell}:{field}:{row[field]}!={old[field]}")
                    break
    return {
        "archived_cells": len(archived),
        "public_annotation_matches": observed,
        "mismatch_count": len(mismatches),
        "first_mismatches": mismatches[:20],
    }


def extract_rows(raw_path: Path, target_cells: list[str], genes: list[str]) -> tuple[np.ndarray, list[str], int]:
    wanted = set(genes)
    values: dict[str, np.ndarray] = {}
    total_gene_rows = 0
    with gzip.open(raw_path, "rb") as handle:
        header = handle.readline().rstrip(b"\r\n").decode("utf-8").split("\t")
        matrix_cells = header[1:]
        positions = {cell: i for i, cell in enumerate(matrix_cells)}
        missing_cells = [cell for cell in target_cells if cell not in positions]
        if missing_cells:
            raise RuntimeError(f"{len(missing_cells)} archived cells absent from GEO matrix")
        indices = np.asarray([positions[cell] for cell in target_cells], dtype=np.int64)
        expected = len(matrix_cells)
        for raw_line in handle:
            total_gene_rows += 1
            tab = raw_line.find(b"\t")
            if tab < 0:
                continue
            gene = raw_line[:tab].decode("utf-8")
            if gene not in wanted:
                continue
            numeric = np.fromstring(raw_line[tab + 1 :].decode("ascii"), sep="\t", dtype=np.float64)
            if numeric.size != expected:
                raise RuntimeError(f"Unexpected value count for {gene}: {numeric.size} != {expected}")
            values[gene] = numeric[indices]
    represented = [gene for gene in genes if gene in values]
    if len(represented) != len(genes):
        raise RuntimeError(f"Requested genes absent from raw matrix: {sorted(set(genes) - set(represented))}")
    matrix = np.vstack([values[gene] for gene in represented])
    return matrix, represented, total_gene_rows


def zmean(log_expression: np.ndarray, row_selector: np.ndarray) -> tuple[np.ndarray, int]:
    block = log_expression[row_selector, :]
    means = block.mean(axis=1)
    sds = block.std(axis=1, ddof=1)
    keep = np.isfinite(sds) & (sds > 0)
    scores = ((block[keep, :] - means[keep, None]) / sds[keep, None]).mean(axis=0)
    return scores, int(keep.sum())


def exact_signed_rank(values: np.ndarray) -> tuple[float, float]:
    """Two-sided exact Wilcoxon signed-rank test for distinct nonzero magnitudes."""
    x = values[np.isfinite(values) & (values != 0)]
    magnitudes = np.abs(x)
    if len(np.unique(magnitudes)) != len(magnitudes):
        return math.nan, math.nan
    order = np.argsort(magnitudes)
    ranks = np.empty(len(x), dtype=int)
    ranks[order] = np.arange(1, len(x) + 1)
    observed = int(ranks[x > 0].sum())
    total_rank = int(ranks.sum())
    counts = np.zeros(total_rank + 1, dtype=np.int64)
    counts[0] = 1
    upper = 0
    for rank in range(1, len(x) + 1):
        for partial_sum in range(upper, -1, -1):
            counts[partial_sum + rank] += counts[partial_sum]
        upper += rank
    denominator = float(2 ** len(x))
    lower_tail = float(counts[: observed + 1].sum()) / denominator
    upper_tail = float(counts[observed:].sum()) / denominator
    return float(min(observed, total_rank - observed)), min(1.0, 2.0 * min(lower_tail, upper_tail))


def main() -> int:
    if len(sys.argv) != 8:
        raise SystemExit(
            "Usage: 05a_rescore_gse132_public_common_core.py <raw UMI.gz> <public annotation.gz> "
            "<archived cell_scores.csv> <archived gene_availability.csv> <common-core genes.txt> "
            "<output directory> <seed>"
        )
    raw_path, annotation_path, archived_path, represented_path, core_path = map(Path, sys.argv[1:6])
    out_dir = Path(sys.argv[6])
    seed = int(sys.argv[7])
    out_dir.mkdir(parents=True, exist_ok=True)

    archived, _ = read_archived_cells(archived_path)
    cell_ids = [row["cell"] for row in archived]
    library_size = np.asarray([float(row["nCount_RNA"]) for row in archived])
    old_score = np.asarray([float(row["SenMayo_zmean"]) for row in archived])
    represented = read_archived_represented(represented_path)
    core = read_gene_set(core_path)
    union = represented + [gene for gene in core if gene not in represented]

    annotation_check = verify_public_annotation(annotation_path, archived)
    if annotation_check["public_annotation_matches"] != len(archived) or annotation_check["mismatch_count"]:
        raise RuntimeError(f"Public annotation validation failed: {annotation_check}")

    print(f"Streaming {raw_path.name} for {len(union)} requested genes", flush=True)
    counts, extracted_genes, total_gene_rows = extract_rows(raw_path, cell_ids, union)
    if np.any(library_size <= 0):
        raise RuntimeError("Non-positive archived library size")

    rep_mask = np.asarray([gene in set(represented) for gene in extracted_genes])
    core_mask = np.asarray([gene in set(core) for gene in extracted_genes])
    candidates: dict[str, np.ndarray] = {
        "log1p_CPM_1e6": np.log1p(counts / library_size[None, :] * 1e6),
        "log1p_CPTT_1e4": np.log1p(counts / library_size[None, :] * 1e4),
        "log1p_raw_UMI": np.log1p(counts),
    }
    calibration: list[dict[str, object]] = []
    candidate_scores: dict[str, np.ndarray] = {}
    for label, expression in candidates.items():
        score, retained = zmean(expression, rep_mask)
        candidate_scores[label] = score
        calibration.append(
            {
                "candidate": label,
                "n_cells": len(cell_ids),
                "represented_genes": int(rep_mask.sum()),
                "retained_nonzero_variance": retained,
                "pearson_vs_archived_SenMayo_zmean": float(np.corrcoef(score, old_score)[0, 1]),
                "rmse_vs_archived_SenMayo_zmean": float(np.sqrt(np.mean((score - old_score) ** 2))),
                "max_abs_difference": float(np.max(np.abs(score - old_score))),
            }
        )
    calibration.sort(key=lambda row: (-float(row["pearson_vs_archived_SenMayo_zmean"]), float(row["rmse_vs_archived_SenMayo_zmean"])))
    chosen = str(calibration[0]["candidate"])
    if float(calibration[0]["pearson_vs_archived_SenMayo_zmean"]) < 0.99:
        raise RuntimeError(f"No public-raw normalization reproduced archived score adequately: {calibration}")

    common_score, retained_core = zmean(candidates[chosen], core_mask)
    groups: dict[tuple[str, str], list[float]] = {}
    for row, value in zip(archived, common_score, strict=True):
        groups.setdefault((row["Patient"], row["compartment"]), []).append(float(value))
    patient_rows: list[dict[str, object]] = []
    for (patient, compartment), values in sorted(groups.items()):
        patient_rows.append(
            {
                "cohort": "GSE132465",
                "patient": patient,
                "compartment": compartment,
                "score": float(np.mean(values)),
                "cell_count": len(values),
                "common_core_n": len(core),
                "retained_nonzero_variance_n": retained_core,
            }
        )
    patient_lookup = {(str(row["patient"]), str(row["compartment"])): float(row["score"]) for row in patient_rows}
    patients = sorted({str(row["patient"]) for row in patient_rows})
    differences = np.asarray([
        patient_lookup[(patient, "Fibroblast")] - patient_lookup[(patient, "Epithelial")]
        for patient in patients
    ])
    rng = np.random.default_rng(seed)
    boot = differences[rng.integers(0, len(differences), size=(10000, len(differences)))].mean(axis=1)
    wilcoxon_statistic, wilcoxon_p = exact_signed_rank(differences)
    summary = {
        "cohort": "GSE132465",
        "normalization": chosen,
        "common_core_n": len(core),
        "retained_nonzero_variance_n": retained_core,
        "patients": len(patients),
        "paired_mean_difference": float(differences.mean()),
        "paired_median_difference": float(np.median(differences)),
        "mean_difference_ci_low": float(np.quantile(boot, 0.025, method="weibull")),
        "mean_difference_ci_high": float(np.quantile(boot, 0.975, method="weibull")),
        "positive_differences": int((differences > 0).sum()),
        "negative_differences": int((differences < 0).sum()),
        "zero_differences": int((differences == 0).sum()),
        "wilcoxon_signed_rank_statistic": wilcoxon_statistic,
        "wilcoxon_signed_rank_p": wilcoxon_p,
    }

    with (out_dir / "gse132_public_normalization_calibration.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(calibration[0]))
        writer.writeheader()
        writer.writerows(calibration)
    with (out_dir / "gse132_common_core_patient_compartment.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(patient_rows[0]))
        writer.writeheader()
        writer.writerows(patient_rows)
    paired_rows = [
        {
            "cohort": "GSE132465",
            "patient": patient,
            "fibroblast_score": patient_lookup[(patient, "Fibroblast")],
            "epithelial_score": patient_lookup[(patient, "Epithelial")],
            "difference": float(difference),
        }
        for patient, difference in zip(patients, differences, strict=True)
    ]
    with (out_dir / "gse132_common_core_paired_values.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(paired_rows[0]))
        writer.writeheader()
        writer.writerows(paired_rows)
    (out_dir / "gse132_common_core_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    provenance = {
        "raw_umi_sha256": sha256(raw_path),
        "public_annotation_sha256": sha256(annotation_path),
        "archived_cell_ledger_sha256": sha256(archived_path),
        "annotation_check": annotation_check,
        "matrix_gene_rows_streamed": total_gene_rows,
        "matrix_cells": 63689,
        "target_cells": len(cell_ids),
        "genes_extracted": extracted_genes,
        "chosen_normalization": chosen,
        "seed": seed,
    }
    (out_dir / "gse132_public_rescore_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"calibration": calibration, "summary": summary}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
