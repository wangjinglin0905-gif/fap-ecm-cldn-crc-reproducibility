from __future__ import annotations

import csv
import gzip
import json
import shutil
import sys
from pathlib import Path

import numpy as np


def load_annotations(path: Path) -> dict[str, dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return {row["Index"]: row for row in csv.DictReader(handle, delimiter="\t")}


def parse_values(raw_line: bytes, expected: int) -> tuple[str, np.ndarray]:
    gene_raw, values_raw = raw_line.rstrip(b"\r\n").split(b"\t", 1)
    values = np.fromstring(values_raw, dtype=np.int32, sep="\t")
    if values.size != expected:
        raise ValueError(f"{gene_raw!r}: expected {expected} values, found {values.size}")
    return gene_raw.decode("utf-8"), values


def main() -> None:
    matrix_path = Path(sys.argv[1])
    annotation_path = Path(sys.argv[2])
    gene_list_path = Path(sys.argv[3])
    output_dir = Path(sys.argv[4])
    output_dir.mkdir(parents=True, exist_ok=True)

    annotations = load_annotations(annotation_path)
    wanted_genes = set(gene_list_path.read_text(encoding="utf-8").splitlines())

    with gzip.open(matrix_path, "rb") as handle:
        header = handle.readline().decode("utf-8").rstrip("\r\n").split("\t")
    cell_ids = header[1:]
    matrix_position = {cell_id: index for index, cell_id in enumerate(cell_ids)}
    if len(cell_ids) != len(annotations):
        raise ValueError(f"Matrix cells={len(cell_ids)}; annotation cells={len(annotations)}")

    candidate_indices = np.array(
        [
            index
            for index, cell_id in enumerate(cell_ids)
            if annotations[cell_id]["Class"] == "Tumor"
            and (
                annotations[cell_id]["Cell_subtype"] == "Myofibroblasts"
                or annotations[cell_id]["Cell_type"] == "Epithelial cells"
            )
        ],
        dtype=np.int32,
    )
    candidate_ids = [cell_ids[index] for index in candidate_indices]
    candidate_position = {cell_id: position for position, cell_id in enumerate(candidate_ids)}
    library_sizes = np.zeros(candidate_indices.size, dtype=np.int64)
    fap_counts = np.zeros(candidate_indices.size, dtype=np.int32)
    col1a1_counts = np.zeros(candidate_indices.size, dtype=np.int32)

    with gzip.open(matrix_path, "rb") as handle:
        handle.readline()
        for line_number, raw_line in enumerate(handle, start=2):
            gene, values = parse_values(raw_line, len(cell_ids))
            selected = values[candidate_indices]
            library_sizes += selected
            if gene == "FAP":
                fap_counts = selected.copy()
            elif gene == "COL1A1":
                col1a1_counts = selected.copy()
            if line_number % 5000 == 0:
                print(f"library-size pass: {line_number - 1} genes", flush=True)

    myofib_positions = np.array(
        [
            position
            for position, cell_id in enumerate(candidate_ids)
            if annotations[cell_id]["Cell_subtype"] == "Myofibroblasts"
        ],
        dtype=np.int32,
    )
    epithelial_positions = np.array(
        [
            position
            for position, cell_id in enumerate(candidate_ids)
            if annotations[cell_id]["Cell_type"] == "Epithelial cells"
        ],
        dtype=np.int32,
    )
    normalized_fap = np.log1p(
        np.divide(
            fap_counts * 10000.0,
            library_sizes,
            out=np.zeros_like(fap_counts, dtype=np.float64),
            where=library_sizes > 0,
        )
    )
    threshold = float(np.quantile(normalized_fap[myofib_positions], 0.75, method="linear"))
    sender_positions = myofib_positions[
        (normalized_fap[myofib_positions] >= threshold)
        & (col1a1_counts[myofib_positions] > 0)
    ]

    selected_positions = np.concatenate([sender_positions, epithelial_positions])
    selected_candidate_ids = [candidate_ids[position] for position in selected_positions]
    selected_matrix_indices = np.array(
        [matrix_position[cell_id] for cell_id in selected_candidate_ids], dtype=np.int32
    )

    entries_path = output_dir / "matrix_entries.tmp"
    genes_found: list[str] = []
    nonzero_count = 0
    with entries_path.open("w", encoding="ascii", newline="\n") as entries:
        with gzip.open(matrix_path, "rb") as handle:
            handle.readline()
            for line_number, raw_line in enumerate(handle, start=2):
                gene, values = parse_values(raw_line, len(cell_ids))
                if gene not in wanted_genes:
                    continue
                genes_found.append(gene)
                selected = values[selected_matrix_indices]
                nonzero_columns = np.flatnonzero(selected)
                row_number = len(genes_found)
                for column in nonzero_columns:
                    entries.write(f"{row_number} {column + 1} {int(selected[column])}\n")
                nonzero_count += nonzero_columns.size
                if len(genes_found) % 250 == 0:
                    print(f"matrix pass: {len(genes_found)} CellChat genes", flush=True)

    matrix_output = output_dir / "counts.mtx"
    with matrix_output.open("w", encoding="ascii", newline="\n") as output:
        output.write("%%MatrixMarket matrix coordinate integer general\n")
        output.write("% strict FAP-high myofibroblast and tumor epithelial cells\n")
        output.write(f"{len(genes_found)} {len(selected_candidate_ids)} {nonzero_count}\n")
        with entries_path.open("r", encoding="ascii") as entries:
            shutil.copyfileobj(entries, output)
    entries_path.unlink()

    (output_dir / "genes.tsv").write_text("\n".join(genes_found) + "\n", encoding="utf-8")
    with (output_dir / "metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["cell_id", "group", "patient", "class", "cell_type", "cell_subtype", "library_size", "FAP_normalized", "COL1A1_count"]
        )
        sender_ids = {candidate_ids[position] for position in sender_positions}
        for cell_id, candidate_pos in zip(selected_candidate_ids, selected_positions):
            row = annotations[cell_id]
            writer.writerow(
                [
                    cell_id,
                    "FAP_high_myofibroblast" if cell_id in sender_ids else "Tumor_epithelial",
                    row["Patient"],
                    row["Class"],
                    row["Cell_type"],
                    row["Cell_subtype"],
                    int(library_sizes[candidate_pos]),
                    f"{normalized_fap[candidate_pos]:.10g}",
                    int(col1a1_counts[candidate_pos]),
                ]
            )

    qc = {
        "matrix_cells": len(cell_ids),
        "annotation_cells": len(annotations),
        "tumor_myofibroblasts": int(myofib_positions.size),
        "tumor_epithelial": int(epithelial_positions.size),
        "fap_normalized_75th_percentile": threshold,
        "strict_sender_cells": int(sender_positions.size),
        "selected_cells": len(selected_candidate_ids),
        "cellchat_genes_requested": len(wanted_genes),
        "cellchat_genes_found": len(genes_found),
        "matrix_nonzero_entries": int(nonzero_count),
    }
    (output_dir / "input_qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(json.dumps(qc, indent=2), flush=True)


if __name__ == "__main__":
    main()
