from __future__ import annotations

import gzip
import os
import shutil
from pathlib import Path

import numpy as np


TASK_ROOT = Path.cwd() / "work" / "reproducibility"
OUTPUT_DIR = TASK_ROOT / "inputs" / "harmony"
MATRIX_PATH = Path(
    os.environ.get(
        "FAP_GSE132465_MATRIX",
        TASK_ROOT / "inputs" / "GSE132465_raw_UMI_matrix.txt.gz",
    )
)
GENE_LIST = OUTPUT_DIR / "qi_top3000_variable_genes.txt"


def parse_values(raw_line: bytes, expected: int) -> tuple[str, np.ndarray]:
    gene_raw, values_raw = raw_line.rstrip(b"\r\n").split(b"\t", 1)
    values = np.fromstring(values_raw, dtype=np.int32, sep="\t")
    if values.size != expected:
        raise ValueError(f"{gene_raw!r}: expected {expected} values, found {values.size}")
    return gene_raw.decode("utf-8"), values


def main() -> None:
    wanted = set(GENE_LIST.read_text(encoding="utf-8").splitlines())
    with gzip.open(MATRIX_PATH, "rb") as handle:
        header = handle.readline().decode("utf-8").rstrip("\r\n").split("\t")
    cells = header[1:]
    entries_path = OUTPUT_DIR / "gse_harmony_entries.tmp"
    genes: list[str] = []
    nonzero_count = 0
    with entries_path.open("w", encoding="ascii", newline="\n") as entries, gzip.open(MATRIX_PATH, "rb") as handle:
        handle.readline()
        for raw_line in handle:
            gene, values = parse_values(raw_line, len(cells))
            if gene not in wanted:
                continue
            genes.append(gene)
            nonzero_columns = np.flatnonzero(values)
            row_number = len(genes)
            for column in nonzero_columns:
                entries.write(f"{row_number} {column + 1} {int(values[column])}\n")
            nonzero_count += int(nonzero_columns.size)
            if len(genes) % 500 == 0:
                print(f"Extracted {len(genes)} genes", flush=True)

    matrix_path = OUTPUT_DIR / "gse_harmony_counts.mtx"
    with matrix_path.open("w", encoding="ascii", newline="\n") as output, entries_path.open(
        "r", encoding="ascii"
    ) as entries:
        output.write("%%MatrixMarket matrix coordinate integer general\n")
        output.write("% GSE132465 genes selected from Qi variable-gene ranking\n")
        output.write(f"{len(genes)} {len(cells)} {nonzero_count}\n")
        shutil.copyfileobj(entries, output)
    entries_path.unlink()
    (OUTPUT_DIR / "gse_harmony_genes.tsv").write_text("\n".join(genes) + "\n", encoding="utf-8")
    (OUTPUT_DIR / "gse_harmony_cells.tsv").write_text("\n".join(cells) + "\n", encoding="utf-8")
    print(f"Matrix: {len(genes)} genes x {len(cells)} cells; {nonzero_count} nonzero entries", flush=True)


if __name__ == "__main__":
    main()
