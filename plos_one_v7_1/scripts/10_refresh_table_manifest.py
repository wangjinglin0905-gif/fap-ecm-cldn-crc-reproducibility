"""Refresh table hashes after formatting; excludes the manifest itself."""
import csv
import hashlib
import os
from pathlib import Path


def main():
    table_dir = Path(__file__).resolve().parents[1] / "tables"
    if os.name == "nt":
        table_dir = Path("\\\\?\\" + str(table_dir).removeprefix("\\\\?\\"))
    output = table_dir / "Supplementary_table_manifest.csv"
    rows = []
    for path in sorted(table_dir.glob("*.csv")):
        if path == output:
            continue
        with path.open(encoding="utf-8-sig", newline="") as handle:
            data = list(csv.reader(handle))
        content = path.read_bytes()
        rows.append({"file": path.name, "rows": len(data) - 1, "columns": len(data[0]), "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest(), "role": "primary supplementary table" if "_Table_" in path.name else "companion source-data table"})
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "rows", "columns", "bytes", "sha256", "role"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Refreshed {len(rows)} table records")


if __name__ == "__main__":
    main()
