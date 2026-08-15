#!/usr/bin/env python3
"""Build the deterministic SHA-256 manifest for the AJCR v6.4 bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "root",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="AJCR release directory (default: parent of scripts/)",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / "CHECKSUMS_SHA256.csv"
    files = sorted(
        (path for path in root.rglob("*") if path.is_file() and path != output),
        key=lambda path: path.relative_to(root).as_posix(),
    )

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["relative_path", "bytes", "sha256"])
        for path in files:
            writer.writerow([path.relative_to(root).as_posix(), path.stat().st_size, sha256(path)])

    print(f"Wrote {len(files)} records to {output}")


if __name__ == "__main__":
    main()
