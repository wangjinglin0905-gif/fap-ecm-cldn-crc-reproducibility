#!/usr/bin/env python3
"""Validate figures, metadata, paths and checksums in the AJCR v6.4 bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from PIL import Image


TEXT_SUFFIXES = {".cff", ".csv", ".json", ".md", ".py", ".r", ".txt", ".yaml", ".yml"}
PRIVATE_PATH_PATTERNS = (
    re.compile(r"(?i)C:[/\\]+Users[/\\]+"),
    re.compile(r"(?i)E:[/\\]"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_figures(root: Path) -> dict:
    figure_dir = root / "figures"
    records = []
    for path in sorted(list(figure_dir.glob("*.png")) + list(figure_dir.glob("*.tiff"))):
        warnings = []
        with Image.open(path) as image:
            dpi = tuple(float(value) for value in image.info.get("dpi", (0.0, 0.0)))
            if image.mode != "RGB":
                warnings.append(f"mode is {image.mode}, expected RGB")
            if min(dpi) < 599.5:
                warnings.append(f"DPI is {dpi}, expected 600")
            records.append(
                {
                    "path": path.relative_to(root).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "format": image.format,
                    "mode": image.mode,
                    "pixels": list(image.size),
                    "dpi": [round(value, 3) for value in dpi],
                    "frames": getattr(image, "n_frames", 1),
                    "warnings": warnings,
                }
            )
    if len(records) != 10:
        raise RuntimeError(f"Expected 10 PNG/TIFF figure files, found {len(records)}")
    if any(record["warnings"] for record in records):
        raise RuntimeError("One or more figure files failed RGB/600-dpi validation")
    report = {"minimum_dpi": 600.0, "files": records}
    output = root / "qa" / "figure_validation_v6.4.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def validate_json(root: Path) -> int:
    count = 0
    for path in root.rglob("*.json"):
        json.loads(path.read_text(encoding="utf-8-sig"))
        count += 1
    return count


def validate_private_paths(root: Path) -> list[str]:
    hits = []
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8-sig", errors="replace")
        if any(pattern.search(text) for pattern in PRIVATE_PATH_PATTERNS):
            hits.append(path.relative_to(root).as_posix())
    return hits


def validate_checksums(root: Path) -> int:
    manifest = root / "CHECKSUMS_SHA256.csv"
    if not manifest.exists():
        raise RuntimeError("CHECKSUMS_SHA256.csv is missing")
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        path = root / row["relative_path"]
        if not path.is_file():
            raise RuntimeError(f"Manifest file missing: {row['relative_path']}")
        if path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"Size mismatch: {row['relative_path']}")
        if sha256(path) != row["sha256"]:
            raise RuntimeError(f"Hash mismatch: {row['relative_path']}")
    expected = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path != manifest
    }
    recorded = {row["relative_path"] for row in rows}
    if expected != recorded:
        raise RuntimeError(f"Checksum coverage mismatch: missing={sorted(expected-recorded)}, extra={sorted(recorded-expected)}")
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", nargs="?", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--skip-checksums", action="store_true", help="Prepare QA before building the checksum manifest")
    args = parser.parse_args()
    root = args.root.resolve()

    figure_report = validate_figures(root)
    json_count = validate_json(root)
    private_hits = validate_private_paths(root)
    if private_hits:
        raise RuntimeError(f"Personal absolute paths remain in: {private_hits}")
    checksum_count = None if args.skip_checksums else validate_checksums(root)
    print(
        json.dumps(
            {
                "status": "PASS",
                "json_files_parsed": json_count,
                "raster_figures_validated": len(figure_report["files"]),
                "checksum_records_validated": checksum_count,
                "personal_absolute_path_hits": 0,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
