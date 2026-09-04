#!/usr/bin/env python3
"""Validate contiguous HTTP range parts, merge them, and record SHA-256."""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from pathlib import Path


PART_RE = re.compile(r"^(\d+)_(\d+)_(\d+)\.part$")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("Usage: merge_verified_range_parts.py <part directory> <output file> <expected bytes>")
    part_dir = Path(sys.argv[1])
    output = Path(sys.argv[2])
    expected_total = int(sys.argv[3])
    records: list[tuple[int, int, int, Path]] = []
    for path in part_dir.glob("*.part"):
        match = PART_RE.match(path.name)
        if match:
            index, start, end = map(int, match.groups())
            records.append((index, start, end, path))
    records.sort()
    if not records:
        raise RuntimeError("No range parts found")

    cursor = 0
    for index, start, end, path in records:
        if start != cursor:
            raise RuntimeError(f"Non-contiguous part {index}: expected start {cursor}, observed {start}")
        expected_size = end - start + 1
        if path.stat().st_size != expected_size:
            raise RuntimeError(f"Invalid part size for {path}: {path.stat().st_size} != {expected_size}")
        cursor = end + 1
    if cursor != expected_total:
        raise RuntimeError(f"Range coverage {cursor} does not equal expected total {expected_total}")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".merging")
    with temporary.open("wb") as destination:
        for _, _, _, path in records:
            with path.open("rb") as source:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                    destination.write(block)
    if temporary.stat().st_size != expected_total:
        raise RuntimeError("Merged temporary file has unexpected size")
    os.replace(temporary, output)

    manifest = {
        "output": str(output.resolve()),
        "bytes": output.stat().st_size,
        "sha256": digest(output),
        "part_count": len(records),
        "range_start": records[0][1],
        "range_end": records[-1][2],
        "validation": "contiguous ranges and exact part sizes verified before merge",
    }
    manifest_path = output.with_suffix(output.suffix + ".range_manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
