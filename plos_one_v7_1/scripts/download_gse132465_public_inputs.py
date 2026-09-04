#!/usr/bin/env python3
"""Download the two public GSE132465 inputs needed for common-core rescoring."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import urllib.request
from pathlib import Path


BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl"
FILES = (
    "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
    "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def download(url: str, destination: Path) -> dict[str, object]:
    partial = destination.with_suffix(destination.suffix + ".part")
    existing = partial.stat().st_size if partial.exists() else 0
    headers = {"User-Agent": "Codex-GSE132465-reproducibility/1.0"}
    if existing:
        headers["Range"] = f"bytes={existing}-"
    request = urllib.request.Request(url, headers=headers)
    start = time.time()
    with urllib.request.urlopen(request, timeout=120) as response:
        status = getattr(response, "status", 200)
        append = existing > 0 and status == 206
        if existing > 0 and not append:
            existing = 0
        mode = "ab" if append else "wb"
        with partial.open(mode) as output:
            while True:
                block = response.read(4 * 1024 * 1024)
                if not block:
                    break
                output.write(block)
    os.replace(partial, destination)
    return {
        "url": url,
        "local_path": str(destination.resolve()),
        "bytes": destination.stat().st_size,
        "sha256": sha256(destination),
        "elapsed_seconds": round(time.time() - start, 3),
    }


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: download_gse132465_public_inputs.py <output directory>")
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    for name in FILES:
        destination = out_dir / name
        if destination.exists() and destination.stat().st_size > 0:
            record = {
                "url": f"{BASE}/{name}",
                "local_path": str(destination.resolve()),
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
                "elapsed_seconds": 0,
                "status": "reused_existing",
            }
        else:
            print(f"Downloading {name}", flush=True)
            record = download(f"{BASE}/{name}", destination)
            record["status"] = "downloaded"
        manifest.append(record)
        print(json.dumps(record, ensure_ascii=False), flush=True)
    (out_dir / "download_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
