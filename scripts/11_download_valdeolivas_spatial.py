from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


RECORD_API = "https://zenodo.org/api/records/7760264"
SAMPLE_PATTERN = re.compile(r"^SN\d+_.+_Rep\d(?:_X)?\.zip$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("work/reproducibility/inputs/valdeolivas_spatial"),
    )
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "FAP-CLDN-reanalysis/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_file(file_record: dict, output_dir: Path) -> dict:
    key = file_record["key"]
    expected_size = int(file_record["size"])
    expected_md5 = file_record["checksum"].split(":", 1)[1]
    destination = output_dir / key
    partial = destination.with_suffix(destination.suffix + ".part")

    if destination.exists() and destination.stat().st_size == expected_size:
        observed_md5 = md5sum(destination)
        if observed_md5 == expected_md5:
            return {
                "file": key,
                "size": expected_size,
                "expected_md5": expected_md5,
                "observed_md5": observed_md5,
                "status": "verified_existing",
            }

    headers = {"User-Agent": "FAP-CLDN-reanalysis/1.0"}
    mode = "wb"
    offset = 0
    if partial.exists():
        offset = partial.stat().st_size
        if offset < expected_size:
            headers["Range"] = f"bytes={offset}-"
            mode = "ab"
        else:
            partial.unlink()
            offset = 0

    request = urllib.request.Request(file_record["links"]["self"], headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        if offset and response.status != 206:
            mode = "wb"
        with partial.open(mode) as handle:
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                handle.write(block)

    if partial.stat().st_size != expected_size:
        raise RuntimeError(
            f"Size mismatch for {key}: {partial.stat().st_size} != {expected_size}"
        )

    observed_md5 = md5sum(partial)
    if observed_md5 != expected_md5:
        raise RuntimeError(f"MD5 mismatch for {key}: {observed_md5} != {expected_md5}")

    partial.replace(destination)
    return {
        "file": key,
        "size": expected_size,
        "expected_md5": expected_md5,
        "observed_md5": observed_md5,
        "status": "downloaded_verified",
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    record = fetch_json(RECORD_API)
    selected = [
        file_record
        for file_record in record["files"]
        if SAMPLE_PATTERN.match(file_record["key"])
        or file_record["key"] == "Pathology_SpotAnnotations.zip"
    ]
    selected.sort(key=lambda item: item["key"])

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(download_file, file_record, args.output_dir): file_record["key"]
            for file_record in selected
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"{result['status']}: {result['file']}", flush=True)

    results.sort(key=lambda item: item["file"])
    manifest_path = args.output_dir / "download_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    metadata = {
        "record_api": RECORD_API,
        "record_id": record["id"],
        "concept_doi": record["conceptdoi"],
        "version_doi": record["doi"],
        "selected_files": len(results),
        "selected_bytes": sum(item["size"] for item in results),
    }
    (args.output_dir / "record_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
