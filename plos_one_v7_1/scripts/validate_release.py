#!/usr/bin/env python3
"""Validate the PLOS ONE v7.1 reproducibility package."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "CHECKSUMS_SHA256.csv"
REPORT = ROOT / "qa" / "release_validation.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def main() -> None:
    checks: list[dict[str, object]] = []

    required = [
        ROOT / "README.md",
        ROOT / "RELEASE_NOTES.md",
        ROOT / "CHECKSUMS_SHA256.csv",
        ROOT / "ledgers" / "numeric_lock_v01_2026-09-04.csv",
        ROOT / "ledgers" / "claim_lock_v01_2026-09-04.csv",
        ROOT / "ledgers" / "reference_verification_v7.1_2026-09-04.csv",
        ROOT / "results" / "single_cell" / "gse132465_cluster_inference.csv",
        ROOT / "results" / "bulk_composition" / "target_purged_partial_correlations.csv",
        ROOT / "results" / "common_core" / "spatial_common_core_patient_effects.csv",
        ROOT / "results" / "spatial_null" / "spatial_matched_null_corrected.csv",
        ROOT / "tables" / "S7_Table_spatial_results_and_provenance.csv",
    ]
    missing = [path.relative_to(ROOT).as_posix() for path in required if not path.is_file()]
    record(checks, "required files", not missing, f"missing={missing or 'none'}")

    main_png = sorted((ROOT / "figures" / "main").glob("Fig*_preview.png"))
    main_tif = sorted((ROOT / "figures" / "main").glob("Fig*_review_600dpi.tiff"))
    supp_png = sorted((ROOT / "figures" / "supporting").glob("S*_preview.png"))
    supp_tif = sorted((ROOT / "figures" / "supporting").glob("S*_review_600dpi.tiff"))
    record(checks, "four main PNG figures", len(main_png) == 4, f"count={len(main_png)}")
    record(checks, "four main TIFF figures", len(main_tif) == 4, f"count={len(main_tif)}")
    record(checks, "six supporting PNG figures", len(supp_png) == 6, f"count={len(supp_png)}")
    record(checks, "six supporting TIFF figures", len(supp_tif) == 6, f"count={len(supp_tif)}")

    table_names = {path.name for path in (ROOT / "tables").glob("*.csv")}
    expected_tables = {
        "S1_Table_dataset_accessions_eligibility_and_checksums.csv",
        "S2_Table_signature_membership_overlap_and_representation.csv",
        "S3_Table_single_cell_compartment_statistics.csv",
        "S4_Table_FAP_specific_models_and_diagnostics.csv",
        "S5_Table_dependent_patient_level_correlations.csv",
        "S6_Table_bulk_composition_overlap_and_sensitivity.csv",
        "S7_Table_spatial_results_and_provenance.csv",
    }
    record(checks, "S1-S7 tables", expected_tables.issubset(table_names), f"missing={sorted(expected_tables - table_names) or 'none'}")

    oversized = [
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.stat().st_size > 50 * 1024 * 1024
    ]
    record(checks, "no file exceeds 50 MiB", not oversized, f"files={oversized or 'none'}")

    private_patterns = [
        re.compile(r"[A-Za-z]:[\\/]Users[\\/]", re.I),
        re.compile(r"[A-Za-z]:[\\/]Program Files[\\/]", re.I),
        re.compile(r"/Users/", re.I),
        re.compile(r"/home/", re.I),
    ]
    text_suffixes = {".csv", ".json", ".md", ".txt", ".py", ".r", ".cff"}
    leaks: list[str] = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in text_suffixes:
            continue
        if path == Path(__file__).resolve():
            continue
        content = path.read_text(encoding="utf-8", errors="replace")
        if any(pattern.search(content) for pattern in private_patterns):
            leaks.append(path.relative_to(ROOT).as_posix())
    record(checks, "no local absolute-path leakage", not leaks, f"files={leaks or 'none'}")

    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    manifest_paths = [row["path"] for row in rows]
    duplicate_paths = sorted({path for path in manifest_paths if manifest_paths.count(path) > 1})
    missing_manifest_files: list[str] = []
    mismatches: list[str] = []
    for row in rows:
        path = ROOT / Path(row["path"])
        if not path.is_file():
            missing_manifest_files.append(row["path"])
            continue
        if path.stat().st_size != int(row["size_bytes"]) or sha256(path) != row["sha256"]:
            mismatches.append(row["path"])
    record(checks, "manifest path uniqueness", not duplicate_paths, f"duplicates={duplicate_paths or 'none'}")
    record(checks, "manifest files present", not missing_manifest_files, f"missing={missing_manifest_files or 'none'}")
    record(checks, "manifest hashes", not mismatches, f"mismatches={mismatches or 'none'}")

    allowed_unmanifested = {"CHECKSUMS_SHA256.csv", "qa/release_validation.json"}
    actual_static = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts
    }
    unmanifested = sorted(actual_static.difference(manifest_paths).difference(allowed_unmanifested))
    record(checks, "no unmanifested static files", not unmanifested, f"files={unmanifested or 'none'}")

    failed = [check for check in checks if check["status"] == "FAIL"]
    payload = {
        "release": "v2.0.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "root": ".",
        "manifest_records": len(rows),
        "summary": {"pass": len(checks) - len(failed), "fail": len(failed)},
        "checks": checks,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload["summary"]))
    print(REPORT)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
