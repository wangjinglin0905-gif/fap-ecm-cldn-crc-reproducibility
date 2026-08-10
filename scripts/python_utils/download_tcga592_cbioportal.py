#!/usr/bin/env python3
"""Download the targeted TCGA-COADREAD PanCancer Atlas expression matrix.

The script intentionally downloads only genes required for the sensitivity
analysis. Run from the repository root. It writes to
data/TCGA592/TCGA592_target_genes_expression.csv.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import time

import requests


REPO = Path(os.environ.get("FAP_REPO_ROOT", Path.cwd())).resolve()
OUT = REPO / "data" / "TCGA592"
OUT.mkdir(parents=True, exist_ok=True)
PROFILE = "coadread_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
SAMPLE_LIST = "coadread_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
API = f"https://www.cbioportal.org/api/molecular-profiles/{PROFILE}/molecular-data/fetch"

FAP13 = ["FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9",
         "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A"]
MATRIX4 = ["COL1A1", "COL1A2", "COL3A1", "FN1"]
RECEPTOR2 = ["SDC4", "CD44"]
SASP25 = ["IL6", "CXCL8", "IL1A", "IL1B", "CCL2", "CCL5", "CXCL1", "CXCL2",
          "CXCL3", "CXCL10", "MMP1", "MMP3", "MMP9", "MMP10", "MMP13",
          "SERPINE1", "PLAU", "TIMP2", "VEGFA", "GDF15", "IGFBP3", "TNF",
          "CSF2", "HGF", "FAS"]


def read_genes(path: Path) -> list[str]:
    genes: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            genes.extend(line.split())
    return sorted(set(genes))


def load_mapping(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            symbol = row.get("SYMBOL", "").strip()
            entrez = row.get("ENTREZID", "").strip()
            if symbol and entrez.isdigit() and symbol not in result:
                result[symbol] = int(entrez)
    return result


def main() -> None:
    senmayo = read_genes(REPO / "config" / "gene_sets" / "SenMayo_125_genes.txt")
    mapping = load_mapping(REPO / "config" / "gene_sets" / "sym2entrez_r.csv")
    symbols = sorted(set(senmayo + FAP13 + MATRIX4 + RECEPTOR2 + SASP25 +
                         ["MKI67", "CDKN2A", "CDKN2B"]))
    present = {symbol: mapping[symbol] for symbol in symbols if symbol in mapping}
    missing = sorted(set(symbols) - set(present))
    if missing:
        print("Symbols without Entrez mapping:", ", ".join(missing))

    records: list[dict] = []
    entrez_ids = sorted(set(present.values()))
    for start in range(0, len(entrez_ids), 30):
        payload = {"entrezGeneIds": entrez_ids[start:start + 30], "sampleListId": SAMPLE_LIST}
        for attempt in range(4):
            try:
                response = requests.post(API, json=payload, timeout=120,
                                         headers={"Accept": "application/json"})
                response.raise_for_status()
                records.extend(response.json())
                break
            except requests.RequestException:
                if attempt == 3:
                    raise
                time.sleep(5)

    entrez_to_symbols: dict[int, list[str]] = {}
    for symbol, entrez in present.items():
        entrez_to_symbols.setdefault(entrez, []).append(symbol)
    values: dict[str, dict[str, float]] = {}
    samples: set[str] = set()
    for record in records:
        sample = record.get("sampleId")
        entrez = record.get("entrezGeneId")
        if not sample or entrez not in entrez_to_symbols:
            continue
        samples.add(sample)
        for symbol in entrez_to_symbols[entrez]:
            values.setdefault(symbol, {})[sample] = record.get("value")

    ordered_samples = sorted(samples)
    ordered_symbols = [s for s in symbols if s in values]
    destination = OUT / "TCGA592_target_genes_expression.csv"
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", *ordered_symbols])
        for sample in ordered_samples:
            writer.writerow([sample, *[values[s].get(sample, "") for s in ordered_symbols]])

    (OUT / "download_metadata.json").write_text(json.dumps({
        "molecular_profile": PROFILE,
        "sample_list": SAMPLE_LIST,
        "n_samples": len(ordered_samples),
        "n_genes": len(ordered_symbols),
        "missing_symbols": missing,
    }, indent=2), encoding="utf-8")
    print(destination)


if __name__ == "__main__":
    main()
