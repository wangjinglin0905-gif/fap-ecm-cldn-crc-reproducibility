#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_gse39582_series_matrix_v2.py
GSE39582_series_matrix.txt.gz format:
  - phenotype lines: !Sample_characteristics_ch1\t"k: v"\t"k: v"... (all samples in one line)
  - expression table: !series_matrix_table_begin -> header "ID_REF" GSM... rows probe\tval\tval...
Output: GSE39582_expr.csv (sample x probe), GSE39582_pheno.csv
"""
import gzip, csv, os, re
from pathlib import Path

REPO = Path(os.environ.get("FAP_REPO_ROOT", Path.cwd())).resolve()
DATA_ROOT = Path(os.environ.get("FAP_DATA_ROOT", REPO / "data")).resolve()
SM = str(DATA_ROOT / "GSE39582" / "GSE39582_series_matrix.txt.gz")
G2P = str(DATA_ROOT / "GSE39582" / "gene2probe.csv")
OUTD = str(DATA_ROOT / "GSE39582")

probes = set()
with open(G2P, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        probes.add(row["probe"].strip())
print(f"[INFO] target probes: {len(probes)}", flush=True)

pheno_keys = {"os.event", "os.delay (months)", "rfs.event", "rfs.delay",
              "tnm.stage", "age.at.diagnosis (year)", "Sex"}

sample_order = []
pheno_rows = []        # list of dict keyed by sample
char_rows = []         # for each characteristic line: list of (key,value) per sample column
table_cols = None      # list of GSM in expression table

def unquote(s):
    return s.strip().strip('"')

# Pass 1: read phenotype + table header
char_rows = []         # for each characteristic line: list of (key,value) per sample column
with gzip.open(SM, "rt", errors="replace", encoding="latin-1") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("!Sample_characteristics_ch1"):
            parts = [unquote(x) for x in line.split("\t")[1:]]
            # parts: ["organism: Homo sapiens", "organism: ...", ...]
            for ci, cell in enumerate(parts):
                if ":" in cell:
                    k, v = cell.split(":", 1)
                    k = k.strip()
                    if k in pheno_keys:
                        while len(char_rows) <= ci:
                            char_rows.append({})
                        char_rows[ci][k] = v.strip()
        elif line.startswith("!series_matrix_table_begin"):
            hdr = next(f).rstrip("\n")
            sample_order = [unquote(x) for x in hdr.split("\t")[1:]]
            break

print(f"[INFO] samples from table header: {len(sample_order)}", flush=True)

# build pheno table (sample x keys)
pheno_out = {}
for ci, gsm in enumerate(sample_order):
    pheno_out[gsm] = {k: (char_rows[ci].get(k, "") if ci < len(char_rows) else "") for k in pheno_keys}

# Pass 2: read expression table, keep target probes
expr_by_gsm = {g: {} for g in sample_order}
with gzip.open(SM, "rt", errors="replace", encoding="latin-1") as f:
    in_tab = False
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("!series_matrix_table_begin"):
            in_tab = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_tab:
            parts = line.split("\t")
            probe = unquote(parts[0])
            if probe == "ID_REF":
                continue
            if probe in probes:
                for gi, gsm in enumerate(sample_order):
                    if gi + 1 < len(parts):
                        try:
                            expr_by_gsm[gsm][probe] = float(unquote(parts[gi + 1]))
                        except ValueError:
                            pass

found = {p for g in expr_by_gsm.values() for p in g}
print(f"[INFO] probes found in table: {len(found)}/{len(probes)}", flush=True)
missing = probes - found
if missing:
    print(f"[WARN] missing probes: {sorted(missing)[:20]}", flush=True)

# write expr
all_probes = sorted(probes)
with open(os.path.join(OUTD, "GSE39582_expr.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["sample"] + all_probes)
    for gsm in sample_order:
        w.writerow([gsm] + [expr_by_gsm[gsm].get(p, "") for p in all_probes])

# write pheno
with open(os.path.join(OUTD, "GSE39582_pheno.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["sample", "os_event", "os_delay_months", "rfs_event", "rfs_delay",
                "tnm_stage", "age_at_diagnosis_year", "Sex"])
    for gsm in sample_order:
        ph = pheno_out.get(gsm, {})
        w.writerow([gsm, ph.get("os.event", ""), ph.get("os.delay (months)", ""),
                    ph.get("rfs.event", ""), ph.get("rfs.delay", ""),
                    ph.get("tnm.stage", ""), ph.get("age.at.diagnosis (year)", ""),
                    ph.get("Sex", "")])

print(f"[DONE] expr {len(sample_order)}x{len(all_probes)} -> GSE39582_expr.csv", flush=True)
print(f"[DONE] pheno {len(sample_order)} -> GSE39582_pheno.csv", flush=True)
