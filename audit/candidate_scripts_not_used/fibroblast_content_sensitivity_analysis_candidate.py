#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fibroblast-content sensitivity analysis for:
"Cross-platform characterization of a FAP-associated stromal matrix state
 and its senescence-related transcriptional covariation in colorectal cancer"

Purpose
-------
Tests whether the primary FAP13-matrix4 association (and the SenMayo
associations) merely tracks fibroblast content. An independent,
non-overlapping pan-fibroblast lineage score (fib5) is used as a
transcript-based proxy for fibroblast content; FAP13 and matrix4 ranks
are residualised against the fib5 rank and the partial Spearman
correlation of the residuals is bootstrapped (5,000 resamples).

Data source
-----------
UCSC Xena, TCGA-COAD/READ HiSeqV2 (pancan-normalised), accessed 2026-08-10:
  https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/HiSeqV2.gz
  https://tcga.xenahubs.net/download/TCGA.READ.sampleMap/HiSeqV2.gz
Primary analysis set: 380 barcode-reviewed primary tumours (sample code 01;
COAD n = 286, READ n = 94) -- identical to the manuscript's primary set.

SenMayo panel: Saul et al., Nat Commun 2022 (MSigDB M45803; 125 genes).
Four genes overlapping FAP13 (CXCL12, WNT2, MMP2, MMP9) are removed, as in
the manuscript; CXCL8 is annotated as IL8 in the HiSeqV2 matrix; BEX3 and
SPX are absent from this annotation version.

Reproduces manuscript values: FAP13-matrix4 rho = 0.930 (n = 380);
READ-only rho = 0.899 (n = 94).

Environment: Python >= 3.10, pandas, numpy, scipy.
"""

import pandas as pd
import numpy as np
from scipy import stats

COAD_FILE = "COAD_HiSeqV2.gz"   # downloaded from UCSC Xena (see header)
READ_FILE = "READ_HiSeqV2.gz"

FAP13 = ["FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9",
         "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A"]
MATRIX4 = ["COL1A1", "COL1A2", "COL3A1", "FN1"]
RECEPTOR2 = ["SDC4", "CD44"]
# Independent pan-fibroblast lineage score: canonical fibroblast markers
# sharing no genes with FAP13 / matrix4 / receptor2
FIB5 = ["PDGFRA", "PDGFRB", "LUM", "DCN", "COL14A1"]

# SenMayo (125 genes; Saul et al. 2022, MSigDB M45803). Full list in
# senmayo_genes.txt alongside this script; loaded here from file.
SENMAYO_FILE = "senmayo_genes.txt"


def load_cohort():
    coad = pd.read_csv(COAD_FILE, sep="\t", index_col=0, compression="gzip")
    read = pd.read_csv(READ_FILE, sep="\t", index_col=0, compression="gzip")
    coad_t = coad[[c for c in coad.columns if c.endswith("-01")]]
    read_t = read[[c for c in read.columns if c.endswith("-01")]]
    common = coad_t.index.intersection(read_t.index)
    return pd.concat([coad_t.loc[common], read_t.loc[common]], axis=1)


def zscore_mean(df, genes):
    sub = df.loc[genes]
    z = sub.apply(lambda r: (r - r.mean()) / r.std(ddof=1), axis=1)
    return z.mean(axis=0)


def partial_spearman(x, y, z):
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)
    resx = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    resy = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return stats.spearmanr(resx, resy)[0]


def boot_ci(x, y, z, n_boot=5000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(x)
    vals = np.empty(n_boot)
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        vals[i] = partial_spearman(x[idx], y[idx], z[idx])
    return np.percentile(vals, [2.5, 97.5])


def main():
    df = load_cohort()
    print(f"primary tumours: {df.shape[1]}")

    senmayo = [g.strip() for g in open(SENMAYO_FILE) if g.strip()]
    senmayo = [g for g in senmayo if g not in FAP13]          # remove 4 overlaps
    senmayo = ["IL8" if g == "CXCL8" else g for g in senmayo]  # HiSeqV2 alias
    senmayo = [g for g in senmayo if g in df.index]
    print(f"SenMayo usable genes: {len(senmayo)}")

    zFAP = zscore_mean(df, FAP13)
    zMAT = zscore_mean(df, MATRIX4)
    zREC = zscore_mean(df, RECEPTOR2)
    zFIB = zscore_mean(df, FIB5)
    zSEN = zscore_mean(df, senmayo)

    print("\n--- marginal Spearman (n = 380) ---")
    print(f"FAP13-matrix4  : rho = {stats.spearmanr(zFAP, zMAT)[0]:.3f}")
    print(f"FAP13-receptor2: rho = {stats.spearmanr(zFAP, zREC)[0]:.3f}, "
          f"P = {stats.spearmanr(zFAP, zREC)[1]:.3f}")
    print(f"SenMayo-FAP13  : rho = {stats.spearmanr(zSEN, zFAP)[0]:.3f}")
    print(f"SenMayo-matrix4: rho = {stats.spearmanr(zSEN, zMAT)[0]:.3f}")

    print("\n--- partial Spearman adjusted for fib5 (fibroblast content) ---")
    for name, x, y in [("FAP13-matrix4", zFAP, zMAT),
                       ("SenMayo-FAP13", zSEN, zFAP),
                       ("SenMayo-matrix4", zSEN, zMAT)]:
        r = partial_spearman(x, y, zFIB)
        lo, hi = boot_ci(x.values, y.values, zFIB.values)
        print(f"{name:16s}: partial rho = {r:.3f} (95% bootstrap CI {lo:.3f}-{hi:.3f})")


if __name__ == "__main__":
    main()
