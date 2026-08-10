from __future__ import annotations

from pathlib import Path
import csv
import os
import shutil

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path(os.environ.get("FAP_REPO_ROOT", Path.cwd())).resolve()
INPUT = Path(os.environ.get("FAP_FIGURE_INPUT", REPO / "legacy_external" / "submission_panels"))
R13 = REPO / "derived_results" / "senescence_corrected"
OUTPUT = REPO / "figures"
SOURCE_DATA = OUTPUT / "source_data"
OUTPUT.mkdir(parents=True, exist_ok=True)
SOURCE_DATA.mkdir(parents=True, exist_ok=True)


def save_pair(image: Image.Image, stem: str) -> None:
    rgb = image.convert("RGB")
    rgb.save(OUTPUT / f"{stem}.png", dpi=(300, 300), optimize=True)
    rgb.save(OUTPUT / f"{stem}.tiff", dpi=(600, 600), compression="tiff_lzw")


def copy_as_pair(source_stem: str, output_stem: str) -> None:
    image = Image.open(INPUT / f"{source_stem}.png")
    save_pair(image, output_stem)


def panel_label(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, size: int) -> None:
    font_path = font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans", weight="bold"))
    font = ImageFont.truetype(font_path, size)
    draw.text((x, y), label, fill="black", font=font)


def combine_horizontal(left_stem: str, right_stem: str, output_stem: str) -> None:
    left = Image.open(INPUT / f"{left_stem}.png").convert("RGB")
    right = Image.open(INPUT / f"{right_stem}.png").convert("RGB")
    target_h = max(left.height, right.height)
    if left.height != target_h:
        left = left.resize((round(left.width * target_h / left.height), target_h), Image.Resampling.LANCZOS)
    if right.height != target_h:
        right = right.resize((round(right.width * target_h / right.height), target_h), Image.Resampling.LANCZOS)
    gap = max(30, round(target_h * 0.02))
    canvas = Image.new("RGB", (left.width + gap + right.width, target_h), "white")
    canvas.paste(left, (0, 0))
    canvas.paste(right, (left.width + gap, 0))
    draw = ImageDraw.Draw(canvas)
    label_size = max(28, round(target_h * 0.035))
    panel_label(draw, 14, 10, "A", label_size)
    panel_label(draw, left.width + gap + 14, 10, "B", label_size)
    save_pair(canvas, output_stem)


def build_senescence_figure() -> None:
    corr_path = R13 / "P1b_multisignature.csv"
    rows = list(csv.DictReader(corr_path.open(encoding="utf-8-sig")))
    values = {row["comparison"]: float(row["rho"]) for row in rows if row.get("rho")}

    signatures = ["SenMayo_clean", "SASP25_clean", "colonSASP_clean", "CellAge_clean"]
    labels = [
        "SenMayo\n(118/120 genes)",
        "SASP25\n(23/24 genes)",
        "Colon-SASP\n(6/7 genes)",
        "CellAge\n(356/366 genes)",
    ]
    cohorts = ["TCGA380", "GSE39582"]
    colors = {"TCGA380": "#2C7FB8", "GSE39582": "#D95F4A"}

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.6), gridspec_kw={"width_ratios": [1.0, 1.0, 1.35]})
    plt.rcParams.update({"font.family": "Arial", "font.size": 9})

    for ax, outcome, title in zip(
        axes[:2],
        ["FAP13", "matrix4"],
        ["Senescence signatures versus FAP13", "Senescence signatures versus matrix4"],
    ):
        y = np.arange(len(signatures))
        offsets = {"TCGA380": -0.10, "GSE39582": 0.10}
        for cohort in cohorts:
            vals = [values[f"{sig}_{outcome}_{cohort}"] for sig in signatures]
            ax.scatter(vals, y + offsets[cohort], s=52, color=colors[cohort], edgecolor="white", linewidth=0.6,
                       label="TCGA (n = 380)" if cohort == "TCGA380" else "GSE39582 tumors (n = 566)", zorder=3)
            for xx, yy in zip(vals, y + offsets[cohort]):
                ax.text(xx + 0.018, yy, f"{xx:.2f}", va="center", fontsize=7.6, color=colors[cohort])
        ax.axvline(0, color="#666666", linewidth=0.8)
        ax.set_xlim(0, 0.92)
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.set_xlabel("Spearman ρ")
        ax.set_title(title, fontsize=10.5, weight="bold", pad=9)
        ax.grid(axis="x", color="#E6E6E6", linewidth=0.7)
        ax.spines[["top", "right"]].set_visible(False)

    sensitivity = [
        ("Raw", 0.8476265873),
        ("Overlap removed", 0.8308981583),
        ("Adjusted: CAF abundance", 0.7816888949),
        ("Adjusted: tumor fraction", 0.8522294038),
        ("Adjusted: age", 0.8488341329),
        ("Adjusted: MKI67", 0.8483034399),
        ("Adjusted: CAF + purity + MKI67", 0.7583808898),
    ]
    ax = axes[2]
    names = [x[0] for x in sensitivity]
    vals = [x[1] for x in sensitivity]
    y = np.arange(len(names))
    ax.scatter(vals, y, s=54, color="#4C78A8", edgecolor="white", linewidth=0.6, zorder=3)
    for xx, yy in zip(vals, y):
        ax.text(xx + 0.006, yy, f"{xx:.3f}", va="center", fontsize=7.5)
    ax.set_xlim(0.70, 0.88)
    ax.set_yticks(y, names)
    ax.invert_yaxis()
    ax.set_xlabel("Partial or raw Spearman ρ")
    ax.set_title("TCGA SenMayo–FAP13 sensitivity", fontsize=10.5, weight="bold", pad=9)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)

    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, frameon=False, loc="lower center", ncol=2, fontsize=8,
               bbox_to_anchor=(0.37, 0.005))
    for idx, ax in enumerate(axes):
        ax.text(0.01, 0.98, chr(ord("A") + idx), transform=ax.transAxes,
                ha="left", va="top", fontsize=12, fontweight="bold")
    fig.suptitle("Senescence-associated transcriptional covariation after overlap control",
                 fontsize=12.2, fontweight="bold", y=1.02)
    fig.tight_layout(rect=(0, 0.08, 1, 0.96), w_pad=2.2)

    png_path = OUTPUT / "Figure5.png"
    tiff_path = OUTPUT / "Figure5.tiff"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(tiff_path, dpi=600, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    Image.open(png_path).convert("RGB").save(png_path, dpi=(300, 300), optimize=True)
    Image.open(tiff_path).convert("RGB").save(tiff_path, dpi=(600, 600), compression="tiff_lzw")

    source_path = SOURCE_DATA / "Figure5_source_data.csv"
    with source_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["panel", "cohort_or_analysis", "signature_or_adjustment", "outcome", "rho"])
        for outcome in ("FAP13", "matrix4"):
            for cohort in cohorts:
                for sig in signatures:
                    writer.writerow(["A" if outcome == "FAP13" else "B", cohort, sig, outcome,
                                     values[f"{sig}_{outcome}_{cohort}"]])
        for name, rho in sensitivity:
            writer.writerow(["C", "TCGA380", name, "SenMayo_vs_FAP13", rho])


def build_cell_de_figure() -> None:
    rows = [
        ("COMP", 4.81, "ECM/remodelling"),
        ("MMP13", 4.51, "ECM/remodelling"),
        ("EPYC", 4.47, "ECM/remodelling"),
        ("COL1A1", 2.43, "ECM/remodelling"),
        ("COL3A1", 2.10, "ECM/remodelling"),
        ("COL1A2", 1.94, "ECM/remodelling"),
        ("FN1", 1.46, "ECM/remodelling"),
        ("SDC4", 0.39, "Receptor"),
        ("CD44", 0.29, "Receptor"),
    ]
    genes = [x[0] for x in rows][::-1]
    effects = [x[1] for x in rows][::-1]
    groups = [x[2] for x in rows][::-1]
    palette = {"ECM/remodelling": "#8F9AA3", "Receptor": "#4C78A8"}

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    y = np.arange(len(genes))
    ax.barh(y, effects, color=[palette[g] for g in groups], edgecolor="white", linewidth=0.6)
    ax.set_yticks(y, genes)
    ax.set_xlabel("Average log2 fold change (FAP-detected vs FAP-undetected cells)")
    ax.set_title("Selected genes from descriptive cell-level differential expression", fontsize=11, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.7)
    for yy, xx in zip(y, effects):
        ax.text(xx + 0.07, yy, f"{xx:.2f}", va="center", fontsize=8)
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[k]) for k in palette]
    ax.legend(handles, list(palette), frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUTPUT / "SupplementaryFigureS3.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT / "SupplementaryFigureS3.tiff", dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    Image.open(OUTPUT / "SupplementaryFigureS3.png").convert("RGB").save(
        OUTPUT / "SupplementaryFigureS3.png", dpi=(300, 300), optimize=True)
    Image.open(OUTPUT / "SupplementaryFigureS3.tiff").convert("RGB").save(
        OUTPUT / "SupplementaryFigureS3.tiff", dpi=(600, 600), compression="tiff_lzw")

    with (SOURCE_DATA / "SupplementaryFigureS3_source_data.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene", "average_log2_fold_change", "display_group"])
        writer.writerows(rows)


# Five main figures: focused evidence chain. The earlier composite bulk figure is
# excluded because it mixed exploratory inventory analyses with the locked
# 380-primary-tumour inference set and used an overlapping protein score.
copy_as_pair("Figure9", "Figure1")
copy_as_pair("Figure1B", "Figure2")
copy_as_pair("Figure2", "Figure3")
copy_as_pair("Figure3", "Figure4")
build_senescence_figure()

for stale_name in ("Figure6.png", "Figure6.tiff", "source_data/Figure6_source_data.csv"):
    stale = OUTPUT / stale_name
    if stale.exists():
        stale.unlink()

# Supplementary figures. Split panels are merged into one submission file.
combine_horizontal("SupplementaryFigureS1A", "SupplementaryFigureS1B", "SupplementaryFigureS1")
combine_horizontal("SupplementaryFigureS2A", "SupplementaryFigureS2B", "SupplementaryFigureS2")
build_cell_de_figure()
combine_horizontal("SupplementaryFigureS4A", "SupplementaryFigureS4B", "SupplementaryFigureS4")
copy_as_pair("SupplementaryFigureS5", "SupplementaryFigureS5")
combine_horizontal("SupplementaryFigureS6A", "SupplementaryFigureS6B", "SupplementaryFigureS6")
copy_as_pair("SupplementaryFigureS7", "SupplementaryFigureS7")
copy_as_pair("Figure4", "SupplementaryFigureS8")
combine_horizontal("Figure6A", "Figure6B", "SupplementaryFigureS9")
copy_as_pair("Figure8", "SupplementaryFigureS10")

# Preserve machine-readable source tables from the corrected R13 analysis.
for name in ("P0_dual_entry_senescence.csv", "P1a_confounder_adjustment.csv",
             "P1b_multisignature.csv", "P1c_prognosis_cox.csv",
             "signature_coverage_and_overlap.csv"):
    shutil.copy2(R13 / name, SOURCE_DATA / name)

print(f"Built figure package at {OUTPUT}")
