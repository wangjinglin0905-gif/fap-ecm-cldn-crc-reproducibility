# -*- coding: utf-8 -*-
"""Figure 8: Immune/stromal deconvolution and secondary analyses (6-panel).
全部数据来自 JCMM_Revision_2026-07-15_v2/reproducibility_R461（真实计算结果）：
  a) MCP-counter: mcpcounter_high_vs_low.csv
  b) EPIC:        epic_high_vs_low.csv
  c) TGFβ-RPPA:   rppa_tgf_beta_correlations_from_raw_api.csv（95% CI 用 Fisher z 变换计算）
  d) PROGENy:     progeny_high_vs_low_and_burden.csv
  e) DoRothEA:    dorothea_selected_emt_tfs.csv
  f) NicheNet:    L3_NicheNet/nichenet_ligand_activity_dataset_specific.csv（top-15 ligand AUROC）
输出：PNG (300 dpi) + TIFF (600 dpi RGB/LZW) + SVG
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['font.size'] = 8

SD = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/fap_jcmm_submission/JCMM_Revision_2026-07-15_v2/reproducibility_R461'
OUTDIR = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0/图片源'
os.makedirs(OUTDIR, exist_ok=True)

# ---------- 读数据 ----------
mcp = pd.read_csv(os.path.join(SD, 'figures/source_data/mcpcounter_high_vs_low.csv'))
epic = pd.read_csv(os.path.join(SD, 'figures/source_data/epic_high_vs_low.csv'))
rppa = pd.read_csv(os.path.join(SD, 'figures/source_data/rppa_tgf_beta_correlations_from_raw_api.csv'))
prog = pd.read_csv(os.path.join(SD, 'figures/source_data/progeny_high_vs_low_and_burden.csv'))
doro = pd.read_csv(os.path.join(SD, 'figures/source_data/dorothea_selected_emt_tfs.csv'))
nich = pd.read_csv(os.path.join(SD, 'results/L3_NicheNet/nichenet_ligand_activity_dataset_specific.csv'))

def spearman_ci(rho, n, alpha=0.05):
    z = np.arctanh(np.clip(rho, -0.999, 0.999))
    se = 1.0 / np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - alpha / 2)
    lo, hi = np.tanh(z - zc * se), np.tanh(z + zc * se)
    return lo, hi

# ---------- 绘图 ----------
fig, axes = plt.subplots(3, 2, figsize=(11.0, 12.5), dpi=100)
fig.subplots_adjust(left=0.075, right=0.97, top=0.96, bottom=0.045, wspace=0.45, hspace=0.55)

# ---- (a) MCP-counter ----
ax = axes[0, 0]
d = mcp.sort_values('median_difference')
colors = ['#C44E52' if f < 0.05 else '#999999' for f in d['fdr_bh']]
ax.barh(d['cell_type'], d['median_difference'], color=colors, edgecolor='black', linewidth=0.4)
for i, (v, fdr) in enumerate(zip(d['median_difference'], d['fdr_bh'])):
    star = '***' if fdr < 0.001 else ('**' if fdr < 0.01 else ('*' if fdr < 0.05 else ''))
    ax.text(v + 0.02, i, star, va='center', fontsize=8)
ax.axvline(0, color='grey', lw=0.7)
ax.set_title('(a) MCP-counter (FAP-high − FAP-low)', fontsize=9.5)
ax.set_xlabel('Median difference', fontsize=8.5)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ---- (b) EPIC ----
ax = axes[0, 1]
d = epic.sort_values('median_difference')
colors = ['#C44E52' if f < 0.05 else '#999999' for f in d['fdr_bh']]
ax.barh(d['cell_type'], d['median_difference'], color=colors, edgecolor='black', linewidth=0.4)
for i, (v, fdr) in enumerate(zip(d['median_difference'], d['fdr_bh'])):
    star = '***' if fdr < 0.001 else ('**' if fdr < 0.01 else ('*' if fdr < 0.05 else ''))
    ax.text(v + 0.008, i, star, va='center', fontsize=8)
ax.axvline(0, color='grey', lw=0.7)
ax.set_title('(b) EPIC fractions (FAP-high − FAP-low)\n(directionality; see convergence warnings)', fontsize=9.5)
ax.set_xlabel('Median difference (fraction)', fontsize=8.5)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ---- (c) TGFβ-RPPA ----
ax = axes[1, 0]
lab = {'SMAD3_RPPA': 'SMAD3 (RPPA)', 'SMAD4_RPPA': 'SMAD4 (RPPA)', 'CDH1_RPPA': 'E-cadherin (RPPA)'}
d = rppa.set_index('protein').loc[['SMAD3_RPPA', 'SMAD4_RPPA', 'CDH1_RPPA']]
ys = np.arange(len(d))[::-1]
for i, (prot, row) in enumerate(d.iterrows()):
    lo, hi = spearman_ci(row['rho'], int(row['n']))
    ax.errorbar(row['rho'], i, xerr=[[row['rho'] - lo], [hi - row['rho']]], fmt='o',
                color='#4C72B0', ecolor='black', elinewidth=1, capsize=3, ms=5)
    star = '***' if row['p_value'] < 0.001 else ('**' if row['p_value'] < 0.01 else '*')
    ax.text(row['rho'] + (0.05 if row['rho'] > 0 else -0.05), i,
            f"ρ = {row['rho']:.3f}{star}", va='center', ha='left' if row['rho'] > 0 else 'right', fontsize=8)
ax.set_yticks(range(len(d)))
ax.set_yticklabels([lab[p] for p in d.index], fontsize=8.5)
ax.axvline(0, color='grey', lw=0.7)
ax.set_title('(c) TGFβ-RNA score vs RPPA proteins\n(Spearman ρ, 95% CI, n = 316)', fontsize=9.5)
ax.set_xlabel('Spearman ρ', fontsize=8.5)
ax.set_xlim(-0.6, 0.6)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ---- (d) PROGENy ----
ax = axes[1, 1]
d = prog.sort_values('mean_difference')
colors = ['#C44E52' if f < 0.05 else '#999999' for f in d['wilcoxon_fdr_bh']]
ax.barh(d['feature'], d['mean_difference'], color=colors, edgecolor='black', linewidth=0.4)
ax.axvline(0, color='grey', lw=0.7)
ax.set_title('(d) PROGENy pathways (FAP-high − FAP-low)\n(no pathway survives BH-FDR)', fontsize=9.5)
ax.set_xlabel('Mean difference (activity)', fontsize=8.5)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ---- (e) DoRothEA ----
ax = axes[2, 0]
d = doro.sort_values('mean_difference')
colors = ['#C44E52' if f < 0.05 else '#999999' for f in d['wilcoxon_fdr_bh']]
ax.barh(d['feature'], d['mean_difference'], color=colors, edgecolor='black', linewidth=0.4)
ax.axvline(0, color='grey', lw=0.7)
ax.set_title('(e) DoRothEA EMT-related TFs\n(no selected TF survives BH-FDR)', fontsize=9.5)
ax.set_xlabel('Mean difference (TF activity)', fontsize=8.5)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ---- (f) NicheNet ----
ax = axes[2, 1]
d = nich.sort_values('auroc', ascending=False).head(15).iloc[::-1]
ax.barh(d['test_ligand'], d['auroc'], color='#55A868', edgecolor='black', linewidth=0.4)
for i, (lig, au) in enumerate(zip(d['test_ligand'], d['auroc'])):
    ax.text(au + 0.004, i, f"{au:.3f}", va='center', fontsize=7.5)
ax.set_xlim(0.5, 0.72)
ax.set_title('(f) NicheNet ligand activity (top 15)\n(receiver-gene audit: 0 genes at FDR < 0.10)', fontsize=9.5)
ax.set_xlabel('AUROC', fontsize=8.5)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ---------- 保存 ----------
out_png = os.path.join(OUTDIR, 'Fig8_deconvolution.png')
out_tiff = os.path.join(OUTDIR, 'Fig8_deconvolution.tiff')
out_svg = os.path.join(OUTDIR, 'Fig8_deconvolution.svg')
fig.savefig(out_png, dpi=300, bbox_inches='tight')
from PIL import Image
fig.savefig(out_tiff, dpi=600, bbox_inches='tight')
Image.open(out_tiff).convert('RGB').save(out_tiff, dpi=(600, 600), compression='tiff_lzw')
fig.savefig(out_svg, bbox_inches='tight')
print('saved:', out_png)
print('saved:', out_tiff)
print('saved:', out_svg)
