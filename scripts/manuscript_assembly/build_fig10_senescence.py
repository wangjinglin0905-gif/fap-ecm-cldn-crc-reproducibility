# -*- coding: utf-8 -*-
"""Figure 10: The FAP-marked stromal program is enriched for senescence-associated
secretory features. 所有数值来自 R12_senescence_results.csv 与
R12_CPTAC_protein_layer.csv（REVIEW_R12 / R12b 脚本真实输出），不编造数据。
输出：PNG (300 dpi) + TIFF (600 dpi) + SVG（论文插图标准三格式）。
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ---------- 字体 ----------
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['font.size'] = 9

BASE = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0'
FIGDIR = os.path.join(BASE, '图片源')
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# 真实数据（来自 R12 结果 CSV）
# ============================================================
# Panel A：热图 ρ 值（行: 评分 × 队列；列: 目标评分）
rows = [
    ('TCGA', 'SenMayo'), ('TCGA', 'SASP'), ('GSE39582', 'SenMayo'), ('GSE39582', 'SASP'),
]
cols = ['FAP13', 'matrix4', 'receptor2']
rho_A = np.array([
    [0.848, 0.751, 0.042],
    [0.653, 0.570, 0.036],
    [0.770, 0.630, 0.135],
    [0.613, 0.511, 0.154],
])
sig_A = [  # 星号（P 值）：*** <0.001, ** <0.01, * <0.05, ns 不标
    ['***', '***', ''],
    ['***', '***', ''],
    ['***', '***', '**'],
    ['***', '***', '***'],
]

# Panel B：MKI67 校正前后（FAP13 ~ SenMayo）
b_raw = {'TCGA': 0.848, 'GSE39582': 0.770}
b_adj = {'TCGA': 0.848, 'GSE39582': 0.759}

# Panel C：CPTAC 蛋白层（n = 97）
c_labels = ['FAP × ECM', 'FAP × REC', 'ECM × REC']
c_rho = [0.812, 0.099, 0.096]
c_p = ['P < 10⁻²³', 'P = 0.34', 'P = 0.35']

# ============================================================
# 绘图
# ============================================================
fig = plt.figure(figsize=(11.5, 5.6), dpi=100)
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.42,
                      left=0.08, right=0.97, top=0.88, bottom=0.14)

# ---------- Panel A：相关矩阵热图 ----------
axA = fig.add_subplot(gs[0, 0])
vmax = 0.9
cmap = matplotlib.colormaps['RdBu_r']
norm = matplotlib.colors.TwoSlopeNorm(vmin=-0.10, vcenter=0, vmax=vmax)
im = axA.imshow(rho_A, cmap=cmap, norm=norm, aspect='auto')
for i in range(rho_A.shape[0]):
    for j in range(rho_A.shape[1]):
        txt = f"{rho_A[i, j]:.3f}{sig_A[i][j]}"
        c = 'white' if abs(rho_A[i, j]) > 0.55 else 'black'
        axA.text(j, i, txt, ha='center', va='center', fontsize=8.5, color=c, fontweight='bold')
axA.set_xticks(range(3))
axA.set_xticklabels(cols, fontsize=10)
axA.set_yticks(range(4))
axA.set_yticklabels([f"{c} · {s}" for c, s in rows], fontsize=9.5)
axA.set_title('A  Senescence scores vs FAP-stromal program\n(Spearman ρ)', fontsize=10, pad=8)
cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.04)
cb.set_label('Spearman ρ', fontsize=9)

# ---------- Panel B：MKI67 校正 ----------
axB = fig.add_subplot(gs[0, 1])
cohorts = ['TCGA\n(n = 380)', 'GSE39582\n(n = 566)']
x = np.arange(2)
w = 0.30
raw = [b_raw['TCGA'], b_raw['GSE39582']]
adj = [b_adj['TCGA'], b_adj['GSE39582']]
b1 = axB.bar(x - w/2, raw, w, label='Raw', color='#4C72B0', edgecolor='black', linewidth=0.6)
b2 = axB.bar(x + w/2, adj, w, label='MKI67-adjusted', color='#C44E52', edgecolor='black', linewidth=0.6)
for xi, v in zip(x - w/2, raw):
    axB.text(xi, v + 0.02, f"{v:.3f}", ha='center', fontsize=8.5)
for xi, v in zip(x + w/2, adj):
    axB.text(xi, v + 0.02, f"{v:.3f}", ha='center', fontsize=8.5)
axB.axhline(0.5, color='grey', ls='--', lw=0.7)
axB.set_ylim(0, 1.0)
axB.set_xticks(x)
axB.set_xticklabels(cohorts, fontsize=9.5)
axB.set_ylabel('Spearman ρ (FAP13 × SenMayo)', fontsize=9)
axB.set_title('B  Proliferation (MKI67) adjustment\n(partial rank correlation)', fontsize=10, pad=8)
axB.legend(frameon=False, fontsize=8.5, loc='lower right')
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)

# ---------- Panel C：CPTAC 蛋白解耦 ----------
axC = fig.add_subplot(gs[0, 2])
x = np.arange(3)
colors = ['#55A868', '#CCB974', '#CCB974']
bars = axC.bar(x, c_rho, 0.55, color=colors, edgecolor='black', linewidth=0.6)
for xi, v, p in zip(x, c_rho, c_p):
    axC.text(xi, v + 0.03, f"ρ = {v:.3f}\n{p}", ha='center', fontsize=8.3)
axC.axhline(0.5, color='grey', ls='--', lw=0.7)
axC.set_ylim(0, 1.0)
axC.set_xticks(x)
axC.set_xticklabels(c_labels, fontsize=9.5)
axC.set_ylabel('Spearman ρ (CPTAC protein, n = 97)', fontsize=9)
axC.set_title('C  Protein-level SASP/ECM covariation\nvs receptor uncoupling', fontsize=10, pad=8)
axC.spines['top'].set_visible(False)
axC.spines['right'].set_visible(False)

# ---------- 保存 ----------
out_png = os.path.join(FIGDIR, 'Fig10_senescence.png')
out_tiff = os.path.join(FIGDIR, 'Fig10_senescence.tiff')
out_svg = os.path.join(FIGDIR, 'Fig10_senescence.svg')
fig.savefig(out_png, dpi=300, bbox_inches='tight')
fig.savefig(out_tiff, dpi=600, bbox_inches='tight')
fig.savefig(out_svg, bbox_inches='tight')
print('saved:', out_png)
print('saved:', out_tiff)
print('saved:', out_svg)
