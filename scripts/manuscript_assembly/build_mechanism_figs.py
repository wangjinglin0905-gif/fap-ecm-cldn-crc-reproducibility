# -*- coding: utf-8 -*-
"""机制假说图 × 2（插画/示意图风格，供汇报与投稿辅助使用）
图 A：主机制假说图 —— FAP+ CAF 衰老相关分泌状态 → SASP → ECM 重塑 → 上皮 SDC4/CD44 界面（解耦）
图 B：证据与创新概念图 —— 四层数据 → FAP 间质程序 + 衰老创新点 + 方法学堡垒 → 受体解耦结论
输出：PNG (300dpi) + TIFF (600dpi RGB/LZW) + SVG（两图）
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 0.8

# 配色（柔和插画风，色盲友好）
C_EPI = '#F5D5B8'; C_EPI_N = '#D98E4A'   # 上皮浅橙
C_ECM = '#C9E3C0'; C_ECM_F = '#6FA35E'   # ECM 浅绿
C_CAF = '#D6C6EC'; C_CAF_N = '#7D5BA6'   # CAF 浅紫
C_SASP = '#E07A5F'                        # SASP 红
C_REC = '#3D5A80'                         # 受体深蓝
C_GRAY = '#4A4A4A'
C_BOX = '#F2F2F2'
C_BOX_B = '#BBBBBB'

OUTDIR = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0/图片源'

def save_all(fig, name):
    import os
    from PIL import Image
    png = os.path.join(OUTDIR, name + '.png')
    tif = os.path.join(OUTDIR, name + '.tiff')
    svg = os.path.join(OUTDIR, name + '.svg')
    fig.savefig(png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(tif, dpi=600, bbox_inches='tight', facecolor='white')
    Image.open(tif).convert('RGB').save(tif, dpi=(600, 600), compression='tiff_lzw')
    fig.savefig(svg, bbox_inches='tight', facecolor='white')
    print('saved:', name, '| PNG + TIFF(600dpi) + SVG')

def box(ax, x, y, w, h, text, fc=C_BOX, ec=C_BOX_B, fs=9, bold=False, tc='black', r=0.06):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.004,rounding_size={r}",
                       fc=fc, ec=ec, lw=1.0, zorder=3)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fs,
            fontweight='bold' if bold else 'normal', color=tc, zorder=4, linespacing=1.35)

def arrow(ax, x1, y1, x2, y2, color=C_GRAY, lw=1.6, style='-|>', ms=14, ls='-'):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms,
                        color=color, lw=lw, linestyle=ls, zorder=2)
    ax.add_patch(a)

def cell(ax, cx, cy, rx, ry, fc, nc, label=None, lfs=8, lc='black', ldy=-0.28, rot=0):
    e = Ellipse((cx, cy), 2*rx, 2*ry, facecolor=fc, edgecolor=nc, lw=1.2, zorder=3)
    ax.add_patch(e)
    n = Ellipse((cx, cy), 2*rx*0.35, 2*ry*0.35, facecolor=nc, edgecolor='none', alpha=0.55, zorder=4)
    ax.add_patch(n)
    if label:
        ax.text(cx, cy + ldy, label, ha='center', va='center', fontsize=lfs, color=lc,
                zorder=5, rotation=rot)

def collagen(ax, x0, x1, y, amp, period, color=C_ECM_F, lw=2.2, n_waves=3, yoff=0.0):
    for i in range(n_waves):
        xs = np.linspace(x0, x1, 200)
        ys = y + yoff + amp * np.sin(2*np.pi*(xs - x0)/period + i*1.1)
        ax.plot(xs, ys, color=color, lw=lw, alpha=0.85, zorder=2)

def receptor_antenna(ax, x, y, h=0.16, color=C_REC, label=None, lfs=7, ldx=0, ldy=-0.14):
    """Y 形受体天线（SDC4/CD44）"""
    ax.plot([x, x], [y, y+h], color=color, lw=1.6, zorder=4)
    ax.plot([x, x-0.05], [y+h, y+h+0.07], color=color, lw=1.6, zorder=4)
    ax.plot([x, x+0.05], [y+h, y+h+0.07], color=color, lw=1.6, zorder=4)
    if label:
        ax.text(x + ldx, y + ldy, label, ha='center', va='center', fontsize=lfs, color=color, zorder=5)

# ============================================================
# 图 A：主机制假说图
# ============================================================
fig, ax = plt.subplots(figsize=(13.0, 7.2), dpi=100)
ax.set_xlim(0, 13); ax.set_ylim(0, 7.2); ax.axis('off')

# 顶部标题
ax.text(6.5, 6.85, 'A senescence-associated stromal program marked by FAP drives ECM remodeling\n'
        'while leaving epithelial SDC4/CD44 receptors uncoupled (hypothesis)',
        ha='center', va='center', fontsize=12.5, fontweight='bold', color=C_GRAY)

# ---------- 左：肿瘤上皮细胞 ----------
box(ax, 0.3, 4.15, 3.2, 1.9, '', fc=C_EPI, ec='none')
ax.text(1.9, 5.95, 'Tumour epithelium', ha='center', fontsize=9.5, fontweight='bold', color='#8a5a2a')
cell(ax, 1.05, 5.05, 0.42, 0.30, C_EPI, C_EPI_N)
cell(ax, 1.95, 4.95, 0.48, 0.33, C_EPI, C_EPI_N)
cell(ax, 2.80, 5.10, 0.40, 0.28, C_EPI, C_EPI_N)
cell(ax, 1.45, 4.45, 0.38, 0.27, C_EPI, C_EPI_N)
cell(ax, 2.35, 4.50, 0.42, 0.29, C_EPI, C_EPI_N)
# 受体天线（SDC4/CD44 在上皮表面）
for ax_, ay_ in [(1.05, 5.38), (1.95, 5.30), (2.80, 5.40), (1.45, 4.74), (2.35, 4.80)]:
    receptor_antenna(ax, ax_, ay_, h=0.13)
ax.text(1.90, 5.62, 'SDC4 / CD44', ha='center', fontsize=8, color=C_REC, fontweight='bold')
ax.text(1.90, 3.98, '“sensing interface”', ha='center', fontsize=8, color=C_REC, style='italic')

# 解耦标注（受体不随程序上调）
ax.plot([3.75, 3.75], [4.4, 5.6], color='#C44E52', lw=1.4, ls=':', zorder=2)
ax.text(3.82, 5.0, 'no co-induction\n(uncoupled)', ha='left', va='center', fontsize=7.5,
        color='#C44E52', style='italic')

# ---------- 中：ECM ----------
box(ax, 4.35, 4.15, 3.9, 1.9, '', fc=C_ECM, ec='none')
ax.text(6.3, 5.95, 'ECM remodelling', ha='center', fontsize=9.5, fontweight='bold', color='#4a6e3d')
collagen(ax, 4.7, 7.9, 5.35, 0.10, 0.9, n_waves=4, yoff=-0.08)
collagen(ax, 4.7, 7.9, 5.55, 0.09, 1.1, n_waves=4, yoff=0.02, color='#8ab87a')
collagen(ax, 4.7, 7.9, 4.85, 0.08, 0.8, n_waves=4, yoff=-0.05)
ax.text(6.3, 4.42, 'COL1A1 · COL1A2 · COL3A1 · FN1', ha='center', fontsize=8, color='#4a6e3d')

# ---------- 右：FAP+ CAF（衰老分泌状态） ----------
box(ax, 8.95, 4.15, 3.75, 1.9, '', fc=C_CAF, ec='none')
ax.text(10.8, 5.95, 'FAP+ CAF (senescence-associated state)', ha='center', fontsize=9.5,
        fontweight='bold', color='#5a3d8a')
cell(ax, 9.7, 5.0, 0.50, 0.35, C_CAF, C_CAF_N)
cell(ax, 10.8, 5.15, 0.55, 0.38, C_CAF, C_CAF_N)
cell(ax, 11.9, 5.0, 0.45, 0.32, C_CAF, C_CAF_N)
ax.text(10.8, 4.45, 'FAP+  ·  SenMayo-enriched  ·  SASP+', ha='center', fontsize=8,
        color='#5a3d8a', fontweight='bold')

# SASP 分泌因子（红点流）
for sx, sy in [(9.7, 5.55), (10.35, 5.7), (10.8, 5.72), (11.3, 5.65), (11.9, 5.52)]:
    c = Circle((sx, sy), 0.045, facecolor=C_SASP, edgecolor='none', zorder=5)
    ax.add_patch(c)
ax.text(10.8, 6.05, 'SASP factors: IL6 · MMP1/3 · TGFβ · CXCL8 …', ha='center', fontsize=7.5,
        color='#b3503a')

# ---------- 箭头流 ----------
# CAF → ECM：SASP 驱动沉积（红色宽箭头）
arrow(ax, 8.95, 5.35, 8.35, 5.35, color=C_SASP, lw=2.2, ms=18)
ax.text(8.65, 5.62, 'SASP', ha='center', fontsize=8, color=C_SASP, fontweight='bold')
ax.text(8.65, 5.32, 'drives', ha='center', fontsize=7, color='#b3503a')
# ECM → 上皮：基质信号（绿色箭头）
arrow(ax, 4.35, 5.0, 3.65, 5.0, color=C_ECM_F, lw=2.2, ms=18)
ax.text(4.0, 5.30, 'ECM signalling', ha='center', fontsize=8, color='#4a6e3d', fontweight='bold')

# ---------- 底部：临床转化 ----------
box(ax, 0.3, 1.35, 3.1, 1.35, 'FAPI-PET imaging\nstromal-state readout', fc='#EAF2F8', ec=C_REC, fs=9)
box(ax, 3.8, 1.35, 3.1, 1.35, 'Senolytic / senomorphic\n(upstream, unvalidated)', fc='#FDF2F0', ec=C_SASP, fs=9)
box(ax, 7.3, 1.35, 3.1, 1.35, 'ECM / SDC4-CD44\ntargeting (candidate)', fc='#F0F7EE', ec=C_ECM_F, fs=9)
box(ax, 10.8, 1.35, 1.9, 1.35, 'pT1\nendoscopic\ndecision', fc='#F5F0E6', ec='#B8860B', fs=9)
ax.text(6.5, 3.15, 'Translational directions (hypothesis-generating; none validated in CRC)',
        ha='center', fontsize=8.5, style='italic', color='#666666')
# 底部链接线
ax.plot([1.85, 1.85, 5.35, 5.35], [2.70, 3.05, 3.05, 2.70], color='#AAAAAA', lw=1.0, zorder=1)
ax.plot([5.35, 5.35, 8.85, 8.85], [2.70, 3.05, 3.05, 2.70], color='#AAAAAA', lw=1.0, zorder=1)
ax.plot([8.85, 8.85, 11.75, 11.75], [2.70, 3.05, 3.05, 2.70], color='#AAAAAA', lw=1.0, zorder=1)

# 底部数据支撑角标
ax.text(0.35, 0.55, 'Evidence: TCGA (n=380) · GSE39582 (n=566) · CPTAC protein (n=97) · scRNA-seq · spatial',
        fontsize=8, color='#888888', style='italic')
ax.text(0.35, 0.18, 'SenMayo×FAP13 ρ=0.848/0.770 · SASP×ECM ρ=0.570/0.511 · SASP×receptor2 ρ=0.036 (NS) · CPTAC FAP×REC ρ=0.099 (NS)',
        fontsize=8, color='#888888', style='italic')

save_all(fig, 'MechanismFigA_FAP_senescence_ECM_interface')
plt.close(fig)

# ============================================================
# 图 B：证据与创新概念图
# ============================================================
fig, ax = plt.subplots(figsize=(13.0, 7.2), dpi=100)
ax.set_xlim(0, 13); ax.set_ylim(0, 7.2); ax.axis('off')

ax.text(6.5, 6.9, 'Multi-layer evidence converges on a senescence-associated FAP stromal program\n'
        'with an uncoupled epithelial receptor interface',
        ha='center', va='center', fontsize=12.5, fontweight='bold', color=C_GRAY)

# ---------- 顶部：四层数据 ----------
layers = [
    ('Bulk RNA-seq\nTCGA (n=380) · GSE39582 (n=566)', 0.5),
    ('Single-cell\nGSE132465 · GSE166555', 3.55),
    ('Spatial\nQi 2022 · Valdeolivas 2024', 6.6),
    ('Proteomics\nCPTAC-COAD (n=97)', 9.65),
]
for txt, x in layers:
    box(ax, x, 5.45, 2.8, 1.05, txt, fc='#EDF1F7', ec=C_REC, fs=9)

# 汇聚箭头
for x in [1.9, 4.95, 8.0, 11.05]:
    arrow(ax, x, 5.4, x, 4.7, color=C_GRAY, lw=1.6)
ax.plot([1.9, 11.05], [4.75, 4.75], color=C_GRAY, lw=1.2, zorder=1)
arrow(ax, 6.5, 4.75, 6.5, 4.25, color=C_GRAY, lw=1.8, ms=16)

# ---------- 中央：FAP 间质程序 ----------
box(ax, 4.0, 2.85, 5.0, 1.35,
    'FAP-marked stromal matrix program\nFAP13 × matrix4  ρ = 0.932 / 0.913\n(validated across all layers)',
    fc=C_CAF, ec=C_CAF_N, fs=10, bold=True, tc='#3d2a5e')

# 方法学武器（左列）
weapons = ['Non-overlapping scores\n(FAP13/matrix4/receptor2)',
           'Matched empirical nulls\n(10,000 draws, P ≤ 9e-4)',
           'Purity adjustment\n(ABSOLUTE, ρ = 0.906)',
           'MKI67 proliferation\nadjustment (ρ = 0.848/0.759)']
for i, t in enumerate(weapons):
    y = 4.35 - i * 0.72
    box(ax, 0.4, y, 3.1, 0.62, t, fc=C_BOX, ec=C_BOX_B, fs=7.5)
    arrow(ax, 3.55, y + 0.31, 3.95, y + 0.31, color='#999999', lw=1.0, ms=10)

ax.text(1.95, 5.05, 'Methodological armour', ha='center', fontsize=8.5, style='italic', color='#888888')

# ---------- 右侧：衰老创新点 ----------
box(ax, 9.7, 3.6, 3.0, 1.75,
    'Senescence axis (novel)\n\nSenMayo × FAP13\nρ = 0.848 / 0.770\n\nSASP × ECM  ρ = 0.570 / 0.511\nSASP × receptor2  ρ = 0.036 (NS)',
    fc='#FDF2F0', ec=C_SASP, fs=8.5)
arrow(ax, 9.05, 3.9, 9.65, 4.1, color=C_SASP, lw=1.8, ms=14)
ax.text(9.35, 4.45, 'enriched for', ha='center', fontsize=7.5, color=C_SASP, style='italic')

# 蛋白层佐证（右侧下）
box(ax, 9.7, 1.65, 3.0, 1.35,
    'Protein layer support\nCPTAC: FAP×ECM ρ=0.812\nFAP×REC ρ=0.099 (NS) → uncoupled',
    fc='#F0F7EE', ec=C_ECM_F, fs=8)
arrow(ax, 11.2, 3.55, 11.2, 3.05, color='#999999', lw=1.0, ms=10)

# ---------- 底部：结论条带 ----------
box(ax, 0.4, 0.55, 12.2, 0.95,
    'Conclusion: FAP marks a senescence-associated ECM program; epithelial SDC4/CD44 is a candidate sensing interface that is NOT co-induced —\n'
    'receptor uncoupling, reproducible across bulk / single-cell / spatial / protein layers  |  Translation: FAPI imaging · senolytic candidates · pT1 decision (all unvalidated)',
    fc='#FBF6E9', ec='#B8860B', fs=9, bold=False)

save_all(fig, 'MechanismFigB_multi_layer_evidence_senescence')
plt.close(fig)

print('DONE: 2 mechanism figures (PNG + TIFF 600dpi + SVG)')
