# FAP–SDC4/CD44 项目独立生物统计与生信方法学审查

**审查日期：** 2026-08-14  
**审查对象：** Scientific Reports 投稿稿、AJCR v5.5/v5.6/v5.7、AJCR v6.0 candidate，及其可获得的原始/派生分析文件  
**审查性质：** 独立方法学复核，不等同于正式统计学签字  
**独立性声明：** 本审查未读取其他代理生成的审查报告。

## 一、结论先行

**当前 v6.0 不宜按“FAP+ CAF 承载衰老程序、形成 CRC 衰老区室化分工”的标题和摘要直接投稿。** 可重复、可保留的新增结果是：在 GSE132465 中，成纤维样区室的 SenMayo/SASP 表达型评分高于上皮区室，而 MKI67 呈反向分布；这一患者配对差异在限制为肿瘤来源细胞后仍存在。它支持“区室间 senescence-associated transcriptional score 不同”，但不能单独诊断衰老细胞或证明“division of labour”。

v6.0 的核心 B1 结论目前不能成立，原因有三项叠加：

1. **分析人群混入正常组织。** 3,462 个“fibroblast-lineage cells”实际由 1,501 个肿瘤来源细胞和 1,961 个正常来源细胞组成；18,539 个上皮细胞则由 17,469 个肿瘤来源和 1,070 个正常来源细胞组成。3,462/18,539 与 1,501/17,469 的差别不是患者阈值，而是是否限制 `Class == Tumor`。
2. **FAP 状态与组织类别、测序深度严重混杂。** 1,025 个 FAP+ 成纤维样细胞中 925 个来自肿瘤、仅 100 个来自正常；2,437 个 FAP− 细胞中 1,861 个来自正常。FAP+ 细胞中位 UMI 数为 13,082，FAP− 为 6,104。限制肿瘤细胞并校正 UMI/检测复杂度后，SenMayo 与 SASP 的 FAP+ 效应均不再显著。
3. **Methods/legend 声称使用“去除 4 个 FAP13 重叠基因后的 SenMayo”，但 v6.0 报出的 β=0.083、P=6.96×10⁻³⁷、区室差 +0.376、β=0.358/0.361 均只能在保留这 4 个重叠基因时精确复现。** 这不是四舍五入差异，而是评分定义不一致。

此外，v6.0 回退到 fib5 单一校正并采用 partial ρ=0.608，遗漏了 v5.7 已完成的 MCP-counter/EPIC 不利敏感性结果。独立复算显示，FAP13–matrix4 的边际 ρ=0.930 在 MCP-counter、EPIC 和联合校正后分别降至 0.113、0.288 和 0.150。因此，稿件必须沿用 v5.7 的完整组成敏感性框架，不能只保留较有利的 fib5 结果。

## 二、审查范围与独立复算依据

主要核对文件：

- `work/extracted/I02_SR_Submitted_Manuscript.extracted.md`
- `work/extracted/I03_AJCR_v5.6.extracted.md`
- `work/extracted/I04_AJCR_v5.7.extracted.md`
- `work/extracted/I06_AJCR_v5.5.extracted.md`
- `work/extracted/I08_AJCR_v6.0_candidate.extracted.md`
- `work/inputs/analysis/fibroblast_content_sensitivity_analysis.py`
- `work/inputs/analysis/satterthwaite_check*.R`
- `work/inputs/analysis/GSE132465_fibroblast_cells.csv`
- GSE132465 annotation、UMI、25,655×63,689 Seurat 对象，以及 COAD/READ Xena 文件和 MCP-counter/EPIC 派生矩阵。

独立复算脚本：

- `work/analysis/independent_recompute_single_cell.R`
- `work/analysis/independent_recompute_bulk_fibroblast_adjustment.R`
- `work/analysis/inspect_gse132465.py`

复算使用 R 4.6.1、SeuratObject、Matrix、lme4/lmerTest。所有单细胞比较均区分细胞数和患者数；患者是跨个体推断的独立单位。

## 三、主要问题，按严重度排序

| 严重度 | 问题 | 对结论的影响 | 必须处理 |
|---|---|---|---|
| **critical** | B1 将正常与肿瘤成纤维样细胞合并，且 FAP 状态几乎等同于组织类别 | β=0.083 不能解释为“FAP+ CAF 衰老评分升高” | 限制 `Class == Tumor`；把正常组织作为独立敏感性/对照层 |
| **critical** | B1 的 FAP+ 定义基于单基因是否检出，FAP+ 细胞 UMI/基因数显著更高；SenMayo 与 UMI Spearman ρ≈0.56 | 检测深度可同时增加 FAP 阳性概率、SenMayo 基因检出和 CDKN2A/B 检出 | 在模型中校正 log UMI/nFeature；进行患者随机斜率或患者内配对；用 UCell/AUCell/伪 bulk 敏感性 |
| **critical** | v6.0 的 SenMayo 核心数字来自保留 4 个 FAP13 重叠基因的 123-gene score，却写成去重后的 120 genes | Methods、Results、Figure 11 和结论不可复现 | 冻结唯一基因集；单细胞实际去重后为 119 个可用基因 |
| **critical** | v6.0 省略 v5.7 MCP-counter/EPIC 组成校正，仅报告 fib5=0.608 | 形成选择性报告风险，夸大“独立于 fibroblast content” | 恢复 v5.7 Table 8 全部规格及不确定性 |
| **major** | Table 3 的 0.119/−0.005 来自重建 CSV，非文中所述 full-library/tumour-fibroblast 标准化管线；n 写成 1,498 | 分析输入、缩放和样本数不一致 | 冻结 v5.7 full-library 结果：n=1,501，0.253/−0.021 |
| **major** | B1 随机截距模型的自由度约 3,400，假设患者间 FAP 效应完全相同 | 不是完全未处理聚类，但仍不能充分支持对 23 名患者总体的外推 | 增加 `(1+FAP_status\|Patient)` 或患者内效应汇总；报告患者异质性 |
| **major** | B2 的配对区室差可复现，但“低 MKI67 + 高 SenMayo/SASP”并不等于稳定衰老 | 普通静息/细胞身份、SASP 非特异炎症程序均可产生相同模式 | 只写 transcriptional features；加入匹配空基因集和独立队列/蛋白验证 |
| **major** | B3 为 7×16=112 次简单相关；图中仅标 nominal P，无 FDR；不是 CellChat 的患者级验证 | 不能支持 SASP–免疫通讯轴，也不能因无显著而证明“无通讯” | 补完整 P/q 表；主文写“未发现 multiplicity-supported pair”，放补充材料 |
| **major** | AJCR 6.0/work 的 `GSE132465_UMI.txt.gz` 损坏/不完整 | 复现者无法得到一致 library total 和基因覆盖 | 以完整 UMI 或 Seurat counts 为唯一输入，记录 SHA-256，删除/隔离损坏副本 |
| **minor** | 多数核心 β、配对中位差未给 95% CI；Figure 11 c/d 与正文所报评分尺度不清楚 | 读者无法判断效应大小和图文一致性 | 给患者 bootstrap CI、模型 CI、源数据表，并据冻结数据重绘 |

## 四、细胞集关系已完全解析

### 4.1 3,462/18,539 与 1,501/17,469 的真实关系

| `Class` | Fibroblast-lineage：Myofibroblasts + Stromal 1–3 | Epithelial cells |
|---|---:|---:|
| Tumor | 1,501 | 17,469 |
| Normal | 1,961 | 1,070 |
| **合计** | **3,462** | **18,539** |

因此：

- 旧稿 1,501/17,469 是**肿瘤来源细胞集**。
- v6.0 B1/B2 的 3,462/18,539 是**肿瘤与配对正常组织合并的细胞集**。
- v6.0 把 3,462 个细胞全部称为 CAF 不正确；其中 1,961 个来自正常组织，不能称为 cancer-associated fibroblasts。
- Table 4 的第三套人群为全部 `Cell_type == Stromal cells`，共 5,933 个细胞，包含内皮、周细胞、平滑肌、肠神经胶质等，并非同一 fibroblast-lineage 集合；其 FAP+ 1,369/FAP− 4,564 不能与 B1 的 1,025/2,437 混用。

### 4.2 FAP 状态的组织类别混杂

| 来源 | FAP− | FAP+ |
|---|---:|---:|
| Normal | 1,861 | 100 |
| Tumor | 576 | 925 |
| **合计** | **2,437** | **1,025** |

- `P(Tumor | FAP+) = 90.2%`
- `P(Tumor | FAP−) = 23.6%`
- FAP+ 中位 UMI=13,082、nFeature=3,178；FAP− 中位 UMI=6,104、nFeature=1,978。

这意味着原 B1 模型中的 FAP_status 同时编码了肿瘤/正常来源和测序复杂度。仅加患者随机截距和 MKI67 不能消除这两个混杂。

## 五、v5.5/v5.7 口径与 Table 3 的最终判定

### 5.1 bulk fibroblast-content adjustment

v5.5 脚本先将 X/Y/Z 排名并残差化，随后又对残差做 Spearman 相关，得到 0.608。标准 partial Spearman 通常应报告“两个排名变量对排名协变量残差的 Pearson 相关”；按这一标准独立复算得到 0.650。

| 目标 | 校正规格 | partial ρ | 95% bootstrap CI | P |
|---|---|---:|---:|---:|
| FAP13–matrix4 | fib5，标准 partial Spearman | 0.650 | 0.572–0.713 | 5.79×10⁻⁴⁷ |
| FAP13–matrix4 | MCP-counter Fibroblasts | 0.113 | −0.014–0.245 | 0.0278 |
| FAP13–matrix4 | EPIC CAFs | 0.288 | 0.176–0.396 | 1.20×10⁻⁸ |
| FAP13–matrix4 | MCP-counter + EPIC | 0.150 | 0.039–0.255 | 0.00342 |
| overlap-removed SenMayo–FAP13 | fib5 | 0.388 | 0.281–0.484 | 4.92×10⁻¹⁵ |
| overlap-removed SenMayo–FAP13 | MCP-counter | 0.378 | 0.263–0.477 | 2.67×10⁻¹⁴ |
| overlap-removed SenMayo–FAP13 | EPIC | 0.575 | 0.490–0.649 | 8.87×10⁻³⁵ |
| overlap-removed SenMayo–FAP13 | MCP-counter + EPIC | 0.384 | 0.266–0.483 | 1.02×10⁻¹⁴ |
| overlap-removed SenMayo–matrix4 | fib5 | 0.133 | 0.020–0.245 | 0.00963 |
| overlap-removed SenMayo–matrix4 | MCP-counter | −0.098 | −0.220–0.040 | 0.0578 |
| overlap-removed SenMayo–matrix4 | EPIC | 0.082 | −0.014–0.188 | 0.110 |
| overlap-removed SenMayo–matrix4 | MCP-counter + EPIC | −0.039 | −0.130–0.063 | 0.449 |

**允许的解释：** FAP13–matrix4 的高边际相关很大部分由 fibroblast abundance 承载；剩余相关对组成代理高度敏感，不能给出单一“composition-free”效应。SenMayo–FAP13 在各组成规格下仍为中等正相关；SenMayo–matrix4 的直接剩余相关不稳定、总体接近零。

**不允许的写法：** “the association survived fibroblast-content adjustment (partial ρ=0.608), proving an activation state independent of fibroblast abundance.”

### 5.2 Table 3 mixed model

在 1,501 个 `Class == Tumor` fibroblast-lineage cells 上，使用完整 Seurat library normalization、在该分析集内逐基因标准化、患者随机截距和 Satterthwaite 自由度，结果精确复现 v5.7：

| 模型 | n cells / patients | slope | 95% Wald CI | t | P |
|---|---:|---:|---:|---:|---:|
| `matrix4 ~ FAP + (1\|Patient)` | 1,501 / 23 | 0.253 | 0.213–0.294 | 12.21 | 9.17×10⁻³³ |
| `receptor2 ~ FAP + (1\|Patient)` | 1,501 / 23 | −0.021 | −0.062–0.019 | −1.03 | 0.302 |

v6.0 的 0.119/−0.005 可由所给重建 CSV 得到，但该 CSV 的缩放参考集与文中 full-library/tumour-fibroblast 管线不一致；R 模型实际使用 1,501 行而不是稿件写的 1,498。**终稿应冻结 v5.7 的 0.253/−0.021 和 n=1,501。**

## 六、B1：FAP+ versus FAP− 成纤维样细胞

### 6.1 SenMayo 基因集与已报告数字的反向追踪

- SenMayo source list：125 genes。
- 去掉与 FAP13 重叠的 `CXCL12/WNT2/MMP2/MMP9` 后：121 source genes。
- 在 25,655-gene Seurat object 中另缺 `BEX3` 和 `CCL3L1`，故**实际可用为 119 genes，不是 120**。
- 若保留 4 个重叠基因，则可用 123 genes。

复算结果：

| 评分 | B1 β | P | B2 患者配对中位差 | B2 cell-level β | MKI67-adjusted β |
|---|---:|---:|---:|---:|---:|
| **123-gene full SenMayo，保留 4 overlaps** | **0.08334** | **6.96×10⁻³⁷** | **0.37549** | **0.35774** | **0.36119** |
| **119-gene overlap-removed SenMayo** | 0.07705 | 1.90×10⁻³¹ | 0.35098 | 0.32315 | 0.32672 |

第一行与 v6.0 的 β=0.083、P=6.96×10⁻³⁷、+0.376、0.358、0.361 精确吻合。因此 v6.0 的核心数字并非来自 Methods/legend 所称的 overlap-removed score。必须统一 score membership、源数据表、Methods、Results 和 Figure 11。

### 6.2 B1 的正确稳健性结果

未校正的 cell-level 随机截距模型可复现非常小的 P 值，但 SenMayo score 与 nCount/nFeature 的 Spearman 相关约为 0.56/0.55，FAP positivity 又直接依赖单基因检出。加入必要校正后：

| 分析 | SenMayo FAP+ effect | 95% Wald CI | P | SASP effect | 95% Wald CI | P |
|---|---:|---:|---:|---:|---:|---:|
| Tumor-only，119-gene score，校正 MKI67 + log UMI + patient random intercept | −0.0133 | −0.0292–0.0025 | 0.0996 | −0.0200 | −0.0539–0.0140 | 0.249 |
| 全 3,462 cells，校正 Class + MKI67 + log UMI + patient random slope | −0.0101 | −0.0275–0.0074 | 0.271 | −0.0304 | −0.0763–0.0155 | 0.207 |

因此，**B1 的“FAP+ cells carry higher SenMayo/SASP”不具备对组织类别、检测深度和患者间效应异质性的稳健性，不能进入标题、摘要结论或 Discussion 第一段。**

随机截距模型并非完全没有处理伪重复，但其固定效应自由度约为细胞数，隐含患者间相同斜率。随机斜率模型的有效自由度约 21–23，更符合对新患者总体的外推。

### 6.3 CDKN2A/CDKN2B/LMNB1 检出率

未经校正的检出率可复现：

| Marker | FAP− | FAP+ | 方向 |
|---|---:|---:|---|
| CDKN2A | 7.1% | 18.8% | FAP+ 高 |
| CDKN2B | 13.3% | 23.9% | FAP+ 高 |
| LMNB1 | 1.6% | 6.4% | FAP+ 高 |

但是，校正 `Class + log UMI + (1|Patient)` 的混合 logistic model 后：

- CDKN2A：OR=0.86，P=0.349
- CDKN2B：OR=0.83，P=0.190
- LMNB1：OR=1.01，P=0.969

此外，经典衰老证据通常要求 **LMNB1 降低/丢失**；v6.0 未校正结果反而显示 FAP+ 细胞 LMNB1 检出更高。因此这三项不能作为“consistent with senescence”的独立支持。应写成“unadjusted detection-rate differences were depth- and class-sensitive”，并移出摘要。

## 七、B2：上皮–间质区室反差

### 7.1 统计结论

B2 的患者配对设计是新增分析中最可靠的一层。即使限制为肿瘤来源细胞，方向仍然一致：

| 人群与评分 | n patients | 中位 fibroblast−epithelial difference | 95% patient-bootstrap CI | 同方向患者 | Wilcoxon P |
|---|---:|---:|---:|---:|---:|
| 全部 22,001 cells，119-gene fixed score | 23 | 0.351 | 0.311–0.366 | 23/23 | 2.38×10⁻⁷ |
| Tumor-only 18,970 cells，重新在肿瘤细胞内标准化 | 23 | 0.466 | 0.410–0.553 | 23/23 | 2.38×10⁻⁷ |
| Tumor-only，且每区室 ≥20 cells | 15 | 0.466 | 0.410–0.551 | 15/15 | 6.10×10⁻⁵ |
| Tumor-only SASP，重新标准化 | 23 | 0.757 | 0.607–0.863 | 23/23 | 2.38×10⁻⁷ |

Tumor-only MKI67 fibroblast−epithelial difference 为负；23 人分析 P=1.81×10⁻⁴，≥20 cells/compartment 的 15 人分析 P=3.05×10⁻⁴。

### 7.2 允许的声明边界

可以写：

> In GSE132465, tumour fibroblast-lineage cells had higher SenMayo/SASP expression-derived scores and lower MKI67 expression than tumour epithelial cells, consistently across patients.

不应写：

> CRC exhibits a compartmentalized senescence architecture, in which FAP+ CAFs are senescent and the epithelial compartment is assigned a proliferative role.

理由：

- SenMayo/SASP 是表达代理，不是稳定细胞周期停滞的诊断。
- 上皮细胞比成纤维细胞更高 MKI67 是可预期的细胞身份差异，不能专门归因于 senescence。
- SenMayo 包含大量分泌、炎症和基质相关基因；跨细胞类型比较需用大小、表达量、检出率及 fibroblast-specificity 匹配的空基因集，排除“任何 fibroblast-enriched signature 都会显著”的替代解释。
- 患者配对与 cell-level mixed model 使用同一数据，不是两次独立验证。
- 目前只有一个单细胞队列做了 senescence compartment analysis；GSE166555 未完成同构验证。

### 7.3 推荐的 B2 主分析顺序

1. Primary：`Class == Tumor`，患者配对，明确每患者各区室细胞数；主结果同时给 23 人全纳入和 ≥20 cells/compartment 的 15 人敏感性。
2. Primary score：冻结 119-gene overlap-removed SenMayo；SASP 为 secondary。
3. Sensitivity：UCell/AUCell 或患者×区室 pseudobulk；不能在 ssGSEA、z-mean、UCell 中择优报告。
4. Specificity：至少 1,000 个匹配大小/平均表达/检出率/fibroblast-specificity 的空基因集，报告观察到的区室效应在 null 中的分位数。
5. 仅把 cell-level mixed model 作为精度/方向支持，不在摘要用 `P<2×10⁻¹⁶` 作为主要证据。

## 八、B3 与 CellChat 证据层级

v6.0 新的 Supplementary Figure S8 不是 CellChat 模型，而是 7 个 fibroblast SASP ligand × 8 个 receptor × 2 个 immune compartments，共 **112 个患者级 Spearman 相关**。图中星号为 nominal P<0.05；最小名义 P=0.0012，若作 Bonferroni 校正为 0.134。由于缺少完整未四舍五入的源表，无法独立重算精确 BH q 值，也不能冻结任何单一 ligand–receptor pair。

建议结果写法：

> The exploratory 112-pair screen did not yield a multiplicity-supported SASP ligand–immune receptor association; complete nominal and FDR-adjusted results are provided in Supplementary Table X.

不要把散在的正/负相关解释为通讯方向，也不要把“未显著”写成“排除了 SASP–immune communication”。

既有 ECM–SDC4/CD44 CellChat 仍属于 E1 假设生成证据：它基于表达和数据库先验，原 pooled network 不是患者级推断；空间 section 嵌套于患者，不能按 section 当独立样本。该结果没有证明物理结合、方向性或受体依赖，更没有验证 senescent FAP+ CAF 的免疫通讯。

## 九、多重比较、效应量与报告规范

建议冻结以下 hypothesis families：

1. B1 primary family：SenMayo、SASP 两个 FAP-status effects；BH 校正。CDKN2A/B/LMNB1 为 marker sensitivity family，单独 BH。
2. B2 primary family：SenMayo patient-paired contrast；SASP 和 MKI67 为 secondary family，BH 校正。
3. B3：112 pairs，统一 BH；图中只标 q<0.05，不标 nominal 星号，或同时显示 P/q 并明确 exploratory。
4. bulk senescence：SenMayo/SASP/CellAge/colon-SASP × FAP13/matrix4/receptor2 应列完整 family 和 q；不能将同一队列中多个高度相关 signature 视为独立复制。

每个主结果应给：独立患者数、细胞数范围、效应量、95% CI、P/q、模型公式、随机效应、score membership、gene coverage、软件版本及源数据文件。

## 十、可进入终稿的冻结数字表

| 模块 | 可冻结数字 | 终稿定位 |
|---|---|---|
| GSE132465 object | 63,689 cells × 25,655 filtered genes；23 patients | Methods/QC |
| Tumor cell set | 1,501 fibroblast-lineage；17,469 epithelial | 所有肿瘤主分析 |
| Pooled cell set | 3,462 fibroblast-lineage=1,501 tumor+1,961 normal；18,539 epithelial=17,469 tumor+1,070 normal | 仅说明旧 B1/B2 的构成，不作为 CAF 主分析 |
| SenMayo coverage | 125 source；去 4 overlaps 后 121；single-cell represented=119 | Methods、Supplementary gene ledger |
| Bulk marginal | FAP13–matrix4 ρ=0.9296，n=380 | 保留，但紧邻组成敏感性 |
| Bulk composition sensitivity | FAP13–matrix4 partial ρ=0.113 MCP、0.288 EPIC、0.150 joint、0.650 fib5 | 必须完整报告 |
| Bulk senescence after composition | SenMayo–FAP13 partial ρ=0.378–0.575；SenMayo–matrix4 −0.098–0.133 | 支持 FAP-associated covariation；不支持独立 matrix coupling |
| Table 3 full-library | matrix4 slope=0.253，P=9.17×10⁻³³；receptor2 slope=−0.021，P=0.302；n=1,501 | 采用 v5.7 口径 |
| B1 robust analysis | Tumor-only + depth: SenMayo β=−0.013，P=0.100；SASP β=−0.020，P=0.249 | 作为敏感性/否定性结果；不可支持 FAP+ senescence 主结论 |
| B2 tumour-only | SenMayo median paired difference=0.466，95% CI 0.410–0.553，23/23，P=2.38×10⁻⁷ | 可作为单细胞主结果，限于 score difference |
| B2 cell-count sensitivity | ≥20 cells/compartment：n=15，difference=0.466，95% CI 0.410–0.551，P=6.10×10⁻⁵ | 必须并列 |
| B3 | 112 tests；无可冻结的 multiplicity-supported pair | 诚实阴性/补充材料 |

## 十一、不能证实或不得进入主结论的项目

1. **不能证实** “FAP+ CAFs have higher SenMayo/SASP after patient-aware adjustment”。正确限制和深度校正后不显著。
2. **不能证实** “β=0.083 comes from a 120-gene overlap-removed SenMayo score”。该数字精确来自保留 overlaps 的 123-gene score。
3. **不能证实** 0.119/−0.005 是 full-library Table 3 最终口径；完整管线对应 0.253/−0.021。
4. **不能证实** 3,462 与 1,501 的差别来自“患者过滤阈值”；它们分别是 pooled tumour+normal 与 tumour-only。
5. **不能证实** CDKN2A/CDKN2B 检出差异为独立衰老证据；校正 Class/UMI 后均为 null。
6. **不能用** LMNB1 结果支持衰老；未校正方向为 FAP+ 更高，校正后 null。
7. **不能证实** CellChat 或 S8 建立 SASP–immune communication axis。
8. **不能称** CPTAC FAP–ECM ρ=0.812 为“protein-level SASP–ECM covariation”；该分析是 FAP protein 与 FAP-excluded ECM score 的相关，不是 SASP protein assay。
9. **不能证实** senolytic treatment rationale；目前证据最高为 E2 表达关联，B3/CellChat 为 E1。
10. **不能声称** senescence compartment finding 已在独立单细胞队列验证；目前仅 GSE132465。
11. **不能使用** “CRC exhibits a compartmentalized senescence architecture”作为确定性结论；可改为“CRC showed compartment-specific senescence-associated transcriptional scores in one single-cell cohort”。

## 十二、复现性缺口

1. `AJCR 6.0/work/GSE132465_UMI.txt.gz` 在读取 29,044 个基因行、最后到 `SLC6A5` 后触发 gzip EOF；完整 raw UMI 应有 33,694 个基因行、文件约 128,941,455 bytes。25,655 是过滤后的 Seurat object 基因数，不是 raw UMI 行数。
2. B1/B2/B3 的原始执行脚本、每细胞 score table、每患者 paired table、完整 112-pair P/q table未在当前输入包中提供。
3. Figure 11 c/d 的坐标尺度与正文所报 z-score medians 不够透明，且横轴仍为 `FALSE/TRUE`；应从冻结 source table 重新绘制。
4. 终稿复现包需包含：唯一 UMI/Seurat 输入及 SHA-256、annotation checksum、score membership ledger、患者×Class×compartment×FAP-status cell-count table、R sessionInfo、模型日志、全部 null/discordant outputs。

## 十三、终稿重构建议

### 推荐保留

- 以 v5.7 的 full-library Table 3 和 MCP-counter/EPIC 组成敏感性为统一底座。
- 保留 B2，但将其写为单队列、患者配对的 “senescence-associated score compartment contrast”。
- 保留 bulk SenMayo–FAP13 组成校正后的中等相关；明确 SenMayo–matrix4 在组成校正后大多为 null。
- 保留 B3 诚实阴性，放补充材料。

### 必须重做/重写

- B1 主分析改为 tumour-only、119-gene score、log UMI/nFeature 校正、patient random slope；把 null 结果完整报告。
- B2 给 tumour-only 23 人与 ≥20 cells 的 15 人结果、patient-bootstrap CI、匹配 null signature。
- 统一 SenMayo gene coverage；不能再写 120 genes。
- Figure 11 依据冻结数据重绘，图注给患者数、cell count、score scale 和 primary test。

### 建议删除或降级

- 从标题、摘要、结论删除 “FAP+ CAFs carry a senescence-associated program”这一确定性表述。
- 删除“protein-level SASP–ECM covariation”与“senolytic-informed upstream intervention”的现时结论。
- 不把极小 cell-level P 值作为创新性卖点。

### 与现有证据最匹配的题目方向

> **A FAP-associated stromal matrix state in colorectal cancer shows compartment-resolved senescence-associated transcriptional features**

该题目保留“衰老相关”作为新增观察，同时不预支“FAP+ CAF 已被证明衰老”或“区室化分工机制”。若作者坚持原 v6.0 标题，则至少需要独立 CRC 单细胞复现、FAP/p16或p21/低 LMNB1 的多标志组织共定位，以及功能性生长停滞/SASP 验证。

## 十四、最终决策

**当前状态：重大重分析后再投稿。**

- B2 可救，并可成为新增单细胞结果。
- B1 不能按现数字继续支撑 FAP+ senescence 主线。
- bulk 应采用 v5.7 完整组成校正和 full-library 口径。
- 在缺乏独立单细胞和多标志蛋白/功能验证时，论文证据上限为“senescence-associated transcriptional covariation / compartment-specific score difference”，不是 senescent CAF diagnosis、机制或治疗依据。
