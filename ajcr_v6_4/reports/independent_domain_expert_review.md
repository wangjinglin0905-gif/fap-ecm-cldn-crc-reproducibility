# CRC/CAF/FAP/细胞衰老独立领域专家审查

**审查日期：** 2026-08-14  
**审查对象：** `I02_SR_Submitted_Manuscript`、`I03_AJCR_v5.6`、`I04_AJCR_v5.7`、`I06_AJCR_v5.5`、`I08_AJCR_v6.0_candidate`；GSE132465 原始 Seurat 对象；现有单细胞 CSV、Satterthwaite 脚本、fib5 组成校正脚本及 SenMayo 基因表。  
**独立性说明：** 本报告未读取其他代理的审查报告。文献判断基于截至 2026-08-14 可检索到的原始研究和领域指南。  
**裁定：** 对 v6.0 当前“compartmentalized senescence/division of labour”主叙事，建议 **Reject–Rebuild**；对数据本身，仍有可救的 **Major Revision** 路径。

---

## 一、执行结论

v6.0 新增单细胞结果中，真正可支持的是：

1. 在 GSE132465 中，成纤维细胞与上皮细胞的 SenMayo/SASP 分数存在稳定的谱系差异，同时上皮细胞的 MKI67 表达较高；
2. FAP 标记的 CRC 间质与 ECM 程序在 bulk、单细胞和蛋白层面协变，而 SDC4/CD44 未随该程序稳定共诱导；
3. bulk 层面 FAP13 与多种衰老/SASP 相关签名协变，但这仍是转录相关性。

现有证据**不能支持**以下表述：

- “FAP+ CAFs are senescent cells”；
- “CRC exhibits a compartmentalized senescence architecture”；
- “stromal–epithelial division of labour”；
- “FAP+ CAFs carry a senescence-associated matrix program”，若该句被理解为同一细胞状态或因果链；
- “SASP drives ECM deposition”；
- “protein-level SASP–ECM covariation”；
- 已定义的 SASP–免疫细胞通讯轴；
- senolytic-informed treatment rationale。

更重要的是，本次对原始 Seurat 对象的独立核查发现，v6.0 的 B1/B2 主分析使用了**肿瘤与正常细胞混合的 3,462 个成纤维样细胞和 18,539 个上皮细胞**。在肿瘤限定后，B1 的 FAP 检出与 SenMayo/SASP 差异对 `nCount_RNA` 校正高度敏感：未校正时显著，校正测序深度后消失。这个问题不是措辞润色可以解决的，而是主结论的识别条件发生了变化。

因此，若不重新冻结肿瘤细胞宇宙、校正检测深度、转为患者级/亚型内推断并完成独立单细胞验证，不建议以“衰老区室化”为标题和摘要卖点投稿 AJCR。

---

## 二、证据边界：当前结果到底能说到哪一步

| 层级 | 当前证据 | 允许的结论 | 禁止跨越的结论 |
|---|---|---|---|
| 观察 | fibroblast 与 epithelial 的 SenMayo/SASP、MKI67 分数不同 | `compartment-skewed senescence-associated transcriptional features` | `compartmentalized senescence` |
| 关联 | bulk FAP13 与 SenMayo/SASP 协变；FAP 与 ECM 协变 | FAP 标记的组织状态伴随衰老/SASP 样转录 | FAP+ CAF 已被证明衰老 |
| 细胞状态 | p16/p15 原始检测率在 FAP-detected 细胞较高，但深度/亚型校正后消失；LMNB1 原始方向相反 | 最多作为不稳定的辅助观察 | bona fide senescence diagnosis |
| 机制 | SenMayo–matrix 在成纤维含量校正后接近零；B3 无一致通讯模式 | 暂无直接机制证据 | SASP 驱动 ECM；sCAF 调控免疫或上皮增殖 |
| 治疗 | 无功能干预 | 提出待验证假说 | senolytic/senomorphic 治疗依据 |

建议全文统一采用以下总括：

> The data identify a fibroblast-enriched, senescence-associated transcriptional pattern rather than bona fide cellular senescence, a functional stromal–epithelial division of labour, or senescence-driven matrix deposition.

---

## 三、P1：投稿前必须解决的问题

### P1-1. B1/B2 使用了肿瘤与正常混合的细胞宇宙，当前数字不能称为“CRC CAF compartment”

原始 Seurat 元数据直接核查结果如下：

| 细胞集合 | 总数 | Tumor | Normal |
|---|---:|---:|---:|
| fibroblast-lineage（Myofibroblasts + Stromal 1–3） | 3,462 | 1,501 | 1,961 |
| epithelial | 18,539 | 17,469 | 1,070 |
| FAP-detected fibroblast-lineage | 1,025 | 925 | 100 |
| FAP-undetected fibroblast-lineage | 2,437 | 576 | 1,861 |

由此可见：

- v6.0 的 `3,462 + 18,539 = 22,001` 正好是**全对象**细胞，而非 tumor-only；
- FAP-detected 组 90.2% 为肿瘤来源，而 FAP-undetected 组 76.4% 为正常来源；
- 因而 B1 的 FAP 分组与 `Class` 严重混杂；
- 1,961 个正常成纤维样细胞不能称为 CAF；
- 用户提供的 `1,882 FAP+ cells` 是 63,689 个全细胞中 FAP UMI 可检出的细胞数，不等于 FAP+ CAF 数。

**必须修改：**

1. 主分析限定为 tumor-only：1,501 fibroblast-lineage + 17,469 epithelial；
2. Methods 提供从 63,689 到各分析集的纳入流程表，逐项写明 `Class`、`Cell_type`、`Cell_subtype`、患者阈值和缺失处理；
3. 所有 “CAF” 仅用于 tumor-derived fibroblast-lineage cells；若包含正常细胞，只能称 fibroblast-lineage cells；
4. 解释 1,501、1,498、3,462、5,933、17,469、18,539、22,001、1,025、1,369、1,882 之间的关系；当前稿件存在至少三个未闭环的细胞宇宙。

### P1-2. B1 的“FAP-detected fibroblast 衰老富集”主要受 FAP 检出与测序深度共同影响

#### 诊断性重算定义

- 对象：tumor-only fibroblast-lineage cells，`n = 1,501`，23 patients；
- FAP 分组：原始 UMI `FAP > 0` 定义为 FAP-detected；
- SenMayo：125 个源基因，去除与 FAP13 重叠的 4 个后有 121 个候选基因，Seurat 对象实际表示 119 个；对基因作 z 标准化后取均值；
- 主模型：

```text
SenMayo_zmean ~ FAP_status + MKI67_z + log1p(nCount_RNA) + Cell_subtype + (1 | Patient)
```

- SASP 模型同式替换结局；P 值为 Satterthwaite 近似，95% CI 为 Wald CI；
- 这是用于压力测试的独立重算，不是作者锁定 pathB 脚本的逐比特复现，故系数尺度不应直接替换 v6.0 的 β=0.083；但同一诊断管线内的敏感性变化足以识别技术混杂。

#### B1 压力测试结果

| 结局与模型 | FAP-detected β | 95% CI | P |
|---|---:|---:|---:|
| SenMayo；`+ MKI67_z + (1|Patient)` | 0.1051 | 0.0745 to 0.1357 | 2.38 × 10⁻¹¹ |
| SenMayo；再加 `log1p(nCount_RNA)` | −0.0163 | −0.0465 to 0.0138 | 0.289 |
| SenMayo；再加 `Cell_subtype`，不加深度 | 0.1342 | 0.1058 to 0.1626 | 6.64 × 10⁻²⁰ |
| SenMayo；同时加深度与亚型 | 0.0238 | −0.0051 to 0.0527 | 0.106 |
| SASP；同时加深度与亚型 | 0.0323 | −0.0429 to 0.1074 | 0.400 |

结果说明：

- 亚型组成本身不能解释 B1；
- 一旦校正 `nCount_RNA`，FAP-detected 与 SenMayo 的关联消失并轻度反向；
- FAP 是否被 scRNA-seq 检出高度依赖文库深度，而 SenMayo/SASP 多基因均值也随检出基因数增加，因此原模型的极小 P 值主要不能作为生物学稳健性的证据；
- 细胞数带来的高统计功效使非常小的偏差也可产生 `P ≈ 10⁻³⁷`，必须优先报告效应大小、置信区间和患者间一致性。

**必须修改：**

1. 不能继续把原 `β=0.083, P=6.96×10⁻³⁷` 作为摘要第一结果；
2. 重做 tumor-only 分析，并至少加入 `log1p(nCount_RNA)`、`nFeature_RNA` 的预先指定敏感性、fibroblast subtype、患者随机效应；避免在同一模型同时放入高度共线的 nCount 与 nFeature；
3. 首选患者级、亚型内 pseudobulk 或每患者内连续 FAP 表达分析；细胞级模型可作补充；
4. 将 `FAP+` 改为 `FAP-detected`，除非有经验证的阈值或蛋白证据；
5. 在独立 CRC scRNA 队列重复同一分析，且不得把已用于其他结论的 GSE166555 简单称作未使用的独立 senescence 验证。

在这些重算前，B1 应从主结论撤下或如实写成对文库深度敏感的阴性敏感性结果。

### P1-3. p16/p15/p21/LMNB1 不能挽救 senescence diagnosis；校正后均不支持 FAP 特异性

在 tumor-only 1,501 个 fibroblast-lineage cells 中，原始检测率如下：

| Marker | FAP-undetected | FAP-detected | 生物学方向 |
|---|---:|---:|---|
| CDKN2A/p16 | 12.0% | 20.2% | 表面上支持 arrest/stress |
| CDKN2B/p15 | 15.5% | 25.5% | 表面上支持 arrest/stress |
| CDKN1A/p21 | 53.3% | 68.3% | v6.0 未报告，表面上支持 |
| LMNB1 | 3.5% | 7.1% | **与经典 lamin B1 loss 相反** |
| MKI67 | 2.6% | 4.5% | 未显示 FAP-detected 细胞更低增殖 |

对每个 marker 使用以下 logistic mixed model：

```text
marker_detected ~ FAP_status + log1p(nCount_RNA) + Cell_subtype + (1 | Patient)
```

`n = 1,501` tumor fibroblast-lineage cells，23 patients。结果为：

| Marker | 调整后 OR（FAP-detected vs undetected） | 95% CI | P | 备注 |
|---|---:|---:|---:|---|
| CDKN2A | 0.890 | 0.620–1.278 | 0.529 | 不支持独立升高 |
| CDKN2B | 0.981 | 0.703–1.369 | 0.911 | 不支持独立升高 |
| CDKN1A | 0.962 | 0.724–1.278 | 0.788 | 不支持独立升高 |
| LMNB1 | 0.856 | 0.478–1.533 | 0.601 | 模型 singular；原始方向与预期相反 |
| MKI67 | 0.836 | 0.418–1.675 | 0.614 | 不支持稳定增殖抑制 |

按 MICSE 与 SenNet 的原则，组织中 senescence 不能由一个 signature 或一个 marker 诊断。至少需要同一细胞内组合显示稳定 cell-cycle arrest 及两个以上不同类别的辅助特征，例如 DNA damage、lysosomal expansion、nuclear/lamina change 和 SASP。参考：[MICSE guidelines, Cell 2024](https://pubmed.ncbi.nlm.nih.gov/39121846/)；[SenNet recommendations, Nat Rev Mol Cell Biol 2024](https://pubmed.ncbi.nlm.nih.gov/38831121/)。

**必须修改：**

- 删除 “consistent with a senescent state” 一类由原始检测率直接推出状态诊断的句子；
- 完整报告 CDKN1A，并正面报告 LMNB1 方向矛盾；
- 将这些 marker 写成 `descriptive, depth-sensitive detection differences`；
- 结论限于 `senescence-associated transcription`。

### P1-4. B2 的方向较稳，但它首先是谱系差异，不是“区室化衰老”

tumor-only 重算中，fibroblast–epithelial 的 SenMayo 差异仍在 23/23 患者中同向；对 18,970 个 tumor cells 使用：

```text
SenMayo_zmean ~ compartment + MKI67_z + log1p(nCount_RNA) + (1 | Patient)
```

fibroblast compartment 的 β = 0.4835，95% CI 0.4763–0.4907，P < 2 × 10⁻¹⁶。说明 B2 不是单纯由正常细胞或文库深度造成。

但是：

- SenMayo 是跨组织的 senescence/SASP panel，包含大量炎症、分泌和基质相关基因；
- fibroblast 与 epithelial 的基础转录程序不同，直接比较 signature 分数会把 lineage identity、activation、quiescence 与 senescence 混在一起；
- 全部患者同向可证明稳定的谱系差异，不能证明全部患者都有同一 senescent CAF 状态；
- 低 fibroblast MKI67 也可能是静息、终末分化或低增殖，并非 senescence；
- [SenePy](https://pubmed.ncbi.nlm.nih.gov/39987255/) 等工作强调 senescence signature 的细胞类型和组织依赖性。

**必须修改：**

1. 把 `compartmentalized senescence` 改为 `compartment-skewed senescence-associated transcription`；
2. 对 fibroblast–epithelial 对比做表达量/检出率匹配的随机基因集零分布；
3. 使用 cell-type-specific senescence 方法或至少以 SenePy、CRC sCAF/CSPM classifier 作一致性验证；
4. 做 leading-edge 分析，说明差异由 arrest、DDR、SASP、lysosome 还是一般 fibroblast activation 基因驱动；
5. 增加 quiescence、S/G2M、stress/inflammation 和 myofibroblast activation 的竞争解释。

### P1-5. v6.0 回退到了 fib5-only 组成校正，不能声称“survived fibroblast-content adjustment”

v5.7 已给出更保守、也更可信的组成敏感性结果：

| FAP13–matrix4 校正 | partial ρ |
|---|---:|
| fib5 transcript proxy | 0.650 |
| MCP-counter Fibroblasts | 0.113（95% CI −0.015–0.250） |
| EPIC CAF | 0.288（95% CI 0.176–0.398） |
| MCP-counter + EPIC joint | 0.150（95% CI 0.040–0.257） |

这说明高达 0.93 的 bulk FAP13–matrix4 相关大部分追踪 fibroblast abundance，剩余 activation-state 分量大小依赖组成 proxy。v6.0 却恢复 `fib5 partial ρ=0.608` 并写成 “not fibroblast content alone”，属于选择性采用最有利口径。

同样，v5.7 中 SenMayo–FAP13 在多种组成校正后仍有 ρ≈0.378–0.575，但 SenMayo–matrix4 在 EPIC/joint 校正后接近零。因此可以说 SenMayo 相关转录伴随 FAP-marked tissue state，不能说 senescence 与 matrix deposition 形成独立耦合。

**必须修改：**

- 以 v5.7 的 MCP-counter、EPIC、joint 和 fib5 全部结果为最终敏感性框架；
- 摘要和 Discussion 使用 `substantially attenuated and proxy-dependent`；
- 不得只报告 fib5 结果；
- Table 3 中 1,498/1,501 细胞及 matrix4 slope 0.119/0.253 的分析管线差异须另行锁定，不能以“Satterthwaite”一词掩盖归一化和细胞纳入差异。

### P1-6. “protein-level SASP–ECM covariation”是事实性错误

CPTAC 部分测量的是 FAP、COL1A1、COL1A2、FN1、SDC4、CD44 等蛋白。`ρ=0.812` 对应 **FAP 与 FAP-excluded matrix3 protein score**，并没有 SASP 蛋白 panel。

v6.0 在 Results 3.11、Discussion、evidence tier 和 Conclusion 中把该结果写成 `protein-level SASP–ECM covariation`，必须全部删除。正确表述是：

> CPTAC supports FAP–ECM protein covariation and receptor uncoupling; it does not provide a protein-level assay of SASP or cellular senescence.

这项错误若保留，属于审稿人很容易识别的证据错配。

### P1-7. “FAP–senescence–matrix”三联主线没有被同一层级直接证明

目前三个事实来自不同分析：

- FAP 与 matrix 相关；
- FAP13 与 SenMayo/SASP 相关；
- SenMayo 与 matrix 的 bulk 边际相关在组成校正后接近零；
- B1 的单细胞 FAP-detected–SenMayo 关联又对 `nCount_RNA` 敏感。

因此标题中的 `FAP+ CAFs carry a senescence-associated matrix program` 把两个并行轴拼成了同一细胞程序，并隐含机制。若要保留，至少需要：

1. tumor-only、同细胞/同患者的 `matrix4 ~ SenMayo + FAP + nCount + subtype + (1|patient)`；
2. 患者级 pseudobulk 与外部 scRNA 复现；
3. SenMayo–matrix 在 FAP、fibroblast abundance 和 activation state 后仍保留；
4. 实验层面证明 senescence 操作能改变 ECM 生成，且 rescue 可逆转。

在此之前，只能写：

> The bulk data show parallel FAP–matrix and FAP–senescence-associated covariation, but do not establish a unified senescence–matrix program or senescence-driven matrix deposition.

### P1-8. 现有 novelty gap 已被 2019–2026 年 CRC 直接研究占据

v6.0 Introduction 声称 CRC 中 sCAF compartment 尚未系统定义，这一表述截至投稿日已不成立：

- [Guo et al., Aging Cell 2019](https://pubmed.ncbi.nlm.nih.gov/31389184/)：人结肠 senescent fibroblast/SASP，GDF15 促进结直肠上皮增殖、迁移、侵袭和 organoid 生长；
- [Yang et al., Front Oncol 2020](https://pubmed.ncbi.nlm.nih.gov/32318333/)：GALC 诱导 senescent fibroblast 并促进 CRC 生长；
- [Linares et al., Nat Commun 2023](https://doi.org/10.1038/s41467-023-36334-1)：铂类在 CAF 中积累，诱导 senescence/SASP 并促进 CRC 进展和耐药；
- [Hattangady et al., Aging 2024](https://pubmed.ncbi.nlm.nih.gov/38385965/)：构建并验证人原代 colon fibroblast core senescence/SASP profile；
- [Ge et al., J Transl Med 2026](https://pubmed.ncbi.nlm.nih.gov/41535933/)：直接使用 GSE166555 识别 CRC senescent fibroblasts，并以 mIHC、共培养、PDO、小鼠和 senolytic 干预验证 CD36 依赖的脂质转移与 CD8 dysfunction；
- [Yang et al., Cancer Res 2026-07-08](https://pubmed.ncbi.nlm.nih.gov/42418705/)：建立 CSPM 识别 CRC sCAF，并以空间、多重 IHC、PDO/PDOX、原位模型和 fibroblast-specific Il1r1 knockout 证明 IL1B–IL1R1 驱动的 sCAF/SASP 耐药机制。

后两篇在本稿 2026-08-14 投稿准备日前已经公开，且证据强度远高于单一公共 scRNA 队列的 signature 对比。继续写 “CRC sCAF remains unresolved” 会被认为文献检索不完整，也无法回应 Scientific Reports 对 originality 的质疑。

当前可能保留的新意只剩较窄的描述性问题：

> patient-aware quantification of FAP-linked senescence-associated transcription together with a fibroblast–epithelial proliferation contrast and explicit negative boundaries for receptor and immune-communication claims.

但由于 FAP-specific B1 在深度校正后不稳，这个窄新意目前也未完全成立。若不增加独立 scRNA benchmark 和组织/功能验证，建议将本稿定位为保守的 multi-cohort reanalysis，而不是首次定义 CRC sCAF architecture。

### P1-9. SenMayo 基因数和重叠处理不闭环

本地 `senmayo_genes.txt` 有 125 个基因；与 FAP13 重叠 4 个（CXCL12、WNT2、MMP2、MMP9），理论剩余 121 个。GSE132465 对象中实际表示 119 个，缺失 `BEX3` 和 `CCL3L1`。因此 v6.0 的 “120 genes after removal of four” 不正确。

Methods 应写为：

> Of 125 source SenMayo genes, four overlapping FAP13 genes were removed; 119 of the remaining 121 genes were represented in GSE132465 (BEX3 and CCL3L1 absent).

此外，25-gene SASP core 与 FAP13 仍共享 `MMP9`。若用于 bulk FAP13–SASP 相关，必须做 overlap-removed sensitivity；CellAge 和 colon-fibroblast SASP 也应逐一审计重叠，不能只因来源独立即假定统计独立。

---

## 四、P2：显著影响可信度和叙事聚焦的问题

### P2-1. MKI67 调整的解释需降级

MKI67 调整最多说明 SenMayo 分数差异不由**这一个增殖转录本**解释。它不能证明：

- stable growth arrest；
- irreversible senescence；
- epithelial proliferation 与 stromal senescence 构成功能互补；
- 调整后 β 代表与增殖完全无关的 senescence。

Senescence 本身包含 growth arrest，机械地把 MKI67 当作普通混杂变量还可能部分调整掉状态定义。建议增加 S/G2M score、quiescence score、EdU/Ki67 功能证据，并使用：

> MKI67 adjustment indicates that the score contrast is not explained by this single proliferation readout; it does not demonstrate stable cell-cycle arrest.

### P2-2. B3 应作为诚实阴性，不能称 cell communication

B3 是患者级 ligand–receptor 表达相关筛查，不是 CellChat/NicheNet 意义上的经过配体、受体、靶基因和细胞比例建模的通讯证据。结果以负相关和零散 nominal P 为主，且进行了多组相关检验。

建议：

- 主文只保留一句；
- 不在主文列出 CCL2×EGFR、GDF15×CCR2/CXCR2 的 nominal P；
- Supplement S8 报告完整矩阵、测试数和 BH q 值；
- 使用：`The exploratory screen did not identify a reproducible patient-level SASP ligand–immune receptor co-variation pattern.`

阴性结果不能证明没有生物通讯，只能说明该数据和该方法未识别到一致模式。

### P2-3. SDC4–CD44 应降为次要、受限的空间界面假说

SDC4/CD44 的主要证据为：

- bulk 无相关；
- 蛋白无稳定共诱导；
- fibroblast 内患者级相关为负；
- corrected mixed model 为 null；
- 跨区室只有 P=0.06–0.10 的趋势；
- spatial closure 为 null。

因此它们不应与 senescence 主线并列，更不应成为主标题、摘要结论或机制图中心。可保留为 Discussion 一段和 Supplement：

> SDC4/CD44 remain an unvalidated epithelial-interface hypothesis; the present data do not support receptor co-induction or receptor-mediated signaling.

若以衰老为主线，claudin、pT1 switch、TSR nomogram、MMR 和大量 receptor 分支应移入补充材料，否则稿件仍像多个旧分析的堆叠，削弱原创问题。

### P2-4. 不能用 P 值代替效应和生物学独立性

细胞级模型即使有患者随机截距，显著性仍受数千细胞驱动。建议主文同时报告：

- β、95% CI、标准化尺度定义；
- 每患者 FAP 效应及 forest plot；
- 患者级配对/亚型内 pseudobulk；
- random slope 或 cluster-robust sensitivity；
- 每患者每亚型细胞数；
- 预先指定的最小细胞阈值。

### P2-5. 需要与 CRC/colon-fibroblast 特异 signature 对照，而非只堆叠通用 gene sets

SenMayo、CellAge、通用 SASP 的同向并不等于独立生物学复制，三者共享炎症、应激和 secretory biology。建议至少加入：

1. Hattangady colon-fibroblast core senescence profile；
2. 2026 Cancer Research CSPM/sCAF classifier；
3. SenePy fibroblast-specific signature；
4. expression/detection-matched random gene-set null；
5. negative-control signatures：quiescence、myofibroblast activation、TGFβ response、wound healing、general inflammation。

如果 SenMayo 仅在通用分泌程序上显著，而 arrest/DDR/lamina/lysosome 轴不一致，就应称 `SASP-like fibroblast activation`，而非 senescence。

---

## 五、P3：术语、呈现和可复现性改进

1. 单细胞中把 `FAP+`/`FAP−` 改为 `FAP-detected`/`FAP-undetected`；FAP 蛋白阳性与 RNA UMI 检出不是同义词。
2. 把 `epithelial proliferative dominance` 改为 `higher epithelial MKI67 expression` 或 `higher epithelial proliferative activity score`，避免“dominance”暗示功能控制。
3. `division of labour` 只能作为未来模型：`a putative division-of-labour model requiring functional validation`，不宜作结果小标题。
4. Figure 11 应显示 tumor-only、每患者点、效应 CI、细胞数和 score scale；小提琴图轴改为 `FAP-undetected/FAP-detected`。
5. Methods 明确 z-score 是在全部 18,970 tumor cells、各 compartment 内，还是仅 fibroblasts 内计算；不同范围的 z-score 不能混用。
6. B1/B2/B3 的完整脚本、sessionInfo、gene list、每患者汇总表和 frozen object metadata 必须进入 GitHub/Zenodo；占位符投稿前清除。
7. 对所有分析写清主/敏感性/探索性层级，避免同一结果在摘要、Discussion 和 Conclusion 被逐次升级。

---

## 六、可直接用于改稿的主张强度与英文措辞

### 6.1 标题

当前标题不建议使用。若只保留现有稳健证据，推荐：

> **Fibroblast-enriched senescence-associated transcription accompanies higher epithelial proliferative activity in colorectal cancer**

若未来患者级、深度校正和外部队列均重新支持 FAP 特异关联，可使用：

> **FAP-marked colorectal fibroblasts show senescence-associated transcriptional enrichment alongside epithelial proliferative activity**

不建议使用 `compartmentalized senescence programs`、`senescence architecture` 或 `division of labour`。

### 6.2 Abstract Background

> Senescence-associated transcription is heterogeneous across cell types, but whether fibroblast-lineage cells in colorectal cancer preferentially carry such features relative to epithelial cells remains incompletely resolved.

注意：不能写 CRC sCAF 尚未研究，必须接入 2019–2026 年直接 CRC 文献并说明本研究的窄问题。

### 6.3 Abstract Results

在完成 tumor-only 冻结后，推荐按如下结构写：

> In tumor-derived cells from GSE132465, fibroblast-lineage cells had higher SenMayo/SASP scores than epithelial cells across patients, whereas epithelial cells had higher MKI67 expression. The compartment contrast persisted after adjustment for MKI67 and library size. By contrast, the apparent SenMayo/SASP enrichment in FAP-detected fibroblasts was sensitive to library-size adjustment and was not supported by depth- and subtype-adjusted arrest-marker detection models.

如果作者重算后 B1 恢复，应填入患者级效应与 CI，而不是只写极小 P 值。

### 6.4 Abstract/Conclusion 的安全结论

> These data identify a fibroblast-enriched senescence/SASP-associated transcriptional pattern alongside higher epithelial proliferative activity. They do not establish bona fide senescent CAFs, a functional stromal–epithelial division of labour, or senescence-driven matrix deposition.

### 6.5 Results 3.12 的结尾

> These analyses localize a SenMayo/SASP-like transcriptional signal to the fibroblast compartment and show higher epithelial proliferative activity. Because the cross-lineage contrast may reflect cell-type-specific baseline transcription and the FAP-detection contrast was library-depth sensitive, the results do not distinguish senescence from fibroblast activation, quiescence, stress, or lineage identity.

### 6.6 Discussion 4.1 的核心解释

> We therefore interpret the signal as senescence-associated transcription superimposed on fibroblast-lineage identity, rather than evidence that FAP-detected CAFs are bona fide senescent cells. The higher unadjusted detection of CDKN2A, CDKN2B and CDKN1A was not retained after library-depth and subtype adjustment, and LMNB1 was not reduced. These discordant features preclude a cellular-senescence diagnosis.

### 6.7 MKI67 解释

> MKI67 adjustment indicates that the score contrast is not explained by this single proliferation transcript; it does not establish persistent growth arrest or distinguish senescence from quiescence.

### 6.8 FAP–matrix–senescence 边界

> Bulk analyses support FAP–matrix covariation and FAP–senescence-associated covariation, but the SenMayo–matrix relationship was largely composition dependent. The data therefore do not establish a unified senescence-associated matrix program or SASP-driven ECM deposition.

### 6.9 B3 阴性结果

> The exploratory patient-level screen did not identify a reproducible SASP ligand–immune receptor co-variation pattern; no immune-communication inference is made.

### 6.10 SDC4/CD44

> SDC4/CD44 remain an unvalidated epithelial-interface hypothesis. The present bulk, fibroblast-level, protein and spatial data do not support receptor co-induction or receptor-mediated signaling.

---

## 七、建议的稿件重构顺序

在不新增湿实验的情况下，最科学的架构不是把所有旧分析套进“衰老区室化”，而是：

1. **主结果一：** FAP-marked stromal/ECM tissue state；完整采用 v5.7 的组成敏感性，明确高 bulk 相关大部分由 fibroblast abundance 携带。
2. **主结果二：** bulk FAP13 与 senescence/SASP-associated transcription 协变，但 SenMayo–matrix 独立关联弱。
3. **主结果三：** tumor-only fibroblast–epithelial SenMayo/SASP 与 MKI67 对比；只称 lineage/compartment-skewed transcription。
4. **敏感性/阴性结果：** FAP-detected B1 对 library depth 敏感；marker GLMM 阴性；B3 无一致通信模式。
5. **次要假说：** SDC4/CD44 空间界面，限于 Discussion/Supplement。
6. **删除或下沉：** claudin 大分支、pT1 switch、nomogram、MMR/TSR 等不能推进中心问题的分析。

若作者坚持“衰老”为第一卖点，则最低升级条件是：

- tumor-only、depth-adjusted、patient-aware B1 重新成立；
- 独立 CRC scRNA 队列复现；
- CRC/colon-fibroblast-specific classifier 一致；
- 至少一组组织多标志共定位或 ex vivo 功能验证。

否则强行使用 senescence-first 标题，既不会解决 Scientific Reports 的 originality 质疑，还会增加 state-diagnosis 和文献遗漏两个新的 desk-reject 风险。

---

## 八、最低湿实验验证标准

仅做 `FAP + p16/p21` mIHC 不足以证明 senescent CAF。建议按 MICSE/SenNet 构建同一细胞多维证据：

### 8.1 组织定位

- fibroblast identity：FAP + PDGFRA/COL1A1/Vimentin，并以 pan-CK 排除上皮；
- arrest：p16 和/或 p21，同时低 Ki67/PCNA；
- auxiliary markers 至少两类：LMNB1 loss、γH2AX/53BP1、GL13/lipofuscin 或 SA-β-gal/lysosomal marker；
- SASP protein：IL6、CXCL8、GDF15、CXCL12 等，需与同一 FAP+ fibroblast 空间共定位；
- 比较肿瘤、癌旁和分期/治疗背景，避免把正常老化基质当作 CAF senescence。

### 8.2 细胞功能

- 分离 FAP-high 与 FAP-low CAF；
- 洗脱刺激后仍保持 EdU 低、克隆形成低，证明 persistent arrest；
- SA-β-gal/GL13、LMNB1、DDR 和 SASP secretion 联合验证；
- 测量 collagen/fibronectin 生成和 ECM deposition，而不仅是转录。

### 8.3 因果与 rescue

- 诱导或清除 senescence 后测 ECM 和 organoid 表型；
- senomorphic/senolytic 干预必须同时证明 senescent-cell burden 确实下降；
- SASP neutralization 或 conditioned-medium rescue 用于区分 secretome 与 FAP/ECM 本身；
- SDC4/CD44 blockade 只能作为独立的次要界面实验，不应预设为 senescence 下游。

---

## 九、审查产生的可复现文件

以下文件为本次独立压力测试所建，未修改任何源稿件：

- 脚本：`work/analysis/domain_tumor_only_check.R`
- 脚本：`work/analysis/domain_marker_check.R`
- 细胞集合与基因可用性：`work/analysis/single_cell_senescence/cell_counts.csv`
- 基因可用性：`work/analysis/single_cell_senescence/gene_availability.csv`
- B1/B2 诊断性结果：`work/analysis/single_cell_senescence/domain_tumor_only_check.txt`
- marker GLMM 与 depth/subtype-adjusted B1：`work/analysis/single_cell_senescence/domain_tumor_marker_check.txt`

原始作者 pathB 脚本及其 frozen derived matrix 未在最初可用材料中，因此 v6.0 的 β=0.083、β=0.131 和 β=0.361 尚不能视为完成独立逐项复算。上述压力测试已经足以证明：主分析必须显式控制 `Class`、library depth 和 cell subtype，并重新生成所有正文数字与图。

---

## 十、最终建议

**不建议按当前 v6.0 投稿。** 最优策略是先暂停标题和摘要定稿，完成 tumor-only、depth-adjusted、patient-aware 的单细胞重算，并用 v5.7 的多种 fibroblast composition proxies 统一 bulk 口径。

若重算结果仍如本次压力测试，即 B1 在深度校正后为 null，则应放弃“FAP+ senescent CAF”为主卖点，把稿件收缩为：

> FAP-associated matrix activation plus a fibroblast-enriched senescence-associated transcriptional pattern, with explicit negative boundaries for senescent-cell diagnosis, immune communication and SDC4/CD44 signaling.

如果作者希望真正以 senescence 为创新核心，则需要独立单细胞复现和多标志/功能实验；仅通过措辞重构无法补足这一证据缺口。
