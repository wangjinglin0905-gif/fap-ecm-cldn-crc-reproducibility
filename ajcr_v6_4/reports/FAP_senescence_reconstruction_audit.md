# FAP–衰老方向论文重构与投稿前审查报告

审查日期：2026-08-14  
目标期刊：*American Journal of Cancer Research*（AJCR）  
审查对象：Scientific Reports 实际投稿稿、AJCR v5.5/v5.6/v5.7、AJCR v6.0 与候选 v6.1、三版 Cover Letter、可获得的原始单细胞对象、派生结果、脚本及图件。  
审查方式：期刊编辑、CRC/CAF/细胞衰老领域、生物信息学/生物统计三条相互独立的审查线；随后进行原始 Seurat 对象复核、关键模型重算、基因集特异性压力测试和图文交叉检查。

## 一、最终裁定

v6.0 不能按现状投稿。三条独立审查线的共同结论是 **Reject–Rebuild**，问题不是英文润色，而是中心结论的识别条件不成立。

可以保留并形成新卖点的结果是：

> 在 GSE132465 的肿瘤来源细胞中，SenMayo 衰老相关转录评分稳定地偏向成纤维细胞区室；这一差异在 23/23 位患者中同向，并超过表达水平、变异度和检出率匹配的随机基因集预期。与此同时，FAP 检出在测序深度和成纤维细胞亚型校正后并不能独立识别这一状态；FAP 主要追踪 ECM 富集状态，而 SenMayo 与 matrix4 在患者层和严格组成校正后可分离。

因此，新稿可使用“**fibroblast-enriched senescence-associated transcription**”作为主线，但不能再使用以下表述：

- `FAP+ CAFs carry a senescence program`；
- `compartmentalized senescence architecture`；
- `stromal–epithelial division of labour`；
- `senescence-associated matrix program`；
- `SASP drives ECM deposition`；
- `protein-level SASP–ECM validation`；
- 已确立的 SASP–免疫通讯轴或 senolytic 治疗依据。

最终重构标题已改为：

> **Fibroblast-enriched senescence-associated transcription and a distinct FAP-linked matrix state in colorectal cancer**

该标题把结果限定在“衰老相关转录”，以并列的正向表述呈现 SenMayo 区室偏向与 FAP–matrix 两条轴，避免标题同时承担结果宣传和防御性否定。

## 二、为什么 v6.0 的中心结论不成立

### 2.1 3,462/18,539 并非“患者阈值不同”，而是混入了正常组织

对 GSE132465 原始 Seurat 元数据的直接核查结果为：

| 细胞集合 | 总数 | Tumor | Normal |
|---|---:|---:|---:|
| Myofibroblasts + Stromal 1–3 | 3,462 | 1,501 | 1,961 |
| Epithelial | 18,539 | 17,469 | 1,070 |
| FAP-detected fibroblast-lineage | 1,025 | 925 | 100 |
| FAP-undetected fibroblast-lineage | 2,437 | 576 | 1,861 |

因此：

- v6.0 的 3,462 + 18,539 是全对象的肿瘤与正常混合集，而非 tumour-only；
- FAP-detected 组 90.2% 来自肿瘤，而 FAP-undetected 组 76.4% 来自正常组织；
- 原 B1 模型中的 FAP 状态同时编码了肿瘤/正常来源；
- 1,961 个正常来源成纤维细胞不能称为 CAF；
- 用户早期记录的 1,882 是 63,689 个全细胞中 FAP UMI 可检出的细胞数，并非 FAP+ CAF 数。

新稿已冻结 tumour-only 主分析集：1,501 个 fibroblast-lineage 细胞、17,469 个 epithelial 细胞、23 位患者。

### 2.2 原 B1 极小 P 值主要受 FAP 检出与 UMI 深度共同影响

FAP-detected 与 FAP-undetected 成纤维样细胞的中位 UMI 约为 13,082 与 6,104。单基因能否检出本身依赖文库深度，多基因平均评分也会随可检出基因数增加，故两者存在共同技术来源。

在 tumour-only 1,501 个成纤维细胞中，独立压力测试得到：

| 结局与模型 | FAP-detected β | 95% CI | P |
|---|---:|---:|---:|
| SenMayo；未加 log UMI | 0.1051 | 0.0745–0.1357 | 2.38×10^-11 |
| SenMayo；加入 log UMI | -0.0163 | -0.0465–0.0138 | 0.289 |
| SenMayo；log UMI + subtype + MKI67 + patient random intercept | 0.0238 | -0.0051–0.0527 | 0.106 |
| SASP；同一完整调整 | 0.0323 | -0.0429–0.1074 | 0.400 |

个体标志物的深度/亚型校正 logistic GLMM 也全部为阴性：

| 标志物 | OR | 95% CI | P |
|---|---:|---:|---:|
| CDKN2A | 0.89 | 0.62–1.28 | 0.529 |
| CDKN2B | 0.98 | 0.70–1.37 | 0.911 |
| CDKN1A | 0.96 | 0.72–1.28 | 0.788 |
| LMNB1 | 0.86 | 0.48–1.53 | 0.601 |
| MKI67 | 0.84 | 0.42–1.67 | 0.614 |

LMNB1 模型出现 singular fit；其原始检出方向还是 FAP-detected 组升高，而经典衰老常见 lamin B1 丢失。故不能用原始 CDKN2A/B 检出率和 LMNB1 方向作为“FAP+ CAF 衰老”的一致性证据。

### 2.3 v6.0 声称使用去重 SenMayo，实际数字来自未去重评分

SenMayo 源集合 125 基因，与 FAP13 重叠 4 个：CXCL12、MMP2、MMP9、WNT2。理论去重后为 121 个；GSE132465 中 BEX3 与 CCL3L1 未表示，因此最终可用为 119 个。

v6.0 的 β=0.083、P=6.96×10^-37、区室差 +0.376、细胞模型 β=0.358/0.361 均只能由保留四个重叠基因的 123-gene 评分精确复现，与 Methods 所写“去除 FAP13 重叠”不一致。新稿统一采用冻结的 119-gene SenMayo，并在 Methods、结果、表格和图例中使用同一口径。

### 2.4 候选 v6.1 与新版 Cover Letter 仍继承了决定性错误

用户后续提供的 `FAP_SDC4_CD44_论文_v6.1_修订版_2026-08-14.docx` 并非已修复的统计版本。逐段提取与交叉检查显示：

- 标题、摘要、Methods 2.13 与 Results 3.12 仍以 3,462 个成纤维样细胞和 18,539 个上皮细胞的 tumour/normal 混合集为基础；
- 仍保留 β=0.083、P=6.96×10^-37、区室差 +0.376 与调整后 +0.361 等不应进入最终稿的数值；
- 仍写 SenMayo 120 genes，而源集合算术与平台表示应为 119；
- matrix/receptor 又回到 1,498-cell、0.119/-0.005 口径，没有采用 v5.7 的 1,501-cell 全库模型；
- Conclusion 与 Cover Letter 继续把 FAP-detected fibroblasts 写成携带衰老程序的 CAF 群体；
- 正文还存在一处整段重复粘贴，提示该文件仍有典型的多版本拼接痕迹。

候选 `AJCR_v6.0_审查报告_2026-08-14.md` 把 v6.0 评价为“方向正确、证据扎实”，但没有识别 tumour/normal 混杂和 UMI 深度混杂，故不能作为独立科学裁定。候选 Cover Letter 同步传播了 +0.376 与 fib5-only 0.608 等旧口径。三份文件仅用于提取作者信息、候选结构和推荐审稿人，不作为数值或结论来源。其 SHA-256 已记录于工作目录的输入账本。

## 三、可保留的新主结果

### 3.1 tumour-only 患者配对区室差异

| 分析 | 结果 | 证据边界 |
|---|---|---|
| SenMayo119，23 位患者 | 23/23 成纤维细胞更高；中位标准化配对差 0.466，95% bootstrap CI 0.410–0.553；exact P=2.38×10^-7 | 稳健的区室评分差异，不诊断衰老细胞 |
| ≥20 cells/compartment 敏感性 | n=15；中位差 0.466，95% CI 0.410–0.551；P=6.10×10^-5 | 对低细胞数患者不敏感 |
| SASP25 | 22/23 成纤维细胞更高；平均 log-normalized 差 0.236；P=4.77×10^-7 | SASP 非特异性更强 |
| MKI67 | 21/23 成纤维细胞更低；平均差 -0.104；P=1.81×10^-4 | 反向增殖分布，但不等于稳定生长停滞 |

### 3.2 表达匹配空基因集检验

5,000 次随机抽样同时匹配基因平均表达、标准差和检出率：

- SenMayo119：观察到的患者平均 fibroblast–epithelial 差为 0.499；null 中位数 0.219，97.5 百分位 0.340；经验 P=0.00020。
- SASP25：观察值 0.739；null 中位数 0.405，97.5 百分位 0.759；经验 P=0.0336。

这一检验说明：SenMayo 的差异幅度超过普通的 lineage-skewed 基因集；但许多匹配随机集合也可出现 23/23 同向，因此“所有患者一致”本身并不等于衰老特异性。SASP 的特异性明显弱于 SenMayo。

### 3.3 FAP、SenMayo 与 matrix4 在患者层可分离

在 tumour-only fibroblast 患者均值中：

- FAP–matrix4：ρ=0.649，P=0.000803；
- FAP–SenMayo：ρ=0.331，P=0.123；
- SenMayo–matrix4：ρ=-0.048，P=0.826。

全库归一化的 1,501-cell 混合模型应冻结为 v5.7 口径：

- matrix4 ~ FAP：slope=0.253，95% CI 0.213–0.294，P=9.17×10^-33；
- receptor2 ~ FAP：slope=-0.021，95% CI -0.062–0.019，P=0.302。

不再采用 v6.0 的 1,498-cell、0.119/-0.005 结果，因为该数值来自不同且未完整描述的缩放管线。

## 四、bulk 结果的最终冻结口径

v6.0 退回了 v5.5 的 fib5 单一校正，并将 ρ=0.608 写成“survived fibroblast-content adjustment”。这属于选择性回退。新稿恢复 v5.7 的完整组成敏感性框架：

| 校正规格 | FAP13–matrix4 | SenMayo–FAP13 | SenMayo–matrix4 |
|---|---:|---:|---:|
| fib5 | 0.650 (0.574–0.715) | 0.388 (0.284–0.485) | 0.133 (0.020–0.248) |
| MCP-counter fibroblasts | 0.113 (-0.015–0.250) | 0.378 (0.270–0.476) | -0.098 (-0.221–0.043) |
| EPIC CAFs | 0.288 (0.176–0.398) | 0.575 (0.490–0.648) | 0.082 (-0.015–0.188) |
| MCP-counter + EPIC | 0.150 (0.040–0.257) | 0.384 (0.273–0.485) | -0.039 (-0.137–0.059) |

解释必须写为：

1. FAP13–matrix4 的 bulk 边际相关很强（TCGA ρ=0.930；GSE39582 ρ=0.913），但主要追踪 fibroblast abundance，剩余量高度依赖代理变量；
2. SenMayo–FAP13 在多种组成规格下较稳定；
3. SenMayo–matrix4 在严格校正后接近零，不能支持统一的“senescence-associated matrix program”；
4. fib5 的 0.650 不能被单独挑出作为主结论。

## 五、机制与蛋白证据边界

### 5.1 CPTAC

CPTAC 可支持 FAP 与不含 FAP 的 COL1A1/COL1A2/FN1 protein score 共变（ρ=0.812，95% CI 0.706–0.879，FDR=1.3×10^-23）。CPTAC 没有形成相应的 SASP 蛋白评分，因此只能称为 **FAP–ECM protein covariation**，不能称为“protein-level SASP–ECM evidence”。

### 5.2 SDC4/CD44

receptor2 在 bulk、tumour fibroblast mixed model 和蛋白层均未随 FAP–matrix 程序稳定共诱导。SDC4/CD44 可以在 Discussion 中作为待验证受体界面，但不能与主轴并列，也不能写成已参与的受体机制。

### 5.3 SASP–免疫分析

旧 Supplementary Figure S8 将 7 个配体与 16 个受体做 112 个任意组合，其中含 CXCL8–TGFBR1 等非 cognate 组合，不能解释为 ligand–receptor 证据。重分析只保留预先列出的 14 个生物学 cognate 组合，在 tumour-only 患者层检验；minimum q=0.598，无一通过 BH 校正。

可用措辞：`No multiplicity-supported cognate pair was identified.`  
不可用措辞：`CellChat confirmed/failed to confirm communication.`

## 六、原创性重新定位

“CRC 中尚未系统刻画 FAP+ senescent CAF”这一空白已不成立。2026 年已有两项更强的 CRC 机制研究：

1. Ge 等在 *Journal of Translational Medicine* 结合单细胞、多重成像与功能实验，证明 senescent fibroblasts 可通过 CD36 相关脂质转移/过氧化导致 CD8+ T 细胞功能障碍；
2. Yang 等在 *Cancer Research* 结合机器学习、空间/多重 IHC、类器官、移植瘤和成纤维细胞特异性 Il1r1 操作，建立 macrophage IL1B–CAF IL1R1–SASP–chemoresistance 机制链。

新稿不再竞争“首个 sCAF”或“机制发现”。可成立的原创点是：

> 以 tumour-only、患者配对、表达匹配空集和多种组成代理为约束，证明 CRC 的 fibroblast-enriched senescence-associated transcription 与 FAP-associated matrix state 并非同一条可互换的轴。

这是方法学澄清与生物学分层，而不是机制证明。它比“FAP 与 ECM 高相关”更有新意，也比 v6.0 的过强叙事更可信。

### 6.1 参考文献扩充与 2026 年文献的实际整合

参考文献由 18 条增加到 27 条。新增文献不是为增加数量，而是分别支撑以下具体论断：

- Isella 与 Calon：CRC bulk 中基质转录对分型和预后信号的贡献；
- González-Gualda 与 SenNet 共识：组织中判定衰老需要多标志、情境化证据；
- Freund：lamin B1 丢失是候选衰老标志，而非单一诊断标准；
- Zimmerman 与 Squair：单细胞数据中患者层级推断及伪重复控制；
- Avila Cobos：bulk 去卷积受参考组成和预处理影响；
- Armingol：配体–受体转录共表达属于通讯推断输入，不等于直接信号测量。

用户指出的两篇 2026 年 CRC–CAF–衰老研究原先已列为参考文献，但旧稿没有充分承担论证功能。本轮将 Ge 等与 Yang 等同时放入 Introduction 和 Discussion：前者作为 CD36 相关脂质转移/免疫功能机制标杆，后者作为 IL1B–IL1R1–SASP–chemoresistance 的空间、类器官和遗传验证标杆。这样既承认最新领域进展，也清楚说明本研究的贡献是组成感知的转录分层，而非重复宣称发现 sCAF 机制。

全部新增条目已经 PubMed/NCBI E-utilities 或出版社原始记录核验；逐条 PMID、DOI 和句子用途见 `reference_expansion_2026-08-14.md`。

## 七、新稿架构调整

### 7.1 已删除或大幅降级

- nomogram、TSR proxy、MMR、左/右半结肠、PROGENy、DoRothEA、NicheNet 等支线；
- 旧 claudin 阶段梯度与大量成员级探索；
- 把 CellChat/空间共表达写成机制的段落；
- senolytic/senomorphic 治疗外推；
- “locked analysis plan”等可能被误解为前瞻预注册的表述。

### 7.2 新正文证据顺序

1. tumour-restricted population：直接定义 63,689 → tumour-only 1,501/17,469，不在正文叙述返修历史；
2. 23 位患者配对的 SenMayo/SASP/MKI67 区室差异；
3. 5,000 个表达匹配空基因集；
4. 外部单细胞队列方向复现；
5. tumour fibroblast 内 FAP 特异效应的阴性深度/亚型校正；
6. 患者层 FAP–matrix、FAP–SenMayo、SenMayo–matrix 分离；
7. v5.7 多代理 bulk 组成敏感性；
8. CPTAC、receptor 与 cognate immune screen 的机制边界；
9. 多标志物组织与功能验证的最低实验标准。

### 7.3 采用“假说—提纲—证据账本—全文收束”的重构方法

本轮不再以旧版本的段落为修订单位，而先冻结中心问题、主结果、第二条轴和证据上限。随后建立 15 项主张账本，每项数字只能绑定一个数据集、样本/细胞宇宙、模型、协变量集合与脚本输出。v5.5、v5.7、v6.0 和候选 v6.1 仅作为片段库；冲突按证据等级和可复现性裁决，不按版本日期裁决。

完成各节后，再进行一次“归纳—收束—总结”通读：结果段先报告稳健观察，在段末集中说明一次证据边界；Discussion 先说明贡献，再与机制研究对照；方法修正史保留在本审查报告而不写入论文。这一流程显著减少了多版本拼接造成的分母漂移、数值回退、章节残影和反复使用 `not/does not/cannot` 的防御性文风。冻结提纲与证据账本见 `manuscript_claim_outline_v2.md`。

### 7.4 证据等级

| 论断 | 最高证据等级 | 允许措辞 |
|---|---:|---|
| fibroblast SenMayo 高于 epithelial | E2，患者配对单细胞关联 | `fibroblast-enriched senescence-associated transcriptional scores` |
| FAP-detected CAF 更衰老 | 不支持 | 必须报告调整后阴性 |
| FAP 与 ECM 同轴 | E2，但 abundance-linked | `FAP-associated matrix state` |
| SenMayo 驱动 ECM | 不支持 | `the two axes were separable` |
| bona fide senescent CAF | E0 | 需 p16/p21、LMNB1、增殖、SA-β-gal/持续停滞等正交证据 |
| SASP–免疫通讯 | 不支持 | `no multiplicity-supported cognate pair` |
| senolytic 治疗 | E0 | 只可列为未来可证伪实验 |

## 八、独立单细胞复现

已按用户授权从 GEO 官方来源下载 GSE166555 原始 count matrices 与 metadata，并冻结以下规则后再分析：

- 仅 tumour 样本；
- CAF-like：原始 `cell_type_str_custom` 中 FBs/CAFs/MyoFBs；
- epithelial：原始 main cell type；
- 每位患者至少 20 个 CAF-like 和 100 个 epithelial 细胞；
- patient-compartment pseudobulk log1p(CPM)；
- 使用同一 119-gene SenMayo，不在外部队列重新选基因；
- 患者配对 exact Wilcoxon 与患者 bootstrap CI。

官方 supplementary files 已完整下载并记录 SHA-256：

- `GSE166555_RAW.tar`：136,632,320 bytes；`c5f7c1927be41f27778f5550de46c25b55bdabd60066b87ba4a6bd0d2c6e6b76`；
- `GSE166555_meta_data.tsv.gz`：80,927,954 bytes；`4aba042a64a71507db4fad088c8a03510a50a4e9b4202b9bab8d6b07e15303d5`。

10 位患者满足冻结的细胞数阈值；外部队列中 119 个预定 SenMayo 基因有 114 个具备可用变异。结果为：

| 外部分析 | 结果 | exact P |
|---|---|---:|
| SenMayo，fibroblast-like vs epithelial | 9/10 患者更高；中位配对差 0.149，95% bootstrap CI 0.044–0.336 | 0.0195 |
| SASP25 | 10/10 更高；中位配对差 0.475，95% CI 0.170–0.740 | 0.00195 |
| MKI67 | 9/10 在 fibroblast-like 更低；中位差 -2.033，95% CI -2.878 至 -0.874 | 0.00391 |

外部 fibroblast-like pseudobulk 还复现了轴分离：

- FAP13–matrix4：ρ=0.879，P=0.000814；
- FAP13–SenMayo：ρ=-0.067，P=0.855；
- SenMayo–matrix4：ρ=-0.406，P=0.244。

这一外部结果显著增强了新叙事：区室方向可复现，FAP–matrix 可复现，但 FAP–SenMayo 和 SenMayo–matrix 不能复现为同一轴。由于两队列的注释和标准化不同，效应量不应直接合并；结果仍是转录关联，不能升级为细胞衰老诊断。

## 九、图件审查与重绘

已重绘的主图均为 600 dpi、RGB、Arial 图中文字，并同时导出 PNG 与 LZW-compressed TIFF：

| 图 | 结论 | 主要风险控制 |
|---|---|---|
| Figure 1 | tumour-only 患者配对 SenMayo/SASP/MKI67 区室分布 | 患者为单位；明确不是衰老诊断 |
| Figure 2 | FAP 调整后阴性与 FAP/matrix/SenMayo 三轴分离 | 展示效应量、95% CI 与阴性结果 |
| Figure 3 | 表达匹配 signature null + bulk 多代理组成校正 | 不选择有利代理；直接展示不确定性 |
| Figure 4 | GSE166555 外部队列 | 同一冻结基因集与患者配对口径 |
| Figure 5 | tumour-only cognate SASP–immune screen | 只保留生物学 cognate pair；BH 校正 |

旧 Figure 11 的 TRUE/FALSE 横轴和包含正常细胞的结果不再使用；旧 S8 任意 ligand×receptor 热图不再使用。AJCR 不接收 supplementary figures/tables，因此所有保留图表均已并入主稿，其余删除。

## 十、AJCR 格式核查

已按 2026-08-14 官方 Instructions/Submission checklist 执行：

- 标题控制为两行内的简洁标题；
- 提供 running title；
- 摘要为单一、无引文段落，当前 268 词，小于 350 词；
- 正文 Franklin Gothic Book 12 pt、单倍行距；
- 表格为 Word 原生表格，置于参考文献后；
- figure legends 位于图前；
- 图置于文末，正文依序引用；
- 图中文字 Arial；
- TIFF/PNG 均为 600 dpi；
- 不保留 supplementary display items。

仍有两个明确投稿阻断项，已在 DOCX 中黄色标记：

1. 补齐 AJCR 投稿系统所需的所有作者邮箱；
2. 填入公开 GitHub URL 与 Zenodo DOI，并上传冻结脚本、source tables、checksum 和 sessionInfo。

## 十一、R 环境核查

本轮分析固定调用：

```text
C:\Program Files\R\R-4.6.1\bin\Rscript.exe --vanilla
R_LIBS_USER=<R-4.6-user-library>
```

TMP/TEMP/TMPDIR 指向项目内临时目录。R 4.5.2 与 4.6.0 备份未进入 PATH 或 `.libPaths()`；未调用与 R 4.6.1 不兼容的 R 4.5.2-compiled `data.table`。关键包可正常加载：SeuratObject 5.4.0、Matrix 1.7.5、lme4 2.0.6、lmerTest 3.2.1、ggplot2 4.0.3、patchwork 1.3.2。

## 十二、投稿前最低湿实验标准

若希望把标题升级回“FAP+ senescent CAF”，最低需要：

1. 在独立 CRC 队列进行 stromal/epithelial segmentation；
2. FAP 与 p16、p21、多重增殖标志、lamin B1 loss 同片共定位；
3. 至少一种衰老相关组织/细胞功能指标，如 SA-β-gal、持久性生长停滞或 DNA-damage response；
4. 分离或培养 FAP-high/FAP-low fibroblast，并对活性、测序深度无关的蛋白标志和 CAF 亚型进行匹配；
5. 若主张 SASP–免疫或 matrix 机制，需通过 conditioned medium、阻断/敲除和 rescue 建立方向性。

在没有这些数据时，当前重构稿的最佳定位是：**严谨的公共数据多层复核，揭示衰老相关转录的区室偏向及其与 FAP–matrix 状态的可分离性。**

## 十三、文件与可复现性记录

- v6.0 候选稿 SHA-256：`1d66ea09b8b42fdfcb7d5b6cc78902de908dd5c9714352388c00593aeccf8554`；
- 原始 Seurat 对象：25,655 genes × 63,689 cells，23 patients；
- 单细胞 tumour-only 重算、marker GLMM、signature-null、cognate screen 均保留脚本、CSV/TXT 结果及 R sessionInfo；
- 外部 GSE166555 下载文件保留大小、SHA-256、原始 GEO accession、分析 ledger 与 session information；
- 最终 Word SHA-256：`1d27e474b725c1da7169a9d0619dd9acf10236d3cf5719c9375a74e6fe22ff53`；DOCX 压缩包完整性检测通过；
- 最终 Word 已通过文档渲染流程导出为 20 页 PDF（SHA-256：`1a521f65feb8794616577946357b3acbae71a9a87b5abfbb4954d03ab62553cb`），并以 150 dpi 全页渲染和 5 张顺序 contact sheet 检查；未见文字溢出、表格截断、空白图框、错位图例或旧 Figure 11/S8 残留；
- 重构 Cover Letter SHA-256：`e1f0ad4cdfbaea95be0d600257d0818b7e3c1ea06eab60a2e6c6cb659864f8c1`；共 2 页，标题与正文一致，三位推荐审稿人条目未跨页拆分，旧标题与旧效应量均为零命中；
- 自动质控结果为 PASS：摘要 268 词、4 个 Word 原生表格、5 个嵌入图、27 条参考文献；正文引用按首次出现顺序完整覆盖 1–27，旧标题、旧 B1 效应量、旧单细胞样本数、旧受体斜率和 supplementary figure 字样均为零命中；
- 10 个最终 PNG/TIFF 文件均通过专用 figure validator：RGB、600 dpi、单帧；TIFF 为 LZW 压缩，无自动警告；
- 另附 `FAP_manuscript_V5.5_V5.7_V6.0_comparison_2026-08-14.md`，逐项记录三版分析口径、数值和叙事的继承/撤回关系。
