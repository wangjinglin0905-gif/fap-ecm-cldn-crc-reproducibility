# FAP–senescence CRC 论文独立医学期刊编辑与结构审查

**审查对象：** Scientific Reports 实际投稿稿（I02）、AJCR v5.5/v5.6/v5.7（I06/I03/I04）、两封 Cover Letter（I01/I05）、既往质控记录（I07）及 AJCR v6.0 候选稿（I08）  
**v6.0 候选稿 SHA-256：** `1d66ea09b8b42fdfcb7d5b6cc78902de908dd5c9714352388c00593aeccf8554`  
**审查边界：** 本轮为只读的期刊编辑/结构审查；未复跑原始代码，未独立核验文献新颖性，未修改稿件或图片。因此，v6.0 新增精确数值在完成代码、源表和图表交叉核验前仍属于“候选结果”；若复核通过，其证据上限为 CRC 单队列单细胞关联证据（E2），不能升级为机制或治疗证据。

## 一、编辑结论

**当前 v6.0 的处理建议：Reject–Rebuild（不宜按现状投稿），但新的科学主线值得保留。**

新增的患者配对单细胞分析，确实比 Scientific Reports 实际投稿稿更直接地回答了一个生物学问题：衰老相关转录特征在 CRC 的细胞区室中如何分布，以及这种信号是否在 FAP 阳性成纤维细胞中富集。它把原稿中主要受组织组成影响的 bulk 相关性，推进到“同一谱系内部”和“同一患者跨区室”的比较，属于实质性增量，具备成为新主轴的潜力。

但 v6.0 尚未完成真正的深层重构，原因有四点：

1. 新主结果仍排在 Results 3.12、Figure 11；3.12 仅约 385/3,355 个结果部分词（约 11.5%），Discussion 4.1 仅约 202/2,217 个讨论部分词（约 9.1%）。正文主体仍由旧的 FAP–matrix、claudin、SDC4/CD44、pT1、MMR、淋巴结模型和多套探索分析占据。这是“追加新节”，不是“以衰老区室化重建证据顺序”。
2. v6.0 回退到了 v5.5 的若干已被 v5.6/v5.7 修正的口径，包括单一 fib5 校正、旧的单细胞归一化/细胞数分支、对 CPTAC 的 SASP 过度解释，以及过强的临床/治疗外推。
3. 新单细胞结果本身仍存在样本集定义、评分基因数、FAP 阳性定义、文库量/检出率偏倚、LMNB1 方向和“预设分析”真实性等 P1 问题。
4. 标题、摘要和结论把“衰老相关转录特征”升级为“CRC 的衰老架构/分工”，且错误地称该分工得到 bulk 和蛋白层支持。现有数据只能定位转录特征，不能诊断真正的细胞衰老，也不能证明区室间功能分工。

**一句话判断：** 新分析“部分回应”了编辑对 originality/new thinking 的担忧，但当前 v6.0 尚未同时解决 validity；若以现状投稿，审稿人很可能认为这是一个有趣但单队列、受细胞身份和技术检出影响的 signature 观察，并发现旧口径回退和内部数字矛盾。

## 二、对 Scientific Reports desk rejection 的闭环判断

编辑的核心表述是论文未呈现“sufficiently valid or original finding”并不足以“stimulate new thinking”。这不是单纯的语言问题，而是中心问题和证据单位的问题。

| 编辑关切 | Scientific Reports 实投稿的状态 | v6.0 新增内容 | 当前闭环状态 | 必须补齐 |
|---|---|---|---|---|
| Originality | 核心为 FAP 与 ECM 共变；直接成纤维校正后大幅衰减，且 FAP–ECM 在 CAF 生物学中预期性较强 | FAP+ 与 FAP− 的谱系内 SenMayo/SASP 比较；23 患者配对的间质–上皮反差 | **部分完成** | 以谱系内 FAP 相关状态为主发现；独立单细胞队列或至少严格患者级敏感性验证；当前文献空白需独立检索确认 |
| Validity | 实投稿已较诚实地报告组成依赖和阴性结果，但中心效应仍容易被视为组织组成现象 | 混合模型、患者配对、MKI67 调整提高了可解释性 | **未闭环** | 统一细胞集/归一化；控制 UMI/基因检出；解决 LMNB1 反向信号；报告效应量 CI 和患者一致性；避免细胞级超小 P 值主导叙事 |
| New thinking | “FAP 是 matrix tissue-state marker”更像方法学澄清，不足以改变 CRC 生物学理解 | “FAP 相关 CAF 状态与上皮增殖在区室间呈不同分布”可形成可检验模型 | **方向成立、措辞过强** | 写成“compartmentalized distribution of senescence-associated transcription”，而非已建立的 senescence architecture/division of labour |
| Mechanism/translation | SDC4/CD44、senolytics、pT1 风险模型均无直接验证 | B3 反而为阴性 | **不支持** | 把 receptor、免疫通讯、senolytic、临床决策降为边界或未来实验，不再与主轴并列 |

## 三、衰老区室化能否成为主轴

### 3.1 可以成为主轴，但主轴应重新命名

最有价值的观察不是“间质 SenMayo 高、上皮 MKI67 高”本身。不同细胞类型的基础转录程序天然不同，这一跨区室反差容易被审稿人视为细胞身份的预期结果。真正更具辨识度的证据是：

- 在同一成纤维谱系内，FAP+ 细胞的 SenMayo/SASP 分数高于 FAP− 细胞；
- 23 位患者的患者级 FAP–SenMayo 相关为正；
- 患者配对的间质–上皮反差为该谱系内结果提供组织区室背景；
- B3 未得到稳定 SASP–免疫受体轴，限定了机制边界。

因此，建议把中心命题从“CRC 存在衰老区室化分工”收紧为：

> **CRC 单细胞转录组中，FAP 阳性成纤维细胞富集衰老相关转录特征；这种 FAP 相关间质状态与上皮较高的增殖信号形成患者内区室反差。**

这是 E2 级观察性命题。所谓 “division of labour” 只能作为 Discussion 中的概念模型，并明确标为待空间和功能验证的假说。

### 3.2 当前标题过度承诺

现标题：

> *Compartmentalized senescence programs in colorectal cancer: FAP+ cancer-associated fibroblasts carry a senescence-associated matrix program while the epithelial compartment retains proliferative dominance*

主要问题：

- “senescence programs/architecture”容易被理解为已识别真正衰老细胞；数据仅为 signature。
- “carry a … matrix program”把 FAP、衰老和 matrix 写成同一被验证程序，但 bulk SenMayo–matrix 关系在更严格组成校正下大幅衰减；CPTAC也没有检测 SASP。
- “while”构成接近功能分工的因果/目的性暗示。

推荐标题：

> **Single-cell analysis identifies a FAP-associated senescence-related stromal state with epithelial proliferative predominance in colorectal cancer**

更保守的备选：

> **Compartmentalized senescence-associated transcription in colorectal cancer links FAP-positive fibroblasts to a matrix-rich stromal state**

上述标题仍需在新分析复核通过后使用；若没有独立单细胞验证，标题中应保留 “single-cell analysis” 以说明证据来源。

## 四、证据等级与措辞审计

| 论断 | 当前证据 | 当前最高等级 | 可用措辞 | 不可用措辞 |
|---|---|---:|---|---|
| FAP+ 成纤维细胞 SenMayo/SASP 较高 | GSE132465，23 患者；混合模型；候选数值尚未在本轮复算 | E0，复核后 E2 | “was associated with higher senescence-related transcript scores” | “senescent FAP+ CAFs were identified” |
| 间质 SenMayo/SASP 高于上皮，MKI67 方向相反 | 同一单细胞队列患者配对 | E0，复核后 E2 | “patient-paired compartment contrast” | “CRC exhibits a senescence architecture/division of labour” |
| FAP13–matrix4 bulk 共变 | 两 bulk 队列；单细胞/蛋白支持 matrix 背景；组成校正高度依赖代理 | E2 | “abundance-linked stromal matrix state” | “activation-state coupling independent of fibroblast abundance” |
| SenMayo 与 FAP13 组成校正后仍相关 | v5.7 多种代理下 partial ρ 约 0.378–0.575 | E2 | “composition-adjusted covariation” | “senescence drives the matrix phenotype” |
| SenMayo 与 matrix4 独立耦合 | v5.7 中多项校正后 −0.098 至 0.133，部分为零 | 不支持强论断 | “marginal association was largely composition-dependent” | “senescence-associated matrix mechanism” |
| CPTAC 支持 SASP–ECM | CPTAC 仅有 FAP–ECM 和 receptor 蛋白；无 SASP 蛋白测定 | E0 | 删除 | “protein-level support for SASP–ECM covariation” |
| 明确 SASP–免疫通讯轴 | B3 以弱/负相关为主，未建立多重校正阳性 | 不支持 | “screen did not support a coherent axis” | 任何免疫抑制机制结论 |
| SDC4/CD44 为功能性上皮接口 | 表达/推断层证据且多层阴性或不显著 | E1 假说 | “candidate spatial hypothesis” | “engaged interface”, “receptor-dependent mechanism” |
| Senolytic/临床决策价值 | 无干预、无前瞻验证；Cox 和模型多为阴性/未外部验证 | E0 | “future testable hypothesis” | “candidate therapeutic angle”, “risk stratification” |

### 必须删除或降级的过强句

1. I08 P0099：“This division of labour was consistent across bulk, single-cell, and proteomic layers…”——错误。区室分布只由 GSE132465 单细胞支持；bulk 无法定位区室，CPTAC未检测衰老/SASP。
2. I08 P0101：“Three independent single-cell analyses…”——三项分析均来自同一 GSE132465，应改为 “three complementary analyses of the same cohort”。
3. I08 P0093/P0105/P0113：“protein-level support for SASP–ECM covariation”——无对应蛋白数据，必须删除。
4. I08 P0093：“link … to a senescence-related mechanism”——观察性签名不能建立机制，应改为 “provide a senescence-related transcriptional context”。
5. I08 P0119：“CRC exhibits a compartmentalized senescence architecture”——改为 “senescence-associated transcript scores showed a compartmentalized distribution in one CRC single-cell cohort”。
6. I08 P0109/P0119 对 “druggable stromal state”“senolytic-informed strategies”的展开超出证据，应只保留一个具体、可证伪的未来实验。

## 五、P1：分析口径和内部一致性问题

### 5.1 v6.0 未采用 v5.7 的最终组成校正口径

v6.0 摘要、Methods 2.4、Results 3.10、Discussion 和 Cover Letter 逻辑仍以 fib5 partial ρ = 0.608 为主，并写成“survived fibroblast-content adjustment”。v5.7/Scientific Reports 实投稿的更完整口径为：

- fib5：partial ρ = 0.650；
- MCP-counter Fibroblasts：0.113，bootstrap CI −0.015 至 0.250；
- EPIC CAF：0.288；
- MCP-counter + EPIC：0.150。

这组结果支持“bulk 高相关主要追踪成纤维细胞丰度；剩余分量依赖代理”，而不支持“排除成纤维细胞含量后仍稳健”。v6.0 必须以 v5.7 为旧分析的事实基线，不能选择性回退到 fib5。

### 5.2 同一 GSE132465 出现至少四套细胞口径

| 位置 | 成纤维/间质细胞 | 上皮细胞 | FAP 分组 | 问题 |
|---|---:|---:|---:|---|
| I08 P0028 数据源 | 1,501 | 17,469 | — | 被定义为 tumour fibroblast-lineage / tumour epithelial |
| I08 P0073 / Table 3 | 1,498 | — | — | 与 Methods 的 1,501 不一致；缺少 3 个细胞的排除规则 |
| I08 P0075 / Table 4 | 5,933 | — | 1,369 / 4,564 | 远大于 1,501 和 3,462；“stromal”范围未定义 |
| I08 Methods 2.13 / Fig. 11 | 3,462 | 18,539 | 1,025 / 2,437 | 与 2.2 口径不同，未给出肿瘤/正常、样本、QC 和注释差异 |

必须给出一个样本流图或 population ledger，逐项说明：原始 63,689 细胞 → 患者/样本类型 → tumour-only 与否 → 细胞注释 → QC → 基因可用性 → 每项分析的最终细胞数。不能仅在不同 Methods 小节分别报数。

### 5.3 旧单细胞模型存在两个互不兼容的“最终口径”

- v5.7：1,501 细胞，full-library log1p(CP10K)，matrix4 slope = 0.253，receptor2 slope = −0.021。
- v6.0：1,498 细胞，matrix4 slope = 0.119，receptor2 slope = −0.005。

两者方向相同，但效应量相差明显，且 v6.0 Methods 2.5 没有完整写明全细胞 UMI 分母、过滤导致 1,498 的原因或自由度。作者必须指定唯一权威脚本、输入文件和模型输出；正文、Table 3、图、摘要及复现包全部从同一冻结结果表取值。

### 5.4 SenMayo 基因数不自洽

I08 多处称“125 genes; 120 after removal of four overlaps”。125−4=121；v5.7 经平台映射和 CXCL8/IL8 处理后报告 119。单细胞数据可能有不同的实际覆盖数，但必须列出：来源基因数、去重叠数、平台可检出数、最终评分数和具体缺失基因。Figure 11 图注也必须使用相同数字。

### 5.5 新单细胞结果存在技术偏倚替代解释

- FAP+ 定义为 `FAP > 0`，高度依赖 UMI 深度和 dropout；FAP+ 细胞可能只是总转录本/检出基因更多。
- CDKN2A/CDKN2B/LMNB1 的“检出率”同样受 library size 影响。应至少调整总 UMI、nFeature_RNA、样本/患者，并做患者级检出率比较。
- SenMayo/SASP 的 z-score universe 必须明确。若基因在各 compartment 内分别中心化，则 stromal–epithelial 分数不可直接比较；必须使用共同的基因标准化参考或可跨细胞类型比较的评分方案。
- 细胞级混合模型给出 P = 10⁻²³–10⁻³⁷，但 β = 0.083/0.131。主文应优先报告 β、95% CI、患者级效应和一致方向患者数，避免让细胞数驱动的极小 P 值制造生物学强效错觉。
- 建议增加连续 FAP、患者内分位数、患者伪 bulk/患者内配对差异，以及 cluster-robust 或患者 bootstrap 敏感性分析。

### 5.6 LMNB1 是必须正面报告的方向冲突

LMNB1 在 FAP+ 中的检出率为 6.4%，FAP− 为 1.6%，方向为更高，而经典衰老表型通常强调 LMNB1 丢失/降低。当前文字把它写成“low in both groups”，没有解释方向冲突。该结果不能作为衰老支持证据，反而应作为重要不一致：可能来自低检出率、文库量偏倚、细胞状态差异，或说明该群体只是 senescence-adjacent transcriptional state。还应报告 CDKN1A/p21（若可用）及其方向，不能只挑选支持标志物。

### 5.7 “预设/锁定分析计划”表述存在完整性风险

I08 P0021 称所有分析遵循 2026-07-27 锁定计划；但项目说明明确表示 B1–B3 是在 desk rejection 后新增。不能把新增分析追溯性地写成预设分析。建议区分：原 bulk/matrix 分析按既有计划执行；单细胞 compartment analysis 是由 bulk 结果驱动的 exploratory extension。无需在正文叙述拒稿历史，但必须准确说明分析是探索性的还是预设的。

## 六、真正的结构重构方案

### 6.1 摘要

当前摘要约 439 词，信息过载，并同时承载 senescence、matrix、receptor、CPTAC、nulls、临床预后。建议按以下分工压缩：

- **Background：** 只写一个问题——CRC 中衰老相关转录特征的细胞区室归属未知。
- **Methods：** GSE132465 为主分析；bulk/CPTAC 为支持层；明确患者级和混合模型。
- **Results：** 先报患者级 compartment contrast 和 FAP+–FAP− 结果，附效应量/CI；随后一句写 bulk FAP–senescence 和 matrix 背景；最后一句写 B3/receptor/prognosis 阴性边界。
- **Conclusion：** “a FAP-associated senescence-related stromal state” + “requires multi-marker validation”；删除机制、senolytic 和 clinical decision language。

### 6.2 Introduction

建议四段：

1. senescence 在 TME 中的作用及“signature ≠ senescent cell”边界；
2. FAP+ CAF/matrix 背景，以及 bulk 无法区分丰度与状态；
3. CRC 的具体空白：FAP 相关 CAF 内部是否富集衰老相关转录特征，及其与上皮增殖分布的患者内关系；
4. 主假说、次假说和分析路线。

当前 claudin、pT1、SDC4/CD44 和临床手术决策占用四个段落，显著稀释衰老主轴。claudin/receptor 应压成一个次要背景句，pT1 临床段建议删除或移至 Discussion 的限制/未来队列。

“We previously analyzed…”（I08 P0018）像修订史，应改成直接的科学陈述，不在稿件中叙述旧路线。

### 6.3 Results

推荐顺序：

1. **3.1 数据集、细胞流和主要分析设计**；
2. **3.2 患者配对的 stromal–epithelial SenMayo/SASP/MKI67 分布**；
3. **3.3 FAP+ vs FAP− 成纤维细胞内关联及患者级/技术敏感性**；
4. **3.4 bulk FAP–senescence 复制及 v5.7 组成校正**；
5. **3.5 FAP–matrix 背景的单细胞和蛋白支持**；
6. **3.6 证据边界：SASP–immune、receptor、survival 均未建立**。

CMS、MMR、claudin 成员、pT1 梯度、淋巴结 nomogram、NicheNet/PROGENy/DoRothEA 和大部分空间 receptor 分析移入补充材料。否则读者会看到一篇旧 FAP–matrix/claudin 论文，末尾附加一段 senescence。

### 6.4 Discussion

建议五节：

1. 精确解释主要 E2 观察；
2. 与其他瘤种 senescent fibroblast 文献的异同，明确为跨癌种 E1 外推；
3. 替代解释：细胞身份、RNA 深度、FAP dropout、LMNB1 方向、单队列；
4. matrix 组成依赖及 receptor/immune/prognosis 阴性边界；
5. 一个最关键验证实验：FAP + p16/p21/LMNB1 多重蛋白并结合空间位置和增殖标记。

“division of labour”可在 4.1 作为工作模型出现一次，不应用于标题、摘要结果和结论的事实陈述。

### 6.5 Conclusion

结论只承担三件事：主要观察、证据边界、下一项验证。当前结论重新汇总 claudin、receptor、T-stage、immune、survival 和三类治疗策略，恢复了旧稿的多主线问题。建议压缩为 3–4 句，不出现 “therapeutic angle”。

## 七、阴性结果如何保留

阴性结果是本稿可信度的重要组成，但应服务于界定主发现：

- **保留在主文：** B3 无稳定 SASP–immune 轴；matrix 相关在成纤维校正后高度衰减且代理依赖；receptor 未共诱导；Cox 阴性。
- **必须正面解释：** LMNB1 方向与经典衰老预期不一致。
- **移至补充：** 详细 claudin、MMR、nomogram、NicheNet/PROGENy/DoRothEA、逐层 receptor 探索。
- **写法：** “did not support a coherent axis” 或 “defined the boundary of the observed state”，不要写成“虽阴性但仍提示机制”。

## 八、图表与编号审计

### 8.1 叙事编号虽顺延，但主图优先级错误

Figure 11 是新主结果，却排在十张旧图之后。建议把它拆分/重组为新 Figure 1–2：

- Figure 1：细胞流、共同评分框架、23 患者配对 stromal–epithelial SenMayo/SASP/MKI67；
- Figure 2：FAP+ vs FAP− 谱系内效应、患者级 FAP–SenMayo、CDKN2A/B/LMNB1 方向与敏感性；
- Figure 3：bulk replication + MCP/EPIC/fib5 组成校正；
- Figure 4：matrix 的单细胞/蛋白支持与 receptor 阴性边界；
- 其余旧图移至补充材料。

Figure 11c/d 的横轴必须由 `FALSE/TRUE` 改为 `FAP−/FAP+`；主图应显示患者级点/连线或患者级效应，而不只呈现细胞级小提琴图。

Supplementary Figure S8 若继续标 nominal asterisks，应同时显示 BH/FDR 结论；否则零散星号会把阴性筛查视觉化为阳性发现。

### 8.2 候选稿仍有机器可见的表格卫生问题

- 七张表的提取结果均出现 `Field` 占位首行，与“28 项校验通过”的描述不一致。
- Table 7 重新出现既往质控记录明确要求删除的 “Pathological TSR: 34/34 node-positive” 行，当前稿件没有足够来源和分析闭环支持该行。
- v6.0 没有纳入 v5.7 的 Table 8 多代理成纤维丰度敏感性结果。
- Table 3 仍为 1,498 细胞，而 Methods 2.2 为 1,501；Figure 11 为另一套 3,462 细胞。

上述问题均需在 DOCX 重新构建、引用冻结和逐页渲染前解决。

## 九、Cover Letter、声明和投稿风险

### 9.1 Cover Letter

Scientific Reports Cover Letter 主要强调技术有效性和诚实阴性结果，但没有给编辑一个改变 CRC 生物学理解的中心句，这与 desk rejection 一致。AJCR v5.5 Cover Letter 又依赖 partial ρ = 0.608，并称跨多层“independent validation”，已不适合作为 v6.0 基础。

新的 Cover Letter 应：

1. 首句提出 CRC 中 compartment-level senescence-associated transcription 的空白；
2. 用 23 患者配对结果和 FAP+ 谱系内结果说明贡献；
3. 用 v5.7 的多代理组成校正说明严谨性；
4. 一句明确不是 senescent-cell diagnosis、immune mechanism 或 prognostic biomarker；
5. 不写 “first” 或 “unresolved” 的绝对新颖性，除非完成截至投稿日的系统检索。

三位建议审稿人的专业方向总体匹配，但当前文件只有 `[email]` 占位。投稿前必须核验现任机构邮箱、近期共同作者/同机构/基金/导师关系及期刊规定的利益冲突窗口。Cover Letter 的签署人和投稿系统通讯作者也应一致。

### 9.2 Data/code availability

I08 Methods 2.14 回退为 `[GitHub repository to be added]` 和 `[Zenodo DOI to be added]`，而 v5.7 已有 release v1.3.0 和 Zenodo DOI。新分析加入后不能机械沿用旧 release；应建立包含 Path B 输入清单、脚本、环境、患者级源表、模型输出及 Figure 11/S8 源数据的新不可变 release/Zenodo 版本，再回填真实链接。

### 9.3 AI disclosure

当前声明仅写 “AI-assisted tools”，缺少工具/模型及版本，不满足工具、用途、人工核验三要素；且本项目的实际用途可能超过单纯语言编辑。必须据实列出参与写作、代码审阅、分析解释和图表制作的工具及版本，并由作者确认最终范围。

### 9.4 Ethics/funding/authorship

公共去标识数据声明基本可保留，但建议说明原始研究具有各自伦理/同意程序，当前研究未新增受试者或样本。基金号、作者顺序、CRediT 和通讯作者信息必须由全体作者逐项确认，不能从旧稿自动继承。

## 十、优先级矩阵

### P1：完成前不得投稿

1. 以 v5.7 为旧分析事实基线，冻结唯一数值表；复核并统一 0.608/0.650/0.113/0.288/0.150 和两套混合模型结果。
2. 建立 GSE132465 population ledger，解释 1,498、1,501、3,462、5,933、17,469、18,539 及两套 FAP+ 数量。
3. 复核 SenMayo 最终基因数、共同标准化 universe、full-library 分母、FAP dropout/UMI 偏倚，并增加患者级敏感性和 95% CI。
4. 正面报告 LMNB1 反向结果；删除 bona fide senescence、mechanism、protein SASP、division-of-labour 事实化和 senolytic/临床外推。
5. 把新增分析标为探索性扩展，不能追溯性宣称全部预设。
6. 以新单细胞结果重排 Results/Figures；旧的 claudin/receptor/clinical 分支降为补充。
7. 删除 `Field` 占位和无来源的 34/34 行；回填含新分析的 GitHub/Zenodo；完善 AI 声明。

### P2：投稿前应完成

1. 将标题改为 senescence-associated transcription/state，而非 senescence architecture。
2. 压缩约 439 词摘要，并按目标期刊最新格式/字数复核。
3. Figure 11 改为 Figure 1/2，改正 FAP 分组标签，主视觉转向患者级效应；S8 显示 FDR 边界。
4. 进行截至投稿日的 CRC senescence/CAF 单细胞文献检索，验证“空白点”，补直接 CRC 文献；跨癌种文献明确为 E1。
5. 新 Cover Letter 以生物学问题和患者级单细胞结果开头，删除单一 fib5“稳健”叙述。

### P3：语言与版式

1. 统一 tumor/tumour、program/programme 等英美拼写。
2. 删除 “We previously analyzed”“three independent analyses”等修订史或证据膨胀措辞。
3. 事实冻结后再统一交叉引用、参考文献、表图编号，并重新逐页渲染检查。

## 十一、最终投稿判断

- **科学方向：** 值得继续；新增单细胞结果可形成比原 Scientific Reports 稿更清晰的最小可发表故事。
- **原创性：** 有潜力，当前为中等，尚未独立文献核验；仅靠跨区室 SenMayo/MKI67 反差不足，必须突出并稳固 FAP+ 成纤维谱系内结果。
- **统计/生信可信度：** 当前候选稿为 Low–Moderate；主要受口径分叉、细胞数不一致和技术偏倚未闭环影响。
- **当前稿件结论：** Reject–Rebuild。
- **完成 P1 后的目标定位：** 一篇以患者配对单细胞关联为核心、bulk/蛋白作为背景支持、明确不作机制和治疗宣称的 AJCR 原创研究稿；能否达到目标期刊标准仍取决于分析复算、文献新颖性核验和最终图文质控。

**作者最终责任：** 所有精确数字、分析预设状态、基金/作者/声明和复现链接均须由作者基于原始代码、数据和期刊现行要求确认。
