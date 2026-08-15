# FAP–SDC4/CD44 论文 V5.5、V5.7 与 V6.0 差异对比

日期：2026-08-14

## 一、结论先行

三版并不是单纯的文字递进，而是经历了两次实质性的分析口径改变。

- **V5.5** 是“FAP–基质状态”为主、衰老为补充解释的版本。它建立了跨 bulk、单细胞和蛋白组的 FAP–ECM 共变主线，但低估了成纤维细胞丰度对 bulk 相关性的影响，且单细胞受体模型使用了较早的归一化/自由度口径。
- **V5.7** 是方法学上最重要的一次校正。它采用全基因库 UMI 归一化、1,501 个肿瘤成纤维细胞、Satterthwaite 自由度，并加入 MCP-counter/EPIC 组成敏感性分析。其结果把结论收紧为“FAP 标记一个明显受成纤维细胞丰度驱动的基质状态”，同时明确否定 SDC4/CD44 与 FAP 的稳定共诱导。
- **V6.0** 将论文改成“区室化衰老”为首要卖点，但新增的 FAP+ 与 FAP− 比较混入了正常组织细胞，并受到 UMI 深度和 CAF 亚型组成的明显混杂。V6.0 的标题、摘要首要效应量及“FAP+ CAF 携带衰老程序”的结论不能作为最终口径。

因此，本次新稿没有简单沿用某一版，而是采取了以下组合：**保留 V5.7 的 bulk 组成敏感性和全库归一化单细胞模型；撤回 V6.0 的 FAP 特异性衰老主张；保留经肿瘤限定、匹配空模型和独立 GSE166555 队列支持的“成纤维细胞区室富集衰老相关转录”结论；并将其与 FAP–matrix 轴明确分开。**

## 二、核心差异总表

| 比较维度 | V5.5 | V5.7 | V6.0 | 本次重构后的处理 |
|---|---|---|---|---|
| 标题与首要叙事 | *Cross-platform characterization of a FAP-associated stromal matrix state and its senescence-related transcriptional covariation in colorectal cancer*；matrix-first | *A reproducible FAP-associated stromal matrix state with senescence-related transcriptional covariation in colorectal cancer*；仍为 matrix-first，但强调可重复性与证据边界 | *Compartmentalized senescence programs… FAP+ CAFs carry a senescence-associated matrix program…*；senescence-first | *Fibroblast-enriched senescence-associated transcription and a distinct FAP-linked matrix state in colorectal cancer*；以并列正向表述呈现两条可分离的轴 |
| bulk FAP13–matrix4 | TCGA ρ=0.930；GSE39582 ρ=0.913 | 相同边际相关，但明确加入多种组成校正 | 沿用 V5.5 的“稳健存在”叙述 | 保留边际相关，同时把组成依赖性作为核心结果而非次要限制 |
| 成纤维细胞丰度校正 | fib5 残差后再做 Spearman，partial ρ=0.608 | 标准秩残差 partial Spearman：fib5 0.650；MCP-counter 0.113；EPIC 0.288；联合 0.150 | 回退到 0.608，未整合 V5.7 的 MCP/EPIC 敏感性 | 冻结 V5.7 全套数值：0.650/0.113/0.288/0.150；结论为 proxy-dependent、显著受细胞丰度影响 |
| 单细胞 matrix4~FAP | 1,498 cells；slope=0.112，t=7.06，P<0.001 | 1,501 个肿瘤成纤维细胞；全库归一化；slope=0.253，95% CI 0.213–0.294，P=9.17×10⁻³³ | 1,498 cells；slope=0.119，t=12.48；并非 V5.7 所述冻结管线 | 采用 V5.7 的 1,501-cell 全库归一化结果 |
| 单细胞 receptor2~FAP | slope=+0.051，P=0.017，弱阳性且与患者级结果矛盾 | slope=−0.021，95% CI −0.062–0.019，P=0.302 | slope=−0.005，P=0.267，样本数仍写 1,498 | 采用 V5.7 的 −0.021、P=0.302；结论为无受体共诱导 |
| SenMayo 定义 | 来源 125 genes，称去除 4 个重叠后为 120；实际平台可用数未完全核清 | 明确区分全平台表示集及去重敏感性；TCGA 去重后实际 119 genes | 仍写 120 genes，并在单细胞新增分析中保留口径不一致 | GSE132465 冻结为 SenMayo119；GSE166555 可获得其中 114 genes；正文逐一说明 |
| 新增单细胞衰老比较 | 无 FAP+/- 主分析；衰老主要在 bulk 层 | 无“FAP+ CAF 衰老”首要结论 | FAP+ vs FAP−：β=0.083，P=6.96×10⁻³⁷；SASP β=0.131；作为摘要首项结果 | 肿瘤限定并校正 log UMI、亚型和患者后：SenMayo β=0.0238，95% CI −0.0051–0.0527，P=0.106；SASP β=0.0323，P=0.400；按阴性结果报告 |
| 细胞宇宙 | 肿瘤成纤维 1,501；肿瘤上皮 17,469 | 同左，定义清楚 | 3,462 成纤维+18,539 上皮，实际混合 tumour/normal，却按 CRC 肿瘤区室表述 | 主分析严格限定 tumour：1,501 fibroblast-lineage +17,469 epithelial cells |
| 区室对比 | 未作为主结果 | 未作为主结果 | SenMayo +0.376、23/23 同向；MKI67 校正后 β=+0.361 | 肿瘤限定后仍稳健：SenMayo 23/23 同向，标准化配对差中位数 0.466，95% CI 0.410–0.553，exact P=2.38×10⁻⁷；但只称“senescence-associated transcription” |
| 衰老标志基因 | bulk CDKN2A/B 高低组比较 | 仍作为 bulk 支持 | 以 FAP+ 细胞较高的 CDKN2A/B 检出率支持衰老 | 深度+亚型校正的 marker GLMM 均不显著；LMNB1 原始方向也不符合经典 lamin B1 loss，故不再作为 FAP+ 衰老证据 |
| 独立单细胞复现 | GSE166555 主要用于 tumour/normal 与 FAP–CLDN 检查，不是衰老区室复现 | 同样未完成独立衰老复现 | 方法中列入 GSE166555，但新衰老结论仍主要来自 GSE132465 | 新下载并冻结 GSE166555；10 名合格患者中 SenMayo 9/10、SASP 10/10 在 CAF-like 区室更高，MKI67 9/10 更低；FAP–SenMayo 不相关 |
| SASP–免疫通讯 | 以探索性推断为主 | 强调共表达不能建立通讯 | 新增 Figure S8，列出散在 nominal P，但称“无强共变” | 改为 14 个预设 cognate 配体–受体比较；minimum BH q=0.598，明确为阴性筛查 |
| 文章风险 | matrix 结果扎实，但组成控制与受体模型容易被审稿人质疑 | 方法最稳健，创新性相对保守 | 标题新颖，但首要证据存在决定性混杂，审稿风险最高 | 以可证伪、可复现的“双轴分离”取代过度机制化叙事 |

## 三、V6.0 首要结果为何需要撤回

### 1. 3,462/18,539 不是患者筛选阈值造成的差异

对原始 GSE132465 Seurat 元数据复核后：

- 3,462 个成纤维样细胞 = 1,501 tumour + 1,961 normal；
- 18,539 个上皮细胞 = 17,469 tumour + 1,070 normal；
- 1,025 个 FAP-detected 成纤维细胞 = 925 tumour + 100 normal；
- 2,437 个 FAP-undetected 成纤维细胞 = 576 tumour + 1,861 normal。

因此，FAP-detected 组有 90.2% 来自肿瘤，而 FAP-undetected 组有 76.4% 来自正常组织。V6.0 的 FAP+/FAP− 效应同时编码了 FAP 检出、肿瘤/正常来源和细胞状态，不能解释为 FAP 特异的 CAF 衰老。

### 2. UMI 深度是第二个关键混杂

FAP-detected 与 FAP-undetected 成纤维细胞的中位 UMI 约为 13,082 与 6,104。FAP 是否被检出和多基因 SenMayo/SASP 得分均高度依赖测序深度。

在 tumour-only 1,501 个成纤维细胞中：

- 未校正深度时，SenMayo 对 FAP 检出仍为阳性；
- 加入 log UMI 后，方向接近零；
- 同时校正 log UMI、CAF subtype 和患者随机截距后，SenMayo β=0.0238，P=0.106；
- 同一模型下 SASP β=0.0323，P=0.400。

这说明 V6.0 的 β=0.083、P=6.96×10⁻³⁷ 主要不能被视为独立的 FAP–衰老关联。

### 3. 单基因检测率不能挽救该结论

V6.0 报告 CDKN2A/CDKN2B 在 FAP+ 细胞中检出率更高，但在 tumour-only、log UMI+CAF subtype+patient 调整的 logistic GLMM 中：

- CDKN2A OR=0.89，P=0.529；
- CDKN2B OR=0.98，P=0.911；
- CDKN1A OR=0.96，P=0.788；
- LMNB1 OR=0.86，P=0.601，且模型 singular；
- MKI67 OR=0.84，P=0.614。

这些结果不支持“FAP+ 细胞具有独立的经典衰老标志组合”。此外，LMNB1 的未校正检出率在 FAP+ 组反而更高，不符合典型 lamin B1 loss 的方向。

## 四、V6.0 中仍可保留的发现

V6.0 的“区室差异”与“FAP 特异性”必须分开判断。肿瘤限定后，成纤维细胞相对上皮细胞的 SenMayo/SASP 富集仍然存在：

- GSE132465：SenMayo 23/23 患者同向，标准化配对差中位数 0.466，exact P=2.38×10⁻⁷；SASP 22/23 同向；MKI67 21/23 在成纤维细胞更低；
- 表达匹配 5,000 次空模型：SenMayo empirical P=0.00020；SASP empirical P=0.0336；
- GSE166555：10 名合格患者中 SenMayo 9/10 同向（exact P=0.0195）、SASP 10/10 同向（P=0.00195）、MKI67 9/10 反向（P=0.00391）。

这些证据支持“肿瘤成纤维细胞区室富集衰老相关转录特征”，但仍不能诊断细胞衰老。它可能包含谱系基线、慢周期状态、应激和 SASP-like 转录等成分，需以 p16/p21、LMNB1 loss、SA-β-gal、增殖排除和空间共定位等多标志验证。

## 五、三版中应冻结的最终统计口径

### 采用 V5.7

- TCGA FAP13–matrix4 的组成敏感性：fib5 0.650；MCP-counter 0.113；EPIC 0.288；MCP+EPIC 0.150。
- tumour fibroblast matrix4~FAP：slope=0.253，95% CI 0.213–0.294，P=9.17×10⁻³³。
- tumour fibroblast receptor2~FAP：slope=−0.021，95% CI −0.062–0.019，P=0.302。
- overlap-removed SenMayo–FAP13 的四种组成校正：0.388、0.378、0.575、0.384；SenMayo–matrix4 对应为 0.133、−0.098、0.082、−0.039。

### 不再采用

- V5.5：matrix4 slope=0.112；receptor2 slope=+0.051、P=0.017；fib5 partial ρ=0.608 作为唯一组成校正。
- V6.0：matrix4 slope=0.119；receptor2 slope=−0.005；FAP+ SenMayo β=0.083；FAP+ SASP β=0.131；混合 tumour/normal 的 +0.376/+0.361 作为标题级证据；SenMayo“120 genes”的表述。

## 六、叙事结构的演变

### V5.5：从 FAP–CLDN 假说退回 FAP–matrix

其价值在于识别了稳定的 FAP13–matrix4 共变，并明确 SDC4/CD44 更像候选上皮界面，而非与 FAP 同步诱导的成纤维程序。问题在于文章仍承载过多旧的 CLDN、空间、nomogram 和机制分支，核心贡献不够集中。

### V5.7：把“强相关”改写成“组成敏感的组织状态”

V5.7 的最大进步不是效应变大，而是承认 bulk 相关的强度依赖于成纤维细胞丰度代理，并用规范化的单细胞模型消除了受体弱阳性的内部矛盾。它是三版中统计口径最适合继承的一版。

### V6.0：创新性提升，但证据越过边界

V6.0 把单细胞结果提到标题和摘要首位，结构上更能回应“stimulate new thinking”。但 FAP+ CAF 衰老这一最关键桥梁并未在 tumour-only、深度校正后成立，而且“CRC 中尚未系统刻画”的新颖性表述也受到 2026 年 CRC senescent-CAF 文献的限制。问题不是措辞稍强，而是首要结论依赖错误分析宇宙。

### 本次新稿：两个相关但可分离的轴

新稿将证据重组为：

1. **区室轴**：CRC 肿瘤成纤维细胞相对上皮细胞富集 SenMayo/SASP-like 转录，并在独立队列复现；
2. **FAP–matrix 轴**：FAP 与 collagen/fibronectin 程序在肿瘤成纤维细胞内稳定耦联，并在 GSE166555 与 CPTAC 层得到支持；
3. **分离证据**：患者级 FAP–SenMayo、SenMayo–matrix4 在两个单细胞队列均不稳定，而 FAP–matrix4 稳定；
4. **阴性边界**：没有 FAP 特异性衰老、SDC4/CD44 共诱导、SASP–免疫通讯或独立预后效应的可靠证据。

这一结构把创新点从未经控制的 FAP+ senescent-CAF 因果链，转为可复现的区室转录分布与 FAP–matrix 轴分层；内部一致性、可复现性和可检验性均明显提高。

## 七、候选 V6.1 的位置

后续提供的候选 V6.1 在版式和章节衔接上做了局部整理，但统计内核仍属于 V6.0：它保留混合 tumour/normal 的 3,462/18,539 细胞宇宙、β=0.083 与 +0.376/+0.361 等标题级数字，并继续使用 SenMayo 120 genes 和 1,498-cell matrix/receptor 口径。配套 Cover Letter 也沿用上述结论。因此，V6.1 没有成为第四套可冻结结果，只能作为措辞、作者信息和投稿材料的候选片段来源。

本次重构没有逐段修改 V6.1，而采用先冻结科学问题和证据账本、再按提纲重写的方式。这样避免把不同版本的分母、模型和解释层级带入同一段文字，也减少了返修历史在正文中留下的防御性语气。

## 八、投稿层面的版本裁定

- **V5.5：不建议直接投稿。** 可作为 matrix-first 的内容库，但需要替换组成校正和单细胞模型。
- **V5.7：可作为统计结果母版。** 其方法学口径应优先于 V5.5/V6.0，但原稿仍过长、分支过多，创新表达不足。
- **V6.0：Reject–Rebuild。** 不能仅通过降低措辞强度修复；必须替换分析宇宙并撤回 FAP+ CAF 衰老首要主张。
- **本次新稿：建议作为下一轮投稿底稿。** 它已把可保留发现、阴性结果和机制边界统一到同一证据结构中；投稿前仍需补齐作者邮箱及公开代码/Zenodo DOI。

## 九、参考文献与文风变化

- V5.5/V5.7 的参考文献主要服务于 FAP–matrix、bulk 队列与方法支线；
- V6.0 加入衰老主线，但对“转录评分不等于衰老细胞”和 2026 年 CRC sCAF 机制文献的定位仍不充分；
- 本次重构将参考文献由 18 条增加至 27 条，新增 CRC 基质来源、SenNet/多标志判定、单细胞伪重复、bulk 去卷积和通讯推断边界等直接支持；
- Ge 与 Yang 两篇 2026 年 CRC–CAF–衰老论文不再只是文末条目，而在 Introduction 和 Discussion 中分别承担“领域已具备何种机制证据”和“本研究为何保持转录层定位”的作用。

最终文风采用“先报告稳健发现、段末一次性限定证据边界”的结构；方法纠错和版本历史留在审查报告，不反复写入正文。

## 十、核查依据

版本文本：

- `I06_AJCR_v5.5.docx` 及其逐段提取稿；
- `I04_AJCR_v5.7.docx` 及其逐段提取稿；
- `I08_AJCR_v6.0_candidate.docx` 及其逐段提取稿。
- `FAP_SDC4_CD44_论文_v6.1_修订版_2026-08-14.docx`、配套 Cover Letter 与审查报告的逐段提取稿。

本次复算：GSE132465 tumour-only population audit、深度/亚型调整混合模型、marker GLMM、表达匹配空模型、GSE166555 独立队列复现、bulk composition sensitivity 与 cognate SASP–immune screen。所有脚本、结果表和 session information 已纳入随稿复现包。
