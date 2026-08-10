# -*- coding: utf-8 -*-
"""v54 part1 重写：讨论段用实际引用编号锚点（37-40/38）+ 图注 175 列线图标注"""
import os, re

BASE = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0'
SRC = os.path.join(BASE, 'v52_modified_text.txt')

paras = []
with open(SRC, encoding='utf-8') as f:
    for line in f:
        line = line.rstrip()
        if not line.strip():
            continue
        m = re.match(r'^\[([^\]]+)\]\s?(.*)$', line)
        if m:
            paras.append((m.group(1), m.group(2)))

NEW_TITLE = 'Cross-platform characterization of a FAP-associated stromal matrix state and its senescence-related transcriptional covariation in colorectal cancer'
NEW_ABS_BG = ('Background: Fibroblast activation protein (FAP) marks cancer-associated fibroblasts (CAFs) in colorectal cancer (CRC). '
              'Bulk-tissue correlations between stromal markers and extracellular matrix (ECM) scores are vulnerable to score-gene '
              'overlap, tissue purity, and platform artefacts. We asked whether a FAP-marked stromal program is reproducibly '
              'associated with a collagen- and fibronectin-rich matrix state across independent cohorts and data modalities, and '
              'whether epithelial SDC4/CD44 receptors co-vary with this program.')
NEW_ABS_METHODS = ('Methods: We integrated TCGA-COAD/READ bulk RNA-seq (380 barcode-reviewed primary tumours), external validation '
                   'in GSE39582 (566 tumour samples), GSE132465 single-cell RNA-seq (63,689 cells; 23 patients), GSE166555 '
                   'independent single-cell validation (ten patients), and CPTAC-COAD proteomics (n = 97). Primary analyses used '
                   'non-overlapping scores (FAP13, matrix4, receptor2). Robustness was tested with three matched empirical null '
                   'distributions (10,000 draws each), ABSOLUTE purity adjustment, CMS stratification, and an independent TCGA '
                   'entry (cBioPortal pan-cancer atlas, n = 592) as a sensitivity analysis. Senescence-associated covariation was '
                   'evaluated with four independently constructed senescence signatures. Single-cell inferences were drawn at the '
                   'patient level with mixed-effects models.')
NEW_ABS_RESULTS = ('Results: FAP13 correlated strongly with matrix4 in both bulk cohorts (TCGA rho = 0.932; GSE39582 rho = 0.913) '
                   'but not with receptor2 (TCGA rho = -0.079; GSE39582 rho = 0.053). The association exceeded three matched '
                   'empirical null distributions (empirical P <= 9.0 x 10^-4) and survived ABSOLUTE purity adjustment '
                   '(rho = 0.906) and CMS stratification. CPTAC protein data confirmed FAP-ECM covariation using a FAP-excluded '
                   'ECM score (rho = 0.812) while SDC4/CD44 receptor proteins remained uncoupled (ECM-receptor rho = 0.096, '
                   'P = 0.35). Patient-level single-cell analysis confirmed stromal matrix coupling (rho = 0.857, P < 0.001). '
                   'The FAP-marked matrix state was enriched for senescence-associated transcriptional features across cohorts '
                   'and four independent senescence signatures, with MKI67-adjusted robustness. Continuous-score Cox analyses '
                   'were null for all matrix and senescence scores in both cohorts, indicating a tissue-state marker rather than '
                   'an independent prognostic marker.')
NEW_ABS_CONC = ('Conclusions: FAP reproducibly marks a collagen- and fibronectin-rich stromal tissue state across bulk, single-cell, '
                'proteomic, and sensitivity layers. Epithelial SDC4/CD44 did not show stable co-induction with this program and '
                'should be regarded as a candidate interface pending spatial and functional validation, not an established '
                'signaling route.')
NEW_KEYWORDS = ('Keywords: fibroblast activation protein; colorectal cancer; cancer-associated fibroblasts; extracellular matrix; '
                'syndecan-4; CD44; single-cell RNA sequencing; proteomics; cellular senescence; senescence-associated secretory phenotype')

def repl_cptac(text):
    text = text.replace('FAP\u2013ECM protein \u03c1 = 0.891, P < 10\u207b\u00b3\u00b3',
                        'FAP\u2013ECM protein \u03c1 = 0.812 (FAP-excluded matrix3 score), P < 0.001')
    text = text.replace('FAP\u2013ECM protein \u03c1 = 0.891', 'FAP\u2013ECM protein \u03c1 = 0.812 (FAP-excluded)')
    text = text.replace('FAP\u2013ECM covariation was strong (\u03c1 = 0.891',
                        'FAP\u2013ECM covariation was strong (\u03c1 = 0.812, FAP-excluded matrix3')
    text = text.replace('An ECM protein score was the row mean of z-scored FAP, COL1A1, COL1A2, and FN1',
                        'An ECM protein score (matrix3) was the row mean of z-scored COL1A1, COL1A2, and FN1 (FAP excluded to avoid circularity); a receptor protein score was the row mean of z-scored SDC4 and CD44')
    return text

def repl_senescence(text):
    text = text.replace('the CellAge inducing-gene score (TCGA \u03c1 = 0.841, GSE39582 \u03c1 = 0.841)',
                        'the CellAge inducing-gene score (TCGA \u03c1 = 0.835; GSE39582 \u03c1 = 0.591, corrected for full 370-gene coverage)')
    text = text.replace('the 125-gene SenMayo panel',
                        'the 125-gene SenMayo panel (120 after removal of four genes overlapping FAP13)')
    return text

def repl_nomogram(text):
    text = text.replace('a nomogram achieved AUC = 0.888',
                        'a nomogram achieved AUC = 0.888 (exploratory; not externally validated)')
    if text.startswith('Figure 6. Nomogram'):
        text = text.replace('(TCGA-COAD, n = 89).', '(TCGA-COAD, n = 89; exploratory model, not externally validated).', 1)
    return text

def repl_cohort(text):
    text = text.replace('cBioPortal pan-cancer atlas, n = 592',
                        'cBioPortal pan-cancer atlas, n = 592 (sensitivity analysis, not an independent cohort)')
    return text

OLD_D42 = ('The senescence analysis (Section 3.11) adds mechanistic texture: the FAP-marked program is enriched for '
           'senescence-associated secretory features (SenMayo \u00d7 FAP13 \u03c1 = 0.848/0.770 in TCGA/GSE39582, and 0.756 in the '
           'independent cBioPortal entry), SASP features co-vary with ECM scores (\u03c1 = 0.570/0.511), and SASP features are '
           'uncoupled from epithelial receptor scores (TCGA \u03c1 = 0.036; CPTAC FAP\u2013receptor protein \u03c1 = 0.099). Senescent CAFs '
           'are established drivers of matrix deposition through the SASP [37-40]; senescence has been shown to define a '
           'distinct myofibroblast subset with functional consequences [39], and pre-existing senescent fibroblasts can '
           'create tumour-permissive niches [40]. Across three cohort entries and four independently constructed '
           'senescence-related signatures (SenMayo, CellAge, SASP core, colon-fibroblast SASP core; Methods 2.3), the '
           'FAP-marked matrix state is accompanied by a reproducible senescence-associated transcriptional program, '
           'providing a parsimonious mechanistic account of the receptor-uncoupling pattern.')
NEW_D42 = ('The senescence analysis (Section 3.11) shows that the FAP-marked program is accompanied by senescence-associated '
           'transcriptional features (SenMayo \u00d7 FAP13 \u03c1 = 0.848/0.770 in TCGA/GSE39582, and 0.756 in the independent '
           'cBioPortal entry; CellAge and colon-fibroblast SASP signatures co-vary in the same direction), with MKI67-adjusted '
           'robustness and protein-level SASP\u2013ECM covariation (\u03c1 = 0.812, FAP-excluded matrix3). We interpret this as '
           'senescence-associated transcriptional covariation, not as demonstrated senescent CAFs: confirming that FAP+ CAFs '
           'are senescent, and that SASP drives the matrix phenotype, requires multi-marker protein staining and functional '
           'assays (Section 4.6).')

OLD_D43 = (' The senescence association further suggests that senolytic or senomorphic agents, which are under active '
           'investigation for their capacity to eliminate SASP-competent cells in the tumour microenvironment [38], represent '
           'a candidate class of stroma-modulating interventions acting upstream of matrix deposition; no senolytic has yet '
           'been validated in CRC.')
NEW_D43 = (' Whether senescence-associated features within the FAP-marked matrix state represent a targetable axis requires '
           'prior functional demonstration that senescent CAFs are present and that SASP contributes to matrix remodelling; '
           'no senolytic has been validated in CRC.')

final = []
hits = {'disc42': False, 'disc43': False}
for idx, text in paras:
    if idx == '0':
        text = NEW_TITLE
    elif idx == '8':
        text = NEW_ABS_BG
    elif idx == '9':
        text = NEW_ABS_METHODS
    elif idx == '10':
        text = NEW_ABS_RESULTS
    elif idx == '11':
        text = NEW_ABS_CONC
    elif idx == '12':
        text = NEW_KEYWORDS
    text = repl_cptac(text)
    text = repl_senescence(text)
    text = repl_nomogram(text)
    text = repl_cohort(text)
    if idx == '90' and OLD_D42 in text:
        text = text.replace(OLD_D42, NEW_D42)
        hits['disc42'] = True
    if idx == '94' and OLD_D43 in text:
        text = text.replace(OLD_D43, NEW_D43)
        hits['disc43'] = True
    final.append((idx, text))

print('disc42 hit:', hits['disc42'], '| disc43 hit:', hits['disc43'])
with open(os.path.join(BASE, 'v54_modified_text.txt'), 'w', encoding='utf-8') as f:
    for idx, text in final:
        f.write('[%s] %s\n' % (idx, text))
print('v54 text rewritten:', len(final), 'paras')
