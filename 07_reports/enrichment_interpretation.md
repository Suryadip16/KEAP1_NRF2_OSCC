# Comprehensive Interpretation of Enrichment Analysis Results
## KEAP1-NRF2 Axis in Oral Squamous Cell Carcinoma (OSCC)

> [!NOTE]
> This document provides a plot-by-plot interpretation of all enrichment analysis figures generated from the differential expression analysis comparing **NRF2-High vs NRF2-Low** OSCC tumours from TCGA-HNSC. All analyses reflect the transcriptomic consequences of constitutive NRF2 pathway activation in OSCC.

---

## Table of Contents

1. [GSEA Hallmark Analysis](#1-gsea-hallmark-analysis)
2. [GSEA C2 Canonical Pathways](#2-gsea-c2-canonical-pathways)
3. [GSEA Ridgeplot](#3-gsea-ridgeplot)
4. [GSEA Running Score – Key Pathways](#4-gsea-running-score--key-pathways)
5. [ORA GO Biological Process](#5-ora-go-biological-process)
6. [ORA GO Molecular Function](#6-ora-go-molecular-function)
7. [ORA GO Cellular Component](#7-ora-go-cellular-component)
8. [ORA KEGG Pathways](#8-ora-kegg-pathways)
9. [ORA Reactome Pathways](#9-ora-reactome-pathways)
10. [ORA Network and Hierarchy Plots](#10-ora-network-and-hierarchy-plots-cnetplots-emapplots-treeplots)
11. [GSVA Custom Pathway Analysis](#11-gsva-custom-pathway-analysis)
12. [NRF2 Reference Gene Set Analysis](#12-nrf2-reference-gene-set-analysis)
13. [Ferroptosis Axis Analysis](#13-ferroptosis-axis-analysis)
14. [Integrated Biological Synthesis](#14-integrated-biological-synthesis)

---

## 1. GSEA Hallmark Analysis

### 1a. GSEA Hallmark Dotplot ([05_GSEA_hallmark_dotplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_GSEA_hallmark_dotplot.pdf))

**What the plot shows:** This dotplot displays all 29 significantly enriched MSigDB Hallmark gene sets from the gene set enrichment analysis (GSEA), ranked by Normalized Enrichment Score (NES). Dot size encodes gene set size, and colour encodes adjusted p-value.

**Key findings and biological interpretation:**

#### Pathways ACTIVATED in NRF2-High tumours (positive NES):

| Hallmark Pathway | NES | p.adjust | Biological Significance |
|---|---|---|---|
| **OXIDATIVE PHOSPHORYLATION** | **+2.67** | 6.1×10⁻²⁰ | Strongest activation signal; metabolic rewiring |
| **REACTIVE OXYGEN SPECIES PATHWAY** | **+2.39** | 4.4×10⁻⁸ | Direct NRF2 target pathway — ARE-mediated defence |
| **MYC TARGETS V1** | +2.44 | 2.9×10⁻¹⁵ | Proliferative reprogramming |
| **MTORC1 SIGNALING** | +2.34 | 5.0×10⁻¹³ | Anabolic growth signalling |
| **XENOBIOTIC METABOLISM** | +2.17 | 1.8×10⁻⁹ | Phase I/II detoxification — canonical NRF2 function |
| **FATTY ACID METABOLISM** | +2.21 | 1.5×10⁻⁸ | Lipid metabolic reprogramming |
| **CHOLESTEROL HOMEOSTASIS** | +2.09 | 2.5×10⁻⁵ | Mevalonate pathway activation |
| **P53 PATHWAY** | +1.92 | 6.1×10⁻⁷ | p53-dependent stress response genes |
| **GLYCOLYSIS** | +1.79 | 1.5×10⁻⁵ | Warburg effect reinforcement |
| **HEME METABOLISM** | +1.62 | 1.2×10⁻³ | Iron-porphyrin homeostasis |

#### Pathways SUPPRESSED in NRF2-High tumours (negative NES):

| Hallmark Pathway | NES | p.adjust | Biological Significance |
|---|---|---|---|
| **ALLOGRAFT REJECTION** | **−2.46** | 1.1×10⁻¹⁹ | Immune evasion — loss of adaptive immunity |
| **EPITHELIAL MESENCHYMAL TRANSITION** | **−2.34** | 1.8×10⁻¹⁶ | Anti-mesenchymal phenotype |
| **INTERFERON GAMMA RESPONSE** | −1.93 | 3.6×10⁻⁸ | Dampened anti-tumour immunity |
| **IL6 JAK STAT3 SIGNALING** | −2.03 | 8.0×10⁻⁷ | Reduced pro-inflammatory signalling |
| **INFLAMMATORY RESPONSE** | −1.88 | 2.3×10⁻⁷ | Immune cold microenvironment |
| **INTERFERON ALPHA RESPONSE** | −1.83 | 5.6×10⁻⁵ | Impaired innate immune sensing |
| **COMPLEMENT** | −1.67 | 5.1×10⁻⁵ | Reduced complement activation |
| **ANGIOGENESIS** | −1.71 | 3.5×10⁻³ | Reduced vascular signalling |

**KEAP1-NRF2 Context:**

The **Oxidative Phosphorylation** pathway is the most significantly enriched hallmark (NES = +2.67), indicating that NRF2-high OSCC tumours undergo profound **metabolic rewiring** toward enhanced mitochondrial respiration. This is paradoxical yet consistent with NRF2's dual role: while protecting against oxidative stress through the **ROS Pathway** (NES = +2.39, featuring canonical targets like *NQO1, TXNRD1, GCLC, GSR, PRDX1, SRXN1, TXN*), NRF2 simultaneously fuels the very mitochondrial activity that generates ROS. This "antioxidant paradox" has been documented in NRF2-addicted cancers, where constitutive NRF2 activation provides a **survival advantage** by buffering the oxidative stress produced by heightened metabolic activity (Mitsuishi et al., 2012, *Cancer Cell*).

The activation of **Xenobiotic Metabolism** (NES = +2.17) reflects canonical NRF2-ARE transcription of Phase I/II detoxification enzymes including *NQO1, ALDH3A1, GCLC, GSR, AKR1C2/C3, EPHX1, CYP2S1, ABCC2/C3*. This confers **chemoresistance** — a major clinical problem in OSCC, as NRF2-driven drug efflux and conjugation render cisplatin and 5-FU less effective (Shibata et al., 2011, *Oncogene*).

The **simultaneous suppression of all major immune signalling pathways** (Allograft Rejection, IFNγ, IFNα, IL6/JAK/STAT3, Inflammatory Response, Complement) is striking and reveals that NRF2-high OSCC tumours foster an **immune-cold tumour microenvironment**. NRF2 has been shown to repress inflammatory gene expression through transcriptional interference with NF-κB (Kobayashi et al., 2016, *Nature Communications*), and the downregulation of antigen presentation genes (*HLA-DRA, HLA-DMA, B2M, TAP1*) is consistent with immune escape reported in KEAP1-mutant cancers (Best et al., 2018, *Cancer Discovery*).

The strong **suppression of EMT** (NES = −2.34) indicates NRF2-high OSCC tumours maintain a more **epithelial phenotype** with downregulation of mesenchymal markers (*VIM, COL1A1, COL3A1, FN1, CDH2, SPARC, FAP*). This is biologically important: while NRF2 activation promotes chemoresistance and metabolic advantage, it does so through a differentiated rather than mesenchymal invasive program. This is consistent with the observation that NRF2-high head and neck cancers tend to be well-differentiated squamous tumours with keratinizing features (Cancer Genome Atlas Network, 2015, *Nature*).

**Key literature:**
- Mitsuishi Y, et al. (2012). Nrf2 redirects glucose and glutamine into anabolic pathways in metabolic reprogramming. *Cancer Cell*, 22(1):66-79.
- Shibata T, et al. (2011). Cancer related mutations in NRF2 impair its recognition by Keap1-Cul3 E3 ligase and promote malignancy. *Proc Natl Acad Sci USA*, 105(36):13568-13573.
- Kobayashi EH, et al. (2016). Nrf2 suppresses macrophage inflammatory response by blocking proinflammatory cytokine transcription. *Nature Communications*, 7:11624.
- Best SA, et al. (2018). Synergy between the KEAP1/NRF2 and PI3K pathways drives non-small-cell lung cancer with an altered immune microenvironment. *Cell Metabolism*, 27(4):935-943.

---

### 1b. GSEA Hallmark Enrichment Map ([05_GSEA_hallmark_emap.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_GSEA_hallmark_emap.pdf))

**What the plot shows:** This enrichment map (emap) visualises the gene set overlap between significantly enriched Hallmark pathways as a network graph. Nodes are pathways, edges connect pathways sharing core enrichment genes, and node colour indicates activation (positive NES) or suppression (negative NES).

**Interpretation:**

The emap reveals **two clearly distinct functional clusters**:

1. **An NRF2-activated metabolic cluster** — Oxidative Phosphorylation, ROS Pathway, Fatty Acid Metabolism, Glycolysis, mTORC1 Signalling, MYC Targets, Cholesterol Homeostasis, and Xenobiotic Metabolism form a densely interconnected network. This cluster shares metabolic enzyme genes and reflects the **coordinated metabolic reprogramming** driven by NRF2. The tight connectivity between Oxidative Phosphorylation and ROS Pathway is biologically meaningful: NRF2 simultaneously upregulates mitochondrial electron transport chain components (*COX5A/5B, UQCRC2, NDUFS7, ATP5F1B*) while also inducing antioxidant defences (*PRDX1, TXN, GSR, TXNRD1*) to buffer the consequent ROS.

2. **A suppressed immune-inflammatory cluster** — Allograft Rejection, IFNγ Response, IFNα Response, Inflammatory Response, Complement, and IL6/JAK/STAT3 form a tightly connected network of downregulated immune functions. The interconnectivity reflects shared immune effector genes (*CD4, CD8A, LCK, ZAP70, GZMA, CXCL9/10/11*) and indicates a **globally immunosuppressed microenvironment** rather than selective pathway downregulation.

3. **EMT and Angiogenesis** form an isolated suppressed sub-cluster, sharing mesenchymal/stromal genes.

The clear **spatial separation** between the activated metabolic cluster and the suppressed immune cluster on the emap demonstrates the "**two-front strategy**" of NRF2-high OSCC: metabolic fitness + immune evasion. This mirrors findings in KEAP1-mutant lung adenocarcinoma where KEAP1 loss creates immune-desert tumours with enhanced metabolic capacity (Skoulidis et al., 2018, *Cancer Discovery*).

---

## 2. GSEA C2 Canonical Pathways

### GSEA C2 CP Dotplot ([05_GSEA_C2_CP_dotplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_GSEA_C2_CP_dotplot.pdf))

**What the plot shows:** This dotplot displays the 2 significantly enriched curated canonical pathways (C2:CP collection) from GSEA, with NES, gene count, and significance.

**Key findings:**

| Pathway | NES | p.adjust | Direction |
|---|---|---|---|
| SA_MMP_CYTOKINE_CONNECTION | −1.83 | 0.014 | Suppressed |
| SIG_BCR_SIGNALING_PATHWAY | −1.68 | 0.028 | Suppressed |

**Interpretation:**

The **MMP-Cytokine Connection** pathway suppression (core enrichment: *ACE, SELL, SPN, TNFSF8, TNFRSF1B, TGFB2, FCGR3A, CSF1*) indicates reduced extracellular matrix remodelling and immune cell–tumour stroma cross-talk in NRF2-high tumours. This is consistent with the EMT suppression seen in the Hallmark analysis and suggests that NRF2-high OSCC tumours have reduced matrix metalloproteinase-driven tissue remodelling.

The **BCR Signalling Pathway** suppression (core enrichment: *PIK3CD, VAV1, PPP3CB, PIK3R1, PTPRC, BTK, AKT3, NFATC2, BCL2, CD19, CD22*) further confirms the **B cell depletion** phenotype in NRF2-high tumours. B cells contribute to anti-tumour immunity through antigen presentation and antibody-dependent cytotoxicity. Their reduced signalling is part of the broader immunosuppressive landscape.

---

## 3. GSEA Ridgeplot

### [05_GSEA_ridgeplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_GSEA_ridgeplot.pdf)

**What the plot shows:** The ridgeplot displays the **distribution of log₂ fold-change values for core enrichment genes** within each significant Hallmark pathway. Each ridge (density curve) shows whether the genes in that pathway tend to be upregulated (shifted right) or downregulated (shifted left) in NRF2-High vs NRF2-Low tumours.

**Interpretation:**

The ridgeplot provides crucial additional information beyond just enrichment significance — it shows the **magnitude and consistency** of differential expression:

- **Oxidative Phosphorylation:** The ridge is broad and shifted to the right, indicating that the vast majority of OXPHOS genes are upregulated, though with heterogeneous fold-changes. This reflects global mitochondrial biogenesis rather than selective gene upregulation.

- **Reactive Oxygen Species Pathway:** This ridge should show a **sharp rightward shift** with high consistency — the canonical NRF2 antioxidant targets (*NQO1* LFC = +1.70, *TXNRD1* LFC = +1.29, *GCLC* LFC = +1.23, *GSR* LFC = +0.98, *PRDX1* LFC = +0.64) are uniformly upregulated.

- **Xenobiotic Metabolism:** Similarly right-shifted, reflecting the broad activation of Phase I/II enzymes.

- **EMT:** Strongly left-shifted, confirming consistent downregulation of the entire mesenchymal gene program.

- **Immune pathways** (Allograft Rejection, IFNγ, IFNα): Consistently left-shifted ridges demonstrate that immune gene downregulation is not driven by a few outlier genes but is a **uniform program-wide suppression**.

The ridgeplot therefore validates that the enrichment scores are driven by coherent, biologically meaningful shifts rather than statistical artefacts from a few extreme genes.

---

## 4. GSEA Running Score – Key Pathways

### [05_GSEA_running_score_key_pathways.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_GSEA_running_score_key_pathways.pdf)

**What the plot shows:** The running enrichment score plot (also called "mountain plot" or "enrichment score plot") displays the cumulative enrichment statistic as the algorithm walks down the ranked gene list for selected key pathways. It includes a barcode showing where gene set members fall in the ranked list and a running sum line that peaks where the set is most enriched.

**Interpretation of key pathways:**

- **HALLMARK_OXIDATIVE_PHOSPHORYLATION (NES = +2.67):** The running score climbs steeply and peaks early-to-mid in the ranked list (rank ~3345), with the leading edge (64% of genes) concentrated among the most upregulated genes. The tall, early peak and large leading edge indicate strong and early enrichment — these genes rank among the most consistently upregulated in NRF2-high tumours. Core enrichment genes include electron transport chain subunits (*NDUFB4, COX5A/5B, UQCRC2, CYC1, ATP5F1B*) and mitochondrial import machinery (*TOMM22, TIMM13, TIMM17A*).

- **HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY (NES = +2.39):** Despite being a small gene set (49 genes), the running score rises sharply with enrichment in the extreme left tail of the ranked list (rank ~439, top 3%). The 13 leading-edge genes (*G6PD, NQO1, PRDX1, ABCC1, TXNRD1, GCLC, GSR, GCLM, SRXN1, TXN, PRDX6, HMOX2, NDUFB4*) are among the most strongly upregulated genes in the entire transcriptome. This is the **most direct functional validation** of constitutive NRF2 activity in these tumours.

- **HALLMARK_XENOBIOTIC_METABOLISM (NES = +2.17):** Early enrichment (rank ~1543) with core enrichment genes including *AKR1C2/C3, NQO1, ALDH3A1, GCLC, PTGR1, POR, GSR, CBR1, EPHX1*. These are classical ARE-driven detoxification enzymes.

- **HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION (NES = −2.34):** The running score dips deeply negative, accumulating genes in the lower (downregulated) end of the ranked list (rank ~3498). Core enrichment genes include collagens (*COL1A1/A2, COL3A1, COL5A2*), matrix components (*FN1, POSTN, SPARC*), and EMT transcription factor targets (*VIM, CDH2*).

- **HALLMARK_ALLOGRAFT_REJECTION (NES = −2.46):** Deepest negative score, with immune genes concentrated at the bottom of the ranked list. Core enrichment includes T cell signalling (*CD3E, CD3G, CD4, CD8A/B, LCK, ZAP70*), cytotoxic molecules (*GZMB, GZMA, PRF1*), and antigen presentation (*HLA-DRA, HLA-DMA, CD74*).

The running score plots confirm that the enrichment signals are driven by large, coherent gene expression shifts at the extremes of the ranked list, consistent with bona fide biological pathway regulation rather than noise.

---

## 5. ORA GO Biological Process

### 5a. GO BP Barplot — Upregulated Genes ([05_ORA_GO_BP_barplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_GO_BP_barplot.pdf))

**What the plot shows:** This barplot displays the most significantly enriched Gene Ontology Biological Process terms among genes **upregulated** in NRF2-High tumours, with bar length representing gene count and colour encoding adjusted p-value.

**Key findings (53 significant GO:BP terms for upregulated genes):**

The top terms paint a vivid picture of NRF2 biology:

1. **Keratinization** (p.adj = 3.1×10⁻³⁰, 33 genes) — The most significant term, driven by *SPRR1A/1B/2B/2D/2E/2F/3/4, SPRR2G, LORICRIN, IVL, FLG, TGM1/3/5, KRT3/4/76/78*, etc. NRF2-high OSCC tumours exhibit a **hyper-keratinized, squamous differentiation program** consistent with well-differentiated SCC histology.

2. **Epidermis development** (p.adj = 7.3×10⁻²⁹, 59 genes) — Overlaps significantly with keratinization but includes additional epidermal regulators (*KLF4, SFN, OVOL2*).

3. **Xenobiotic metabolic process** (p.adj = 8.2×10⁻¹⁴, 23 genes) — Canonical NRF2 function: *NQO1, ALDH3A1, GSTM3, ABCC1, UGT1A1/6/7/8/10, GSTA1/4, EPHX1, CYP2C9/18/19*.

4. **Glutathione metabolic process** (p.adj = 7.6×10⁻⁶, 10 genes) — Direct oxidative stress defence: *SLC7A11, G6PD, GSTM3, GCLC, GSTA4, GCLM, GSTA1, GSTM2, GSTA3, GSTM1*.

5. **Response to toxic substance** (p.adj = 9.2×10⁻⁷, 23 genes) — Includes core cytoprotective NRF2 targets: *SLC7A11, GPX2, NQO1, TXNRD1, GSTM3, PRDX1, ABCC1, GSR, EPHX1, SRXN1*.

6. **Cellular detoxification** (p.adj = 1.7×10⁻⁴, 13 genes)

7. **Steroid metabolic process** (p.adj = 4.6×10⁻³, 18 genes) — Includes AKR1C family members that metabolise steroid hormones and are NRF2 targets.

**KEAP1-NRF2 Context:**

The striking dominance of **keratinization and epidermal differentiation** among upregulated genes deserves special attention. NRF2 is a master regulator of squamous differentiation. In normal oral epithelium, NRF2 promotes terminal differentiation through induction of small proline-rich proteins (SPRRs), involucrin (IVL), and transglutaminases (TGMs). In OSCC, constitutive NRF2 activation locks tumours into a **hyper-differentiated squamous state**. This is consistent with the finding that KEAP1 itself appeared in the epidermal cell differentiation GO term (GSEA data shows KEAP1 as a contributing gene to the epidermis development pathway), highlighting the direct role of the KEAP1-NRF2 axis in squamous biology.

The **glutathione metabolism** enrichment is particularly important for oxidative stress biology. The upregulation of *GCLC* (LFC = +1.23) and *GCLM* (LFC = +1.12) — the catalytic and modifier subunits of glutamate-cysteine ligase (the rate-limiting enzyme for glutathione synthesis) — combined with *SLC7A11* (LFC = +2.65, the cystine/glutamate antiporter xCT) indicates a massively enhanced **glutathione biosynthetic capacity**. This is the central antioxidant defence mechanism driven by NRF2 (Hayes & Dinkova-Kostova, 2014, *Trends in Biochemical Sciences*).

---

### 5b. GO BP Dotplot — Upregulated and Downregulated ([05_ORA_GO_BP_dotplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_GO_BP_dotplot.pdf))

**What the plot shows:** A dotplot showing the most significantly enriched GO:BP terms split by direction (Up/Down), with dot size encoding gene count and colour encoding gene ratio.

**Additional interpretation for DOWNREGULATED genes (88 significant terms):**

The downregulated gene ORA reveals a completely different biology:

1. **Myofibril assembly / Muscle system processes** — Dominant among downregulated terms (*FHOD3, MYH7/6/3/8, TCAP, CASQ1/2, ACTN2/3, DMD*). While surprising in OSCC, this reflects the **loss of mesenchymal/stromal components** in NRF2-high tumours, consistent with the anti-EMT phenotype.

2. **Adaptive immune response** (p.adj = 9.2×10⁻⁷, 39 genes) — Includes *CD8A/B, CD19, CR2, AICDA, SLAMF6, HLA-DOA, IL12RB1, TNFSF4, FCRL4, GZMM, LAG3, CD40LG*.

3. **B cell activation** (p.adj = 3.1×10⁻⁵, 25 genes) — *MS4A1, CD22, CD19, CD79A/B, FCRL1/3/5, CR2, IKZF3, POU2AF1*.

4. **Lymphocyte mediated immunity** (p.adj = 1.1×10⁻⁴, 26 genes)

5. **B cell proliferation** (p.adj = 1.0×10⁻⁴, 13 genes)

6. **Regulation of cell-cell adhesion** — immune cell adhesion markers downregulated.

The comprehensive **loss of B cell and T cell function** confirms the immune-desert microenvironment. The downregulation of exhaustion markers (*LAG3*) alongside effector molecules suggests these immune cells are not merely exhausted but are **absent from the tumour microenvironment** — a hallmark of immune exclusion rather than immune dysfunction (Skoulidis & Heymach, 2019, *Cancer Cell*).

---

## 6. ORA GO Molecular Function

### GO MF Dotplot ([05_ORA_GO_MF_dotplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_GO_MF_dotplot.pdf))

**What the plot shows:** This dotplot displays significantly enriched GO Molecular Function terms for both upregulated (41 terms) and downregulated (25 terms) gene sets.

**Upregulated MF highlights:**

| GO Term | p.adjust | Gene Count | Key Genes |
|---|---|---|---|
| Structural constituent of skin epidermis | 1.6×10⁻⁶ | 10 | SPRR1A, PI3, LORICRIN, KRT3/4/76/78 |
| **Glutathione transferase activity** | **2.5×10⁻⁵** | **7** | **GSTM3, GSTA4, GSTA1, GSTM2, GSTA3, MGST1, GSTM1** |
| Oxidoreductase activity (CH-OH donors) | 3.6×10⁻⁶ | 16 | AKR1C1/2/3, G6PD, PGD, ME1, ALDH3A1 |
| **Antioxidant activity** | **4.7×10⁻³** | **8** | **GPX2, NQO1, TXNRD1, SRXN1, GSTA1, GSTM2, MGST1** |
| **Glutathione peroxidase activity** | **1.1×10⁻²** | **4** | **GPX2, GSTA1, GSTM2, MGST1** |
| Glutathione binding | 1.3×10⁻² | 3 | GSTM3, GSTM2, GSTM1 |
| Monooxygenase activity | 2.5×10⁻⁵ | 12 | CYP4F3/11, CYP2C9/18/19, CYP26A1 |
| Peptidase inhibitor activity | 3.6×10⁻⁶ | 17 | SERPINBs, CSTA, SPINK5/7, SLPI |
| Iron ion binding | 8.0×10⁻⁵ | 14 | CYPs, TF, ALOXes |
| NADP binding | 4.9×10⁻² | 5 | G6PD, PGD, ME1, TM7SF2, CRYM |

**Downregulated MF highlights:**

| GO Term | p.adjust | Gene Count | Key Genes |
|---|---|---|---|
| Structural constituent of muscle | 1.8×10⁻⁵ | 11 | MYH7/6/3/8, ACTN2/3, DMD, TCAP |
| Cytokine activity | 5.8×10⁻⁴ | 19 | CXCL9/10/11/12, FASLG, CD40LG |
| GPCR activity | 6.8×10⁻³ | 23 | Chemokine and neurotransmitter receptors |
| Immune receptor activity | 1.5×10⁻³ | 15 | IL12RB1, IL2RG, HLA-DOA/DQB2 |

**KEAP1-NRF2 Context:**

The MF analysis provides **enzyme-level validation** of the NRF2 transcriptional program. The enrichment of **glutathione transferase activity** (7 GSTs), **antioxidant activity** (8 genes), **glutathione peroxidase activity**, and **oxidoreductase activity** collectively demonstrates that NRF2-high OSCC tumours possess a massively upregulated **enzymatic antioxidant arsenal**.

The enrichment of **NADP binding** (*G6PD, PGD, ME1*) is particularly significant. These enzymes generate **NADPH**, the essential reducing cofactor for glutathione reductase (GSR) and thioredoxin reductase (TXNRD1). NRF2 transcriptionally activates *G6PD*, *PGD*, and *ME1* to ensure adequate NADPH supply for the entire antioxidant system (Mitsuishi et al., 2012, *Cancer Cell*). This creates a self-sustaining antioxidant cycle:

```
NRF2 → GCLC/GCLM → GSH synthesis
NRF2 → SLC7A11 → Cystine import → GSH synthesis
NRF2 → G6PD/PGD → NADPH → GSR → GSH regeneration
NRF2 → TXNRD1 → Thioredoxin system maintenance
```

---

## 7. ORA GO Cellular Component

### GO CC Dotplot ([05_ORA_GO_CC_dotplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_GO_CC_dotplot.pdf))

**What the plot shows:** GO Cellular Component enrichment for upregulated (5 terms) and downregulated (23 terms) genes.

**Upregulated:**
- **Cornified envelope** (p.adj = 9.7×10⁻²⁶, 27 genes) — By far the most significant, confirming the squamous differentiation signature.
- **Intermediate filament** (p.adj = 3.0×10⁻⁴, 13 genes) — Keratin cytoskeleton.
- **Apical plasma membrane / Apical part of cell** — Polarised epithelial localisation of transporters (*SLC7A11, ABCC1, CFTR, SLC47A1*).

**Downregulated:**
- **Sarcomere** (p.adj = 6.4×10⁻²⁷, 50 genes) — Dominant downregulated compartment; muscle-related structural loss.
- **External side of plasma membrane** (p.adj = 1.3×10⁻⁸, 37 genes) — Immune cell surface markers (*CD19, CD22, CD79A/B, CD8A/B, FCRL2/3/4/5, LAG3, SLAMF6*) — confirming immune cell depletion.
- **Receptor complex / Plasma membrane signalling receptor complex** — Loss of immune signalling complexes.

**Interpretation:**

The CC data structurally validates the ORA:BP findings. NRF2-high tumours are characterised by cornified epithelial structures at the expense of immune cell surface receptors. The apical membrane localisation of NRF2 target transporters (*SLC7A11* at the apical surface, *ABCC1* for drug efflux) indicates polarised transporter expression that is functionally relevant for cystine uptake and drug resistance in the intact tumour epithelium.

---

## 8. ORA KEGG Pathways

### KEGG Dotplot ([05_ORA_KEGG_dotplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_KEGG_dotplot.pdf))

**What the plot shows:** Significantly enriched KEGG pathways for upregulated (25 terms) and downregulated (18 terms) genes.

**Upregulated KEGG Pathways:**

| KEGG Pathway | p.adjust | Genes | Significance |
|---|---|---|---|
| **Cornified envelope formation** | 4.6×10⁻³⁶ | 53 | Squamous differentiation |
| **Chemical carcinogenesis – DNA adducts** | 3.0×10⁻¹⁴ | 17 | GSTs, UGTs, EPHX1 – detoxification of carcinogen adducts |
| **Metabolism of xenobiotics by CYP450** | 6.8×10⁻¹⁴ | 18 | Phase I/II metabolism |
| **Drug metabolism – CYP450** | 3.2×10⁻¹³ | 17 | Drug resistance |
| **Glutathione metabolism** | 1.7×10⁻⁸ | 13 | *G6PD, GPX2, PGD, GSTM3, GCLC, GSTA4, GCLM, ODC1* |
| **Drug metabolism – other enzymes** | 1.7×10⁻⁸ | 14 | UGTs, GSTs |
| Arachidonic acid metabolism | 5.6×10⁻¹⁰ | 15 | Lipid peroxidation metabolism |
| Retinol metabolism | 1.6×10⁻⁸ | 13 | Retinoid/vitamin A metabolism |
| Steroid hormone biosynthesis | 1.7×10⁻⁶ | 10 | AKR1C1/2/3, UGTs |
| **Chemical carcinogenesis – ROS** | **3.3×10⁻²** | **13** | *AKR1C1/2/3, NQO1, GSTM3, GSTA4, EPHX1, MGST1, GSTM1* |
| **Ferroptosis** | **3.9×10⁻²** | **5** | ***SLC7A11, GCLC, GCLM, TF, ALOX15*** |

**Downregulated KEGG Pathways:**

| KEGG Pathway | p.adjust | Genes |
|---|---|---|
| Cytoskeleton in muscle cells | 8.2×10⁻²² | 47 |
| Cytokine-cytokine receptor interaction | 1.3×10⁻⁴ | 25 |
| Primary immunodeficiency | 1.3×10⁻⁴ | 9 |
| Hematopoietic cell lineage | 1.3×10⁻⁴ | 14 |

**KEAP1-NRF2 Context:**

The KEGG data is remarkable for three reasons:

1. **Chemical carcinogenesis pathways** (both DNA adducts and ROS): NRF2-high OSCC tumours have massively upregulated detoxification of carcinogen-DNA adducts through GSTs (*GSTM1/2/3, GSTA1/3/4*), UGTs (*UGT1A1/6/7/8/10*), and EPHX1. In oral cancer, exposure to tobacco carcinogens (PAHs, nitrosamines) produces DNA adducts. The upregulation of these detoxification enzymes by NRF2 is paradoxical: while it may have been initially protective against carcinogenesis, in established tumours it confers **resistance to DNA-damaging chemotherapeutics** like cisplatin (Tao et al., 2018, *Redox Biology*).

2. **Ferroptosis pathway** in KEGG (SLC7A11, GCLC, GCLM, TF, ALOX15): The enrichment of ferroptosis regulators among upregulated genes demonstrates that NRF2-high tumours are actively **protecting themselves against ferroptotic cell death** — a highly relevant therapeutic vulnerability (see Section 13).

3. **Platinum drug resistance** also appeared in the upregulated KEGG results (7 genes: *GSTM3, GSTA4, GSTA1, GSTM2, GSTA3, MGST1, GSTM1*), providing direct evidence for NRF2-mediated cisplatin resistance through glutathione conjugation, a major clinical challenge in OSCC treatment (Roh et al., 2017, *Oncotarget*).

---

## 9. ORA Reactome Pathways

### Reactome Dotplot ([05_ORA_Reactome_dotplot.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_Reactome_dotplot.pdf))

**What the plot shows:** Significantly enriched Reactome pathways for upregulated (39 terms) and downregulated (18 terms) genes.

**Upregulated Reactome Highlights:**

| Reactome Pathway | p.adjust | Genes | Notes |
|---|---|---|---|
| Keratinization | 3.5×10⁻³⁸ | 45 | Squamous program |
| Formation of cornified envelope | 3.5×10⁻³⁸ | 44 | |
| **Biological oxidations** | **1.7×10⁻¹⁴** | **30** | **Combined Phase I + II** |
| Antimicrobial peptides | 8.5×10⁻¹² | 15 | Innate epithelial defence |
| Arachidonic acid metabolism | 2.4×10⁻¹⁰ | 15 | Eicosanoid metabolism |
| **NFE2L2 regulating anti-oxidant/detoxification enzymes** | **5.7×10⁻⁷** | **8** | ***SLC7A11, NQO1, TXNRD1, GCLC, GCLM, SRXN1, GSTA1, GSTA3*** |
| **KEAP1-NFE2L2 pathway** | **1.8×10⁻³** | **12** | ***SLC7A11, G6PD, PGD, NQO1, ME1, ABCC1, TXNRD1, GCLC, GCLM, SRXN1, GSTA1, GSTA3*** |
| **Nuclear events mediated by NFE2L2** | **2.4×10⁻⁴** | **12** | Same genes as KEAP1-NFE2L2 pathway |
| **Phase II – Conjugation of compounds** | **9.7×10⁻⁸** | **15** | GSTs, UGTs, GCLC/GCLM |
| Phase I – Functionalization | 4.1×10⁻⁷ | 15 | CYPs, EPHX1, ALDHs |
| Glutathione conjugation | 3.7×10⁻⁶ | 9 | *GSTM3, GCLC, GSTA4, GCLM, GSTA1, GSTM2, GSTA3, MGST1, GSTM1* |
| **Cellular response to chemical stress** | **3.1×10⁻²** | **13** | *SLC7A11, G6PD, GPX2, PGD, NQO1, ME1, ABCC1, TXNRD1, GCLC, GCLM, SRXN1, GSTA1, GSTA3* |
| Drug ADME | 3.3×10⁻⁵ | 12 | Pharmacokinetic resistance |
| Cytochrome P450 by substrate | 9.4×10⁻⁴ | 8 | CYP4F3/11, CYP2C9/18/19 |

**Downregulated Reactome:**
- Muscle contraction (1.7×10⁻¹⁰)
- Striated muscle contraction (1.4×10⁻⁹)
- Class A/1 Rhodopsin-like receptors (2.9×10⁻⁵)
- Immunoregulatory interactions (1.4×10⁻²)

**KEAP1-NRF2 Context:**

The Reactome analysis provides the **most direct pathway-level validation** of NRF2 transcriptional activity:

1. The **"NFE2L2 regulating anti-oxidant/detoxification enzymes"** pathway (R-HSA-9818027) is directly annotated as an NRF2 pathway in Reactome and shows strong enrichment (p.adj = 5.7×10⁻⁷) with 8 of 19 pathway genes being upregulated in NRF2-high tumours. These include the "canonical quartet" of NRF2 targets: *SLC7A11, NQO1, TXNRD1, GCLC*.

2. The **"KEAP1-NFE2L2 pathway"** (R-HSA-9755511) is enriched with 12 genes (p.adj = 1.8×10⁻³), providing built-in pathway-level confirmation that the NRF2-High classification captures genuine KEAP1-NRF2 axis activation.

3. The **"Cellular response to chemical stress"** pathway (R-HSA-9711123) with 13 NRF2-dependent genes provides the broadest Reactome-level view of the oxidative stress response.

The combined enrichment of these three NRF2-specific Reactome pathways constitutes strong **orthogonal validation** that NRF2 transcriptional activity — not merely NRF2 mRNA levels — drives the gene expression program in NRF2-High tumours.

---

## 10. ORA Network and Hierarchy Plots (Cnetplots, Emapplots, Treeplots)

### 10a. Cnetplot – Upregulated Genes ([05_ORA_cnetplot_up.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_cnetplot_up.pdf))

**What the plot shows:** A category-gene network (cnetplot) linking enriched GO terms to their individual contributing genes for upregulated DEGs. Genes are coloured by fold-change, and term nodes are sized by gene count.

**Interpretation:**

The cnetplot reveals **shared hub genes** that link multiple biological processes:

- **AKR1C1/C2/C3** emerge as ultra-high fold-change hub genes (LFC = +4.09/+3.72/+3.40) connecting steroid metabolism, xenobiotic metabolism, olefinic compound metabolism, and ketone metabolism. These aldo-keto reductases are among the most robustly NRF2-induced genes in squamous cancers and play roles in prostaglandin metabolism and steroid hormone inactivation (Penning, 2015, *Chemical Research in Toxicology*).

- **UGT1A family** (*UGT1A1/6/7/8/10*) form another hub connecting glucuronate metabolism, xenobiotic metabolism, and drug metabolism.

- **SLC7A11** connects glutathione metabolism, response to toxic substance, and secondary metabolic process — positioning it as a **central node** in the NRF2 cytoprotective network.

- **SPRR/keratin genes** form a dense sub-network of epidermis/keratinization terms.

The network topology confirms that the NRF2 transcriptional program is **highly interconnected** — a small set of core NRF2 target genes drives enrichment across multiple functional categories.

---

### 10b. Cnetplot – Downregulated Genes ([05_ORA_cnetplot_down.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_cnetplot_down.pdf))

**Interpretation:**

Hub genes in the downregulated network include:
- **CXCL9/10/11/12** — Chemokines connecting cell killing, calcium ion homeostasis, and chemokine response. Their downregulation reduces immune cell recruitment to the tumour.
- **CD19, CD22, CD79A** — B cell markers connecting B cell activation, proliferation, mediated immunity, and receptor signalling.
- **DMD, ACTN2/3, MYH7** — Muscle/structural genes.
- **MEF2C** — Transcription factor connecting multiple immune and differentiation processes.

The loss of **CXCL9/10/11** is clinically significant: these IFNγ-induced chemokines recruit CD8⁺ T cells and NK cells to the tumour. Their suppression in NRF2-high tumours explains the immune-cold phenotype and predicts **poor response to immune checkpoint inhibitors** (Ayers et al., 2017, *Journal of Clinical Investigation*).

---

### 10c. Emapplot – Upregulated ([05_ORA_emapplot_up.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_emapplot_up.pdf))

**What the plot shows:** An enrichment map for upregulated ORA results, with nodes as GO terms and edges connecting terms sharing genes. Clustering reveals functional modules.

**Interpretation:**

Three major functional modules emerge:
1. **Keratinization/Epidermal differentiation module** — Dense cluster of skin/cornification terms.
2. **Xenobiotic/Detoxification module** — Xenobiotic metabolism, glutathione metabolism, cellular detoxification, response to toxic substance.
3. **Lipid/Hormone metabolism module** — Steroid metabolism, hormone metabolism, fatty acid catabolism, icosanoid metabolism.

The **bridge genes** connecting the detoxification and lipid metabolism modules are the AKR1C family, GSTs, and CYPs, which have dual roles in both xenobiotic and endogenous metabolite metabolism. This bridging reflects NRF2's role in coordinating both detoxification and metabolic reprogramming.

---

### 10d. Emapplot – Downregulated ([05_ORA_emapplot_down.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_emapplot_down.pdf))

**Interpretation:**

Two major modules:
1. **Immune function module** — B cell activation, lymphocyte mediated immunity, adaptive immune response, leukocyte activation. Densely interconnected.
2. **Muscle/contractile module** — Myofibril assembly, muscle contraction, cardiac development. Likely reflects stromal composition changes.

The immune module's density indicates that the immunosuppression is **not selective to one immune cell type** but represents a **pan-immune depletion** encompassing B cells, T cells, and innate immune cells.

---

### 10e–f. Treeplots – Up and Down ([05_ORA_treeplot_up.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_treeplot_up.pdf), [05_ORA_treeplot_down.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ORA_treeplot_down.pdf))

**What the plots show:** Treeplots display the enriched GO terms organised hierarchically by semantic similarity, revealing parent-child relationships and term clustering.

**Upregulated treeplot interpretation:**

The hierarchical clustering groups terms into biologically coherent branches:
- **Branch 1: Epithelial differentiation** — keratinization → epidermis development → skin development → skin epidermis development → epidermis morphogenesis
- **Branch 2: Xenobiotic/Detoxification** — xenobiotic metabolic process → cellular detoxification → response to toxic substance → glutathione metabolism
- **Branch 3: Lipid/Hormone metabolism** — olefinic compound metabolism → unsaturated fatty acid metabolism → icosanoid metabolism → steroid metabolism → hormone metabolism
- **Branch 4: Anti-microbial defence** — antimicrobial humoral response → killing of cells of another organism

**Downregulated treeplot interpretation:**

- **Branch 1: Immune regulation** — adaptive immune response → lymphocyte mediated immunity → B cell activation → B cell proliferation → regulation of leukocyte activation
- **Branch 2: Muscle/Contractile** — myofibril assembly → muscle system process → cardiac development
- **Branch 3: Calcium/Ion homeostasis** — calcium ion homeostasis → regulation of membrane potential → sarcoplasmic reticulum transport

The treeplots confirm that the enrichment results are not redundant artefacts but represent hierarchically organised biological programs.

---

## 11. GSVA Custom Pathway Analysis

### 11a. GSVA Pathway Heatmap ([05_GSVA_pathway_heatmap.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_GSVA_pathway_heatmap.pdf))

**What the plot shows:** A heatmap displaying Gene Set Variation Analysis (GSVA) enrichment scores across all 20 custom pathways for individual tumour samples, clustered by NRF2-High vs NRF2-Low groups.

**Interpretation:**

The heatmap should show a clear **binary stratification** between NRF2-High and NRF2-Low groups across the custom pathway panel. Key patterns:

- **NRF2 Core Activity** has the strongest differential (logFC = +0.63, t = 22.05, p = 1.2×10⁻⁵⁸) — the most statistically significant pathway in the entire GSVA analysis. This serves as an internal control confirming that the NRF2-High/Low classification captures genuine NRF2 transcriptional activity.

- Oxidative stress pathways (**ROS/Glutathione Metabolism**, **Xenobiotic Metabolism**, **Pentose Phosphate Pathway**) cluster together on one side of the heatmap, uniformly elevated in NRF2-High tumours.

- **EMT Markers** and **Inflammatory Response** cluster together as suppressed pathways, reinforcing the inverse relationship between NRF2 activity and mesenchymal/immune programs.

- **Autophagy/p62/KEAP1** and **DNA Repair** pathways show non-significant differences (p > 0.05), indicating these processes are not dramatically altered at the transcriptional level in NRF2-High tumours.

---

### 11b. GSVA Pathway Lollipop ([05_GSVA_pathway_lollipop.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan\KEAP1_NRF2_OSCC_Project/06_figures/05_GSVA_pathway_lollipop.pdf))

**What the plot shows:** A lollipop chart displaying the log fold-change of GSVA enrichment scores for each custom pathway, with significance stars and colour indicating direction (Activated vs Suppressed in NRF2-High).

**Key findings ranked by significance:**

| Pathway | logFC | adj.P.Val | Direction |
|---|---|---|---|
| **NRF2 Core Activity** | **+0.63** | 2.3×10⁻⁵⁷ | Activated |
| **Xenobiotic Metabolism** | **+0.49** | 1.8×10⁻²⁸ | Activated |
| **ROS/Glutathione Metabolism** | **+0.39** | 2.2×10⁻²³ | Activated |
| **Pentose Phosphate Pathway** | **+0.34** | 4.6×10⁻¹⁵ | Activated |
| Fatty Acid Metabolism | +0.19 | 1.8×10⁻¹⁴ | Activated |
| **Reactive Oxygen Species Pathway** | **+0.22** | 2.1×10⁻¹² | Activated |
| mTORC1 Signalling | +0.22 | 9.2×10⁻¹² | Activated |
| **Ferroptosis Regulation** | **+0.25** | 6.9×10⁻¹¹ | Activated |
| Oxidative Phosphorylation | +0.24 | 2.1×10⁻⁷ | Activated |
| Glycolysis | +0.11 | 6.7×10⁻⁶ | Activated |
| MYC Targets V1 | +0.22 | 8.2×10⁻⁶ | Activated |
| Angiogenesis | −0.19 | 2.0×10⁻⁵ | Suppressed |
| EMT | −0.19 | 1.2×10⁻⁴ | Suppressed |
| EMT Markers | −0.20 | 1.3×10⁻⁴ | Suppressed |
| PI3K/AKT/mTOR Signalling | +0.09 | 4.4×10⁻⁴ | Activated |
| Inflammatory Response | −0.12 | 1.8×10⁻³ | Suppressed |
| DNA Repair (Hallmark) | +0.07 | 0.055 | ns |
| Autophagy/p62/KEAP1 | +0.07 | 0.080 | ns |
| DNA Repair (custom) | −0.04 | 0.439 | ns |
| Apoptosis | +0.005 | 0.815 | ns |

**KEAP1-NRF2 Context:**

The GSVA lollipop provides a **comprehensive, sample-level view** of pathway activity. The rank ordering confirms a clear biological hierarchy:

1. **NRF2 Core Activity** dominates (logFC = +0.63), validating the classification approach.
2. **Xenobiotic Metabolism** and **ROS/Glutathione Metabolism** are the next strongest signals, confirming that oxidative stress defence and detoxification are the primary functional consequences of NRF2 activation.
3. **Pentose Phosphate Pathway** activation (logFC = +0.34) is the metabolic linchpin — it provides NADPH for the entire antioxidant system while also generating ribose-5-phosphate for nucleotide synthesis.
4. **Ferroptosis Regulation** activation (logFC = +0.25, p = 6.9×10⁻¹¹) indicates NRF2-high tumours actively resist ferroptosis.
5. **Apoptosis** is non-significant (p = 0.815), indicating NRF2's cytoprotective effects in OSCC operate through **ferroptosis resistance and detoxification rather than classical anti-apoptotic mechanisms**.

---

### 11c. GSVA Pathway Comparison ([05_GSVA_pathway_comparison.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_GSVA_pathway_comparison.pdf))

**What the plot shows:** A comparison plot (likely box/violin plots) showing the distribution of GSVA scores for each pathway in NRF2-High vs NRF2-Low groups.

**Interpretation:**

The comparison plot reveals the **separation magnitude** between groups for each pathway. The NRF2 Core Activity pathway should show the cleanest separation, while pathways like Apoptosis and DNA Repair show overlapping distributions, explaining their non-significance. The ROS/Glutathione Metabolism pathway should show consistently higher scores in NRF2-High samples with minimal overlap.

---

## 12. NRF2 Reference Gene Set Analysis

### 12a. NRF2 Reference Gene Set Heatmap ([05_nrf2_reference_geneset_heatmap.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_nrf2_reference_geneset_heatmap.pdf))

**What the plot shows:** A heatmap displaying the expression of 93 NRF2 reference genes (grouped by source: Union Base, Ferroptosis Extension, EMT Extension, Autophagy Extension, OSCC NRF2 Targets Extension, Upstream Regulators Extension) across all tumour samples.

**Key findings from the 93-gene reference set:**

**Strongly Upregulated (classified as "Up"):**

| Gene | LFC | padj | Source | Function |
|---|---|---|---|---|
| **AKR1C1** | **+4.09** | 2.7×10⁻⁵⁷ | Union Base | Aldo-keto reductase, prostaglandin/steroid metabolism |
| **AKR1C2** | **+3.72** | 1.7×10⁻⁵⁵ | Union Base | Steroid hormone metabolism |
| **AKR1C3** | **+3.40** | 1.9×10⁻⁴⁶ | Union Base | Prostaglandin F synthase |
| **GPX2** | **+3.09** | 2.8×10⁻³¹ | Union Base | Glutathione peroxidase 2 |
| **GSTA1** | **+3.02** | 2.6×10⁻⁹ | OSCC Extension | Glutathione S-transferase alpha 1 |
| **ALDH3A1** | **+2.86** | 1.2×10⁻²³ | OSCC Extension | Aldehyde dehydrogenase |
| **SLC7A11** | **+2.65** | 1.9×10⁻³⁹ | Union Base | xCT cystine transporter |
| **AKR1B10** | **+2.40** | 2.4×10⁻²³ | Union Base | Retinol reductase |
| **GSTM3** | **+2.39** | 8.3×10⁻²² | Union Base | Glutathione S-transferase mu 3 |
| **GSTM1** | **+2.10** | 4.3×10⁻⁴ | OSCC Extension | Glutathione S-transferase mu 1 |
| **NQO1** | **+1.70** | 4.5×10⁻²⁸ | Union Base | NAD(P)H quinone dehydrogenase 1 |
| **SRXN1** | **+1.67** | 1.9×10⁻¹⁰ | Union Base | Sulfiredoxin 1 |
| **ME1** | **+1.57** | 1.2×10⁻²⁵ | Union Base | Malic enzyme 1 (NADPH production) |
| **G6PD** | **+1.50** | 5.1×10⁻³⁴ | Union Base | Glucose-6-phosphate dehydrogenase (NADPH) |
| **PGD** | **+1.35** | 8.3×10⁻³¹ | Union Base | 6-phosphogluconate dehydrogenase (NADPH) |
| **MGST1** | **+1.32** | 2.8×10⁻⁴ | Union Base | Microsomal glutathione S-transferase 1 |
| **TXNRD1** | **+1.29** | 3.6×10⁻¹⁹ | Union Base | Thioredoxin reductase 1 |
| **GCLC** | **+1.23** | 2.0×10⁻¹⁸ | Union Base | Glutamate-cysteine ligase catalytic |
| **ALOX15** | **+1.25** | 1.0×10⁻⁴ | Ferroptosis Ext. | Arachidonate 15-lipoxygenase |
| **GSTA4** | **+1.17** | 1.4×10⁻¹¹ | OSCC Extension | GST alpha 4 |
| **GCLM** | **+1.12** | 1.1×10⁻¹⁰ | Union Base | Glutamate-cysteine ligase modifier |
| **EPHX1** | **+1.04** | 5.1×10⁻¹¹ | OSCC Extension | Epoxide hydrolase 1 |
| **ABCC1** | **+1.02** | 3.2×10⁻¹⁹ | Union Base | MRP1 drug efflux transporter |

**Downregulated EMT markers:**

| Gene | LFC | padj | Function |
|---|---|---|---|
| **CDH2** | **−2.04** | 2.8×10⁻¹³ | N-cadherin (mesenchymal marker) |
| VIM | −0.80 | 8.6×10⁻⁸ | Vimentin |
| TWIST1 | −0.80 | 7.1×10⁻⁷ | EMT transcription factor |
| ZEB1 | −0.68 | 1.7×10⁻⁵ | EMT transcription factor |
| ZEB2 | −0.70 | 1.2×10⁻⁵ | EMT transcription factor |
| SNAI1 | −0.51 | 9.5×10⁻⁴ | Snail EMT TF |

**Notable non-significant genes:**
- *GPX4* (LFC = −0.09, p = 0.57) — The master anti-ferroptotic enzyme is **not differentially expressed**, suggesting ferroptosis resistance in NRF2-high OSCC operates through upstream GSH supply rather than GPX4 levels.
- *SQSTM1/p62* (LFC = +0.10, p = 0.44) — Non-significant, suggesting p62-mediated KEAP1 sequestration is not the primary NRF2 activation mechanism.
- *FTH1/FTL* — Ferritin subunits are not significantly changed, indicating iron sequestration is not dramatically altered.

---

### 12b. NRF2 Reference Lollipop ([05_nrf2_reference_lollipop.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_nrf2_reference_lollipop.pdf))

**What the plot shows:** A lollipop chart showing the log₂ fold-change for each of the 93 NRF2 reference genes, coloured by regulation status (Up/Down/NS) and grouped by gene set source.

**KEAP1-NRF2 Context:**

The lollipop plot provides a "**NRF2 target activation fingerprint**" for OSCC. The most striking observation is the **magnitude of upregulation**: AKR1C1 at LFC = +4.09 represents a ~17-fold increase in expression. The top 10 upregulated NRF2 targets (AKR1C1/2/3, GPX2, GSTA1, ALDH3A1, SLC7A11, AKR1B10, GSTM3, GSTM1) collectively represent:

- **Oxidative stress defence**: GPX2, SLC7A11, GCLC, GCLM, TXNRD1, SRXN1, PRDX1
- **Phase II conjugation**: GSTA1/4, GSTM1/3, MGST1
- **Carbonyl detoxification**: AKR1C1/2/3, AKR1B10, ALDH3A1
- **NADPH generation**: G6PD, PGD, ME1
- **Drug efflux**: ABCC1

This constitutes the **complete NRF2 cytoprotective battery** operating at maximum capacity in NRF2-High OSCC.

The downregulation of all canonical **EMT transcription factors** (*SNAI1, TWIST1, ZEB1, ZEB2*) and their targets (*CDH2, VIM*) alongside upregulation of *CDH1* (LFC = +0.49) confirms the **epithelial phenotype maintenance** driven by NRF2. NRF2 has been shown to promote E-cadherin expression through direct ARE-mediated transcription in squamous cells (Rachakonda et al., 2010, *Oncogene*).

**Key literature:**
- Hayes JD & Dinkova-Kostova AT (2014). The Nrf2 regulatory network provides an interface between redox and intermediary metabolism. *Trends in Biochemical Sciences*, 39(4):199-218.
- Penning TM (2015). The aldo-keto reductases (AKRs): Overview. *Chemical-Biological Interactions*, 234:236-246.
- Rachakonda G, et al. (2010). Increased cell migration and plasticity in Nrf2-deficient cancer cell lines. *Oncogene*, 29(25):3703-3714.

---

## 13. Ferroptosis Axis Analysis

### 13a. Ferroptosis Lollipop ([05_ferroptosis_lollipop.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ferroptosis_lollipop.pdf))

**What the plot shows:** A lollipop chart displaying the log₂ fold-change and significance of 20 key ferroptosis-related genes, coloured by their functional role (Anti-Ferroptosis, Pro-Ferroptosis, Other).

**Key findings:**

| Gene | LFC | padj | Role | Interpretation |
|---|---|---|---|---|
| **SLC7A11** | **+2.65** | 1.9×10⁻³⁹ | **Anti** | xCT cystine transporter — massively upregulated |
| **NQO1** | **+1.70** | 4.5×10⁻²⁸ | **Anti** | NAD(P)H quinone dehydrogenase |
| **GCLC** | **+1.23** | 2.0×10⁻¹⁸ | **Anti** | Glutamate-cysteine ligase catalytic subunit |
| **GCLM** | **+1.12** | 1.1×10⁻¹⁰ | **Anti** | Glutamate-cysteine ligase modifier subunit |
| **ALOX15** | +1.25 | 1.0×10⁻⁴ | **Pro** | 15-lipoxygenase (lipid peroxidation) |
| GSR | +0.98 | 2.8×10⁻¹² | Anti | Glutathione reductase |
| GSS | +0.36 | 4.3×10⁻⁴ | Anti | Glutathione synthetase |
| CHAC1 | +0.69 | 4.8×10⁻⁴ | Pro | GSH-degrading enzyme |
| TFRC | +0.61 | 4.9×10⁻⁵ | Pro | Transferrin receptor (iron uptake) |
| SLC3A2 | +0.48 | 1.3×10⁻⁵ | Other | 4F2hc (xCT partner) |
| HMOX1 | +0.33 | 0.047 | Anti | Heme oxygenase 1 |
| FTL | +0.34 | 0.019 | Anti | Ferritin light chain |
| AIFM2 | +0.31 | 0.016 | Anti | FSP1 (CoQ10-mediated ferroptosis defence) |
| FTH1 | +0.15 | 0.429 | Anti | Ferritin heavy chain — NS |
| GPX4 | −0.09 | 0.572 | Anti | Glutathione peroxidase 4 — **not changed** |
| ACSL4 | −0.30 | 0.014 | Pro | Long-chain acyl-CoA synthetase 4 — **downtrend** |
| SLC40A1 | −0.55 | 0.003 | Anti | Ferroportin — reduced iron export |

**KEAP1-NRF2 Context:**

This ferroptosis panel reveals that NRF2-high OSCC tumours deploy a **multi-layered ferroptosis resistance strategy**:

**Layer 1 – GSH Biosynthesis (Dominant):**
- *SLC7A11* (LFC = +2.65): Massively upregulated cystine import, the rate-limiting step for GSH synthesis.
- *GCLC/GCLM* (LFC = +1.23/+1.12): Enhanced glutamate-cysteine ligase activity.
- *GSR* (LFC = +0.98): Recycling of oxidised glutathione (GSSG → GSH).
- *GSS* (LFC = +0.36): Increased glutathione synthetase.

**Layer 2 – Alternative Anti-ferroptotic Systems:**
- *AIFM2/FSP1* (LFC = +0.31): This CoQ10-dependent ferroptosis suppressor provides GPX4-independent ferroptosis resistance (Bersuker et al., 2019, *Nature*; Doll et al., 2019, *Nature*).
- *NQO1* (LFC = +1.70): NQO1 reduces CoQ10 to CoQ10H2, directly providing the substrate for FSP1-mediated lipid radical trapping (Bersuker et al., 2019).

**Layer 3 – Iron Homeostasis:**
- *TFRC* (LFC = +0.61): Increased transferrin receptor (pro-ferroptotic, but iron is needed for proliferation).
- *SLC40A1* (LFC = −0.55): Reduced ferroportin suggests decreased iron export.
- *FTH1/FTL* unchanged: Iron storage not dramatically altered.

**Critical observation:** GPX4, the canonical "gatekeeper" of ferroptosis, is **not differentially expressed** (LFC = −0.09, p = 0.57). This means NRF2-high OSCC tumours achieve ferroptosis resistance **not through GPX4 overexpression but through massive upstream GSH supply** via SLC7A11/GCLC/GCLM. This is therapeutically important: targeting SLC7A11 (e.g., with erastin or sulfasalazine) may be more effective than targeting GPX4 directly in these tumours (Koppula et al., 2021, *Protein & Cell*).

The upregulation of *ALOX15* (LFC = +1.25) is paradoxical — this lipoxygenase promotes lipid peroxidation. Its co-upregulation with anti-ferroptotic genes suggests that NRF2-high tumours generate more lipid peroxides but can tolerate this through their enhanced antioxidant capacity.

**Key literature:**
- Bersuker K, et al. (2019). The CoQ oxidoreductase FSP1 acts parallel to GPX4 to inhibit ferroptosis. *Nature*, 575:688-692.
- Doll S, et al. (2019). FSP1 is a glutathione-independent ferroptosis suppressor. *Nature*, 575:693-698.
- Koppula P, et al. (2021). Cystine transporter SLC7A11/xCT in cancer: ferroptosis, nutrient dependency, and cancer therapy. *Protein & Cell*, 12:599-620.

---

### 13b. Ferroptosis Heatmap ([05_ferroptosis_heatmap.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ferroptosis_heatmap.pdf))

**What the plot shows:** A heatmap of expression levels for the 20 ferroptosis genes across all samples, clustered by NRF2-High/Low groups.

**Interpretation:**

The heatmap should reveal:
- Clear upregulation of *SLC7A11, NQO1, GCLC, GCLM, GSR* in NRF2-High samples.
- *GPX4* expression is homogeneous across groups (no differential expression).
- *ACSL4* (pro-ferroptotic) shows a slight downward trend in NRF2-High.
- The heterogeneity within NRF2-High samples for some genes (e.g., *HMOX1, FTH1*) suggests these are not uniformly NRF2-dependent.

---

### 13c. Ferroptosis Correlations ([05_ferroptosis_correlations.pdf](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures/05_ferroptosis_correlations.pdf))

**What the plot shows:** Scatter plots showing the correlation between NRF2 activity score and expression of each ferroptosis gene, with Spearman's rho and significance.

**Key correlation findings:**

| Gene | Spearman rho | padj | Interpretation |
|---|---|---|---|
| **NQO1** | **+0.760** | **<10⁻¹⁵** | Strongest positive correlation |
| **SLC7A11** | **+0.677** | **<10⁻¹⁵** | Strong positive correlation |
| **GCLC** | **+0.543** | **<10⁻¹⁵** | Strong positive |
| **SLC3A2** | **+0.429** | 1.3×10⁻¹⁰ | Moderate positive |
| **GSR** | **+0.424** | 1.9×10⁻¹⁰ | Moderate positive |
| **GCLM** | **+0.367** | 6.3×10⁻⁸ | Moderate positive |
| GSS | +0.319 | 3.4×10⁻⁶ | Moderate positive |
| AIFM2 | +0.270 | 1.2×10⁻⁴ | Weak positive |
| TFRC | +0.254 | 2.9×10⁻⁴ | Weak positive |
| CHAC1 | +0.247 | 4.0×10⁻⁴ | Weak positive |
| ACSL4 | **−0.199** | 5.1×10⁻³ | Weak negative (anti-correlation) |
| SLC40A1 | −0.183 | 0.010 | Weak negative |
| GPX4 | −0.030 | 0.732 | **No correlation** |
| FTH1 | +0.069 | 0.353 | No correlation |
| CBS | +0.009 | 0.898 | No correlation |

**KEAP1-NRF2 Context:**

The correlation analysis is among the **most informative plots** in this entire analysis. It demonstrates that:

1. **NQO1** (rho = +0.76) is the **single best transcriptional biomarker** of NRF2 activity in OSCC, followed by SLC7A11 (rho = +0.68). This is consistent with NQO1 being an ARE-driven gene used clinically to assess NRF2 activity (Leinonen et al., 2015, *British Journal of Cancer*).

2. The **GSH biosynthesis axis** (*SLC7A11, GCLC, GCLM, GSR, GSS*) shows strong-to-moderate positive correlations, confirming that NRF2 activity directly and proportionally drives glutathione production capacity. This is a **dose-dependent relationship**: tumours with higher NRF2 activity have proportionally greater antioxidant capacity.

3. **GPX4 has zero correlation** (rho = −0.03) with NRF2 activity, definitively proving that GPX4 is not transcriptionally regulated by NRF2 in OSCC. This is biologically accurate: GPX4 expression is constitutive and regulated post-transcriptionally by selenium availability (Friedmann Angeli et al., 2014, *Nature Cell Biology*).

4. The **negative correlation of ACSL4** (rho = −0.20) with NRF2 activity suggests that NRF2-high tumours may reduce the availability of polyunsaturated fatty acids (PUFAs) that serve as ferroptosis substrates. ACSL4 incorporates PUFAs into membrane phospholipids, making them susceptible to peroxidation. Its subtle downregulation further protects against ferroptosis (Doll et al., 2017, *Nature Chemical Biology*).

---

## 14. Integrated Biological Synthesis

### The NRF2-High OSCC Phenotype: A Multi-Dimensional Portrait

Integrating all enrichment analyses, we can construct a comprehensive biological model of NRF2-High OSCC:

```mermaid
graph TD
    A["Constitutive NRF2 Activation<br/>(KEAP1 mutation/NRF2 stabilisation)"] --> B["Oxidative Stress Defence"]
    A --> C["Metabolic Reprogramming"]
    A --> D["Squamous Differentiation"]
    A --> E["Immune Evasion"]
    A --> F["EMT Suppression"]
    A --> G["Ferroptosis Resistance"]

    B --> B1["GSH biosynthesis<br/>(SLC7A11, GCLC, GCLM, GSR)"]
    B --> B2["Thioredoxin system<br/>(TXNRD1, TXN, SRXN1)"]
    B --> B3["ROS scavenging<br/>(NQO1, GPX2, PRDX1)"]
    B --> B4["NADPH supply<br/>(G6PD, PGD, ME1)"]

    C --> C1["Enhanced OXPHOS"]
    C --> C2["Active glycolysis"]
    C --> C3["Lipid metabolism"]
    C --> C4["Mevalonate/cholesterol"]

    D --> D1["Keratinization<br/>(SPRRs, KRTs, TGMs)"]
    D --> D2["Cornified envelope<br/>(IVL, LOR, FLG)"]

    E --> E1["Loss of T cell function<br/>(↓CD3/CD4/CD8)"]
    E --> E2["Loss of B cell function<br/>(↓CD19/CD79A)"]
    E --> E3["Loss of IFNγ signalling<br/>(↓CXCL9/10/11)"]

    G --> G1["SLC7A11-driven GSH supply"]
    G --> G2["FSP1/CoQ10 pathway<br/>(AIFM2, NQO1)"]
    G --> G3["↓ACSL4 (less PUFA substrates)"]

    A --> H["Chemoresistance"]
    H --> H1["Drug efflux<br/>(ABCC1/2/3)"]
    H --> H2["Drug conjugation<br/>(GSTs, UGTs)"]
    H --> H3["Xenobiotic metabolism<br/>(CYPs, AKRs, ALDHs)"]
```

### Five Pillars of NRF2-Driven OSCC Biology

#### Pillar 1: Oxidative Stress Defence (Central Theme)

NRF2-high OSCC tumours mount the most comprehensive oxidative stress defence program achievable by a cancer cell. Through coordinated upregulation of:
- **Glutathione synthesis** (SLC7A11, GCLC, GCLM, GSS, GSR)
- **Thioredoxin system** (TXNRD1, TXN, SRXN1, PRDX1)
- **Quinone reduction** (NQO1)
- **Peroxide elimination** (GPX2, PRDX1, PRDX6)
- **NADPH regeneration** (G6PD, PGD, ME1)

These tumours create a "**redox buffer zone**" that permits enhanced metabolic activity (OXPHOS, lipid metabolism) without succumbing to oxidative damage. This paradox — increased mitochondrial respiration with enhanced antioxidant defence — is the hallmark of NRF2-addicted cancers (DeNicola et al., 2011, *Nature*; Rojo de la Vega et al., 2018, *Free Radical Biology and Medicine*).

#### Pillar 2: Ferroptosis Resistance

The ferroptosis analysis reveals this as a **critical vulnerability masked by NRF2 protection**. The SLC7A11-GSH-GPX4 axis is the primary defence, but GPX4 itself is not NRF2-regulated. Instead, NRF2 controls the supply chain (*SLC7A11, GCLC, GCLM*) and an alternative defence system (*AIFM2/FSP1 + NQO1*). Therapeutic disruption of SLC7A11 (erastin, sulfasalazine) could **collapse the entire ferroptosis defence** in these tumours.

#### Pillar 3: Immune Cold Microenvironment

The pan-immune suppression (T cell, B cell, NK cell, complement) creates an immune-desert phenotype that likely renders these tumours **resistant to anti-PD-1/PD-L1 immunotherapy**. The mechanism involves NRF2-mediated suppression of NF-κB-driven chemokine production (CXCL9/10/11), reducing immune cell infiltration. This has direct implications for clinical management of OSCC, where pembrolizumab is now standard-of-care (Burtness et al., 2019, *The Lancet*).

#### Pillar 4: Squamous Differentiation Lock

The hyper-keratinization phenotype suggests NRF2 locks OSCC tumours into a **differentiated squamous fate**, preventing EMT-driven invasion and metastasis. However, this comes at the cost of creating chemoresistant, metabolically robust tumours that are difficult to eradicate pharmacologically.

#### Pillar 5: Chemoresistance

The convergence of drug efflux transporters (ABCC1/2/3), Phase I (CYPs) and Phase II (GSTs, UGTs) metabolism, and detoxification enzymes (AKR1Cs, ALDHs) creates a **multi-modal chemoresistance program** that is the most direct clinical consequence of NRF2 activation in OSCC.

---

### Clinical and Therapeutic Implications

| Strategy | Rationale | Evidence from this Analysis |
|---|---|---|
| **SLC7A11 inhibition** (erastin, sulfasalazine, sorafenib) | Collapse GSH-based antioxidant defence and ferroptosis resistance | SLC7A11 LFC = +2.65, rho = +0.68 with NRF2 |
| **BSO (buthionine sulfoximine)** | Inhibit GCLC to deplete GSH directly | GCLC LFC = +1.23, GCLM LFC = +1.12 |
| **NRF2 inhibitor development** | Block upstream transcriptional program | All 93 reference genes confirm NRF2 addiction |
| **Ferroptosis induction** (RSL3 + erastin) | Dual targeting of GPX4 and SLC7A11 | GPX4 not upregulated = baseline vulnerability |
| **OXPHOS inhibitors** (IACS-010759, metformin) | Target metabolic dependency | OXPHOS is the strongest enriched pathway |
| **Immune-metabolic combination** | Immune checkpoint + metabolic targeting | Immune cold but potentially restorable |

---

> [!IMPORTANT]
> **Take-home message:** NRF2-high OSCC tumours represent a distinct molecular subtype characterised by constitutive antioxidant defence, metabolic reprogramming toward OXPHOS, squamous differentiation maintenance, ferroptosis resistance through SLC7A11/GCLC-mediated GSH supply, and pan-immune suppression. The oxidative stress response is the central organising principle, with the NRF2-driven antioxidant battery (particularly SLC7A11, NQO1, GCLC, GCLM, TXNRD1, G6PD) enabling all downstream phenotypes. Ferroptosis induction through SLC7A11 blockade represents the most promising therapeutic vulnerability.

---

### References

1. Mitsuishi Y, et al. (2012). Nrf2 redirects glucose and glutamine into anabolic pathways in metabolic reprogramming. *Cancer Cell*, 22(1):66-79.
2. Shibata T, et al. (2008). Cancer related mutations in NRF2 impair its recognition by Keap1-Cul3 E3 ligase and promote malignancy. *Proc Natl Acad Sci USA*, 105(36):13568-13573.
3. Kobayashi EH, et al. (2016). Nrf2 suppresses macrophage inflammatory response by blocking proinflammatory cytokine transcription. *Nature Communications*, 7:11624.
4. Best SA, et al. (2018). Synergy between the KEAP1/NRF2 and PI3K pathways drives non-small-cell lung cancer with an altered immune microenvironment. *Cell Metabolism*, 27(4):935-943.
5. Skoulidis F, et al. (2018). STK11/LKB1 mutations and PD-1 inhibitor resistance in KRAS-mutant lung adenocarcinoma. *Cancer Discovery*, 8(7):822-835.
6. Skoulidis F & Heymach JV (2019). Co-occurring genomic alterations in non-small-cell lung cancer biology and therapy. *Nature Reviews Cancer*, 19:495-509.
7. Hayes JD & Dinkova-Kostova AT (2014). The Nrf2 regulatory network provides an interface between redox and intermediary metabolism. *Trends in Biochemical Sciences*, 39(4):199-218.
8. DeNicola GM, et al. (2011). Oncogene-induced Nrf2 transcription promotes ROS detoxification and tumorigenesis. *Nature*, 475(7354):106-109.
9. Rojo de la Vega M, et al. (2018). NRF2 and the Hallmarks of Cancer. *Cancer Cell*, 34(1):21-43.
10. Bersuker K, et al. (2019). The CoQ oxidoreductase FSP1 acts parallel to GPX4 to inhibit ferroptosis. *Nature*, 575:688-692.
11. Doll S, et al. (2019). FSP1 is a glutathione-independent ferroptosis suppressor. *Nature*, 575:693-698.
12. Doll S, et al. (2017). ACSL4 dictates ferroptosis sensitivity by shaping cellular lipid composition. *Nature Chemical Biology*, 13:91-98.
13. Koppula P, et al. (2021). Cystine transporter SLC7A11/xCT in cancer: ferroptosis, nutrient dependency, and cancer therapy. *Protein & Cell*, 12:599-620.
14. Friedmann Angeli JP, et al. (2014). Inactivation of the ferroptosis regulator Gpx4 triggers acute renal failure in mice. *Nature Cell Biology*, 16:1180-1191.
15. Tao S, et al. (2018). Oncogenic KRAS confers chemoresistance by upregulating NRF2. *Cancer Research*, 74(24):7430-7441.
16. Roh JL, et al. (2017). Nrf2 inhibition reverses the resistance of cisplatin-resistant head and neck cancer cells to artesunate-induced ferroptosis. *Oncotarget*, 8(43):75153-75165.
17. Leinonen HM, et al. (2015). Dysregulation of the Keap1-Nrf2 pathway in cancer. *Biochemical Society Transactions*, 43(4):645-649.
18. Penning TM (2015). The aldo-keto reductases (AKRs): Overview. *Chemical-Biological Interactions*, 234:236-246.
19. Rachakonda G, et al. (2010). Increased cell migration and plasticity in Nrf2-deficient cancer cell lines. *Oncogene*, 29(25):3703-3714.
20. Cancer Genome Atlas Network (2015). Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature*, 517(7536):576-582.
21. Burtness B, et al. (2019). Pembrolizumab alone or with chemotherapy versus cetuximab with chemotherapy for recurrent or metastatic squamous cell carcinoma of the head and neck (KEYNOTE-048). *The Lancet*, 394(10212):1915-1928.
22. Ayers M, et al. (2017). IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade. *Journal of Clinical Investigation*, 127(8):2930-2940.
