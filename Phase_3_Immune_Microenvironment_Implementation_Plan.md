# Comprehensive Implementation Plan: Phase 3 (REVISED)
## Immune Microenvironment, Checkpoint Architecture & Immunotherapy Relevance in the KEAP1/NRF2 OSCC Cohort

**Project Title**: KEAP1/NRF2 Axis in Oral Squamous Cell Carcinoma (OSCC) — TCGA-HNSC Cohort  
**Author**: Lead Bioinformatics Researcher & Cancer Systems-Biology Specialist  
**Cohort Scope**: TCGA Oral Cavity Squamous Cell Carcinoma (n = 224 primary tumors; 112 NRF2-High vs. 112 NRF2-Low; 16 matched normal tissues)  
**Governing Document**: [scope_of_work.txt](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/scope_of_work.txt) — Phase 3 (Sections 3.1 & 3.2)  
**Status**: Revised for Study Coherence — KEAP1-NRF2-Oxidative Stress Axis as Central Organising Principle  

---

## 0. Coherence Review Notes

> [!IMPORTANT]
> This version has been revised following a critical coherence audit of the prior plan against the study scope and Phases 1 & 2 findings. Key changes are:
> 1. **Central axis restored**: Every module is now explicitly anchored to the KEAP1-NRF2-oxidative stress axis as the biological driver; Phase 3 is framed as interrogating the *immunological consequences* of constitutive antioxidant activation, not immune biology in isolation.
> 2. **Antioxidant-Immune Bridge Module added** (Module 3.2): A dedicated analytical module directly correlates NRF2 antioxidant targets (NQO1, SLC7A11, GCLC/GCLM, TXNRD1, G6PD) with immune cell infiltration — this is the mechanistic spine linking Phases 1/2 to Phase 3 and is absent from the prior version.
> 3. **HPV data unavailable**: HPV status is not present in the TCGA metadata for this cohort and is therefore omitted from all analyses. Smoking status remains the primary clinico-environmental covariate.
> 4. **ESTIMATE input files clarified**: `estimate_input_matrix.gct`, `estimate_input_matrix.tsv`, and `estimate_filtered.gct` are **archival provenance files** documenting what went in to the ESTIMATE algorithm. The primary analytical input is `estimate_scores.tsv`. These files are NOT new expression matrices to be re-analysed.
> 5. **Smoking framed mechanistically**: The tobacco-ROS-NRF2-immune axis is now explicitly modelled, not treated as a generic clinical covariate.
> 6. **Ferroptosis-immune interface added**: The Phase 2 finding that SLC7A11 and glutathione synthesis are central NRF2 targets is now connected to immune biology (SLC7A11 / xCT suppresses T cell activity via glutamate deprivation and cystine competition).
> 7. **ROS-checkpoint coupling made explicit**: IDO1 is now correctly framed as a ROS-regulated tryptophan catabolism gene, not merely a checkpoint gene.

---

## 1. Executive Summary & Biological Rationales

### 1.1 The Central Thesis

This study established in Phases 1 & 2 that NRF2-High OSCC tumours are characterised by:
- Constitutive antioxidant defence (SLC7A11/xCT, NQO1, GCLC/GCLM, TXNRD1, G6PD) creating a **glutathione-rich, ROS-suppressed intracellular environment**.
- Simultaneous suppression of all major immune signalling pathways (GSEA: Allograft Rejection NES = −2.46; IFNγ NES = −1.93; Inflammatory Response NES = −1.88).
- Significant EMT suppression and squamous differentiation lock mediated by a **novel multi-transglutaminase signature** (TGM1/3/5/6/7 all upregulated).
- Ferroptosis resistance through SLC7A11-GSH supply chain, with GPX4 paradoxically uncorrelated with NRF2 activity ($\rho = -0.03$) — resistance operates upstream through cystine import and GSH synthesis rather than at the GPX4 effector level.
- A **"primed instability" redox state**: three lipoxygenase family members (ALOX12, ALOX12B, ALOX15) are co-upregulated alongside the antioxidant programme, meaning NRF2-High tumours generate abnormally high lipid peroxide loads while simultaneously neutralising them.
- **NQO1** as the sole independently prognostic canonical NRF2 target gene (aHR = 1.16, p = 0.035; C-index 0.635) and **CACNA1A** loss as the single best discriminating non-canonical biomarker (HR = 0.58, C-index 0.674).

**Phase 3 tests the hypothesis** that the redox microenvironment created by constitutive NRF2 activity is mechanistically responsible for immune exclusion. The plan below interrogates the *transcriptional and deconvolution evidence* for six causal mechanisms:

1. **Redox-Driven Chemokine Suppression → Immune Exclusion**: NRF2 suppresses NF-κB transcription of CXCL9, CXCL10, CXCL11 — confirmed by our Phase 2 GSEA (Kobayashi et al., 2016). These chemokines are required for CD8+ T-cell homing. Their absence in NRF2-High tumours explains the immune-desert ESTIMATE findings ($P = 3.49 \times 10^{-8}$).

2. **SLC7A11 / xCT Antioxidant Exporter → T-Cell Suppression**: The NRF2 target SLC7A11 (xCT), the most strongly upregulated NRF2 gene in our cohort (LFC = +2.65), exports glutamate while importing cystine. In the TME, this depletes extracellular cysteine, which T cells require for activation and proliferation. This is a direct, mechanistically demonstrable link between the NRF2 antioxidant programme and immune suppression (Jiang et al., 2021, *Nature*).

3. **ROS-Suppression → Loss of Immunogenic Cell Death Signalling**: Endogenous ROS — suppressed by NRF2 — is required for activation of inflammasomes, HMGB1 release, calreticulin exposure, and immunogenic cell death (ICD). NRF2-mediated ROS quenching may ablate these danger signals, preventing dendritic cell maturation and antigen cross-presentation.

4. **TNFSF18/GITRL Paradox → Treg Expansion in an Immune-Cold TME** *(Novel Phase 2 Finding)*: TNFSF18 (GITRL, the GITR co-stimulatory ligand) is among the top 30 most significantly upregulated DEGs in NRF2-High tumours (LFC = +2.46, padj = 7.2 × 10⁻²⁰) — paradoxically, this T-cell co-stimulatory molecule is massively elevated in the most immunologically cold tumours. This has not been reported in any cancer type and may represent tumour-intrinsic GITRL expression that preferentially expands Tregs (not effector T cells), reinforcing immunosuppression while maintaining the appearance of immune-active signalling.

5. **ALOX12/12B/15 "Primed Instability" → Conditional Immune-Activating Signal** *(Novel Phase 2 Finding)*: ALOX12 (LFC = +1.44), ALOX12B (LFC = +1.57), and ALOX15 (LFC = +1.25) are co-upregulated alongside the NRF2 antioxidant programme, creating a metabolic tightrope where the tumour generates elevated lipid peroxides but counterbalances them with massive antioxidant capacity. Lipid peroxidation products (4-HNE, 15-HETE, PGE2) are known to suppress NK cell function and polarise macrophages toward M2 — this "primed instability" state may actively shape the immunosuppressive TME.

6. **Multi-TGM Differentiation Lock → Stromal Immune Interface** *(Novel Phase 2 Finding)*: Five transglutaminase family members (TGM1/3/5/6/7) are co-upregulated in NRF2-High tumours alongside the antioxidant programme, establishing a cornification-differentiation lock. TGM2 (which is independently prognostic, HR = 1.20, p = 0.005) catalyses extracellular matrix crosslinking that stiffens the stromal barrier — which may physically impede immune cell infiltration and contribute to the T-cell excluded phenotype.

### 1.2 Analytic Strategy

Phase 3 proceeds in logical layered order:
- **Layer 1** (Module 3.1): Characterise *what* the immune microenvironment looks like in NRF2-High vs Low (CIBERSORT + ESTIMATE).
- **Layer 2** (Module 3.2): Establish *why* — the mechanistic antioxidant-immune bridge correlating NRF2 targets with immune infiltrates.
- **Layer 3** (Module 3.3): Profile immune checkpoints and immunomodulators in the redox context.
- **Layer 4** (Module 3.4): Quantify immunotherapy resistance through TIDE and emerging biomarkers.
- **Layer 5** (Module 3.5): Integrate immune metrics with survival to determine clinical relevance.
- **Layer 6** (Module 3.6): Master synthesis and consensus model.

---

## 2. Ingestion & Validation of Existing Data Assets

### 2.1 Primary Analytical Inputs

The workspace directory `05_results/cibsersort_data/` contains the following verified pre-processed assets:

| File | Role | Format | Samples |
| :--- | :--- | :--- | :--- |
| `cibersortx_results.csv` | **Primary immune deconvolution input** | CSV | 224 tumors, 22 LM22 columns + QC metrics |
| `estimate_scores.tsv` | **Primary stromal/immune burden input** | TSV | 224 tumors: StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity |
| `01_data/processed/vst_normalised.tsv` | **Gene expression for all molecular analyses** | TSV | 240 samples (224 tumor + 16 normal), all genes |
| `01_data/processed/metadata_nrf2_classified.tsv` | **Clinical + NRF2 classification metadata** | TSV | 240 samples with NRF2_group, NRF2_GSVA_score, clinical covariates |

### 2.2 Archival Provenance Files (Not Re-analysed)

The following files document the inputs used to run the ESTIMATE algorithm for the previous project. They are preserved for reproducibility but are **not used as analytical inputs** in this project — the scores have already been derived:
- `estimate_filtered.gct` (9,370 genes × 224 samples GCT 1.2): Input expression matrix sent to ESTIMATE.
- `estimate_input_matrix.gct` and `estimate_input_matrix.tsv`: Pre-filtered expression matrices used to construct the ESTIMATE GCT file.
- `cibersortx_input_matrix.tsv`: Expression matrix submitted to CIBERSORTx.

### 2.3 Sample Concordance Audit
- **CIBERSORTx**: 224/224 tumor `File_ID`s matched. **112 NRF2-High, 112 NRF2-Low. Zero missing.**
- **ESTIMATE**: 224/224 tumor `File_ID`s matched. **Zero missing.**
- **HPV Status**: Confirmed absent from `metadata_nrf2_classified.tsv`. HPV will not be used as a covariate. The primary clinico-environmental covariates available are: Age, Stage, Grade, and Smoking status.

### 2.4 Verified ESTIMATE Preliminary Findings

| Metric | NRF2-High Median | NRF2-Low Median | Wilcoxon $P$ | Spearman $\rho$ | Interpretation |
| :--- | :---: | :---: | :---: | :---: | :--- |
| `StromalScore` | -624.2 | -213.6 | $5.75 \times 10^{-6}$ | −0.316 | Stroma-depleted in NRF2-High |
| `ImmuneScore` | -65.8 | +419.5 | $1.02 \times 10^{-6}$ | −0.354 | Immune-desert in NRF2-High |
| `ESTIMATEScore` | -568.0 | +214.6 | $3.49 \times 10^{-8}$ | −0.393 | Global TME desertification |
| `TumorPurity` | 0.867 | 0.804 | $3.49 \times 10^{-8}$ | +0.393 | Higher purity in NRF2-High |

*These findings independently validate that NRF2-High OSCC represents an immune-cold, stroma-depleted entity — consistent with Phase 2 GSEA showing suppression of Allograft Rejection (NES = −2.46), IFNγ (NES = −1.93), and Inflammatory Response (NES = −1.88) pathways.*

### 2.5 CAF Coverage Plan

CIBERSORTx LM22 is hematopoietic-only and does not deconvolve CAFs. Coverage is achieved through three complementary approaches:
1. **ESTIMATE StromalScore** as a validated bulk proxy for stromal content.
2. **Canonical CAF marker expression** from VST matrix: *ACTA2*, *FAP*, *COL1A1*, *COL1A2*, *PDGFRA*, *PDGFRB*, *TGFB1*, *POSTN*. Note: Phase 2 GSEA showed *COL1A1*, *COL3A1*, *FAP*, *POSTN* are among core EMT-suppressed genes — this is a key contextual finding for CAF interpretation.
3. **TIDE CAF exclusion score** from Module 3.4.

---

## 3. Statistical Testing & Multiple Testing Philosophy

**Governing principle**: "Wherever it is not necessary, do not do unnecessary FDR correction." Pre-specified, biologically motivated hypotheses are tested without multiple testing correction. Exploratory screens receive BH-FDR.

### 3.1 Targeted Hypotheses — Raw P-values, Effect Sizes & 95% CI

**Covariates in all group comparisons**: NRF2 group (High vs Low). Available clinical covariates are Age, Stage, Grade, and Smoking status. Smoking status receives a dedicated targeted interaction analysis (Module 3.1, analysis 5). HPV data are unavailable and are excluded.

| Hypothesis Domain | Variables Tested | Statistical Test |
| :--- | :--- | :--- |
| **Primary immune infiltration** | CD8+ T cells, M1, M2, M1/M2 ratio, Tregs, CD8/Treg ratio, NK cells (activated + resting), Dendritic cells | Two-sided Wilcoxon rank-sum, Cliff's $\delta$, HL median diff 95% CI |
| **ESTIMATE microenvironmental burden** | StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity | Wilcoxon, Cliff's $\delta$, Spearman $\rho$ with NRF2 GSVA |
| **Antioxidant-immune bridge** | Spearman $\rho$ between each NRF2 antioxidant target (NQO1, SLC7A11, GCLC, GCLM, TXNRD1, G6PD, GSR) and each primary immune population | Spearman correlation, raw $P$ |
| **Canonical checkpoints** | PDCD1, CD274, PDCD1LG2, CTLA4, LAG3, HAVCR2, TIGIT, IDO1 | Wilcoxon, Spearman $\rho$ |
| **TIDE core metrics** | TIDE Score, T-cell Dysfunction, T-cell Exclusion | Wilcoxon, Cliff's $\delta$ |

### 3.2 Exploratory Wide Screens — BH-FDR Required

| Screen | n Tests (approx.) | Rationale for FDR |
| :--- | :--- | :--- |
| Full 22-cell LM22 screen | 22 | All 22 simultaneously; not all pre-specified |
| Extended co-inhibitory / co-stimulatory receptor screen | ~40 | Broad family-level |
| Cytokine & chemokine repertoire | ~50 | Exploratory |
| APM gene screen | ~12 | Broad family |
| Extended NRF2 target – immune correlation matrix | ~140 pairwise | Genome-wide correlations |

---

## 4. Modular Analytical Architecture

```
Phase 3 Pipeline — Anchored to KEAP1-NRF2-Oxidative Stress Axis:
│
├── Module 3.1: Immune Microenvironment Deconvolution
│     CIBERSORT LM22 + ESTIMATE + CAF markers
│     Grade-stratified immune landscape (Grade Paradox)
│     → What does the immune landscape look like?
│
├── Module 3.2: Antioxidant-Immune Bridge (KEY MODULE)
│     NRF2 antioxidant targets × immune infiltration correlations
│     SLC7A11/xCT-T-cell suppression axis
│     ROS-suppression & immunogenic cell death signalling
│     TNFSF18/GITRL paradox — Treg expansion mechanism
│     ALOX12/15 primed instability — lipid peroxide immune axis
│     Multi-TGM lock — ECM-mediated immune exclusion
│     → WHY is the landscape immunologically cold?
│
├── Module 3.3: Immune Checkpoint & Immunomodulatory Architecture
│     Redox-checkpoint coupling (IDO1 as tryptophan-ROS nexus)
│     Extended immunomodulator screen (TNFSF18 specifically flagged)
│     Chemokine trafficking (CXCL9/10/11 as NF-κB/NRF2 targets)
│     Antigen processing & presentation machinery
│
├── Module 3.4: TIDE & Immunotherapy Resistance Prediction
│     TIDE Score / Dysfunction / Exclusion
│     CYT, IFNγ GEP, Immunophenoscore
│     ICB response prediction
│
├── Module 3.5: Immune-Prognostic Integration
│     Bivariate NRF2 × immune survival stratification
│     Non-canonical prognostic DEG × immune correlates (CACNA1A, PITX2, GATA3)
│     Multivariable Cox models with immune + clinical + purity covariates
│
├── Module 3.6: Master TME Synthesis & Redox-Immune Model
│     Integrative consensus figure
│     Mechanistic model: ROS suppression → immune desert
│
└── Module 3.7: Novel Angles — Cross-Phase Integration (NEW)
      TNFSF18/GITRL-Treg paradox quantification
      ALOX-primed instability & immune interface
      Grade Paradox immune stratification
      CACNA1A / PITX2 / GATA3 immune landscape correlates
      Smoking-NRF2 immune dichotomy
```

### Module 3.1: Immune Microenvironment Deconvolution

1. **Cross-Platform Concordance**: Spearman correlation between ESTIMATE `ImmuneScore` and total CIBERSORTx leukocyte fraction across all 224 tumors (validates both methods).
2. **Lineage Aggregation**: Aggregate 22 LM22 subsets into macro-lineages; calculate M1/M2 and CD8/Treg ratios.
3. **NRF2-Stratified Comparison**: NRF2-High vs NRF2-Low for all primary immune candidates (raw $P$, Cliff's $\delta$).
4. **Tertile Dose-Response Analysis**: Compare across NRF2 Low/Mid/High tertiles to assess dose-dependent immune remodeling.
5. **Smoking-Stratified Immune Interaction** (targeted, no FDR): Evaluate immune infiltration (CD8+ T cells, M1/M2 ratio, Tregs) across smoking status (Current / Reformed / Never) within NRF2 groups — framed around the tobacco-ROS-NRF2-immune axis.
6. **ESTIMATE Analysis**: Ingest `estimate_scores.tsv`; compare StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity across NRF2 groups.
7. **CAF Profiling**: Correlate StromalScore and NRF2 GSVA score with canonical CAF marker expression from `vst_normalised.tsv`. *Contextual note*: Phase 2 showed *FAP*, *COL1A1*, *POSTN*, *VIM* are downregulated in NRF2-High tumors (EMT suppression), suggesting CAF depletion rather than CAF enrichment — this hypothesis will be tested.
8. **Histologic Grade Paradox — Immune Stratification** *(Novel angle)*: Phase 1 survival analysis revealed that NRF2 activity is *inversely* correlated with histologic grade (G1 median NRF2 GSVA = +0.165; G3/G4 median = −0.324; Kruskal-Wallis p < 0.001). This means well-differentiated NRF2-High tumours are both antioxidant-armed and immune-cold. Test whether the immune-desert phenotype (ImmuneScore, CD8+ T cells) follows the same grade-NRF2 coupling — i.e., G1/NRF2-High tumours should have the coldest TME despite appearing clinically "indolent" histologically.

### Module 3.2: Antioxidant-Immune Bridge (NEW — Critical Module)

This module is the mechanistic heart of Phase 3 — it provides the causal link between the oxidative stress biology established in Phases 1/2 and the immune landscape in Phase 3.

1. **NRF2 Antioxidant Battery × Immune Cell Correlation Matrix**:
   - For each of the 7 core antioxidant targets (NQO1, SLC7A11, GCLC, GCLM, TXNRD1, G6PD, GSR), compute Spearman $\rho$ with all 22 LM22 cell fractions and 4 ESTIMATE scores.
   - These genes have pre-established Spearman correlations with NRF2 GSVA score from Phase 2 (NQO1 $\rho = +0.76$; SLC7A11 $\rho = +0.68$), making this a controlled, sequentially logical analysis.

2. **SLC7A11 / xCT T-Cell Suppression Axis**:
   - SLC7A11 exports glutamate and imports cystine, depleting extracellular cysteine and glutamate in the TME.
   - Cysteine is essential for T-cell activation (Jiang et al., 2021, *Nature*).
   - Analysis: Scatter plots of *SLC7A11* expression vs CD8+ T cells, Tregs, and cytolytic activity (CYT: GZMA/PRF1); test whether *SLC7A11*-high tumors have lower CD8+ and higher Treg infiltration.

3. **ROS-Suppression & Immunogenic Cell Death (ICD) Signalling**:
   - ROS is required for HMGB1 release, calreticulin surface exposure, and ATP secretion — the three core ICD danger signals.
   - NRF2 suppresses ROS via TXNRD1, GPX2, NQO1, and GSR.
   - Analysis: Correlate NRF2 GSVA score with *HMGB1*, *CALR*, and pro-inflammatory DAMP-sensing receptors (*TLR4*, *AGER*/RAGE, *NLRP3*).

4. **Glutathione-Depletion Immunotherapy Hypothesis**:
   - Assess whether tumors with highest *GCLC* + *GCLM* expression (greatest GSH synthesis capacity) show lowest immune infiltration.

5. **TNFSF18 / GITRL Paradox — Treg Expansion in an Immune-Cold TME** *(Novel — first in any cancer type)*:
   - TNFSF18 (GITRL, the ligand for the GITR Treg co-stimulatory receptor) is the top significantly upregulated co-stimulatory DEG in NRF2-High tumours (LFC = +2.46, padj = 7.2 × 10⁻²⁰), yet these tumours are profoundly immune-cold.
   - **Hypothesis**: Tumour-intrinsic GITRL preferentially activates and expands Tregs (not effector T cells) through GITR ligation, actively enforcing immunosuppression — a previously undescribed NRF2-GITRL-Treg axis.
   - Analysis: Correlate *TNFSF18* expression with CIBERSORTx Treg fraction and CD8+ T cells; compare NRF2-High vs NRF2-Low; test whether *TNFSF18* high tumours have highest Treg/CD8 ratios. Context: this was listed as a specifically targeted gene in Module 3.3 extended screen (*TNFRSF18*) — flag TNFSF18 as a priority result.

6. **ALOX12/12B/15 "Primed Instability" — Lipid Peroxide Immune Axis** *(Novel — multi-ALOX pattern not reported in OSCC)*:
   - Phase 2 DEG analysis identified simultaneous upregulation of *ALOX15* (LFC = +1.25), *ALOX12* (LFC = +1.44), and *ALOX12B* (LFC = +1.57) in NRF2-High tumours — pro-ferroptotic lipoxygenases co-elevated with the antioxidant programme.
   - **Immune relevance**: ALOX15-derived lipid mediators (15-HETE, lipoxins) and ALOX12 products are known to suppress NK cell cytotoxicity, promote M2 macrophage polarization, and generate PGE2-related immunosuppressive metabolites.
   - Analysis: Correlate the composite ALOX12+ALOX12B+ALOX15 expression score with NK cell (resting + activated), M1/M2 ratio, and ESTIMATE ImmuneScore. Test whether high-ALOX / high-NRF2 tumours show the deepest immune suppression, consistent with lipid peroxide-driven immunosuppression rather than pure ROS quenching.

7. **Multi-TGM Differentiation Lock — ECM Crosslinking and Immune Exclusion** *(Novel — 5-TGM signature not reported in any cancer)*:
   - Phase 2 identified 5 transglutaminase family members upregulated in NRF2-High tumours (TGM1/3/5/6/7). TGM2 is independently prognostic (Phase 1: HR = 1.20, p = 0.005; multivariable aHR = 1.17, p = 0.024) and catalyses ECM crosslinking and matrix stiffening.
   - **Immune relevance**: ECM stiffening via TGM2-mediated crosslinking is a known physical barrier to T-cell infiltration (Thomas & Bhatt, 2022). Stiffened collagen matrices trap T cells in peritumoral regions, preventing intratumoural penetration.
   - Analysis: Correlate *TGM2* expression with ESTIMATE StromalScore, TumorPurity, and CD8+ T-cell infiltration. Test the composite TGM score (average of TGM1/2/3/5/7) vs TIDE Exclusion Score — does the multi-TGM signature predict mechanical immune exclusion?

### Module 3.3: Immune Checkpoint & Immunomodulatory Architecture

1. **Targeted Checkpoint Profiling** (no FDR):
   - Extract VST expression for: *PDCD1* (PD-1), *CD274* (PD-L1), *PDCD1LG2* (PD-L2), *CTLA4*, *LAG3*, *HAVCR2* (TIM-3), *TIGIT*, *IDO1*.
   - **Special note on IDO1**: IDO1 encodes indoleamine 2,3-dioxygenase, which converts tryptophan to kynurenine. IDO1 is transcriptionally activated by IFNγ and is partially regulated by ROS/oxidative stress. Given that NRF2-High tumors suppress IFNγ signaling (Phase 2 GSEA: NES = −1.93) and suppress ROS, IDO1 should be interpreted in this context rather than treated as a generic checkpoint.
   - Wilcoxon comparisons and Spearman $\rho$ vs NRF2 GSVA score.

2. **Extended Immunomodulator Screen** (BH-FDR):
   - Co-inhibitory: *BTLA*, *CD96*, *CD160*, *CD244*, *VSIR*, *SIGLEC15*, *LAIR1*.
   - Co-stimulatory: *CD28*, *ICOS*, *TNFRSF4*, *TNFRSF9*, *TNFRSF18*, *CD27*, *CD40LG*, *CD40*, *CD80*, *CD86*.

3. **Chemokine Trafficking — NF-κB / NRF2 Suppression Axis** (targeted, no FDR):
   - T-cell recruiting / NF-κB-driven: *CXCL9*, *CXCL10*, *CXCL11* (downregulated in NRF2-High per Phase 2; confirmed by GSEA IFNγ suppression).
   - Myeloid / MDSC attracting: *CCL2*, *CXCL8* (IL-8), *CCL5*, *CXCL1*, *CXCL2*.
   - Immunosuppressive: *TGFB1*, *IL10*, *VEGFA*, *IL6*.
   - *Mechanistic framing*: CXCL9/10/11 downregulation in NRF2-High is mechanistically attributed to NRF2-mediated NF-κB transcriptional suppression (Kobayashi et al., 2016). This should be explicitly stated in figure legends and the report.

4. **Antigen Processing & Presentation Machinery (APM)** (BH-FDR):
   - MHC Class I: *HLA-A*, *HLA-B*, *HLA-C*, *B2M*.
   - Peptide loading: *TAP1*, *TAP2*, *TAPBP*, *ERAP1*, *ERAP2*.
   - MHC Class II: *HLA-DRA*, *HLA-DRB1*, *HLA-DPA1*, *CIITA*.
   - *Context*: Phase 2 GSEA Allograft Rejection core enrichment confirmed downregulation of *HLA-DRA*, *HLA-DMA*, *B2M*, *TAP1* in NRF2-High tumors — APM analysis validates this at the individual gene level.

### Module 3.4: TIDE & Immunotherapy Resistance Prediction

1. **TIDE Input Matrix Preparation**: Mean-centered log2 expression from `vst_normalised.tsv` (tumor samples only).
2. **TIDE Execution**: TIDE Score, Dysfunction, Exclusion, MDSC, CAF, TAM M2, MSI, Predicted Responder.
3. **Concordance Biomarkers**:
   - **Cytolytic Activity Score (CYT)**: Geometric mean of *GZMA* and *PRF1*.
   - **IFNγ T cell-inflamed GEP**: Ayers et al. 18-gene signature.
   - **Immunophenoscore (IPS)**.
4. **ICB Resistance Hypothesis**: Test whether NRF2-High tumors show significantly higher TIDE Scores (= predicted ICB resistance), consistent with Phase 2 findings of pan-immune suppression.

### Module 3.5: Immune-Prognostic Integration

1. **Bivariate Kaplan-Meier**: NRF2 × CD8+ T-cell combined strata; NRF2 × TIDE Responder status.
2. **Non-canonical Prognostic DEG × Immune Landscape Correlates** *(Novel cross-phase integration)*:
   - Phase 1 identified three non-canonical biomarkers with strong prognostic power in this exact cohort:
     - *CACNA1A* (HR = 0.58, C-index = 0.674; top protective biomarker — loss marks high-risk OSCC)
     - *PITX2* (HR = 1.29, aHR = 1.25, p = 0.007; top adverse transcription factor co-activated with NRF2)
     - *GATA3* (HR = 0.77, aHR = 0.81; epithelial lineage enforcer lost in NRF2-High tumours)
   - These DEGs were identified in the NRF2-High vs NRF2-Low transcriptomic landscape. Their correlation with immune infiltration has never been tested. Analysis: Correlate each with CD8+ T cells, Treg fraction, ImmuneScore, and TIDE Score. Test whether *CACNA1A* low tumours (the high-risk group) also have the coldest immune TME — this would establish a novel multi-level biomarker hierarchy (redox → immune → prognosis).
3. **Smoking-NRF2 Immune Dichotomy** *(Novel — informed by Phase 1 survival finding)*:
   - Phase 1 showed NRF2 is prognostically lethal in **non-smokers** (median OS 4.14 vs >7.5 years; NRF2_High vs Low) but the survival curves completely overlap in smokers (p = 0.569). This implies that NRF2 in non-smokers reflects genuine oncogene addiction (genomic mutations in *NFE2L2/KEAP1*), whereas in smokers NRF2 is a physiological reaction to tobacco ROS.
   - **Immune hypothesis**: Intrinsic NRF2 activation (non-smokers) may create a more deeply immune-excluded TME than exogenous NRF2 induction (smokers), because the antioxidant-immune suppression machinery is genuinely tumour-cell-autonomous rather than environmentally reactive.
   - Analysis: Within non-smokers and within smokers separately, compare ImmuneScore, CD8+ T cells, and Treg fraction between NRF2-High and NRF2-Low. Test whether the immune exclusion magnitude (Cliff's $\delta$) is larger in non-smokers than in smokers.
4. **Multivariable Cox**: NRF2 GSVA score + CD8+ T cells + M2 macrophages + TumorPurity + TIDE Score + Age + Stage + Grade + Smoking. TumorPurity is included because NRF2-High tumors have significantly higher purity ($P = 3.49 \times 10^{-8}$) — a potential confounder for apparent immune differences. HPV data are unavailable and excluded.

### Module 3.6: Master Synthesis

Integrative consensus figure summarizing the complete mechanistic model: **KEAP1 loss / NRF2 stabilization → constitutive antioxidant programme (SLC7A11, NQO1, GCLC/GCLM, TXNRD1) → ROS suppression + cysteine depletion → CXCL9/10/11 downregulation + ICD signal ablation + T-cell starvation → immune-desert TME → ICB resistance + poor immunotherapy response**.

---

## 5. Comprehensive Catalog of Publication-Grade Figures (24 Figures)

**Formats**: Vector PDF + 300 DPI PNG for all figures.  
**Palette**: `#E63946` NRF2-High, `#457B9D` NRF2-Low, `#2A9D8F` Normal; diverging RdBu for correlation heatmaps; viridis for density.

### Batch A: Immune Landscape Characterisation (Module 3.1)

| Fig | File | Plot Type | Question & Key Features |
| :--- | :--- | :--- | :--- |
| **1** | `01_cibersort_landscape_complex_heatmap` | ComplexHeatmap | **Global TME Portrait**: 22 LM22 fractions × 224 tumors ordered by NRF2 GSVA score. Annotations: NRF2 group, NRF2 GSVA score, Smoking status, Stage, Grade, Age, TumorPurity, ESTIMATE ImmuneScore. |
| **2** | `02_targeted_immune_cell_comparisons_grid` | 8-panel Violin + Jitter | **Targeted Infiltration**: CD8+ T, M1, M2, M1/M2 ratio, Tregs, CD8/Treg ratio, NK activated, DC activated. Raw Wilcoxon $P$, Cliff's $\delta$. |
| **3** | `03_cibersort_composition_stacked_barplots` | Stacked Percent Bar | **TME Compositional Shifts**: NRF2-High vs Low and NRF2 tertiles; shows lineage expansion/contraction. |
| **4** | `04_nrf2_immune_continuous_correlation_heatmap` | Correlation Heatmap | **Redox-Immune Dose Response**: Spearman $\rho$ between continuous NRF2 GSVA score and 22 LM22 subsets + 4 ESTIMATE scores. |
| **5** | `05_immune_effect_size_lollipop_plot` | Lollipop / Dot Plot | **Prioritised Effect Sizes**: Cliff's $\delta$ for all 22 LM22 populations + functional ratios, ordered by effect size, colored by raw $P$-value. |
| **6** | `06_estimate_stroma_purity_violin_grid` | 4-Panel Violin + Orthogonal Scatter | **ESTIMATE Desertification**: StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity; raw Wilcoxon $P$, Cliff's $\delta$. + scatter: ESTIMATE ImmuneScore vs CIBERSORTx total leukocyte fraction (cross-platform validation). |
| **7** | `07_caf_markers_stromal_scatter_grid` | Multi-panel Scatter + Regression | **CAF Profiling**: *ACTA2*, *FAP*, *COL1A1*, *PDGFRA*, *TGFB1* expression vs ESTIMATE StromalScore and NRF2 GSVA score; Spearman $\rho$, 95% CI. Note: Phase 2 showed these genes are among EMT-suppressed genes — figures to reflect this directionality. |
| **8** | `08_smoking_immune_ros_interaction_grid` | Grouped Boxplot Grid | **Tobacco-ROS-NRF2-Immune Axis**: CD8+, M1/M2 ratio, Tregs by Smoking History (Current / Reformed / Never) within NRF2 groups. Mechanistic framing: tobacco ROS activates NRF2 → immune exclusion. |

### Batch B: Antioxidant-Immune Bridge (Module 3.2 — NEW)

| Fig | File | Plot Type | Question & Key Features |
| :--- | :--- | :--- | :--- |
| **9** | `09_antioxidant_immune_correlation_matrix` | Annotated Heatmap | **Antioxidant Battery × Immune Infiltration**: 7 antioxidant targets (NQO1, SLC7A11, GCLC, GCLM, TXNRD1, G6PD, GSR) vs 22 LM22 fractions + ESTIMATE metrics. Spearman $\rho$ with significance annotation. |
| **10** | `10_slc7a11_xct_tcel_suppression_scatter` | 4-Panel Scatter + Smooth | **SLC7A11/xCT → T-Cell Deprivation**: SLC7A11 expression vs CD8+ T cells, Tregs, M2 macrophages, CYT (cytolytic activity). Spearman $\rho$, raw $P$. |
| **11** | `11_ros_suppression_icd_damp_signalling` | Multi-panel Violin + Scatter | **ROS Quenching & ICD Loss**: TXNRD1/NQO1 expression vs HMGB1, CALR, TLR4, NLRP3; NRF2 GSVA score vs each DAMP. Tests whether redox suppression ablates danger signalling. |
| **12** | `12_gsh_biosynthesis_immune_desert` | 2-Panel Scatter | **GSH Synthesis Capacity vs Immune Burden**: Composite GCLC+GCLM expression vs ESTIMATE ImmuneScore and CD8+ T cells; tests dose-dependent GSH-immune suppression. |

### Batch C: Checkpoint & Immunomodulatory Architecture (Module 3.3)

| Fig | File | Plot Type | Question & Key Features |
| :--- | :--- | :--- | :--- |
| **13** | `13_canonical_checkpoints_violin_grid` | 8-panel Violin + Box | **Core Immune Checkpoints**: PDCD1, CD274, PDCD1LG2, CTLA4, LAG3, HAVCR2, TIGIT, IDO1. Raw Wilcoxon $P$, Cliff's $\delta$. IDO1 interpreted in IFNγ-suppression context. |
| **14** | `14_checkpoint_nrf2_continuous_scatter_grid` | Multi-panel Scatter + Smooth | **Redox-Checkpoint Coupling**: NRF2 GSVA score vs each canonical checkpoint; Spearman $\rho$, $P$. |
| **15** | `15_extended_immunomodulators_heatmap` | ComplexHeatmap | **Extended Receptor/Ligand Screen**: 35+ co-inhibitory and co-stimulatory genes; BH-FDR significance tags. |
| **16** | `16_chemokine_nfkb_nrf2_suppression_grid` | Bi-directional Lollipop | **NF-κB / NRF2 Chemokine Conflict**: CXCL9/10/11 (NRF2-suppressed T-cell recruiters) vs CCL2/IL8/VEGFA (myeloid recruiters). $\log_2$FC and $-\log_{10}P$. Mechanistically labelled. |
| **17** | `17_antigen_presentation_apm_machinery` | Multi-panel Violin | **APM Downregulation**: HLA-A/B/C, B2M, TAP1/2, TAPBP, HLA-DRA, CIITA. Context: Phase 2 GSEA confirmed HLA-DRA and B2M in Allograft Rejection core enrichment. |

### Batch D: TIDE & ICB Resistance Prediction (Module 3.4)

| Fig | File | Plot Type | Question & Key Features |
| :--- | :--- | :--- | :--- |
| **18** | `18_tide_primary_scores_comparison` | 3-Panel Violin | **TIDE Core Metrics**: TIDE Score, Dysfunction, Exclusion between NRF2-High and Low; Wilcoxon $P$. |
| **19** | `19_tide_exclusion_mediators_breakdown` | 3-Panel Violin | **Exclusion Mechanisms**: MDSC, CAF, TAM M2 exclusion scores; links to ESTIMATE StromalScore (CAF) and CIBERSORT (M2). |
| **20** | `20_tide_icb_responder_contingency` | Stacked Bar + Forest | **ICB Resistance Prediction**: Proportions of Predicted Responders vs Non-Responders; Fisher's exact OR and 95% CI. |
| **21** | `21_cytolytic_gep_immunophenoscore` | 3-Panel Violin + Scatter | **Effector Immune Signatures**: CYT (GZMA/PRF1), Ayers 18-gene IFNγ GEP, IPS across NRF2 groups. |

### Batch E: Prognostic Integration & Master Synthesis (Modules 3.5–3.6)

| Fig | File | Plot Type | Question & Key Features |
| :--- | :--- | :--- | :--- |
| **22** | `22_bivariate_nrf2_immune_km_grid` | Multi-panel Kaplan-Meier | **Immune-Prognostic Synergy**: KM curves for NRF2 × CD8+ T-cell strata and NRF2 × TIDE Responder strata; log-rank $P$, risk tables. |
| **23** | `23_tme_prognostic_multivariate_forest` | Forest Plot | **Multivariable Cox**: NRF2 GSVA + CD8+ T cells + M2 + TumorPurity + TIDE + Age + Stage + Grade + Smoking. HR and 95% CI. |
| **24** | `24_redox_immune_mechanistic_consensus_model` | Integrative Multi-track Panel | **Master Synthesis**: Per-sample multi-track track plot: NRF2 GSVA score → antioxidant targets (SLC7A11, NQO1, GCLC) → CXCL9/10/11 → CD8+ T cells → ESTIMATE ImmuneScore → TIDE Score → Survival risk. Visualises the causal chain across all 224 tumors. |

### Batch F: Novel Cross-Phase Integration Angles (Module 3.7 — NEW)

> [!NOTE]
> These 6 figures arise directly from **novel and significant findings in Phases 1 & 2** not previously present in Phase 3. They represent the primary discovery-oriented contribution of Phase 3.

| Fig | File | Plot Type | Question & Key Features |
| :--- | :--- | :--- | :--- |
| **25** | `25_tnfsf18_gitrl_treg_paradox` | 3-Panel Violin + Scatter + Bar | **TNFSF18/GITRL-Treg Immunosuppression Paradox**: (i) TNFSF18 expression NRF2-High vs Low (violin + Wilcoxon $P$); (ii) TNFSF18 expression vs Treg fraction scatter with Spearman $\rho$; (iii) TNFSF18 vs CD8+ T cells. Tests the novel hypothesis that tumour-intrinsic GITRL enforces Treg expansion despite an immune-cold TME. |
| **26** | `26_alox_primed_instability_immune_axis` | Multi-panel Scatter + Heatmap | **ALOX "Primed Instability" Immune Interface**: (i) Composite ALOX score (ALOX12+ALOX12B+ALOX15) vs NK activated cells, M1/M2 ratio, ESTIMATE ImmuneScore; (ii) ALOX score heatmap overlaid on NRF2 GSVA ordering; (iii) ALOX score vs NRF2 GSVA correlation. Tests whether elevated lipoxygenase products in NRF2-High tumours contribute to NK cell suppression and M2 polarisation. |
| **27** | `27_grade_paradox_immune_landscape` | 3-Panel Violin Grid + Bubble Plot | **Histologic Grade Paradox — Immune Stratification**: Compare ImmuneScore, CD8+ T cells, and TIDE Score across Grade × NRF2 strata (G1/NRF2-High vs G1/NRF2-Low vs G3/NRF2-High etc.). Tests whether G1/NRF2-High tumours — clinically appearing "well-differentiated and indolent" — are actually the most immune-excluded. |
| **28** | `28_noncanonicaldeg_immune_correlates` | Multi-panel Scatter Grid | **Non-Canonical Prognostic DEG × Immune TME**: CACNA1A (HR=0.58), PITX2 (HR=1.29), and GATA3 (HR=0.77) expression vs CD8+ T cells, Treg fraction, ImmuneScore, and TIDE Score. Spearman $\rho$, raw $P$. Tests whether these top survival biomarkers also stratify immune landscape — potential multi-level marker hierarchy (redox → immune → outcome). |
| **29** | `29_tgm_ecm_lock_immune_exclusion` | 3-Panel Scatter + Heatmap | **Multi-TGM ECM Lock and T-Cell Exclusion**: TGM2 expression (independently prognostic, HR=1.20) vs CD8+ T cells, ESTIMATE StromalScore, TIDE Exclusion Score; composite TGM score (TGM1+2+3+5+7) vs TIDE Exclusion; per-sample TGM heatmap ordered by NRF2 GSVA score. Tests whether ECM crosslinking by the multi-TGM programme physically contributes to immune exclusion. |
| **30** | `30_smoking_nrf2_immune_dichotomy` | 2×3 Faceted Boxplot | **Smoking-NRF2 Immune Dichotomy**: Compare ImmuneScore, CD8+ T cells, and Treg fraction between NRF2-High vs NRF2-Low *separately* within non-smokers and smokers. Tests whether intrinsic NRF2 activation (non-smokers) produces deeper immune exclusion than exogenous tobacco-driven NRF2 induction (smokers) — directly testing the Phase 1 survival dichotomy at the immune level. |

---

## 6. Deliverables Pipeline

### 6.1 Scripts (`08_scripts/06_immune_analysis/`)

| Script | Figures Generated | Key Inputs |
| :--- | :--- | :--- |
| `00_tide_prediction.py` | Computational output | `vst_normalised.tsv`, `metadata_nrf2_classified.tsv` (Runs local `tidepy`) |
| `01_cibersort_tme_characterization.R` | Figs 1–8 | `cibersortx_results.csv`, `estimate_scores.tsv`, `metadata_nrf2_classified.tsv` |
| `02_antioxidant_immune_bridge.R` | Figs 9–12 | `vst_normalised.tsv`, `metadata_nrf2_classified.tsv`, `cibersortx_results.csv`, `estimate_scores.tsv` |
| `03_checkpoint_immunomodulators.R` | Figs 13–17 | `vst_normalised.tsv`, `metadata_nrf2_classified.tsv` |
| `04_tide_immunotherapy_prediction.R` | Figs 18–21 | `vst_normalised.tsv`, `metadata_nrf2_classified.tsv`, `tide_predictions_scores.tsv` |
| `05_immune_survival_integration.R` | Figs 22–23 | All above outputs + survival data from `metadata_nrf2_classified.tsv` |
| `06_master_phase3_synthesis.R` | Fig 24 | All above outputs |
| `07_novel_crossphase_angles.R` | Figs 25–30 | `vst_normalised.tsv`, `cibersortx_results.csv`, `estimate_scores.tsv`, `metadata_nrf2_classified.tsv` |

### 6.2 Results (`05_results/06_immune_analysis/`)

- `cibersort_lineage_ratios_summary.tsv`: Per-sample LM22 fractions, aggregate lineages, M1/M2, CD8/Treg ratios.
- `cibersort_statistical_comparison_table.tsv`: All 22 cells + ratios: Median High/Low, Cliff's $\delta$, raw Wilcoxon $P$, BH $q$.
- `estimate_integration_table.tsv`: ESTIMATE scores, cross-platform validation metrics.
- `antioxidant_immune_correlation_table.tsv`: Full 7 × 26 Spearman correlation matrix (antioxidant genes × immune populations/ESTIMATE scores).
- `checkpoint_expression_statistics.tsv`: Canonical and extended checkpoint expression statistics.
- `tide_predictions_scores.tsv`: Complete TIDE output.
- `immune_survival_cox_models.tsv`: Uni- and multivariate Cox results.
- `novel_angles_correlation_table.tsv`: TNFSF18/GITRL-Treg; ALOX composite × immune; TGM composite × TIDE; CACNA1A/PITX2/GATA3 × immune; smoking-stratified immune effect sizes.
- `master_phase3_immune_summary_table.tsv`: Full per-patient integration table.

### 6.2 Results (`05_results/06_immune_analysis/`)

- `cibersort_lineage_ratios_summary.tsv`: Per-sample LM22 fractions, aggregate lineages, M1/M2, CD8/Treg ratios.
- `cibersort_statistical_comparison_table.tsv`: All 22 cells + ratios: Median High/Low, Cliff's $\delta$, raw Wilcoxon $P$, BH $q$.
- `estimate_integration_table.tsv`: ESTIMATE scores, cross-platform validation metrics.
- `antioxidant_immune_correlation_table.tsv`: Full 7 × 26 Spearman correlation matrix (antioxidant genes × immune populations/ESTIMATE scores).
- `checkpoint_expression_statistics.tsv`: Canonical and extended checkpoint expression statistics.
- `tide_predictions_scores.tsv`: Complete TIDE output.
- `immune_survival_cox_models.tsv`: Uni- and multivariate Cox results.
- `master_phase3_immune_summary_table.tsv`: Full per-patient integration table.

### 6.3 Figures (`06_figures/06_immune_analysis/`)
All 30 figures as `.pdf` (vector) and `.png` (300 DPI).

### 6.4 Report (`07_reports/`)
`oscc_keap1_nrf2_phase3_immune_report.md`: Figure-by-figure narrative explicitly connecting oxidative stress biology (Phases 1/2) to immune landscape findings (Phase 3), with mechanistic interpretation and therapeutic implications.

---

## 7. Execution Sequence

1. **Step 0** — Script 00 (Python): Module 3.4 TIDE prediction via `00_tide_prediction.py`. Runs local `TIDEpy` package on mean-centered VST expression matrix to compute TIDE scores, dysfunction, exclusion, MDSC, CAF, and TAM M2 metrics (`05_results/06_immune_analysis/tide_predictions_scores.tsv`).
2. **Step 1** — Script 01 (R): Module 3.1 deconvolution characterisation (Figs 1–8).
3. **Step 2** — Script 02 (R): Module 3.2 antioxidant-immune bridge including novel TNFSF18, ALOX, and TGM analyses (Figs 9–12). *This is the mechanistic core — must precede checkpoint analysis to inform interpretation.*
4. **Step 3** — Script 03 (R): Module 3.3 checkpoint & immunomodulators (Figs 13–17).
5. **Step 4** — Script 04 (R): Module 3.4 TIDE & ICB resistance (Figs 18–21). Ingests TIDE output generated in Step 0.
6. **Step 5** — Script 05 (R): Module 3.5 survival integration including CACNA1A/PITX2/GATA3 immune correlates and smoking-NRF2 immune dichotomy (Figs 22–23).
7. **Step 6** — Script 06 (R): Module 3.6 master synthesis (Fig 24).
8. **Step 7** — Script 07 (R): Module 3.7 novel cross-phase angles (Figs 25–30). *Run after all prior modules so TNFSF18, ALOX, TGM, DEG-immune correlates can be contextualised against the full TME landscape.*
9. **Step 8** — Author comprehensive Phase 3 report.
## Immune Microenvironment, Checkpoint Architecture & Immunotherapy Relevance in the KEAP1/NRF2 OSCC Cohort

**Project Title**: KEAP1/NRF2 Axis in Oral Squamous Cell Carcinoma (OSCC) — TCGA-HNSC Cohort  
**Author**: Lead Bioinformatics Researcher & Cancer Systems-Biology Specialist  
**Cohort Scope**: TCGA Oral Cavity Squamous Cell Carcinoma (n = 224 primary tumors; 112 NRF2-High vs. 112 NRF2-Low; 16 matched normal tissues)  
**Governing Document**: [scope_of_work.txt](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/scope_of_work.txt) — Phase 3 (Sections 3.1 & 3.2)  
**Status**: Updated with Ingested CIBERSORTx & ESTIMATE Assets (Ready for Implementation)  

---

## 1. Executive Summary & Biological Rationales

Phase 3 transitions the project from intrinsic tumor cell genomics, redox pathway biology, and survival stratification (Phases 1 & 2) into the **extrinsic tumor microenvironment (TME)**. In oral squamous cell carcinoma (OSCC), sustained NRF2 hyperactivation not only orchestrates metabolic rewiring and ferroptosis resistance, but also profoundly shapes the immune and stromal landscapes:

1. **Redox-Driven Immune Exclusion & "Immune-Desert" Phenotype**: Hyperactive NRF2 suppresses endogenous reactive oxygen species (ROS) and downregulates pro-inflammatory chemokine cascades (e.g., *CXCL9*, *CXCL10*), fostering an immunologically "cold" or "excluded" phenotype that impedes cytotoxic CD8+ T-cell infiltration into the tumor core.
2. **Stromal Remodeling & Myeloid Co-option**: NRF2 activation fosters recruitment of immunosuppressive myeloid-derived suppressor cells (MDSCs) and polarizes tumor-associated macrophages (TAMs) toward the pro-tumorigenic M2 phenotype, accompanied by cancer-associated fibroblast (CAF) remodeling.
3. **Immune Checkpoint Plasticity & Primary Immunotherapy Failure**: Dysregulated antioxidant signaling alters cell-surface checkpoint ligand density (*CD274* / PD-L1, *PDCD1LG2* / PD-L2) and impairs antigen processing/presentation machinery (APM / MHC-I), directly driving primary resistance to immune checkpoint blockade (ICB).

This implementation plan provides a rigorous, publication-grade analytical framework to dissect these phenomena across **224 well-characterized OSCC tumors**, integrating the user's pre-processed CIBERSORTx deconvolution data and pre-computed ESTIMATE scores, expanding to stromal/CAF compartments, profiling immune checkpoints, computing TIDE scores, and delivering a 21-figure visual suite.

---

## 2. Ingestion & Validation of Existing Data Assets

The workspace directory `05_results/cibsersort_data/` contains high-value pre-processed assets for this exact cohort:
*   `cibersortx_results.csv`: Complete deconvolution for **all 224 primary tumors** across 22 leukocyte subsets (LM22), plus empirical deconvolution metrics (`P-value`, `Correlation`, `RMSE`).
*   `cibersortx_input_matrix.tsv`: Filtered expression matrix utilized for CIBERSORTx.
*   `estimate_scores.tsv`: Pre-computed ESTIMATE scores for **all 224 primary tumors** (`StromalScore`, `ImmuneScore`, `ESTIMATEScore`, `TumorPurity`).
*   `estimate_filtered.gct`: GCT 1.2 formatted filtered expression matrix (9,370 genes x 224 tumor samples).
*   `estimate_input_matrix.gct` & `estimate_input_matrix.tsv`: Normalized input matrices utilized for the ESTIMATE pipeline.

### 2.1 Sample Concordance & Data Integrity Audit
*   **Sample ID Matching**: 100% of the 224 sample identifiers in both `cibersortx_results.csv` (`Mixture`) and `estimate_scores.tsv` (`SampleID`) match the 224 primary tumor `File_ID`s in `01_data/processed/metadata_nrf2_classified.tsv`.
*   **Cohort Balance**: Exactly **112 NRF2-High** and **112 NRF2-Low** samples. Zero missing cases.
*   **Initial Biological Validation of ESTIMATE Scores**:
    *   **StromalScore**: Strongly depleted in NRF2-High tumors (Median -624.2 vs. -213.6; Wilcoxon $P = 5.75 \times 10^{-6}$; Spearman $\rho = -0.316$).
    *   **ImmuneScore**: Strongly depleted in NRF2-High tumors (Median -65.8 vs. 419.5; Wilcoxon $P = 1.02 \times 10^{-6}$; Spearman $\rho = -0.354$).
    *   **ESTIMATEScore**: Combined microenvironment score significantly lower in NRF2-High (Wilcoxon $P = 3.49 \times 10^{-8}$; Spearman $\rho = -0.393$).
    *   **TumorPurity**: Significantly higher in NRF2-High tumors (Median 0.867 vs. 0.804; Wilcoxon $P = 3.49 \times 10^{-8}$; Spearman $\rho = +0.393$).
    *   *Interpretation*: These preliminary findings provide striking independent confirmation that NRF2-High OSCC tumors represent an **immune-desert and stroma-depleted ("immune-cold") entity** with higher cellular tumor purity.

### 2.2 Bridging the CAF & Stromal Dimension
*   *Methodological Resolution*: CIBERSORT LM22 is strictly hematopoietic (22 immune subsets) and does not capture **Cancer-Associated Fibroblasts (CAFs)**, which are explicitly mandated by [scope_of_work.txt](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/scope_of_work.txt) Section 3.1.
*   *Solution*: 
    1. Directly ingest the verified `estimate_scores.tsv` to quantify `StromalScore` and `TumorPurity`.
    2. Cross-correlate `StromalScore` with canonical CAF markers (*ACTA2*, *FAP*, *COL1A1*, *COL1A2*, *PDGFRA*, *PDGFRB*, *TGFB1*, *POSTN*) derived from `vst_normalised.tsv`.
    3. Extract the **TIDE CAF exclusion score** during Section 3.2 TIDE modeling to determine whether CAF-mediated exclusion differs between NRF2 classes.
    4. Perform an **orthogonal cross-platform concordance check**: Correlate ESTIMATE `ImmuneScore` against the total sum of CIBERSORTx leukocyte fractions across all 224 tumors.

---

## 3. Statistical Testing & Multiple Testing Philosophy (FDR Guidance)

Per lead researcher best practices and explicit project instructions: **"Wherever it is not necessary, do not do unnecessary FDR correction."** Reflexive genome-wide or family-wise corrections across pre-specified candidate hypotheses drastically inflate Type II error (false negatives), obscuring genuine biological interactions.

### 3.1 Targeted Hypotheses (NO Multiple Testing Correction Applied)
The following pre-planned, hypothesis-driven tests evaluate established biological axes and will be reported with **exact, two-sided raw p-values**, effect sizes, and 95% confidence intervals:
1. **Primary Immune Infiltration Candidates (LM22)**:
   *   CD8+ T cells (cytotoxic antitumor effectors)
   *   M1 Macrophages (pro-inflammatory) vs. M2 Macrophages (immunosuppressive)
   *   **M1/M2 Ratio** (macrophage polarization balance)
   *   Regulatory T cells (Tregs)
   *   **CD8+ / Treg Ratio** (immune balance / cytotoxic capacity index)
   *   Natural Killer (NK) cells (activated vs. resting)
   *   Dendritic cells (activated vs. resting)
2. **ESTIMATE Bulk Microenvironment Metrics**:
   *   `StromalScore`, `ImmuneScore`, `ESTIMATEScore`, `TumorPurity`.
3. **Canonical Immune Checkpoints (A Priori Defined Targets)**:
   *   *PDCD1* (PD-1), *CD274* (PD-L1), *PDCD1LG2* (PD-L2), *CTLA4*, *LAG3*, *HAVCR2* (TIM-3), *TIGIT*, *IDO1*.
4. **Core TIDE Biomarkers**:
   *   Overall TIDE Score, T-cell Dysfunction Score, T-cell Exclusion Score.
5. **Bivariate Associations**:
   *   Direct Spearman correlation ($\rho$) between continuous NRF2 GSVA score and the primary candidate cell types, ESTIMATE scores, or checkpoints.

*Statistical tests for targeted hypotheses*: Two-sided Wilcoxon rank-sum test (Mann-Whitney U), Hodges-Lehmann median difference with 95% CI, and non-parametric effect sizes (**Cliff's delta** $\delta$ and **Cohen's d**).

### 3.2 Exploratory Wide Screens (BH-FDR Correction Rigorously Applied)
Multiple testing correction via the **Benjamini-Hochberg (BH) false discovery rate (q-value)** is strictly reserved for broad, exploratory screens:
1. **Full 22-Cell LM22 Profile Screen**: Screen across all 22 sub-populations (where raw $P$ and BH $q$ are presented side-by-side in data tables).
2. **Extended Checkpoint & Immunomodulator Panel**: Comprehensive screen of 40+ extended receptors/ligands (co-stimulatory, co-inhibitory, and B7 family members).
3. **Cytokine & Chemokine Repertoire**: Screen of 50+ immune-active signaling molecules.
4. **Antigen Processing & Presentation Machinery (APM)**: Screen of MHC-I/II and peptide loading genes.

---

## 4. Modular Analytical Architecture

```
Phase 3 Modular Pipeline:
├── Module 3.1: Immune Microenvironment Deconvolution (CIBERSORT LM22 + ESTIMATE + CAFs)
├── Module 3.2: Immune Checkpoint & Immunomodulatory Architecture
├── Module 3.3: TIDE & Computational Immunotherapy Response Prediction
├── Module 3.4: Immune-Prognostic Integration (Survival x TME Stratification)
└── Module 3.5: Master Synthesis & Multi-Omics TME Consensus Model
```

### Module 3.1: Immune Microenvironment Deconvolution
1. **Lineage Aggregation & Cross-Platform Concordance**:
   *   Aggregate fine-grained subsets into macro-lineages (Total T, Total Macrophages, Total B, Total Myeloid).
   *   Calculate functional ratios: $\text{M1/M2 Ratio}$, $\text{CD8/Treg Ratio}$.
   *   Cross-platform validation: Spearman correlation between ESTIMATE `ImmuneScore` and CIBERSORTx total leukocyte fraction.
2. **Comparative Infiltration Analysis**:
   *   Compare NRF2-High vs. NRF2-Low tumors across all 22 populations and aggregate lineages.
   *   Compare across NRF2 tertiles (Low vs. Mid vs. High) to assess dose-dependent microenvironmental remodeling.
3. **Continuous Redox-Immune Coupling**:
   *   Compute Spearman correlation matrix between continuous NRF2 GSVA score and all 22 cell subsets + macro-lineages.
4. **Stromal & CAF Quantification**:
   *   Ingest `estimate_scores.tsv` (`StromalScore`, `ImmuneScore`, `ESTIMATEScore`, `TumorPurity`).
   *   Correlate NRF2 scores and `StromalScore` with canonical CAF markers (*ACTA2*, *FAP*, *COL1A1*, *PDGFRA*, *TGFB1*).
5. **Tobacco Smoking-Stratified Immune Interaction**:
   *   Evaluate immune infiltration across smoking status (Smoker vs. Reformed vs. Non-Smoker) within NRF2 classes, testing whether tobacco-induced ROS enhances immune exclusion.

### Module 3.2: Immune Checkpoint & Immunomodulatory Architecture
1. **Targeted Checkpoint Profiling**:
   *   Extract VST-normalized expression for canonical checkpoint genes: *PDCD1*, *CD274*, *PDCD1LG2*, *CTLA4*, *LAG3*, *HAVCR2*, *TIGIT*, *IDO1*.
   *   Perform Wilcoxon tests and linear regression against NRF2 GSVA scores.
2. **Extended Immunomodulator Super-Family Screen**:
   *   Co-inhibitory: *BTLA*, *CD96*, *CD160*, *CD244*, *VSIR* (VISTA), *SIGLEC15*, *LAIR1*.
   *   Co-stimulatory: *CD28*, *ICOS*, *TNFRSF4* (OX40), *TNFRSF9* (4-1BB), *TNFRSF18* (GITR), *CD27*, *CD40LG*, *CD40*, *CD80*, *CD86*.
3. **Chemokine & Cytokine Trafficking Axis**:
   *   T-cell homing / CXCR3 axis: *CXCL9*, *CXCL10*, *CXCL11*, *CXCL12*.
   *   Myeloid / MDSC / Treg recruitment: *CCL2*, *CCL5*, *CCL20*, *CXCL8* (IL-8), *CXCL1*, *CXCL2*.
   *   Immunosuppressive cytokines: *TGFB1*, *IL10*, *VEGFA*, *IL6*.
4. **Antigen Processing & Presentation Machinery (APM)**:
   *   MHC Class I: *HLA-A*, *HLA-B*, *HLA-C*, *B2M*.
   *   Peptide loading: *TAP1*, *TAP2*, *TAPBP*, *ERAP1*, *ERAP2*.
   *   MHC Class II & master transactivator: *HLA-DRA*, *HLA-DRB1*, *HLA-DPA1*, *CIITA*.

### Module 3.3: TIDE & Computational Immunotherapy Response Prediction
1. **TIDE Input Matrix Preparation**:
   *   Prepare expression matrix formatted for TIDE (mean-centered log2 expression).
2. **TIDE Execution**:
   *   Derive primary output metrics:
       *   **TIDE Score**: Integrative signature of immune evasion (higher = resistant).
       *   **Dysfunction Score**: T-cell dysfunction in tumors with high cytotoxic T-lymphocyte (CTL) infiltration.
       *   **Exclusion Score**: Exclusion of CTLs mediated by immunosuppressive cells.
       *   Component exclusion scores: **MDSC**, **CAF**, **TAM M2**.
       *   **MSI Score**: Microsatellite instability biomarker.
       *   **Predicted Responder**: Categorical ICB response classification (True / False).
3. **Concordance with Emerging Biomarkers**:
   *   **Cytolytic Activity Score (CYT)**: Geometric mean of *GZMA* and *PRF1*.
   *   **IFN-$\gamma$ / T cell-inflamed GEP**: Ayers et al. 18-gene expanded immune signature.
   *   **Immunophenoscore (IPS)**: Quantifying antigen presentation, effector cells, checkpoints, and suppressor cells.

### Module 3.4: Immune-Prognostic Integration (Survival x TME Stratification)
1. **Bivariate Survival Stratification**:
   *   Kaplan-Meier survival analyses evaluating combined NRF2 and immune infiltration states:
       *   *CD8+ T-cell High / NRF2-Low* (Presumed best prognosis, inflamed)
       *   *CD8+ T-cell Low / NRF2-High* (Presumed worst prognosis, immune-excluded)
       *   Intermediate discordance groups.
2. **TIDE Prognostic Validation**:
   *   Kaplan-Meier curves for TIDE predicted Responders vs. Non-Responders across NRF2 strata.
3. **Multivariable Cox Modeling**:
   *   Integrate immune cell fractions, ESTIMATE scores, TumorPurity, TIDE score, NRF2 GSVA score, and clinical covariates (Age, Stage, Grade, Smoking) to determine if immune metrics explain or independently modify NRF2-driven survival risk. HPV data are unavailable and are excluded from this model.

---

## 5. Comprehensive Catalog of Publication-Grade Figures

All figures will be generated in two formats:
*   **Vector PDF**: Formatted at exact publication dimensions, ready for vector graphics manipulation (Adobe Illustrator / Inkscape).
*   **Raster PNG**: 300 DPI high-resolution for manuscripts, presentations, and Markdown embeds.
*   **Cohesive Palette**: Consistent with Phases 1 & 2 (`#E63946` for NRF2-High, `#457B9D` for NRF2-Low, `#2A9D8F` for Normal/Control, with diverging RdBu for expression and viridis for density).

| Figure ID | File Name | Plot Type | Biological Question & Key Features |
| :--- | :--- | :--- | :--- |
| **Fig 1** | `01_cibersort_landscape_complex_heatmap.pdf/png` | ComplexHeatmap (Annotated) | **Global Landscape**: Heatmap of all 22 LM22 cell fractions across 224 tumors ordered by NRF2 GSVA score. Top annotations: NRF2 group, GSVA score, Smoking status, Stage, Grade, Age, Sex. Bottom annotations: Total T, Total Macrophage, ESTIMATE ImmuneScore, TumorPurity. |
| **Fig 2** | `02_targeted_immune_cell_comparisons_grid.pdf/png` | Multi-panel Box/Violin + Jitter | **Targeted Infiltration**: 8-panel grid comparing CD8+ T cells, M1, M2, M1/M2 ratio, Tregs, CD8/Treg ratio, NK cells, and Dendritic cells between NRF2-High vs. Low. Shows median, IQR, jittered points, exact raw Wilcoxon p-value, Cliff's delta. |
| **Fig 3** | `03_cibersort_composition_stacked_barplots.pdf/png` | Stacked Percent Barplot | **TME Compositional Shifts**: Cumulative fractional composition of immune infiltrates comparing NRF2-High vs. NRF2-Low and across NRF2 tertiles. Demonstrates relative expansion/contraction of lineages. |
| **Fig 4** | `04_nrf2_immune_correlation_matrix.pdf/png` | Correlogram / Heatmap | **Continuous Redox Coupling**: Spearman correlation between continuous NRF2 GSVA score and all 22 LM22 cell types + summary ratios, with significance markers (*, **, ***). |
| **Fig 5** | `05_immune_infiltration_effect_size_volcano.pdf/png` | Volcano / Lollipop Plot | **Infiltration Effect Sizes**: Cliff's delta (or median difference) on x-axis vs $-\log_{10}(P\text{-value})$ on y-axis, highlighting significantly depleted (CD8, M1) vs enriched (M2, Tregs) populations in NRF2-High tumors. |
| **Fig 6** | `06_immune_infiltration_pca_umap.pdf/png` | 2D Dimension Reduction (PCA/UMAP) | **Immune Space Segregation**: Unsupervised clustering of tumors based solely on immune cell infiltration, colored by NRF2 group and GSVA score gradient, evaluating global microenvironmental divergence. |
| **Fig 7** | `07_estimate_stroma_tumor_purity_analysis.pdf/png` | 4-Panel Violin Grid + Cross-Validation | **ESTIMATE Microenvironmental Desertification**: Comparison of `StromalScore`, `ImmuneScore`, `ESTIMATEScore`, and `TumorPurity` between NRF2-High and Low (raw Wilcoxon $P$ and Cliff's $\delta$), plus orthogonal scatter correlation between ESTIMATE `ImmuneScore` and CIBERSORTx total leukocyte fraction. |
| **Fig 8** | `08_caf_markers_expression_correlation.pdf/png` | Multi-panel Scatter + Regression | **Cancer-Associated Fibroblasts (CAFs)**: Expression of canonical CAF genes (*ACTA2*, *FAP*, *COL1A1*, *PDGFRA*, *TGFB1*) vs. NRF2 GSVA score and ESTIMATE `StromalScore` with linear fit, 95% CI band, and Spearman $\rho$. |
| **Fig 9** | `09_smoking_stratified_immune_infiltration.pdf/png` | Grouped Boxplot Grid | **Smoking Interaction**: Infiltration of CD8+ T cells, M1/M2 ratio, and Tregs stratified by Smoking History (Current Smoker vs Reformed vs Non-Smoker) across NRF2 groups. |
| **Fig 10** | `10_canonical_checkpoints_expression_grid.pdf/png` | Multi-panel Violin + Box | **Core Immune Checkpoints**: Comparison of *PDCD1* (PD-1), *CD274* (PD-L1), *PDCD1LG2* (PD-L2), *CTLA4*, *LAG3*, *HAVCR2* (TIM-3), *TIGIT*, *IDO1* between NRF2-High and Low. Raw Wilcoxon $P$, effect size. |
| **Fig 11** | `11_checkpoint_nrf2_continuous_scatter_grid.pdf/png` | Multi-panel Scatter + Smooth | **Redox-Checkpoint Coupling**: Continuous correlation of NRF2 GSVA score with core checkpoint expression levels, annotated with Spearman $\rho$ and p-values. |
| **Fig 12** | `12_extended_immunomodulators_heatmap.pdf/png` | ComplexHeatmap | **Broad Immunomodulator Screen**: Clustered heatmap of 35+ co-inhibitory and co-stimulatory receptors/ligands, annotated with BH-FDR significance tags. |
| **Fig 13** | `13_chemokine_cytokine_trafficking_grid.pdf/png` | Bi-directional Lollipop / Forest | **T-cell Homing vs Myeloid Recruitment**: Differential expression ($\log_2\text{FC}$ and $-\log_{10}P$) of T-cell recruiting chemokines (*CXCL9/10/11*) vs myeloid recruiting factors (*CCL2*, *IL8*, *VEGFA*, *TGFB1*). |
| **Fig 14** | `14_antigen_presentation_apm_machinery.pdf/png` | Multi-panel Violin Grid | **Antigen Processing & Presentation**: Expression of MHC-I (*HLA-A/B/C*, *B2M*), MHC-II (*HLA-DRA*, *CIITA*), and peptide loaders (*TAP1/2*, *TAPBP*) comparing NRF2 groups. |
| **Fig 15** | `15_tide_benchmark_scores_comparison.pdf/png` | Multi-panel Violin/Box | **TIDE Primary Scores**: Distribution of TIDE Score, T-cell Dysfunction Score, and T-cell Exclusion Score between NRF2-High and NRF2-Low. |
| **Fig 16** | `16_tide_exclusion_mediators_breakdown.pdf/png` | 3-Panel Violin Grid | **Exclusion Mechanisms**: Breakdown of specific cell-type exclusion scores: MDSC score, CAF score, and TAM M2 score between NRF2 groups. |
| **Fig 17** | `17_tide_responder_prediction_contingency.pdf/png` | Stacked Proportion Bar + Odds Ratio | **ICB Response Prediction**: Proportion of predicted Responders vs Non-Responders in NRF2-High vs Low tumors; Fisher's exact test, Odds Ratio and 95% CI. |
| **Fig 18** | `18_cytolytic_activity_and_gep_scores.pdf/png` | 2-Panel Violin + Scatter | **Cytolytic & Inflamed Scores**: Cytolytic activity index (CYT: *GZMA* + *PRF1*) and Ayers 18-gene T-cell inflamed GEP score across NRF2 groups. |
| **Fig 19** | `19_immune_nrf2_bivariate_km_survival.pdf/png` | Multi-panel Kaplan-Meier Grid | **Immune-Prognostic Synergy**: KM curves for CD8+ T-cell / NRF2 combined strata and TIDE Responder / NRF2 combined strata with log-rank tests and risk tables. |
| **Fig 20** | `20_tme_prognostic_multivariate_forest.pdf/png` | Publication Forest Plot | **Multivariable Independence**: Forest plot of Hazard Ratios (95% CI) for NRF2, CD8+ T cells, M2 macrophages, TumorPurity, TIDE score, and clinical covariates. |
| **Fig 21** | `21_consensus_tme_schematic_overview.pdf/png` | Integrative Multi-metric Matrix | **Master Synthesis**: Integrative multi-track summary aligning NRF2 status, CIBERSORT infiltration, ESTIMATE scores, Checkpoint levels, TIDE evasion scores, and survival risk. |

---

## 6. Deliverables Pipeline & Data Artifacts

Following the established project structure:

### 6.1 Scripts Directory (`08_scripts/06_immune_analysis/`)
*   `01_cibersort_tme_characterization.R`: Ingests `05_results/cibsersort_data/cibersortx_results.csv`, calculates aggregate lineages, ratios, Wilcoxon tests (raw $P$ and Cliff's $\delta$), correlations, and generates Figures 1–6.
*   `02_stroma_caf_estimate_analysis.R`: Ingests `05_results/cibsersort_data/estimate_scores.tsv`, performs cross-platform validation with CIBERSORTx, profiles CAF markers and smoking interactions, generates Figures 7–9.
*   `03_checkpoint_immunomodulators.R`: Profiles canonical and extended checkpoints, APM machinery, and chemokine cascades, generates Figures 10–14.
*   `04_tide_immunotherapy_prediction.R`: Formats expression matrix, runs/integrates TIDE scores and CYT/GEP signatures, generates Figures 15–18.
*   `05_immune_survival_integration.R`: Bivariate KM curves, Cox models, and multivariate forest plots (incorporating TumorPurity), generates Figures 19–20.
*   `06_master_phase3_synthesis_plots.R`: Assembles the master TME matrix and consensus synthesis (Figure 21).

### 6.2 Results Directory (`05_results/06_immune_analysis/`)
*   `cibersort_lineage_and_ratios_summary.tsv`: Per-sample table with 22 LM22 fractions + aggregate lineages + M1/M2 and CD8/Treg ratios.
*   `cibersort_statistical_comparison_table.tsv`: Summary table for all 22 cell types and ratios (Median High, Median Low, Diff, Cliff's delta, raw Wilcoxon $P$, BH $q$-value).
*   `estimate_scores_and_caf_markers.tsv`: ESTIMATE scores, cross-validation metrics, and CAF gene expression per sample.
*   `checkpoint_expression_and_statistics.tsv`: Expression and statistical comparisons for canonical and extended checkpoints.
*   `tide_predictions_and_scores.tsv`: Complete TIDE output (TIDE, Dysfunction, Exclusion, MDSC, CAF, TAM M2, MSI, Responder).
*   `immune_survival_cox_models.tsv`: Univariate and multivariate survival models incorporating immune metrics and tumor purity.
*   `master_phase3_immune_summary_table.tsv`: Comprehensive integration table merging all patient-level metrics.

### 6.3 Figures Directory (`06_figures/06_immune_analysis/`)
*   All 21 figure files output as `.pdf` and `.png` (300 DPI).

### 6.4 Comprehensive Scientific Report (`07_reports/`)
*   `oscc_keap1_nrf2_phase3_immune_report.md`: Exhaustive, figure-by-figure analytical narrative linking redox biology, ferroptosis evasion, myeloid polarization, T-cell exclusion, and immunotherapy failure in OSCC.

---

## 7. Execution Timeline & Next Steps

Upon your approval:
1. **Step 1**: Initialize `08_scripts/06_immune_analysis/` and implement Script 01 & 02 (CIBERSORT Ingestion, ESTIMATE Analysis, CAFs, Figures 1–9).
2. **Step 2**: Implement Script 03 (Checkpoint Architecture, APM, Chemokines, Figures 10–14).
3. **Step 3**: Implement Script 04 (TIDE Execution, CYT/GEP Scoring, Figures 15–18).
4. **Step 4**: Implement Script 05 & 06 (Survival Integration, Master Synthesis, Figures 19–21).
5. **Step 5**: Compile the master tables and author the comprehensive scientific report.
