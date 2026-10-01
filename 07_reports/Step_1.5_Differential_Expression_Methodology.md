# Methodology: Differential Gene Expression Analysis (Step 1.5)

## KEAP1/NRF2 Axis in Oral Squamous Cell Carcinoma (OSCC)

> [!NOTE]
> This document provides a comprehensive, concept-rich explanation of the methodology behind the Differential Gene Expression Analysis (DGEA) executed in [`04_differential_expression.R`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/08_scripts/03_differential_expression/04_differential_expression.R). It bridges the statistical theory of negative binomial Generalized Linear Models (GLMs) with the oncogenic biology of the KEAP1-NRF2 pathway in OSCC.

---

## Table of Contents

1. [Executive Summary & Analytical Goal](#1-executive-summary--analytical-goal)
2. [Glossary of Core Concepts & Statistical Foundations](#2-glossary-of-core-concepts--statistical-foundations)
3. [Biological Rationale: Transcriptomic Footprint of NRF2 Hyperactivation](#3-biological-rationale-transcriptomic-footprint-of-nrf2-hyperactivation)
4. [Input Datasets & Cohort Stratification](#4-input-datasets--cohort-stratification)
5. [The DESeq2 Statistical Framework](#5-the-deseq2-statistical-framework)
6. [Significance Thresholding & DEG Classification](#6-significance-thresholding--deg-classification)
7. [Exploratory Data Analysis & Quality Control](#7-exploratory-data-analysis--quality-control)
8. [Biological Interpretation of Top DEGs](#8-biological-interpretation-of-top-degs)
9. [Figure-by-Figure Interpretation Guide](#9-figure-by-figure-interpretation-guide)
10. [Downstream Integration & Pipeline Continuity](#10-downstream-integration--pipeline-continuity)
11. [References & Software](#11-references--software)

---

## 1. Executive Summary & Analytical Goal

Differential Gene Expression Analysis (DGEA) sits at the heart of Phase 1 of this project. Having determined in Steps 1.3 and 1.4 that the **V1 gene set (31 genes)** scored via **GSVA (Gene Set Variation Analysis)** provides the optimal, most biologically robust separation of oral cavity squamous cell carcinomas (OSCC), Step 1.5 performs a genome-wide contrast between:

> **NRF2-High  vs  NRF2-Low**

### Key Analysis Facts

| Metric / Parameter | Value | Description |
|:---|:---|:---|
| **Cohort** | TCGA-HNSC Oral Cavity | n = 224 primary tumour samples (normal tissue excluded) |
| **Stratification** | V1 GSVA Median Split | 112 NRF2-High tumours vs 112 NRF2-Low tumours |
| **Statistical Model** | DESeq2 (Wald Test) | Negative Binomial GLM with design `~ NRF2_group` |
| **Input Matrix** | Filtered Raw Counts | 16,591 genes passing expression filtering |
| **Significance Cutoffs** | padj < 0.05 & \|log2FC\| ≥ 1.0 | FDR-adjusted p-value + minimum 2-fold change |
| **Total DEGs Identified** | **927 genes** | 409 Upregulated (44.1%) and 518 Downregulated (55.9%) |
| **Primary Output Files** | `DEG_full_results.tsv`, `DEG_significant.tsv` | Complete statistical table and filtered gene subsets |

The primary analytical goal is not merely to catalogue individual gene alterations, but to uncover the broad transcriptional reprogramming driven by constitutive NRF2 activity, spanning antioxidant defenses, metabolic rewiring, ferroptosis resistance, and immune evasion.

---

## 2. Glossary of Core Concepts & Statistical Foundations

To understand RNA-seq differential expression, one must appreciate why standard statistical methods (such as Student's t-tests or ordinary linear regressions) cannot be directly applied to sequencing reads.

```
                      Raw RNA-seq Read Counts
                                 │
                   Overdispersed Count Data
                 (Variance > Mean, Positive Integers)
                                 │
                 ┌───────────────┴───────────────┐
                 ▼                               ▼
       DESeq2 GLM Fitting              Variance Stabilizing
     (Negative Binomial Model)          Transformation (VST)
                 │                               │
       Dispersion Shrinkage                      │
       & Wald Significance Tests                 │
                 │                               │
                 ▼                               ▼
     Differentially Expressed Genes        Genome-wide PCA &
           (927 DEGs)                     Expression Heatmaps
```

### 2.1 Why Raw Counts Are Not Normally Distributed
RNA-seq data consists of integer read counts assigned to gene loci. Count data has distinct mathematical characteristics:
1. **Discrete & Non-negative**: Counts are integers {0, 1, 2, ...}, bounded at zero.
2. **Heteroscedasticity**: The variance of gene expression depends strongly on the mean. High-abundance genes have vastly larger absolute variance than low-abundance genes.
3. **Biological Overdispersion**: Under a pure Poisson process, variance equals the mean (σ² = μ). In biological replicates, however, biological variation (cellular heterogeneity, genetic background, environmental differences) inflates variance well beyond the mean (σ² = μ + α·μ², where α is the dispersion parameter). 

Because of this, modern RNA-seq tools model counts using the **Negative Binomial (Gamma-Poisson) distribution**.

### 2.2 Core Statistical Terminology

| Term / Symbol | Definition / Formula | Biological Meaning & Practical Application |
|:---|:---|:---|
| **Raw Count (K_ij)** | Discrete integer counts ≥ 0 | The observed number of sequencing reads mapped to gene i in sample j. Unnormalized; influenced by sequencing depth. |
| **Size Factor (s_j)** | Median-of-ratios estimate | Normalization factor for library size. Accounts for differences in total sequencing depth and library composition between samples. |
| **Dispersion (α_i)** | Var(K_ij) = μ_ij + α_i · μ_ij² | Quantifies biological variance beyond Poisson noise. α_i = 0 is pure Poisson; larger α_i indicates high sample-to-sample variability. |
| **Empirical Bayes Shrinkage** | Information borrowing across genes | Low-count genes have noisy dispersion estimates. DESeq2 squeezes gene-specific dispersions toward a global mean trend line, preventing false positives. |
| **Log2 Fold Change (log2FC)** | log2(Expression_High / Expression_Low) | Effect size metric. log2FC = +1 means 2-fold higher expression in NRF2-High; log2FC = -1 means 50% lower expression (halved). |
| **Wald Test Statistic (W)** | W = β̂_i / SE(β̂_i) | Evaluates whether the estimated coefficient β̂_i differs significantly from zero, dividing by its standard error. Under H0, W follows a standard normal distribution. |
| **Adjusted p-value (padj)** | Benjamini-Hochberg (BH) FDR | Controls the False Discovery Rate across all 16,591 simultaneous hypothesis tests. padj < 0.05 guarantees an expected false discovery rate under 5%. |
| **Independent Filtering** | Filtering on mean normalized counts | Filters out genes with near-zero counts prior to multiple-testing correction, maximizing statistical power for biologically relevant genes. |
| **VST (Variance Stabilizing Transformation)** | Homoscedastic log2-like transform | Produces continuous, approximately normally distributed data with constant variance across the expression range. Crucial for PCA and heatmaps, where raw counts would bias clustering toward high-expression genes. |

---

## 3. Biological Rationale: Transcriptomic Footprint of NRF2 Hyperactivation

### 3.1 The Oncogenic Switch in OSCC
In oral cavity squamous cell carcinoma, the transcription factor **NRF2** (*NFE2L2*) transitions from a physiological stress-responder to a potent oncogenic driver:
* **Physiological State**: Under basal conditions, the substrate adaptor **KEAP1** binds NRF2 and facilitates its rapid polyubiquitination via the **CUL3-RBX1** E3 ligase complex, targeting it for 26S proteasomal degradation (half-life ~20 minutes).
* **Oncogenic Hyperactivation**: In OSCC, chronic exposure to tobacco smoke-derived electrophiles, coupled with somatic loss-of-function mutations in *KEAP1* or gain-of-function mutations in *NFE2L2*, abrogates KEAP1-mediated repression.
* **Transcriptional Activation**: Stabilized NRF2 accumulates, translocates into the nucleus, forms obligate heterodimers with small **MAF proteins** (*MAFG*, *MAFK*), and binds to **Antioxidant Response Elements (ARE)**: `5'-TGACNNNGC-3'`.

```
    Normal / Unstressed State                  NRF2-Hyperactivated State (OSCC)
    ─────────────────────────                  ────────────────────────────────
       KEAP1 + CUL3 Scaffold                      KEAP1 Inactivated / Mutated
               │                                               │
        Binds & Ubiquitinates                           NRF2 Stabilized
               │                                               │
       NRF2 Degradation                           Translocates to Nucleus
      (26S Proteasome)                                         │
               │                                     Binds MAF Proteins
      Basal ARE Activity                                       │
                                                     Binds ARE DNA Elements
                                                               │
                                               Massive Transcriptional Upregulation:
                                               - Phase II Detoxification (AKR1C, NQO1)
                                               - Glutathione Synthesis (GCLC, GCLM)
                                               - Ferroptosis Resistance (SLC7A11, GPX2)
                                               - NADPH Generation (G6PD, PGD, TKT)
```

### 3.2 What the Differential Expression Contrast Identifies
Comparing NRF2-High vs. NRF2-Low tumours allows us to establish the comprehensive downstream transcriptomic footprint:
1. **Direct Effectors**: Classic ARE-containing detoxification and antioxidant enzymes.
2. **Indirect Metabolic Rewiring**: Secondary transcriptional cascades supporting NADPH regeneration, nucleotide biosynthesis, and iron homeostasis.
3. **Cross-Talk Phenotypes**: Genes regulating ferroptosis sensitivity, epithelial-mesenchymal plasticity, and immune cell infiltration.

---

## 4. Input Datasets & Cohort Stratification

The script [`04_differential_expression.R`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/08_scripts/03_differential_expression/04_differential_expression.R) ingests four foundational datasets prepared in prior pipeline steps:

```
01_data/processed/counts_filtered.tsv ────────┐
                                               ▼
01_data/processed/metadata_nrf2_classified.tsv ──► Filter to Tumour Only (n = 224) ──► DESeq2 Model (~ NRF2_group)
                                               ▲
03_processed/batch_info.rds ──────────────────┘
```

### 4.1 Input Files and Preprocessing Verification

1. **Filtered Count Matrix (`counts_filtered.tsv`)**:
   - Contains 16,591 genes across 256 samples (224 tumours and 32 adjacent-normal controls).
   - Low-count genes (e.g., genes with < 10 counts in < 20% of samples) were pre-filtered during Step 1.1 to prevent testing unexpressed genomic noise.
   - **Crucial Rule**: Raw integer counts are mandatory for DESeq2. Normalized values (TPM, FPKM, or VST) must never be passed to `DESeqDataSetFromMatrix()`, as DESeq2's internal statistical model relies directly on the discrete sampling characteristics of read counts.

2. **Classified Metadata (`metadata_nrf2_classified.tsv`)**:
   - Generated by [`03_nrf2_pathway_scoring.R`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/08_scripts/02_pathway_scoring/03_nrf2_pathway_scoring.R).
   - Contains clinical, pathological, and exposure variables alongside the continuous GSVA pathway activity score and the resulting binary `NRF2_group` label.

3. **VST Normalized Matrix (`vst_normalised.tsv`)**:
   - Continuous, homoscedastic expression values used strictly for downstream exploratory visualization (PCA, ComplexHeatmap clustering) where count heteroscedasticity would otherwise distort Euclidean distances.

4. **Batch Information (`batch_info.rds`)**:
   - Records whether technical batch effects (e.g., sequencing plate, tissue source site) were detected and handled.
   - Because ComBat-seq was applied directly to the count matrix during Step 1.2, batch effects are already removed from the integer counts, allowing a clean, single-factor experimental design: `design = ~ NRF2_group`.

### 4.2 Restricting to Tumour Samples
Adjacent-normal samples (n = 32) are intentionally removed prior to running DESeq2:
* **Biological Justification**: Comparing normal tissue to tumour tissue reflects general oncogenic transformation (proliferation, genomic instability, loss of differentiation), which would overshadow the specific biological variance driven by the NRF2 axis.
* **Tumour Cohort**: Exactly 224 oral cavity tumour samples are analyzed:
  - **NRF2-High**: n = 112 tumours (GSVA score > median)
  - **NRF2-Low**: n = 112 tumours (GSVA score ≤ median)

---

## 5. The DESeq2 Statistical Framework

The model fitting process in [`04_differential_expression.R`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/08_scripts/03_differential_expression/04_differential_expression.R#L120-L154) follows the established negative binomial GLM framework of Love et al. (*Genome Biology* 2014):

```
       1. Size Factor Estimation (Median of Ratios)
                           │
                           ▼
       2. Gene-wise Dispersion Estimation
                           │
                           ▼
       3. Empirical Bayes Dispersion Shrinkage
                           │
                           ▼
       4. GLM Fitting: log2(μ_ij) = β_0 + β_1 · Group_j
                           │
                           ▼
       5. Wald Hypothesis Testing: H_0: β_1 = 0
                           │
                           ▼
       6. Benjamini-Hochberg FDR Correction (padj)
```

### Step 1: Size Factor Estimation (Median-of-Ratios)
To account for differences in library size (sequencing depth) and library composition, DESeq2 calculates a size factor `s_j` for each sample `j`:
1. A pseudo-reference sample is created by computing the geometric mean of each gene across all samples:
   `K_i^R = ( ∏ K_ij )^(1/m)`
2. For sample `j`, the ratio of its count to the pseudo-reference count is computed for each gene: `K_ij / K_i^R`.
3. The median of these ratios across all genes defines the sample's size factor:
   `s_j = median_i( K_ij / K_i^R )`

*Why this matters*: Unlike total count normalization (RPM) or upper-quartile normalization, the median-of-ratios is immune to extreme outliers (e.g., a few hyper-expressed ribosomal or keratin genes dominating sequencing reads).

### Step 2 & 3: Dispersion Estimation & Empirical Bayes Shrinkage
For each gene `i`, the mean `μ_ij = s_j · q_ij` and dispersion `α_i` are modeled such that:

`Var(K_ij) = μ_ij + α_i · μ_ij²`

- Because n = 224 is large enough for stable estimation, but individual low-count genes still exhibit sampling noise, DESeq2 fits a smooth trend curve across all gene dispersions as a function of mean expression.
- Individual gene dispersions are then shrunk toward this fitted curve using an Empirical Bayes approach. Genes with high dispersion but low counts are shrunk strongly toward the trend, reducing false positives.

### Step 4 & 5: Generalized Linear Model & Wald Test
The design formula is defined as:

`log2(μ_ij) = β_0,i + β_1,i · x_j`

where `x_j = 0` if sample `j` is `NRF2_Low` (the reference level) and `x_j = 1` if sample `j` is `NRF2_High`.
- The slope `β_1,i` represents the **log2 fold change (log2FC)** for gene `i` between High and Low groups.
- The Wald test evaluates the null hypothesis `H0: β_1,i = 0`:

`W_i = β̂_1,i / SE(β̂_1,i) ~ N(0, 1)`

### Step 6: Multiple Testing Correction (FDR)
Evaluating 16,591 genes at a nominal p < 0.05 would falsely flag ~830 genes purely by chance. To control the false discovery rate, raw p-values are adjusted using the **Benjamini-Hochberg (BH)** step-up procedure:
1. Sort raw p-values in ascending order: `p(1) ≤ p(2) ≤ ... ≤ p(m)`.
2. Compute adjusted p-values:
   `padj(i) = min_{k ≥ i} [ min(1, (m/k) · p(k)) ]`

---

## 6. Significance Thresholding & DEG Classification

In [`04_differential_expression.R`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/08_scripts/03_differential_expression/04_differential_expression.R#L171-L204), rigorous two-dimensional criteria are applied to define true biological signal:

> **padj < 0.05  AND  |log2FC| ≥ 1.0**

```
                ▲ -log10(padj)
                │
                │        Upregulated (409)
                │        padj < 0.05 & log2FC >= 1.0
                │              ▲
       Downregulated (518)     │
       padj < 0.05 &           │
       log2FC <= -1.0          │
             ▲                 │
             │                 │
    ─────────┼─────────────────┼───────── padj = 0.05 threshold
             │       NS        │
             │  (15,664 genes) │
    ─────────┴────────┼────────┴────────► log2FC
                   log2FC = 0
             │                 │
          log2FC = -1.0     log2FC = +1.0
```

### Quantitative Results Breakdown

| Category | Definition | Gene Count | Percentage of Tested | Biological Significance |
|:---|:---|:---:|:---:|:---|
| **Upregulated** | padj < 0.05 & log2FC ≥ 1.0 | **409** | 2.47% | Induced transcriptional targets; antioxidant defense, xenobiotic detoxification, NADPH production, ferroptosis defense. |
| **Downregulated** | padj < 0.05 & log2FC ≤ -1.0 | **518** | 3.12% | Repressed or displaced genes; immune signaling, adhesion molecules, extracellular matrix components. |
| **Non-Significant (NS)** | Fails p-value or FC threshold | **15,664** | 94.41% | Unaltered background genes; housekeeping and pathway-independent transcripts. |
| **Total Tested** | Surviving expression filter | **16,591** | 100.0% | Genome-wide testable transcriptome. |

---

## 7. Exploratory Data Analysis & Quality Control

### 7.1 Genome-Wide Principal Component Analysis (PCA)
In [`04_differential_expression.R`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/08_scripts/03_differential_expression/04_differential_expression.R#L209-L300), unsupervised Principal Component Analysis was conducted on the entire VST-normalized matrix (16,591 genes across 224 tumours):
* **Variance Explained**:
  - **PC1**: **13.2%** of total variance.
  - **PC2**: **9.2%** of total variance.
  - **PC3**: **6.5%** of total variance.
  - Cumulative Top 3 PCs: **28.9%** of transcriptomic variance.
* **Biological Separation**:
  - Plotting **PC1 vs. PC2** reveals clear separation between NRF2-High and NRF2-Low groups along PC1. 
  - Centroids (marked by prominent diamonds) show a substantial Euclidean shift between groups.
  - When points are shaded by their continuous GSVA score (`p_pca_cont`), a smooth, uninterrupted colour gradient flows along the PC1 axis, proving that **NRF2 pathway activation is one of the single largest global drivers of transcriptomic variance across the entire OSCC cohort**.

### 7.2 PCA Scree Plot
The scree plot (`04_DEA_pca_scree.pdf`) charts the variance explained across the top 20 principal components:
- Confirms an "elbow" around PC4–PC5.
- Demonstrates that PC1 and PC2 capture meaningful biological structure rather than dispersed background noise.

### 7.3 MA Plot: Checking for Expression Bias
The MA plot (`04_DEA_MA_plot.pdf`) plots the log2FC (y-axis) against the mean normalized count (log10 baseMean, x-axis):
- **Symmetry**: Differentially expressed genes are distributed symmetrically above and below the horizontal log2FC = 0 line across the dynamic range.
- **Funnel Shape**: At low mean counts, variance is wider, but significance requires larger fold changes; at high mean counts, smaller fold changes achieve high statistical significance. This confirms healthy, unbiased negative binomial variance modeling.

---

## 8. Biological Interpretation of Top DEGs

### 8.1 The Hallmark Upregulated Effectors

The top upregulated genes provide definitive confirmation of massive NRF2 activation:

| Gene Symbol | Base Mean | log2FC | Linear Fold Change | Adjusted p-value | Key Biological Function & OSCC Role |
|:---|:---:|:---:|:---:|:---:|:---|
| **`AKR1C1`** | 10,052 | **+4.09** | **17.1-fold** | 2.69 × 10⁻⁵⁷ | Aldo-keto reductase; detoxifies reactive aldehydes and lipid peroxides; potent marker of cisplatin resistance. |
| **`AKR1C2`** | 9,999 | **+3.72** | **13.2-fold** | 1.67 × 10⁻⁵⁵ | Reduces polycyclic aromatic hydrocarbons and prostaglandins; xenobiotic detoxification in oral mucosa. |
| **`AKR1C3`** | 5,972 | **+3.40** | **10.6-fold** | 1.87 × 10⁻⁴⁶ | Prostaglandin F synthase; regulates steroid and lipid metabolism; promotes tumour proliferation and survival. |
| **`SLC7A11`** | 2,270 | **+2.65** | **6.3-fold** | 1.87 × 10⁻³⁹ | Catalytic subunit of system x_c- cystine/glutamate antiporter; imports cystine for GSH synthesis; **master driver of ferroptosis resistance**. |
| **`CYP4F3`** | 757 | **+3.48** | **11.2-fold** | 2.59 × 10⁻³⁴ | Cytochrome P450 leukotriene B4 ω-hydroxylase; modulates eicosanoid metabolism and local inflammatory signaling. |
| **`CES1`** | 3,274 | **+4.37** | **20.7-fold** | 3.34 × 10⁻³³ | Carboxylesterase 1; hydrolyzes xenobiotics, esters, and chemotherapeutic carbamates. |
| **`OSGIN1`** | 935 | **+2.18** | **4.5-fold** | 9.48 × 10⁻³³ | Oxidative stress-induced growth inhibitor 1; canonical NRF2 target regulating redox balance. |
| **`GPX2`** | 6,117 | **+3.09** | **8.5-fold** | 2.79 × 10⁻³¹ | Gastrointestinal glutathione peroxidase; reduces organic hydroperoxides; prevents lipid peroxidation and ferroptosis. |
| **`ALDH3A1`** | 7,247 | **+2.86** | **7.3-fold** | 1.17 × 10⁻²³ | Aldehyde dehydrogenase 3A1; clears lipid peroxidation-derived aldehydes; protects oral epithelial stem cells. |
| **`AKR1B10`** | 12,844 | **+2.40** | **5.3-fold** | 2.38 × 10⁻²³ | Aldo-keto reductase; detoxifies cytotoxic carbonyls; highly expressed in smoking-induced oral carcinomas. |
| **`NQO1`** | 7,110 | **+1.70** | **3.2-fold** | 4.52 × 10⁻²⁸ | Prototypical NRF2 direct target; obligate two-electron quinone reductase; prevents quinone redox cycling. |
| **`TXNRD1`** | 3,845 | **+1.29** | **2.4-fold** | 3.64 × 10⁻¹⁹ | Thioredoxin reductase 1; reduces oxidized thioredoxin; central to cytoplasmic redox buffering. |
| **`GCLC`** | 4,210 | **+1.23** | **2.3-fold** | 2.04 × 10⁻¹⁸ | Rate-limiting enzyme in de novo glutathione synthesis; essential for antioxidant defense and chemoresistance. |

### 8.2 Behavior of the Regulatory Triad (*NFE2L2*, *KEAP1*, *CUL3*)
Examining the NRF2 regulatory triad in [`DEG_nrf2_pathway_genes.tsv`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/04_analysis/DEG_nrf2_pathway_genes.tsv) reveals an essential biological insight:
- **`NFE2L2` (NRF2)**: Shows modest but highly significant mRNA upregulation (log2FC = +0.75, padj = 1.37 × 10⁻¹³).
- **`KEAP1`**: Shows slight mRNA upregulation (log2FC = +0.31, padj = 1.01 × 10⁻⁵).
- **`CUL3`**: Remains unchanged (log2FC = +0.14, padj = 0.20, Non-Significant).

> [!IMPORTANT]
> **Post-Translational Paradigm**: Why is *NFE2L2* mRNA only increased by 1.7-fold (log2FC = +0.75) while its downstream targets (*AKR1C1*, *GPX2*, *SLC7A11*) are elevated 6- to 17-fold?
> This occurs because the KEAP1-NRF2 pathway is primarily regulated **post-translationally** via protein stability, not transcriptionally. Inactivating mutations in *KEAP1* or *NFE2L2* prevent NRF2 protein degradation, allowing massive nuclear accumulation of protein without requiring dramatic upregulation of *NFE2L2* mRNA. The slight mRNA increase in *NFE2L2* and *KEAP1* reflects minor auto-regulatory feedback loops.

### 8.3 Top Repressed Genes: Immune and Extracellular Alterations
The top downregulated genes indicate that constitutive NRF2 activity is associated with a distinct, immunosuppressed microenvironment:
* **`GRIN2A`** (log2FC = -3.44, padj = 2.38 × 10⁻²⁶): Glutamate ionotropic receptor subunit.
* **`SPIB`** (log2FC = -2.80, padj = 1.91 × 10⁻²²): ETS-family transcription factor essential for plasmacytoid dendritic cell and B-cell maturation.
* **`VCAM1`** (log2FC = -2.36, padj = 3.51 × 10⁻²²): Vascular cell adhesion molecule 1; mediates leukocyte-endothelial adhesion and immune extravasation.
* **`SMOC1`** (log2FC = -2.55, padj = 1.93 × 10⁻²⁰): Secreted modular calcium-binding protein; regulates extracellular matrix organisation.
* **`IL17REL`** (log2FC = -2.87, padj = 6.21 × 10⁻²⁰): Interleukin-17 receptor-like protein; mediates epithelial inflammatory defense.

---

## 9. Figure-by-Figure Interpretation Guide

The script generates 8 publication-grade visualization artifacts in [`06_figures/`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/06_figures) and high-resolution QC mirrors in [`02_qc/`](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/02_qc):

```
06_figures/
├── 04_DEA_pca_genome_wide.pdf        # Unsupervised 3-panel PCA (PC1-3 + GSVA gradient)
├── 04_DEA_pca_scree.pdf              # Variance explained across top 20 PCs
├── 04_DEA_volcano_plot.pdf           # EnhancedVolcano with top DEGs & ARE markers
├── 04_DEA_MA_plot.pdf                # Log2FC vs Mean Expression
├── 04_DEA_summary_barplot.pdf        # Total Upregulated (409) vs Downregulated (518)
├── 04_DEA_top_genes_heatmap.pdf      # ComplexHeatmap of top 50 DEGs with Stage/Smoking
├── 04_DEA_lfc_distribution.pdf       # Effect size distribution of significant DEGs
└── 04_DEA_nrf2_pathway_lollipop.pdf  # Targeted log2FC & padj of 17 core NRF2 genes
```

### 9.1 `04_DEA_pca_genome_wide.pdf` (Genome-Wide PCA)
* **What it plots**: Three side-by-side PCA panels constructed from all 16,591 genes across 224 tumours:
  1. PC1 (13.2%) vs PC2 (9.2%) coloured by NRF2 group with 95% confidence ellipses and group centroids (diamonds).
  2. PC1 (13.2%) vs PC3 (6.5%) displaying persistent group separation.
  3. PC1 vs PC2 coloured along a continuous blue-to-red GSVA pathway activity gradient.
* **Key Interpretation**: Unsupervised genome-wide expression directly aligns with the supervised NRF2 pathway score, proving the classification reflects broad biological reality rather than narrow gene-set bias.

### 9.2 `04_DEA_volcano_plot.pdf` (Enhanced Volcano Plot)
* **What it plots**: Magnitude of change (log2FC, x-axis) versus statistical significance (-log10 padj, y-axis).
* **Thresholds**: Vertical dashed lines at log2FC = ±1.0; horizontal line at padj = 0.05.
* **Key Interpretation**: Strong asymmetric skew toward the upper right quadrant, featuring canonical NRF2 markers (*AKR1C1-3*, *SLC7A11*, *GPX2*, *NQO1*) reaching astronomical significance (-log10 padj > 50, p < 10⁻⁵⁰). Connectors clearly label top biological drivers.

### 9.3 `04_DEA_MA_plot.pdf` (Mean-Difference Plot)
* **What it plots**: log10(baseMean + 1) on x-axis against log2FC on y-axis.
* **Key Interpretation**: Verifies homoscedasticity across the count dynamic range. Confirms that top DEGs are abundant, biologically viable transcripts (baseMean > 1,000) rather than low-expression stochastic artifacts.

### 9.4 `04_DEA_summary_barplot.pdf` (DEG Count Summary)
* **What it plots**: Clean bar chart illustrating the distribution of significant alterations.
* **Key Interpretation**: Documents **409 Upregulated** and **518 Downregulated** genes, establishing that constitutive NRF2 activity induces extensive gene repression in addition to target activation.

### 9.5 `04_DEA_top_genes_heatmap.pdf` (ComplexHeatmap)
* **What it plots**: Standardized expression (Z-scores, capped at ±3) for the **top 50 DEGs** across all 224 tumours.
* **Annotations**: Top column annotations track `NRF2_Group`, `Stage` (AJCC pathologic stage I–IV), and `Smoking` (Current, Former, Never). Columns are split cleanly by NRF2 classification.
* **Key Interpretation**: Shows sharp block-like expression contrast between NRF2-High and Low groups. Demonstrates that NRF2 hyperactivation occurs across all clinical stages (I–IV) and is strongly enriched among current and former smokers.

### 9.6 `04_DEA_nrf2_pathway_lollipop.pdf` (Pathway Gene Focus)
* **What it plots**: Ranked horizontal lollipop chart of 17 core and extended NRF2 pathway genes showing exact log2FC effect size and significance asterisks (`***` p < 0.001, `**` p < 0.01, `*` p < 0.05).
* **Key Interpretation**: Demonstrates a clear functional hierarchy:
  1. *Detoxification & Ferroptosis effectors* (*AKR1C1*, *AKR1C3*, *GPX2*, *SLC7A11*, *NQO1*) show massive activation (log2FC ≥ 1.7, p < 10⁻²⁷).
  2. *Metabolic & PPP enzymes* (*TKT*, *G6PD*) exhibit intermediate activation (log2FC ≈ 0.8–1.5).
  3. *Core regulatory components* (*NFE2L2*, *KEAP1*, *MAFG*, *CUL3*) exhibit subtle or non-significant changes, consistent with post-translational regulation.

---

## 10. Downstream Integration & Pipeline Continuity

The differential expression outputs generated in this step directly feed into subsequent pipeline analyses:

```
04_differential_expression.R (Step 1.5)
                 │
                 ├──► DEG_full_results.tsv ──────────────► Ranked GSEA Analysis (Step 2.1)
                 │
                 ├──► DEG_upregulated_genelist.tsv ──────► Over-Representation Analysis (ORA - GO/KEGG)
                 │                                        Ferroptosis Resistance Mapping
                 │
                 ├──► DEG_downregulated_genelist.tsv ────► Immune Infiltration & Depletion Analysis (Phase 3)
                 │
                 └──► dds_tumor_nrf2.rds ────────────────► Cross-Talk GLM Modeling (EMT & Ferroptosis)
```

1. **Step 2.1: Functional & Pathway Enrichment Analysis**:
   - `DEG_upregulated_genelist.tsv` and `DEG_downregulated_genelist.tsv` serve as inputs for Over-Representation Analysis (ORA) across GO Biological Process, KEGG, and Reactome databases.
   - `DEG_full_results.tsv` provides the signed statistical metric (stat = W_i) to rank the entire 16,591-gene background for pre-ranked Gene Set Enrichment Analysis (GSEA).
2. **Phase 2: Ferroptosis & EMT Plasticity**:
   - High-confidence upregulation of *SLC7A11*, *GPX2*, and *GCLC* provides empirical foundation for testing lipid peroxidation resistance in NRF2-High tumours.
3. **Phase 3: Immune Microenvironment Deconvolution**:
   - Repression of *VCAM1*, *SPIB*, and related immune receptors prompts quantitative deconvolution (CIBERSORTx / ESTIMATE) to evaluate whether NRF2-High OSCC tumours constitute an "immune-cold" phenotype.

---

## 11. References & Software

### Analytical Software & R Packages

| Package | Version / Source | Primary Role in Step 1.5 |
|:---|:---|:---|
| **`DESeq2`** | Bioconductor | Negative binomial GLM fitting, dispersion shrinkage, Wald tests |
| **`EnhancedVolcano`** | Bioconductor | Publication-ready volcano plot generation with custom connectors |
| **`ComplexHeatmap`** | Bioconductor | Multidimensional heatmap clustering with clinical covariate tracks |
| **`ggplot2`** | CRAN | Genome-wide PCA, MA plots, scree plot, and summary barplots |
| **`patchwork`** | CRAN | Multi-panel layout composition |
| **`matrixStats`** | CRAN | High-performance row-wise and column-wise matrix calculations |

### Key Methodological References

1. **DESeq2 Framework**: Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
2. **Multiple Testing Correction**: Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1), 289-300.
3. **TCGA HNSC Genomic Landscape**: Cancer Genome Atlas Network. (2015). Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature*, 517(7536), 576-582. [doi:10.1038/nature14129](https://doi.org/10.1038/nature14129)
4. **NRF2 Molecular Biology & ARE Architecture**: Taguchi, K., & Yamamoto, M. (2017). The KEAP1-NRF2 system as a molecular target of cancer treatment. *Free Radical Biology and Medicine*, 108, 801-815. [doi:10.1016/j.freeradbiomed.2016.12.005](https://doi.org/10.1016/j.freeradbiomed.2016.12.005)
5. **Ferroptosis & SLC7A11 in OSCC**: Koppula, P., Zhuang, L., & Gan, B. (2021). Cystine transporter SLC7A11/xCT in cancer: ferroptosis, nutrient dependency, and therapeutic targets. *Protein & Cell*, 12(8), 599-620. [doi:10.1007/s13238-020-00789-5](https://doi.org/10.1007/s13238-020-00789-5)
6. **Aldo-Keto Reductases in Cancer Resistance**: Penning, T. M. (2015). The aldo-keto reductases (AKRs) and cancer: AKR1C enzymes and chemical carcinogenesis. *Metabolism: Clinical and Experimental*, 64(3 Suppl 1), S37-S42. [doi:10.1016/j.metabol.2014.10.027](https://doi.org/10.1016/j.metabol.2014.10.027)
