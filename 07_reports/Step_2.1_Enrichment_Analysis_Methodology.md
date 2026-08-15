# Methodology: Functional & Pathway Enrichment Analysis (Step 2.1)

## KEAP1/NRF2 Axis in Oral Squamous Cell Carcinoma (OSCC)

> [!NOTE]
> This document provides a detailed, beginner-friendly explanation of the methodology behind the enrichment analysis performed in [05_functional_enrichment.R](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/08_scripts/04_enrichment/05_functional_enrichment.R). It is written so that someone with minimal exposure to pathway enrichment analysis can understand both the *computational logic* and the *biological reasoning* underpinning every step.

---

## Table of Contents

1. [Why Enrichment Analysis?](#1-why-enrichment-analysis)
2. [Glossary of Key Terms](#2-glossary-of-key-terms)
3. [Biological Context: The KEAP1-NRF2 Axis in OSCC](#3-biological-context-the-keap1-nrf2-axis-in-oscc)
4. [Input Data — What Goes In](#4-input-data--what-goes-in)
5. [Gene ID Mapping — Translating Between Identifier Systems](#5-gene-id-mapping--translating-between-identifier-systems)
6. [Section 2: Over-Representation Analysis (ORA)](#6-section-2-over-representation-analysis-ora)
7. [Section 3: Gene Set Enrichment Analysis (GSEA)](#7-section-3-gene-set-enrichment-analysis-gsea)
8. [Section 4: Gene Set Variation Analysis (GSVA)](#8-section-4-gene-set-variation-analysis-gsva)
9. [Section 5: Ferroptosis Focus Analysis](#9-section-5-ferroptosis-focus-analysis)
10. [Section 6: NRF2 Extended Reference Gene Set Analysis](#10-section-6-nrf2-extended-reference-gene-set-analysis)
11. [Statistical Corrections Applied](#11-statistical-corrections-applied)
12. [Visualisations — How to Read Each Plot](#12-visualisations--how-to-read-each-plot)
13. [Output Summary](#13-output-summary)
14. [References & Software](#14-references--software)

---

## 1. Why Enrichment Analysis?

In Phase 1 of this project, we performed **differential gene expression analysis (DGEA)** and identified **927 genes** that are significantly different in expression between NRF2-High and NRF2-Low OSCC tumours. These are called **differentially expressed genes (DEGs)**.

However, a raw list of 927 individual gene names is difficult to interpret biologically. A single gene rarely acts alone — genes work together in **pathways** (coordinated chains of molecular interactions that accomplish a biological function). Enrichment analysis addresses a fundamental question:

> **"Among the hundreds of genes that changed, are there recognisable biological themes or pathways that are being turned on or shut down?"**

For example, if 30 out of your 409 upregulated genes happen to belong to the glutathione metabolism pathway, that is far more than you would expect by random chance from ~16,000 tested genes. This tells you that the NRF2-High tumours are specifically boosting their antioxidant defences — not just randomly changing gene expression.

This script uses **three complementary enrichment approaches** to answer this question from different angles:

| Approach | What it tests | Input | Analogy |
|----------|--------------|-------|---------|
| **ORA** | Are specific pathways over-represented in our DEG list? | Binary list: DEG yes/no | "Did the classroom disproportionately elect students from the chess club?" |
| **GSEA** | Do pathway genes tend to cluster at the top or bottom of the full ranked gene list? | All genes, ranked by significance and direction | "If you line up all students by exam score, do chess club members cluster near the top?" |
| **GSVA** | How active is each pathway in each individual tumour sample? | Per-sample expression matrix | "Give each student a 'chess skill score' and compare across groups" |

Using all three together gives robust, cross-validated results. If a pathway is flagged by ORA, confirmed by GSEA, and also shows differential GSVA scores, we can be highly confident it is biologically meaningful.

---

## 2. Glossary of Key Terms

| Term / Abbreviation | Definition |
|---------------------|------------|
| **DEG** | Differentially Expressed Gene — a gene whose expression level differs significantly between two conditions (here: NRF2-High vs NRF2-Low). |
| **log2FC (log2 Fold Change)** | A measure of how much a gene's expression differs between groups. log2FC = 1 means the gene is expressed 2x higher; log2FC = -2 means it is 4x lower. |
| **padj (adjusted p-value)** | A p-value corrected for multiple testing. Because we test thousands of genes simultaneously, some will appear significant by chance. Adjustment (typically by the **Benjamini-Hochberg / BH method**) controls the proportion of false discoveries. |
| **FDR (False Discovery Rate)** | The expected proportion of "significant" results that are actually false positives. A padj < 0.05 means we accept that up to 5% of our "significant" results might be false. |
| **Gene Set** | A predefined group of genes that share a common biological function, pathway membership, or regulatory relationship. Example: "Glutathione metabolism" = {GCLC, GCLM, GSR, GSS, GPX1, GPX2, ...}. |
| **GO (Gene Ontology)** | A standardised classification system that describes gene functions across three domains: **BP** (Biological Process — what the gene does), **MF** (Molecular Function — how the gene product works at the molecular level), and **CC** (Cellular Component — where in the cell the gene product is located). |
| **KEGG** | Kyoto Encyclopaedia of Genes and Genomes — a database of manually curated pathway maps linking genes to metabolic and signalling pathways. |
| **Reactome** | A peer-reviewed, open-source database of biological pathways and reactions in humans. |
| **MSigDB** | The Molecular Signatures Database — the largest collection of annotated gene sets, maintained by the Broad Institute. Contains multiple collections (H = Hallmark, C2 = Curated, C5 = GO-based, etc.). |
| **Hallmark Gene Sets** | A refined collection of 50 gene sets from MSigDB that represent well-defined biological states or processes (e.g., "EMT", "MYC Targets", "Inflammatory Response"). These are considered the "gold standard" summary of biology. |
| **NES (Normalised Enrichment Score)** | In GSEA, the enrichment score normalised for gene set size, allowing comparison across different gene sets. Positive NES = enriched at the top of the ranked list (upregulated); Negative NES = enriched at the bottom (downregulated). |
| **Entrez ID** | A unique numerical identifier for each gene, assigned by the NCBI (National Center for Biotechnology Information). Many enrichment databases use Entrez IDs rather than gene symbols. |
| **VST (Variance Stabilising Transformation)** | A mathematical transformation applied to RNA-seq count data (by DESeq2) that stabilises the variance across the expression range, making the data approximately normally distributed and suitable for distance-based analyses. |
| **limma** | A widely-used R/Bioconductor package for differential analysis of genomics data, using linear models and empirical Bayes moderation. |
| **Spearman correlation (rho)** | A non-parametric measure of how monotonically two variables are related. Unlike Pearson correlation, it does not assume linearity — it asks "do they consistently go up together?" Range: -1 (perfect inverse) to +1 (perfect agreement). |

---

## 3. Biological Context: The KEAP1-NRF2 Axis in OSCC

### What is the KEAP1-NRF2 pathway?

**NRF2** (encoded by the gene *NFE2L2*) is a master transcription factor — a protein that switches on dozens of other genes simultaneously. Specifically, NRF2 activates genes that protect cells against **oxidative stress** (damage caused by reactive oxygen species, or ROS — highly reactive molecules that can damage DNA, proteins, and lipids).

Under normal conditions, NRF2 is kept at low levels by **KEAP1** (Kelch-like ECH-associated Protein 1), which acts as its "brake." KEAP1 binds NRF2 in the cytoplasm and targets it for degradation via the **CUL3-RBX1** ubiquitin ligase complex. This is like a recycling truck that constantly picks up NRF2 and destroys it.

When the cell encounters oxidative stress (e.g., from tobacco smoke, radiation, or tumour metabolism), KEAP1's ability to bind NRF2 is disrupted. NRF2 escapes degradation, enters the nucleus, pairs with small **MAF proteins** (MAFG, MAFK), and binds to **ARE (Antioxidant Response Element)** sequences in DNA to activate its target genes.

### Why does this matter in OSCC?

In **oral squamous cell carcinoma (OSCC)**, the KEAP1-NRF2 axis is frequently hijacked by cancer cells:

- **Tobacco-induced ROS** — Smoking directly produces ROS in oral tissue, chronically activating NRF2. This is why smoking history is a prioritised covariate in this study.
- **KEAP1 mutations or silencing** — Cancer cells may lose KEAP1 function through mutations or epigenetic silencing, leaving NRF2 permanently "on."
- **Constitutive NRF2 activation** — When NRF2 is always active, it provides cancer cells with a survival advantage: enhanced detoxification of chemotherapy drugs, resistance to cell death (including ferroptosis), altered metabolism, and immune evasion.

### What this enrichment analysis asks

Having already identified 927 DEGs between NRF2-High and NRF2-Low OSCC tumours, this analysis asks:

1. **What biological processes are these DEGs involved in?** (ORA)
2. **Which pathways show coordinated expression shifts across the whole genome?** (GSEA)
3. **How do specific pathway activities differ per-sample between tumour groups?** (GSVA)
4. **Is there evidence that NRF2-High tumours resist ferroptosis?** (Ferroptosis focus)
5. **How do known KEAP1-NRF2 axis genes behave in our OSCC data?** (Reference gene set)

---

## 4. Input Data — What Goes In

The enrichment analysis begins with outputs from Phase 1:

### 4.1 DEG Results

| File | Contents | Use |
|------|----------|-----|
| [DEG_full_results.tsv](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/04_analysis/DEG_full_results.tsv) | All ~16,591 tested genes with log2FC, p-value, adjusted p-value | **GSEA** (needs all genes, ranked) |
| [DEG_significant.tsv](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/04_analysis/DEG_significant.tsv) | 927 significant DEGs (padj < 0.05, \|log2FC\| >= 1.0) | **ORA** (needs binary hit list) |

### 4.2 Expression Matrix

The **VST-normalised expression matrix** ([vst_normalised.tsv](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/01_data/processed/vst_normalised.tsv)) contains stabilised expression values for every gene across all 224 tumour samples. This is used for:
- **GSVA** (needs per-sample expression to calculate pathway scores)
- **Heatmaps** (visualising gene expression patterns)
- **Correlation analyses** (relating NRF2 scores to individual gene expression)

### 4.3 Sample Metadata

The [metadata_nrf2_classified.rds](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/03_processed/metadata_nrf2_classified.rds) file contains clinical information and NRF2 pathway classification for each sample:
- **NRF2_group**: `NRF2_High` (112 samples) or `NRF2_Low` (112 samples)
- **NRF2_GSVA_score**: Continuous pathway activity score (from Step 1.4)
- Clinical covariates: smoking status, tumour stage, etc.

### 4.4 NRF2 Reference Gene Set

The starting point is the **Union column** from [gsva_gene_lists.xlsx](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/10_references/gsva_gene_lists.xlsx), containing 43 curated KEAP1-NRF2 axis genes. This is then extended with additional biologically relevant genes (see Section 10).

---

## 5. Gene ID Mapping — Translating Between Identifier Systems

### Why mapping is necessary

Genes have multiple naming systems. Our expression data uses **HGNC gene symbols** (human-readable names like `NQO1`, `HMOX1`), but many pathway databases (GO, KEGG) use **NCBI Entrez IDs** (numerical codes like `1728`, `3162`).

Before running enrichment, we must translate our gene list from symbols to Entrez IDs using the `bitr()` function from clusterProfiler, which queries the **org.Hs.eg.db** annotation database (a comprehensive mapping of human gene identifiers).

### What happens during mapping

```
Gene Symbol -> org.Hs.eg.db lookup -> Entrez ID
  NQO1      ->                      -> 1728
  HMOX1     ->                      -> 3162
  SLC7A11   ->                      -> 23657
```

Not every gene symbol maps successfully — some may be obsolete aliases, non-coding transcripts, or novel genes not yet in the database. The script reports the mapping success rate (typically >85%).

### Universe definition

For ORA, we also define a **universe** — the full set of genes that were tested in the differential expression analysis. This is crucial because ORA uses a statistical test (hypergeometric / Fisher's exact) that needs to know the total "pool" of genes from which our DEGs were drawn. Using the wrong universe can inflate or deflate results.

---

## 6. Section 2: Over-Representation Analysis (ORA)

### 6.1 The Core Question

> "Given my list of, say, 409 upregulated genes, are certain biological categories represented more often than expected by chance?"

### 6.2 The Statistical Logic

ORA uses the **hypergeometric test** (or equivalently, Fisher's exact test). Imagine an urn containing all ~16,000 genes we tested. Some of those genes (say 200) belong to "glutathione metabolism." If we draw 409 genes randomly from this urn, how many glutathione genes would we expect to get? If our DEG list contains significantly more glutathione genes than expected, that pathway is "enriched."

Formally, for each pathway:

```
                    Pathway genes in DEGs     Expected by chance
Fisher's test:     ---------------------- vs ---------------------
                    Total DEGs                Total genes tested
```

The test produces a **p-value** — the probability of seeing this many (or more) pathway genes in our DEG list purely by chance. A small p-value (< 0.05 after correction) means the enrichment is statistically significant.

### 6.3 Separate Testing for Up and Down

We run ORA separately for **upregulated** and **downregulated** DEGs. This is important because a pathway could contain both up- and down-regulated genes, and mixing them could cancel out the signal. Separating them tells us:
- What biological processes are **activated** in NRF2-High tumours (upregulated genes)
- What biological processes are **suppressed** in NRF2-High tumours (downregulated genes)

### 6.4 Databases Tested

The script tests DEGs against three complementary pathway databases:

#### Gene Ontology (GO) — Three Sub-ontologies

| Ontology | Code | What it describes | Example |
|----------|------|-------------------|---------|
| Biological Process | **BP** | What biological function does the gene contribute to? | "glutathione metabolic process", "response to oxidative stress" |
| Molecular Function | **MF** | What molecular activity does the gene product perform? | "oxidoreductase activity", "glutathione transferase activity" |
| Cellular Component | **CC** | Where in the cell is the gene product located? | "mitochondrial matrix", "extracellular space" |

GO terms are organised hierarchically (like a tree). For instance, "glutathione metabolic process" is a child of "sulfur compound metabolic process," which is a child of "organic substance metabolic process." This hierarchy means results can be highly redundant.

**Simplification**: To reduce redundancy, the script applies `simplify()` with a semantic similarity cutoff of 0.7. This groups highly similar GO terms and keeps only the most significant representative from each cluster.

#### KEGG (Kyoto Encyclopaedia of Genes and Genomes)

KEGG provides hand-drawn pathway maps showing how genes interact in metabolic and signalling cascades. Examples relevant to KEAP1-NRF2 in OSCC:
- `hsa04216`: Ferroptosis
- `hsa00480`: Glutathione metabolism
- `hsa00980`: Metabolism of xenobiotics by cytochrome P450
- `hsa04151`: PI3K-Akt signalling pathway

#### Reactome

Reactome provides highly detailed, peer-reviewed reaction-level pathway annotations. It often captures finer mechanistic detail than KEGG. Examples:
- "Detoxification of Reactive Oxygen Species"
- "Biological oxidations"
- "Phase II conjugation of compounds"

### 6.5 Significance Thresholds

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `pvalueCutoff` | 0.05 | Only report terms with BH-adjusted p-value < 0.05 |
| `qvalueCutoff` | 0.1 | Additional FDR filter using the q-value method |
| `pAdjustMethod` | "BH" | Benjamini-Hochberg correction for multiple testing |

### 6.6 Enrichment Metrics Reported

For each enriched term, the script outputs:

| Metric | Meaning |
|--------|---------|
| **GeneRatio** | Fraction of DEGs in this pathway (e.g., 15/409 = 3.7% of upregulated DEGs) |
| **BgRatio** | Fraction of background genes in this pathway (e.g., 200/16000 = 1.25% of all tested genes) |
| **p.adjust** | BH-corrected p-value |
| **Count** | Number of DEGs in this pathway |
| **geneID** | The actual gene symbols that hit this pathway |

---

## 7. Section 3: Gene Set Enrichment Analysis (GSEA)

### 7.1 Why GSEA in Addition to ORA?

ORA has a fundamental limitation: it requires a binary cutoff. Genes are either "significant" or not. This means:
- Genes with padj = 0.051 (just barely non-significant) are completely ignored
- The magnitude of change is lost — a gene with log2FC = 5 is treated the same as one with log2FC = 1

**GSEA** eliminates both problems. It uses **all genes**, ranked by a continuous metric, and asks:

> "Do the genes in a particular pathway tend to be concentrated at the top (upregulated) or bottom (downregulated) of the ranked list, rather than scattered randomly?"

### 7.2 Building the Ranked Gene List

The script creates a single ranking metric for every gene:

```
Rank metric = sign(log2FC) x -log10(p-value)
```

This formula captures both:
- **Direction**: `sign(log2FC)` is +1 for upregulated genes and -1 for downregulated genes
- **Significance**: `-log10(p-value)` gives higher scores to more statistically significant genes

The result is a ranked vector where:
- The top has highly significant, strongly upregulated genes (large positive values)
- The bottom has highly significant, strongly downregulated genes (large negative values)
- The middle has genes with no significant change (near zero)

### 7.3 The GSEA Algorithm

For each gene set (pathway), GSEA performs these steps:

1. **Walk down the ranked list** from top to bottom
2. Every time you encounter a gene that IS in the pathway, take a step **up** (weighted by its rank metric)
3. Every time you encounter a gene that IS NOT in the pathway, take a step **down**
4. Track the cumulative score — this trace is the **running enrichment score**
5. The **Enrichment Score (ES)** is the maximum deviation from zero

If pathway genes cluster at the top of the list -> ES is large and positive -> pathway is "activated"
If pathway genes cluster at the bottom -> ES is large and negative -> pathway is "suppressed"
If pathway genes are scattered randomly -> ES stays near zero -> no enrichment

The ES is then **normalised** for gene set size to produce the **NES (Normalised Enrichment Score)**, allowing comparison across gene sets of different sizes.

### 7.4 Statistical Significance in GSEA

GSEA determines significance via **permutation testing**: it shuffles the gene labels thousands of times, re-computes ES each time, and builds a null distribution. The p-value is the fraction of permuted ES values that are as extreme as (or more extreme than) the observed ES.

The `eps = 1e-300` parameter sets the minimum reportable p-value, allowing extremely significant results to retain their precision.

### 7.5 MSigDB Collections Used

The script runs GSEA against three MSigDB collections:

#### Hallmark Gene Sets (Collection H)

The **50 Hallmark gene sets** are the most interpretable and widely used collection. They were derived by clustering thousands of gene sets from other collections and extracting their most coherent, non-redundant representatives. Each Hallmark set captures a well-defined biological state.

Key Hallmarks relevant to KEAP1-NRF2 in OSCC:

| Hallmark | Biological Relevance to NRF2 |
|----------|------------------------------|
| REACTIVE_OXYGEN_SPECIES_PATHWAY | Direct NRF2 target domain — NRF2 is the master regulator of ROS defence |
| XENOBIOTIC_METABOLISM | NRF2 transcriptionally activates Phase I/II drug-metabolising enzymes (CYPs, UGTs, GSTs) |
| EPITHELIAL_MESENCHYMAL_TRANSITION | NRF2 has known crosstalk with EMT regulators (SNAI1, ZEB1); relevant to OSCC invasion |
| MYC_TARGETS_V1 | MYC and NRF2 cooperate in metabolic reprogramming |
| PI3K_AKT_MTOR_SIGNALING | PI3K-AKT stabilises NRF2 by inhibiting GSK3B-mediated degradation |
| GLYCOLYSIS | NRF2 promotes the pentose phosphate pathway (G6PD, PGD, TKT), diverting glucose metabolism |
| INFLAMMATORY_RESPONSE | NRF2 has immunomodulatory effects relevant to tumour immune evasion |
| APOPTOSIS | NRF2 can promote resistance to apoptosis |
| DNA_REPAIR | Relevant to NRF2's role in therapy resistance |

#### Canonical Pathways (Collection C2:CP)

These are curated pathway gene sets from KEGG, Reactome, WikiPathways, PID, and BioCarta. They provide finer-grained pathway resolution than Hallmark sets. The top 30 by |NES| are reported.

#### GO Biological Process (Collection C5:GO:BP)

The MSigDB version of GO:BP gene sets, filtered to sets of 15-500 genes (to avoid both very small and very broad terms).

---

## 8. Section 4: Gene Set Variation Analysis (GSVA)

### 8.1 How GSVA Differs from ORA and GSEA

| Feature | ORA | GSEA | GSVA |
|---------|-----|------|------|
| Resolution | Group-level | Group-level | **Per-sample** |
| Input | DEG list | Ranked gene list | Expression matrix |
| Output | Enriched terms | Enriched terms + NES | **Pathway score per sample** |
| Statistical test | Hypergeometric | Permutation | Post-hoc (e.g., limma) |

GSVA is unique because it assigns each individual tumour sample a pathway activity score. This allows us to:
- Compare pathway distributions between NRF2-High and NRF2-Low groups
- Correlate pathway scores with clinical variables
- Visualise pathway activity across all samples in a heatmap

### 8.2 The GSVA Algorithm

For each gene set in each sample, GSVA:

1. Ranks all genes by their expression in that sample
2. Computes a **Kolmogorov-Smirnov-like statistic** — essentially, how much the pathway genes deviate from uniform distribution in the ranked expression list
3. This produces an "enrichment score" that indicates whether the pathway genes tend to be among the highly-expressed or lowly-expressed genes in that particular sample

The `kcdf = "Gaussian"` parameter tells GSVA to assume the VST-normalised expression values follow an approximately Gaussian (normal) distribution, which is appropriate for VST-transformed data.

### 8.3 Custom Pathway Gene Sets

The script constructs **8 custom, biologically curated gene sets** plus **12 MSigDB Hallmark sets**, covering the key pathway domains specified in the scope of work:

| Custom Gene Set | Genes | Biological Rationale |
|----------------|-------|----------------------|
| **ROS_Glutathione_Metabolism** | 27 genes (GSR, GCLC, GCLM, GPX1-4, SOD1/2, TXN, TXNRD1/2, HMOX1, NQO1, ...) | Core NRF2-driven antioxidant programme — these are the direct effectors of ROS neutralisation |
| **Xenobiotic_Metabolism** | 27 genes (CYPs, AKRs, ALDHs, ABCs, UGTs, ...) | Phase I/II detoxification enzymes — NRF2 activates these to metabolise and export carcinogens and chemotherapy drugs |
| **Ferroptosis_Regulation** | 21 genes (SLC7A11, GPX4, GCLC/M, FTH1/FTL, ACSL4, ...) | Genes that control ferroptosis sensitivity — a key NRF2-regulated cell death pathway |
| **EMT_Markers** | 17 genes (VIM, CDH1/2, SNAI1/2, ZEB1/2, TWIST1, ...) | EMT signature genes — NRF2 has bidirectional crosstalk with EMT in cancer invasion |
| **NRF2_Core_Activity** | 31 genes | The benchmarked V1 gene set used for NRF2 pathway scoring in Phase 1 |
| **Pentose_Phosphate_Pathway** | 10 genes (G6PD, PGD, TKT, TALDO1, ...) | NRF2 promotes PPP for NADPH production (essential for glutathione regeneration) |
| **Autophagy_p62_KEAP1** | 14 genes (SQSTM1, ATG5/7, ULK1, ...) | The p62/SQSTM1-KEAP1 feedback loop — autophagy defects can activate NRF2 |
| **DNA_Repair** | 19 genes (BRCA1/2, ATM, ATR, PARP1, ...) | NRF2 activation has been linked to altered DNA repair capacity and therapy resistance |

### 8.4 Statistical Testing of GSVA Scores

After GSVA generates per-sample pathway scores, the script uses **limma** to test whether scores differ between NRF2-High and NRF2-Low groups:

1. A **linear model** is fit: `GSVA_score ~ NRF2_group` for each pathway
2. **Empirical Bayes moderation** (via `eBayes()`) borrows information across pathways to produce more stable variance estimates
3. Results include **logFC** (difference in mean GSVA score), **t-statistic**, **p-value**, and **BH-adjusted p-value**

This is similar to doing a t-test for each pathway, but with better statistical properties due to the empirical Bayes shrinkage.

---

## 9. Section 5: Ferroptosis Focus Analysis

### 9.1 What is Ferroptosis?

**Ferroptosis** is a form of regulated cell death driven by iron-dependent lipid peroxidation. Unlike apoptosis (programmed cell death that cancer cells often evade), ferroptosis occurs when:

1. **Lipid peroxides** accumulate in cell membranes (toxic oxidised fats)
2. The cell's **antioxidant defence** (specifically, the cystine-glutathione-GPX4 axis) fails to neutralise these peroxides
3. The resulting membrane damage kills the cell

### 9.2 Why NRF2 and Ferroptosis?

NRF2 is arguably the most potent transcriptional shield against ferroptosis. It directly activates:

- **SLC7A11** (xCT) — imports cystine into the cell, the rate-limiting step for glutathione synthesis
- **GCLC / GCLM** — catalyse the first step of glutathione synthesis
- **GSS / GSR** — complete glutathione synthesis and regeneration
- **GPX4** — the enzyme that directly neutralises lipid peroxides (the last line of defence)
- **FTH1 / FTL** — ferritin subunits that sequester free iron (less free iron = less lipid peroxidation)
- **HMOX1** — degrades heme (a source of free iron), though this has complex pro/anti-ferroptotic effects
- **SLC40A1** (ferroportin) — exports iron from the cell

This means NRF2-High tumours are predicted to be **resistant to ferroptosis** — which has major therapeutic implications, since inducing ferroptosis is being explored as a cancer treatment strategy.

### 9.3 What this Section Does

The ferroptosis focus analysis performs three sub-analyses:

#### A. Expression Heatmap
Visualises the expression of 20 core ferroptosis genes across all 224 tumour samples, split by NRF2 group. Genes are annotated by their functional role:
- **Anti-ferroptosis** (green): SLC7A11, GPX4, GCLC, GCLM, GSS, GSR, FTH1, FTL, SLC40A1, HMOX1, NQO1, AIFM2, CBS, CTH
- **Pro-ferroptosis** (orange): ACSL4, LPCAT3, ALOX15, TFRC, CHAC1

#### B. Spearman Correlation Analysis
For each ferroptosis gene, computes the **Spearman correlation (rho)** between the gene's expression and the NRF2 GSVA pathway activity score across all 224 tumour samples.

A positive rho means the gene's expression increases as NRF2 activity increases — expected for NRF2 target genes like SLC7A11. A negative rho would indicate suppression with rising NRF2 activity. P-values are BH-corrected.

#### C. Lollipop Chart
Displays the log2FC of each ferroptosis gene from the DESeq2 analysis, colour-coded by role (anti- vs pro-ferroptosis), with significance stars.

---

## 10. Section 6: NRF2 Extended Reference Gene Set Analysis

### 10.1 Gene Set Construction

The reference gene set is built in layers:

| Layer | Source | Genes | Purpose |
|-------|--------|-------|---------|
| **Union (Base)** | [gsva_gene_lists.xlsx](file:///c:/Users/dipak/Desktop/Riku/NyberMan/KEAP1_NRF2_OSCC_Project/10_references/gsva_gene_lists.xlsx) — Union column | 43 genes | Core KEAP1-NRF2 genes validated across multiple gene set versions |
| **Ferroptosis** | Literature-curated | ~17 genes | SLC7A11, GPX4, ACSL4, LPCAT3, FTH1, FTL, etc. |
| **EMT** | Literature-curated | 10 genes | VIM, CDH1/2, SNAI1/2, TWIST1, ZEB1/2, etc. |
| **Autophagy** | p62-KEAP1 axis literature | 6 genes | SQSTM1, ATG5, ATG7, ULK1, MAP1LC3B, GABARAPL1 |
| **OSCC NRF2 Targets** | OSCC-specific literature | ~23 genes | Additional canonical targets relevant in oral cancer context |
| **Upstream Regulators** | Pathway biology | 9 genes | KEAP1, CUL3, NFE2L2, RBX1, GSK3B, AKT1, PIK3CA, etc. |

After deduplication, the extended set contains approximately 80-90 unique genes.

### 10.2 What this Section Examines

1. **How many reference genes are differentially expressed?** — Validates that the NRF2-High vs NRF2-Low classification captures real pathway biology
2. **Expression heatmap** — Shows the expression patterns of all reference genes across tumour samples, with row annotations for DEG status and gene category
3. **Lollipop chart** — Displays the fold change and significance of every reference gene, coloured by category

---

## 11. Statistical Corrections Applied

### Multiple Testing Problem

When testing thousands of gene sets or pathways simultaneously, some will appear significant by chance. For example, testing 10,000 GO terms at p < 0.05 would produce ~500 false positives even if no real enrichment exists.

### Benjamini-Hochberg (BH) Correction

All p-values in this analysis are corrected using the **Benjamini-Hochberg** method, which controls the **False Discovery Rate (FDR)** — the expected proportion of false positives among all rejected hypotheses.

The procedure:
1. Rank all p-values from smallest to largest
2. Multiply each p-value by (total number of tests) / (rank)
3. The resulting "adjusted p-values" (padj) control FDR

A padj < 0.05 threshold means we accept that up to 5% of our "significant" results may be false discoveries.

### Where BH Correction is Applied

| Analysis | What is corrected |
|----------|-------------------|
| ORA (GO, KEGG, Reactome) | p-values across all tested terms |
| GSEA (Hallmark, C2, C5) | p-values across all tested gene sets |
| GSVA (limma) | p-values across all tested pathways |
| Ferroptosis correlations | p-values across all tested genes |

---

## 12. Visualisations — How to Read Each Plot

### Dotplot (ORA / GSEA)
- **Y-axis**: Pathway / GO term name
- **X-axis**: Gene ratio (ORA) or NES (GSEA)
- **Dot size**: Number of genes or gene set size
- **Dot colour**: Adjusted p-value (darker = more significant)
- **Interpretation**: Large, dark dots positioned far from zero indicate the most significant and biologically impactful pathways

### Barplot (ORA)
- **Y-axis**: GO term description
- **X-axis**: -log10(adjusted p-value)
- **Bar colour**: Red = upregulated, Blue = downregulated
- **Interpretation**: Longer bars = more significant enrichment

### Gene-Concept Network (cnetplot)
- **Large nodes**: Enriched pathways/terms
- **Small nodes**: Individual genes
- **Edges**: Connect genes to the pathways they belong to
- **Interpretation**: Reveals which genes are shared across multiple pathways (hub genes) and which pathways overlap

### Enrichment Map (emapplot)
- **Nodes**: Enriched terms (sized by gene count)
- **Edges**: Semantic similarity between terms (thicker = more similar)
- **Clusters**: Groups of closely related terms
- **Interpretation**: Reduces redundancy by showing how enriched terms relate to each other thematically

### Treeplot
- Organises enriched GO terms into a hierarchical dendrogram (tree structure)
- Groups similar terms into labelled clusters
- **Interpretation**: Provides a high-level "executive summary" of what biological themes dominate

### Running Enrichment Score Plot (GSEA)
- **Top panel**: The running enrichment score as GSEA walks down the ranked list
- **Middle panel**: Vertical bars showing where pathway genes fall in the ranked list
- **Bottom panel**: The ranking metric values
- **Interpretation**: A peak at the left = pathway genes are enriched among upregulated genes; a trough at the right = enriched among downregulated genes

### Ridge Plot (GSEA)
- Shows the **distribution of fold changes** for genes in each enriched gene set
- **Interpretation**: Shifted to the right = most genes in the set are upregulated; shifted to the left = most are downregulated

### Heatmap (GSVA / Ferroptosis / Reference)
- **Rows**: Genes or pathways
- **Columns**: Individual tumour samples (split by NRF2 group)
- **Colour**: Z-score (blue = below average, white = average, red = above average)
- **Interpretation**: Blocks of red in NRF2-High samples indicate pathway genes that are specifically upregulated in that group

### Violin Plot (GSVA)
- Shows the **full distribution** of GSVA pathway scores for NRF2-High vs NRF2-Low
- Internal lines show 25th, 50th (median), and 75th percentiles
- **Interpretation**: Separated violins with non-overlapping medians indicate differential pathway activity

### Lollipop Chart (GSVA / Ferroptosis / Reference)
- Each "lollipop" represents a pathway or gene
- **Stem length**: Magnitude of log2FC or logFC
- **Dot size**: Statistical significance (-log10 padj)
- **Colour**: Direction (red = activated/upregulated, blue = suppressed/downregulated)
- **Stars**: Significance markers (* p < 0.05, ** p < 0.01, *** p < 0.001)

### Scatter Plot (Ferroptosis Correlations)
- Each dot is one tumour sample
- **X-axis**: NRF2 GSVA score
- **Y-axis**: VST expression of a ferroptosis gene
- **Trend line**: Linear regression fit with 95% confidence band
- **rho value**: Spearman correlation coefficient
- **Interpretation**: A positive slope with significant rho means the gene's expression tracks with NRF2 activity

---

## 13. Output Summary

### Result Tables (saved in `04_analysis/enrichment/`)

| File | Description |
|------|-------------|
| `ORA_GO_BP_up.tsv` / `ORA_GO_BP_down.tsv` | GO Biological Process ORA results |
| `ORA_GO_MF_up.tsv` / `ORA_GO_MF_down.tsv` | GO Molecular Function ORA results |
| `ORA_GO_CC_up.tsv` / `ORA_GO_CC_down.tsv` | GO Cellular Component ORA results |
| `ORA_KEGG_up.tsv` / `ORA_KEGG_down.tsv` | KEGG pathway ORA results |
| `ORA_Reactome_up.tsv` / `ORA_Reactome_down.tsv` | Reactome pathway ORA results |
| `GSEA_Hallmark_results.tsv` | GSEA results for MSigDB Hallmark gene sets |
| `GSEA_C2_CP_results.tsv` | GSEA results for Curated Canonical Pathways |
| `GSEA_C5_GO_BP_results.tsv` | GSEA results for GO Biological Process |
| `GSVA_pathway_scores.tsv` | Per-sample GSVA scores for all pathways |
| `GSVA_differential_pathways.tsv` | limma differential pathway activity results |
| `ferroptosis_analysis.tsv` | Ferroptosis gene DEG status + NRF2 correlations |
| `nrf2_extended_reference_geneset.tsv` | Extended KEAP1-NRF2 reference gene set with DEG annotations |
| `enrichment_summary.tsv` | Summary counts for all analyses |

### Publication Figures (saved in `06_figures/`)

| Figure | Analysis | Content |
|--------|----------|---------|
| `05_ORA_GO_BP_dotplot.pdf` | ORA | GO:BP enrichment (up + down) |
| `05_ORA_GO_MF_dotplot.pdf` | ORA | GO:MF enrichment |
| `05_ORA_GO_CC_dotplot.pdf` | ORA | GO:CC enrichment |
| `05_ORA_KEGG_dotplot.pdf` | ORA | KEGG pathway enrichment |
| `05_ORA_Reactome_dotplot.pdf` | ORA | Reactome pathway enrichment |
| `05_ORA_cnetplot_up.pdf` | ORA | Gene-concept network (upregulated GO:BP) |
| `05_ORA_cnetplot_down.pdf` | ORA | Gene-concept network (downregulated GO:BP) |
| `05_ORA_emapplot_up.pdf` | ORA | Enrichment map (upregulated GO:BP) |
| `05_ORA_emapplot_down.pdf` | ORA | Enrichment map (downregulated GO:BP) |
| `05_ORA_treeplot_up.pdf` | ORA | Hierarchical clustering (upregulated GO:BP) |
| `05_ORA_treeplot_down.pdf` | ORA | Hierarchical clustering (downregulated GO:BP) |
| `05_ORA_GO_BP_barplot.pdf` | ORA | Combined barplot (up + down GO:BP) |
| `05_GSEA_hallmark_dotplot.pdf` | GSEA | Hallmark NES dotplot |
| `05_GSEA_C2_CP_dotplot.pdf` | GSEA | C2:CP top 30 NES dotplot |
| `05_GSEA_running_score_key_pathways.pdf` | GSEA | Running score plots for key NRF2-relevant pathways |
| `05_GSEA_ridgeplot.pdf` | GSEA | Ridge plot of fold change distributions |
| `05_GSEA_hallmark_emap.pdf` | GSEA | Hallmark enrichment map |
| `05_GSVA_pathway_heatmap.pdf` | GSVA | Pathway activity heatmap across all samples |
| `05_GSVA_pathway_comparison.pdf` | GSVA | Violin plots for significant pathways |
| `05_GSVA_pathway_lollipop.pdf` | GSVA | Lollipop chart of differential pathway scores |
| `05_ferroptosis_heatmap.pdf` | Ferroptosis | Ferroptosis gene expression heatmap |
| `05_ferroptosis_correlations.pdf` | Ferroptosis | NRF2 vs ferroptosis gene scatter plots |
| `05_ferroptosis_lollipop.pdf` | Ferroptosis | Ferroptosis gene DEG lollipop |
| `05_nrf2_reference_geneset_heatmap.pdf` | Reference | Extended reference gene set heatmap |
| `05_nrf2_reference_lollipop.pdf` | Reference | Reference gene set fold change lollipop |

---

## 14. References & Software

### R Packages

| Package | Version | Citation |
|---------|---------|----------|
| **clusterProfiler** | >= 4.0 | Wu T, et al. (2021). *Innovation*. DOI: 10.1016/j.xinn.2021.100141 |
| **enrichplot** | >= 1.16 | Yu G (2023). enrichplot: Visualization of Functional Enrichment Result |
| **ReactomePA** | >= 1.40 | Yu G, He QY (2016). *Mol BioSyst*. DOI: 10.1039/c5mb00663e |
| **GSVA** | >= 1.44 | Hanzelmann S, et al. (2013). *BMC Bioinformatics*. DOI: 10.1186/1471-2105-14-7 |
| **fgsea** | >= 1.22 | Korotkevich G, et al. (2021). *bioRxiv*. DOI: 10.1101/060012 |
| **msigdbr** | >= 7.5 | Dolgalev I (2022). msigdbr: MSigDB Gene Sets for Multiple Organisms |
| **limma** | >= 3.52 | Ritchie ME, et al. (2015). *Nucleic Acids Res*. DOI: 10.1093/nar/gkv007 |
| **org.Hs.eg.db** | >= 3.17 | Carlson M. org.Hs.eg.db: Genome wide annotation for Human |
| **DESeq2** | >= 1.36 | Love MI, et al. (2014). *Genome Biol*. DOI: 10.1186/s13059-014-0550-8 |
| **ComplexHeatmap** | >= 2.12 | Gu Z, et al. (2016). *Bioinformatics*. DOI: 10.1093/bioinformatics/btw313 |

### Databases

| Database | URL | Description |
|----------|-----|-------------|
| Gene Ontology (GO) | http://geneontology.org | Standardised gene function annotations |
| KEGG | https://www.kegg.jp | Pathway maps for metabolism and signalling |
| Reactome | https://reactome.org | Peer-reviewed biological pathway database |
| MSigDB | https://www.gsea-msigdb.org | Curated gene set collections |

### Key Methodological References

1. Subramanian A, et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43):15545-15550.
2. Hanzelmann S, et al. (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*, 14:7.
3. Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS-B*, 57(1):289-300.
4. Rojo de la Vega M, et al. (2018). NRF2 and the Hallmarks of Cancer. *Cancer Cell*, 34(1):21-43.
5. Dodson M, et al. (2019). NRF2 plays a critical role in mitigating lipid peroxidation and ferroptosis. *Redox Biology*, 23:101107.
