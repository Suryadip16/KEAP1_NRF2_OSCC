# Novel Insights from NRF2-High vs NRF2-Low OSCC Enrichment Analysis

## Observations That Are Understudied or Uncharacterised in the Published OSCC/HNSCC Literature

> [!IMPORTANT]
> This report identifies observations from the DEG and enrichment analysis that, based on a targeted literature review, appear to represent either (a) mechanistic connections not previously described in OSCC/HNSCC specifically, or (b) findings that have been characterised in other cancer types but remain unvalidated in the oral cavity context. Each is framed as a **testable hypothesis**. The absence of a finding from the published literature was confirmed through web searches conducted on 2026-08-13; it is possible that very recent publications not yet indexed may have addressed some of these points.

---

## Insight 1: Ferroptosis Resistance in OSCC Operates Through the GSH Supply Chain, Not GPX4

**Observation from our data:**

The canonical "gatekeeper" of ferroptosis, **GPX4**, is **not differentially expressed** between NRF2-High and NRF2-Low OSCC tumours (log₂FC = −0.09, padj = 0.57, ρ = −0.03 with NRF2 activity). Instead, ferroptosis resistance is achieved through massive upregulation of the upstream **glutathione supply chain**: SLC7A11 (log₂FC = +2.65), GCLC (+1.23), GCLM (+1.12), GSR (+0.98), and GSS (+0.36). In parallel, the GPX4-independent defence system—**AIFM2/FSP1** (+0.31) and its substrate regenerator **NQO1** (+1.70)—is activated.

**What is known vs. what is novel:**

- The NRF2-ferroptosis axis is documented in general cancer biology (Dodson et al., 2019).
- The FSP1/AIFM2 pathway is known generically (Bersuker et al., 2019; Doll et al., 2019).
- **What has NOT been characterised in OSCC specifically** is the observation that GPX4 has zero correlation with NRF2 activity (ρ = −0.03) while the upstream supply chain shows strong dose-dependent correlation (NQO1 ρ = +0.76, SLC7A11 ρ = +0.68, GCLC ρ = +0.54). This dissociation — supply chain up, effector enzyme flat — has not been quantified in OSCC patient tumour data. It reframes the therapeutic strategy: targeting SLC7A11 (erastin, sulfasalazine) should be prioritised over targeting GPX4 (RSL3) in NRF2-high oral cancers.

**Testable Hypothesis:**

> **H1:** In NRF2-high OSCC cell lines (e.g., HSC-4, SAS), SLC7A11 inhibition (erastin/sulfasalazine) will induce ferroptosis more effectively than GPX4 inhibition (RSL3), and the combination of SLC7A11 blockade + FSP1 inhibition (iFSP1) will produce synergistic cell death, reflecting the dual-layer protection revealed by this transcriptomic analysis.

---

## Insight 2: TNFSF18 (GITRL) Is Upregulated in the Immune-Cold NRF2-High Microenvironment — A Paradoxical Co-Stimulatory Signal

**Observation from our data:**

**TNFSF18** (GITRL, the ligand for the GITR co-stimulatory receptor) is among the top 30 most significantly upregulated DEGs in NRF2-High tumours (log₂FC = +2.46, padj = 7.2×10⁻²⁰). This is paradoxical: TNFSF18/GITRL is a **T-cell co-stimulatory molecule** typically associated with immune activation, yet it is strongly upregulated in tumours that simultaneously show profound pan-immune suppression (downregulated CD3E, CD4, CD8A/B, CXCL9/10/11, GZMB, PRF1).

**What is known vs. what is novel:**

- GITR/GITRL biology is studied in immunotherapy contexts generally.
- NRF2-driven immune exclusion is documented in lung adenocarcinoma (Best et al., 2018; Skoulidis et al., 2018).
- **What has NOT been described** is the specific co-occurrence of TNFSF18 upregulation with NRF2-driven immune exclusion in any cancer type. Literature searches return no direct NRF2→TNFSF18 regulatory link. The upregulation could reflect: (a) a tumour-intrinsic GITRL expression that paradoxically promotes Treg expansion (maintaining immunosuppression), (b) an NRF2-driven epithelial differentiation program that incidentally includes TNFSF18 (since GITRL is expressed by epithelial cells), or (c) a compensatory signal in response to absent immune infiltration. This has not been mechanistically dissected in OSCC.

**Testable Hypothesis:**

> **H2:** TNFSF18/GITRL upregulation in NRF2-high OSCC tumours preferentially expands regulatory T cells (Tregs) rather than effector T cells, contributing to the immune-cold phenotype. This can be tested by co-culturing NRF2-high OSCC cell lines with mixed T-cell populations and measuring Treg vs. effector T-cell proportions ± GITRL-blocking antibody.

---

## Insight 3: The ALOX15 Paradox — Pro-Ferroptotic Enzyme Upregulated Alongside Anti-Ferroptotic NRF2 Targets Creates a "Primed Instability" State

**Observation from our data:**

**ALOX15** (15-lipoxygenase), a pro-ferroptotic enzyme that oxidises membrane phospholipids, is significantly **upregulated** (log₂FC = +1.25, padj = 1.0×10⁻⁴) alongside the anti-ferroptotic NRF2 programme. Additionally, **ALOX12** (+1.44) and **ALOX12B** (+1.57) are also upregulated. This creates a state where the cell simultaneously produces more lipid peroxides AND more antioxidant defence.

**What is known vs. what is novel:**

- ALOX15's role in ferroptosis is established generically.
- NRF2-ALOX15 co-upregulation has been noted as a theoretical concept in cancer biology.
- **What has NOT been characterised in OSCC** is the simultaneous upregulation of three lipoxygenase family members (ALOX15, ALOX12, ALOX12B) creating what we term a "primed instability" state. This means NRF2-high OSCC tumours exist in a delicate redox equilibrium — they generate abnormally high levels of lipid peroxides but counterbalance them with massive antioxidant capacity. Therapeutic disruption of the NRF2 antioxidant programme would unleash pre-existing lipid peroxidation, potentially triggering catastrophic ferroptosis far more effectively than in tumours without elevated lipoxygenase expression.

**Testable Hypothesis:**

> **H3:** NRF2-high OSCC cells with elevated ALOX12/15 expression will undergo ferroptosis more rapidly and completely upon NRF2 inhibition (e.g., brusatol, ML385) compared to NRF2-low cells, because the pre-existing lipid peroxidation load is unmasked. Baseline lipid peroxidation levels (measured by C11-BODIPY or MDA) will be paradoxically higher in NRF2-high cells despite their viability, confirming the "primed instability" state.

---

## Insight 4: Neurotensin (NTS) Is Among the Most Strongly Upregulated Genes — An Unstudied NRF2-Neuropeptide Axis in OSCC

**Observation from our data:**

**NTS** (neurotensin) is the 9th most significant DEG in the entire analysis (log₂FC = +4.34, padj = 1.4×10⁻²¹). This is a 20-fold increase in expression. NTS encodes a 13-amino acid neuropeptide known to promote cell proliferation, survival, and migration in pancreatic and colon cancer through NTSR1-mediated signalling (activating MAPK, PI3K/AKT, and NF-κB).

**What is known vs. what is novel:**

- NTS/NTSR1 signalling is studied in pancreatic, breast, lung, and colon cancers.
- NRF2 has NO established regulatory relationship with NTS in any published study, in any cancer type.
- **NTS has not been studied in the context of OSCC biology at all**, and its association with the NRF2-high molecular subtype is entirely novel. Possible mechanistic explanations include: (a) NTS may be a direct or indirect NRF2 transcriptional target (ARE motifs in the NTS promoter have not been investigated), (b) NTS upregulation may reflect the squamous differentiation programme (neuroendocrine-like features in well-differentiated SCC), or (c) NTS-NTSR1 autocrine signalling may contribute to the proliferative advantage of NRF2-high tumours.

**Testable Hypothesis:**

> **H4a:** NTS expression is directly regulated by NRF2 in OSCC. This can be tested by NRF2 knockdown/overexpression in OSCC cell lines followed by NTS qPCR/Western blot, and ChIP-seq for NRF2 binding at the NTS locus.
>
> **H4b:** NTS-NTSR1 autocrine signalling contributes to NRF2-high OSCC proliferation or migration. This can be tested by NTSR1 antagonist (SR-48692) treatment of NRF2-high vs. NRF2-low OSCC cell lines and measuring proliferation, migration, and downstream MAPK/AKT activation.

---

## Insight 5: CES1 (Carboxylesterase 1) Is Massively Upregulated — A Potential NRF2-Driven Prodrug Resistance Mechanism

**Observation from our data:**

**CES1** is the single most highly upregulated gene by fold change in the entire DEG list (log₂FC = +4.37, padj = 3.3×10⁻³³). CES1 is a serine esterase that hydrolyses ester- and amide-bond-containing drugs, including the chemotherapeutic prodrugs capecitabine and irinotecan. It also metabolises endogenous lipids (cholesterol esters, triacylglycerols).

**What is known vs. what is novel:**

- CES1's role in drug metabolism is established in pharmacology.
- CES1 is known to be regulated by nuclear receptors (PXR, CAR) in hepatocytes.
- **CES1 has NOT been identified as an NRF2 target gene**, and its massive upregulation (>20-fold) in NRF2-high OSCC has not been reported. This is notable because CES1 upregulation could confer resistance to ester-bond-containing chemotherapeutics and may also participate in lipid remodelling relevant to ferroptosis (by altering cholesterol ester pools). The co-upregulation of CES1 with the canonical NRF2 detoxification programme (GSTs, UGTs, CYPs, AKRs) suggests CES1 may be part of a broader, previously unrecognised NRF2-driven Phase I metabolic programme in squamous epithelium.

**Testable Hypothesis:**

> **H5:** CES1 is transcriptionally regulated by NRF2 through ARE motifs, and its upregulation in NRF2-high OSCC confers resistance to ester-bond-containing prodrugs. This can be tested by: (a) NRF2 ChIP-qPCR at the CES1 promoter, (b) CES1 expression changes upon NRF2 silencing, and (c) CES1 inhibitor (bis-para-nitrophenyl phosphate) sensitisation of NRF2-high OSCC cells to capecitabine.

---

## Insight 6: NRF2-Driven Squamous Differentiation Lock With Simultaneous EMT Suppression — A Mechanism Specific to OSCC

**Observation from our data:**

NRF2-high tumours show the most significant upregulation in **keratinization** (padj = 3.1×10⁻³⁰, 33 genes) and simultaneously the strongest suppression of **EMT** (NES = −2.34). Every canonical EMT transcription factor is downregulated: SNAI1 (−0.51), TWIST1 (−0.80), ZEB1 (−0.68), ZEB2 (−0.70), along with mesenchymal markers CDH2 (−2.04) and VIM (−0.80). Meanwhile, multiple transglutaminases are massively upregulated (TGM3 = +2.89, TGM7 = +2.48, TGM6 = +2.37, TGM1 = +1.51, TGM5 = +1.06) — enzymes that catalyse protein cross-linking in the cornified envelope.

**What is known vs. what is novel:**

- NRF2's role in keratinocyte differentiation is established in skin biology (Keap1-null mice develop hyperkeratosis).
- NRF2-EMT crosstalk is studied, but often in the context of NRF2 *promoting* EMT (in pancreatic and breast cancers).
- **What has NOT been described in OSCC** is the coordinated upregulation of **five transglutaminase family members** (TGM1/3/5/6/7), establishing a "differentiation lock." This multi-TGM signature has not been reported as an NRF2-associated feature in any cancer. The clinical implication is that NRF2-high OSCC tumours may be well-differentiated and less prone to distant metastasis, but are simultaneously chemoresistant and immune-evasive — a distinct clinical profile from the poorly-differentiated, EMT-driven aggressive OSCC subtype.

**Testable Hypothesis:**

> **H6:** NRF2 directly activates transglutaminase gene transcription (particularly TGM1 and TGM3) through ARE sequences, and this locks OSCC cells into a differentiated squamous state that is resistant to EMT induction by TGFβ. This can be tested by: (a) NRF2 ChIP-seq at TGM promoters, (b) TGFβ treatment of NRF2-high vs. NRF2-low cell lines measuring EMT markers, and (c) correlation of TGM expression with histological differentiation grade in the TCGA-HNSC cohort.

---

## Summary Table

| # | Insight | Novelty Level | Key Data Point | Hypothesis |
|---|---------|---------------|----------------|------------|
| 1 | GPX4-independent ferroptosis resistance via SLC7A11-GCLC supply chain + AIFM2/NQO1 backup | **OSCC-specific quantification is novel** | GPX4 ρ = −0.03 vs NQO1 ρ = +0.76 | SLC7A11 inhibition > GPX4 inhibition in NRF2-high OSCC |
| 2 | TNFSF18/GITRL upregulation in immune-cold NRF2-high tumours | **Novel in any cancer type** | TNFSF18 log₂FC = +2.46 amid pan-immune suppression | GITRL expands Tregs, reinforcing immunosuppression |
| 3 | ALOX12/12B/15 co-upregulation: "primed instability" | **OSCC-specific multi-ALOX pattern is novel** | ALOX15 +1.25, ALOX12 +1.44, ALOX12B +1.57 | NRF2 inhibition triggers catastrophic ferroptosis |
| 4 | NTS (neurotensin) massively upregulated | **Entirely novel — no NRF2-NTS link in literature** | NTS log₂FC = +4.34 (20-fold) | NTS-NTSR1 autocrine loop in NRF2-high OSCC |
| 5 | CES1 as putative NRF2 target | **Novel — CES1 not known as NRF2 target** | CES1 log₂FC = +4.37 (highest FC) | CES1 confers ester-prodrug resistance |
| 6 | Multi-TGM differentiation lock | **Novel multi-TGM signature in any cancer** | 5 TGMs upregulated + all EMT TFs downregulated | NRF2→TGMs locks squamous differentiation |

> [!NOTE]
> **Insights NOT included because they are well-established:**
> - NRF2-driven chemoresistance via GST/CYP/ABCC drug metabolism (extensively published in OSCC)
> - NRF2-driven immune suppression via NF-κB interference (published in lung cancer, emerging in HNSCC)
> - NRF2 as a driver of OXPHOS + antioxidant defence simultaneously (published in multiple cancer types, including recent OSCC metastasis studies)
> - NRF2 regulation of the pentose phosphate pathway for NADPH supply (Mitsuishi et al., 2012, *Cancer Cell*)
> - SLC7A11 as a ferroptosis vulnerability in HNSCC (studied in HNSCC cell lines)

---

*Report generated from TCGA-HNSC oral cavity subset (n = 224 tumour samples, 112 NRF2-High, 112 NRF2-Low). DEG criteria: |log₂FC| ≥ 1.0, padj < 0.05 (BH-corrected). Enrichment analyses: ORA (GO/KEGG/Reactome), GSEA (Hallmark/C2/C5), GSVA (20 custom + Hallmark sets).*
