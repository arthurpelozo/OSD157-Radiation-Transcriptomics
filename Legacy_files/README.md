# OSD-157: Transcriptomic Analysis of Ionizing Radiation Response

**Author:** Arthur Pelozo   
**Date:** December 2025   
**Platform:** Agilent Whole Human Genome Microarray (Single-Channel)

## Project Overview

This repository contains the bioinformatics analysis pipeline and results for the NASA GeneLab dataset **OSD-157**. The study evaluates the transcriptomic response of human peripheral blood mononuclear cells exposed *ex vivo* to varying doses of Gamma radiation (0 to 8 Gy) after 48 hours.

* **Primary Goal:** To validate the radiation dose-response using canonical pathways (p53, Inflammation).
* **Secondary Goal:** To investigate potential overlaps between acute radiation response and **Parkinson's Disease (PD)** mechanisms, specifically focusing on **Mitochondrial Quality Control**.

## Key Findings

### 1. Validation
* Strong activation of DNA damage response (p53) and cell cycle arrest (G2/M) pathways, confirming the biological efficacy of the radiation (**NES = +2.32**, p < 10⁻¹³).

### 2. Broad PD Hypothesis 
* **Zero Enrichment:** No statistically significant enrichment was found for the generic *KEGG Parkinson's Disease* pathway (**NES = 0.64**, p = 0.99), ruling out a broad disease signature (e.g., neurodegeneration) in blood.

### 3. Targeted Mitophagy Discovery 
Upon decomposing the pathway into custom functional modules, a specific, dose-dependent **"Push-Pull" Mitophagy Mechanism** was identified:
* **Activation (+):** Upregulation of autophagic execution markers (*SQSTM1*, *MAP1LC3B*; **logFC +0.064**, p = 0.0001) and core sensors (*PINK1*, *PRKN*).
* **Disinhibition (-):** Significant downregulation of negative regulators (Deubiquitinases *USP30*, *USP15*; **logFC -0.086**, p = 0.0001).

**Conclusion:** Radiation actively promotes mitochondrial clearance by simultaneously "pressing the gas" (machinery) and "releasing the brakes" (regulators).

### 4. Systemic Risk Analysis (The "Dual Signature")
To validate whether mitochondrial stress drives broader PD molecular features, a "Clean PD Risk Score" (excluding mitophagy genes) was modeled against Radiation Dose and Mitophagy activity:
* **Parallel Hits:** Modeling confirmed that Radiation **independently** drives both the metabolic cleanup (Mitophagy, p=0.020) and the transcriptomic vulnerability (Risk Genes, p=0.0026).
* **Dual Signature:** Unlike simple mediation, the radiation effect on risk remained robust (**p = 0.016**) even when controlling for mitophagy. This suggests the cell is under a "Double Attack": a metabolic cost to clean damage AND a direct upregulation of susceptibility genes.
## Dataset Structure

* **Source:** [NASA GeneLab (OSD-157)](https://osdr.nasa.gov/bio/repo/data/studies/OSD-157)
* **Platform:** Agilent Whole Human Genome Microarray 4x44K
* **Design:** 5 Donors (Biological Replicates) x 5 Doses (0, 0.5, 2, 5, 8 Gy).
* **Sample Count:** 25 samples.
* **Files:**
    * `GLDS-157_micoarray_E-GEOD-44201.raw.1.zip`
    * `GLDS-157_micoarray_E-GEOD-44201.raw.2.zip`
    * `GLDS-157_micoarray_E-GEOD-44201.raw.3.zip`

## Analysis Pipeline

The analysis was performed in R, evolving from standard enrichment to custom module interrogation and causal modeling.

### 1. Preprocessing
* Annotation recovery directly from raw data objects.
* Background correction (`normexp`) + Quantile Normalization.
* **Log2 Transformation:** Critical step applied to stabilize variance in raw linear data.
* Probe collapsing via `avereps`.

### 2. Differential Expression & Enrichment
* **Paired Linear Model:** `Y ~ Dose + Donor` using `limma` to isolate radiation effects from donor heterogeneity.
* **GSEA/GSVA:** Multilevel enrichment analysis (`fgsea`) and per-sample pathway scoring (`GSVA`) using Gaussian kernels.

### 3. Causal Inference Framework (The "Dual Signature" Test)
To distinguish between **Mediation** (Radiation → Mitophagy → Risk) and **Parallel Effects** (Radiation → Risk & Radiation → Mitophagy), we employed a three-step regression framework using a non-overlapping "Clean PD Risk Score":

* **Model A (Total Effect):** `Risk_Score ~ Dose + Donor`
    * *Hypothesis:* Does radiation directly increase the expression of susceptibility genes?
* **Model B (Mediator Activation):** `Mitophagy_Score ~ Dose + Donor`
    * *Hypothesis:* Does radiation actively trigger the mitochondrial cleanup machinery?
* **Model C (Independence Test):** `Risk_Score ~ Dose + Mitophagy_Score + Donor`
    * *Hypothesis:* Does the radiation effect on risk disappear when controlling for mitophagy? (If yes = Mediation; If no = Dual/Parallel Signature).

## Repository Contents

This repository is organized into the following core analysis files:

* **`Clean_Full_Code.R`**: 
  The complete, reproducible R script containing the entire pipeline (Data Loading, Log2 Normalization, Limma Modeling, GSEA/GSVA, and Mediation Analysis).

* **`Full_Research_Log.md`**: 
  The compiled computational notebook. Click this file to view the full analysis, statistical outputs, and visualizations directly on GitHub.

* **`Full_Research_Log_files/`**: 
  Folder automatically generated containing the high-resolution plots (Heatmaps, Barplots, Regression plots) linked within the `.md` log.

* **`OSD157_Report.txt`**: 
  A plain text executive summary listing the study's specific objectives, methodology, and the final scientific conclusions regarding the "Dual Signature".

* **`README.md`**: 
  Project overview and documentation (this file).

## Abstract

Ionizing radiation is a potent environmental stressor, yet its potential to induce molecular signatures associated with neurodegenerative vulnerability, specifically Parkinson’s Disease (PD), remains understudied in non-neuronal tissues. This study investigated the transcriptomic response of human peripheral blood mononuclear cells exposed *ex vivo* to gamma radiation (0–8 Gy) to determine if acute exposure triggers specific PD-related biological mechanisms.

Using the NASA GeneLab OSD-157 dataset, we applied a rigorous bioinformatics pipeline incorporating background correction, **Log2 transformation**, and quantile normalization. Data were analyzed using paired linear modeling (`limma`), Gene Set Enrichment Analysis (GSEA), and Gene Set Variation Analysis (GSVA). A "Clean PD Risk Score" (comprising 15 non-mitophagy susceptibility genes) was constructed to test the causal relationship between mitochondrial stress and systemic molecular vulnerability.

The efficacy of the radiation exposure was validated by robust activation of the p53 signaling pathway (**NES = +2.42**, *p* < 10⁻¹²) and a massive induction of lysosomal machinery (**NES = +2.82**, *p* < 10⁻¹⁷). While broad PD modules associated with bioenergetic collapse (Complex I) and synaptic degeneration showed no significant enrichment, a targeted analysis revealed the specific activation of the mitochondrial quality control system (**NES = +1.72**, *p* = 0.019). This response followed a coordinated **"Push-Pull" mechanism**: the transcriptional upregulation of autophagic execution markers (*SQSTM1*, *MAP1LC3B*; **logFC +0.064**, *p* = 0.0001) occurred concurrently with the significant repression of negative regulators (*USP30*, *USP15*; **logFC -0.086**, *p* = 0.0001), effectively lowering the threshold for mitophagy.

Statistical modeling revealed a **"Dual Vulnerability Signature."** Radiation dose linearly drives the expression of independent PD risk genes (**Coefficient = +0.012**, *p* = 0.0026), establishing a "transcriptomic scar." Multivariate analysis demonstrated that this risk accumulation occurs in parallel to, and independent of, the mitophagic response (*p* = 0.016 in the combined model). We conclude that ionizing radiation acts as a systemic stressor that simultaneously imposes a metabolic cleanup cost while independently upregulating a latent transcriptomic profile of genetic vulnerability with a baseline increase in genetic vulnerability of ~0.85% per Gray.

