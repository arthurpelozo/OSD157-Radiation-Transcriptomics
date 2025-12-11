# OSD-157: Transcriptomic Analysis of Ionizing Radiation Response
**Author:** Arthur Pelozo
**Date:** December 2025 
**Platform:** Agilent Whole Human Genome Microarray (Single-Channel)

## Project Overview
This repository contains the bioinformatics analysis pipeline and results for the NASA GeneLab dataset **OSD-157**. The study evaluates the transcriptomic response of human peripheral blood mononuclear cells exposed *ex vivo* to varying doses of Gamma radiation (0 to 8 Gy) after 48 hours.

**Primary Goal:** To validate the radiation dose-response using canonical pathways (p53, Inflammation).
**Secondary Goal:** To test the hypothesis that acute radiation induces a molecular signature similar to **Parkinson's Disease (PD)** in leukocytes.

## Key Findings
1. **Validation:** Strong activation of DNA damage response (p53) and cell cycle arrest (G2/M) pathways, confirming the biological efficacy of the radiation.
2. **Hypothesis Test:** No statistically significant enrichment was found for the KEGG Parkinson's Disease pathway (NES = 0.79, FDR = 0.95), suggesting no direct molecular link in this tissue context.

## Dataset Structure
* **Platform:** Agilent Whole Human Genome Microarray 4x44K
* **Design:** 5 Donors (Biological Replicates) x 5 Doses (0, 0.5, 2, 5, 8 Gy).
* **Sample Count:** 25 samples.

## Analysis Pipeline
The analysis was performed in R using `limma` for linear modeling and `fgsea` for functional enrichment.
1. **Normalization:** Background correction (`normexp`) + Quantile Normalization.
2. **Modeling:** Paired design (`~ Dose + Donor`) to account for high inter-individual variability.
3. **GSEA:** Multilevel enrichment analysis using Hallmark and KEGG collections.

## Repository Contents
* `scripts/`: R source code for the analysis pipeline.
* `results/`: Processed tables of Differentially Expressed Genes (DEGs) and GSEA results.
* `docs/`: Full scientific report.

Data set: https://osdr.nasa.gov/bio/repo/data/studies/OSD-157

GLDS-157_micoarray_E-GEOD-44201.raw.3.zip
92.22 MB
Tue Dec 05 2017

GLDS-157_micoarray_E-GEOD-44201.raw.1.zip
178.25 MB
Tue Dec 05 2017

GLDS-157_micoarray_E-GEOD-44201.raw.2.zip
180.83 MB
Tue Dec 05 2017


Abstract

This study evaluated global gene expression changes in human peripheral blood samples exposed ex vivo to increasing doses of gamma radiation (0 to 8 Gy) following a 48-hour incubation period. Using the NASA GeneLab OSD-157 dataset, a robust bioinformatics pipeline was applied to perform normalization, linear modeling, and Gene Set Enrichment Analysis (GSEA). The results demonstrate a severe transcriptomic response characterized by the activation of the p53 signaling pathway and cell cycle arrest. Furthermore, the hypothesis that ionizing radiation induces molecular signatures associated with Parkinson’s Disease (PD) in leukocytes was specifically tested. No statistically significant enrichment was observed for the PD pathway (NES = 0.79, FDR = 0.95), suggesting an absence of molecular overlap in this acute experimental model.
