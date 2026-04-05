# Radiation Dose-Responsive Transcriptomic Remodeling in Human PBMC Reveals Reproducible Neurodegeneration-Relevant System Signatures Without Evidence of a Broad Parkinson Program

## Abstract

Ionizing radiation is a strong systems-level stressor, but its peripheral blood transcriptomic consequences are often interpreted through narrow disease-centric hypotheses. Here, we used an exploratory systems-biology framework to determine which molecular systems are dose-responsive in human PBMC, which changes are neurodegeneration-relevant, and how reproducible these signals are across analytic methods and donors.

We analyzed 25 Agilent microarrays (5 donors, 5 radiation doses: 0, 0.5, 2, 5, and 8 Gy; 48 h). Preprocessing included normexp background correction, log2 transformation, quantile normalization, and probe-to-gene collapsing. Dose effects were modeled with donor-aware linear modeling. Pathway-level evidence was triangulated using GSVA (sample-level pathway activity), FGSEA (rank-based enrichment), and CAMERA (competitive enrichment accounting for inter-gene correlation). Multiple testing was controlled with false discovery rate (FDR). Reproducibility was evaluated by cross-method concordance and leave-one-donor-out (LODO) stability analyses.

Across the global module scan, 281 modules were parsed, 278 mapped, and 236 tested. Significant modules were identified by GSVA (75), FGSEA (63), and CAMERA (32), with 28 modules significant across all three methods. Using explicit GSVA-based classification, 36 modules were enriched (positive dose slope, FDR < 0.05), 39 were impoverished (negative dose slope, FDR < 0.05), and 57 showed low-effect/non-informative behavior in this dataset. Robust enriched examples included ISR death-factor signaling (GSVA logFC/Gy = 0.0836, FDR = 0.0019) and antigen presentation (0.0641, FDR = 0.0020). Robust impoverished examples included cyto-ribosome assembly (-0.1067, FDR = 0.0013) and mtDNA replication (-0.0906, FDR = 0.0019). PD-focused analysis remained selective: Familial PD genes were significant (0.04548, FDR = 0.01056), while broad PD-all was not. Cross-method concordance was high (Spearman rho = 0.908 globally; rho = 0.904 in neurodegeneration-relevant subset). LODO analyses identified modules with perfect sign consistency across donor removals.

Radiation induces reproducible, dose-structured transcriptomic reprogramming in PBMC characterized by concurrent activation of stress/inflammatory pathways and suppression of maintenance/proteostasis pathways. Neurodegeneration-relevant signals are present but selective, supporting a molecular vulnerability interpretation rather than a broad disease-state claim in peripheral blood.

## Introduction
Radiation biology in peripheral blood is frequently reduced to binary narratives, such as whether a specific disease signature is or is not present. That framing can obscure the core systems question: how does dose-dependent stress reorganize cellular programs, and which of those programs map onto neurodegeneration-relevant biology [1,2]? PBMC provide a useful peripheral readout of immune-metabolic stress adaptation, but they are not a proxy for direct neuronal pathology. Consequently, the strongest scientific contribution from this dataset is expected to be mechanistic patterning, not diagnosis-level inference [3].

Previous work has shown that stress systems involving inflammation, mitochondrial control, and protein homeostasis can converge on neurodegeneration-relevant molecular states across tissues [4,5]. However, in small donor cohorts, robust interpretation requires explicit control of donor heterogeneity, conservative multiple-testing procedures, and independent pathway inference strategies rather than reliance on a single method [6]. These considerations are central for avoiding overinterpretation.

This study was designed as an exploratory but statistically disciplined systems analysis of OSD-157/GLDS-157 PBMC data. The primary question was not whether radiation causes Parkinson disease, but rather which molecular systems are dose-responsive, which changes are neurodegeneration-relevant, and how reproducible those findings are. We therefore integrated donor-aware dose modeling with three complementary pathway frameworks (GSVA, FGSEA, CAMERA), followed by reproducibility analyses including cross-method concordance and leave-one-donor-out stability [7,8].

The resulting framework allows separation of robust findings from exploratory ones and supports biologically meaningful claims without clinical overreach.

## Methods
### Study design and data
The dataset comprised 25 Agilent single-channel PBMC microarrays (5 donors: a, b, X, Y, Z; 5 doses: 0, 0.5, 2, 5, 8 Gy; 48 h). Module definitions were derived from the PD2_Mm resource, and analyses were performed on both PD-focused subsets and global module coverage.

### Preprocessing and expression matrix construction
Raw intensities were processed using a standard microarray workflow: normexp background correction with offset, log2 transformation, quantile normalization, and probe-level collapsing to gene-level expression estimates. This sequence was used to reduce technical variance and ensure comparability across arrays before dose modeling.

### Donor-aware dose modeling
Dose was modeled as a continuous variable, with donor adjustment to account for baseline inter-individual differences. Donor-aware models were used throughout the pipeline to reduce confounding between dose effects and donor structure. In the gene-level and pathway-level trend frameworks, effect size was represented as slope per Gy (reported as GSVA logFC/Gy in pathway analyses).

### Pathway inference: rationale for triangulation
Three complementary methods were used because they answer different statistical questions.

1. GSVA estimates pathway activity per sample and then tests dose trends at pathway level.
2. FGSEA evaluates whether pathway genes are enriched among strongest ranked gene-level effects.
3. CAMERA provides a competitive test that accounts for inter-gene correlation within pathways.

Agreement among these methods was interpreted as stronger evidence than significance from any single method.

### Multiple testing and evidence classes
False discovery rate (FDR) control was applied to pathway-level analyses. For transparent classification in the global scan:

1. Enriched pathways: GSVA slope > 0 and GSVA FDR < 0.05.
2. Impoverished pathways: GSVA slope < 0 and GSVA FDR < 0.05.
3. Low-effect/non-informative pathways: absolute GSVA slope < 0.01, GSVA FDR > 0.25, FGSEA FDR > 0.25, and CAMERA FDR > 0.25.

This conservative low-effect class was included to explicitly report what the dataset does not strongly inform.

### PD-focused and neurodegeneration-relevant analyses
A PD-focused module set was analyzed separately to avoid dilution by global effects. A clean PD risk score (excluding mitophagy-overlap genes) was modeled against dose with donor adjustment. Neurodegeneration-relevant subsets in global analyses were defined using pathway-name filters consistent with PD, AD, ALS, HD, prion, and related mechanistic domains (mitochondrial, proteostasis, DNA repair, inflammatory signaling).

### Reproducibility analyses
Two reproducibility layers were performed.

1. Cross-method concordance: Spearman correlation between GSVA effect slopes and FGSEA NES across modules.
2. Leave-one-donor-out (LODO): pathway slope recalculated after removing each donor in turn; stability summarized by sign consistency and slope dispersion.

Confidence-interval-oriented visualization was generated for top neurodegeneration-relevant modules to pair effect size with uncertainty.

### Reporting conventions
All major claims were tied to quantitative outputs from workspace evidence tables and figures. Where data were unavailable (for example, external clinical outcomes), no inference was made.

## Results
### Dataset coverage and global signal burden
Global scan coverage was high: 281 modules parsed, 278 mapped, and 236 testable modules. Significant burden differed by method but remained substantial: GSVA identified 75 significant modules (FDR < 0.05), FGSEA 63, and CAMERA 32. Notably, 28 modules were significant across all three methods, indicating a robust high-confidence core.

Using explicit GSVA-centered operational criteria, modules were classified as follows: 36 enriched, 39 impoverished, and 57 low-effect/non-informative. This distribution supports broad systems remodeling rather than isolated pathway perturbation.

### Enriched pathways: robust upshifts in stress and inflammatory architecture
Among enriched modules, the largest positive GSVA slopes included ISR target-gene death factors (0.0836/Gy, FDR = 0.0019), adaptive immune antigen presentation (0.0641/Gy, FDR = 0.0020), and inflammasome-linked innate signatures (0.0638/Gy, FDR = 0.0097). Additional enriched modules included RAAS/PANapoptosis and tissue-damage signaling blocks.

These convergent findings indicate dose-linked amplification of stress signaling, inflammatory coordination, and cell-fate decision pathways. Because these pathways are central to injury response and chronic stress adaptation, the pattern is biologically plausible under radiation exposure and relevant to neurodegeneration-linked stress biology [2,4].

### Impoverished pathways: suppression of maintenance and homeostatic programs
The strongest impoverished modules included cyto-ribosome assembly (-0.1067/Gy, FDR = 0.0013), mtDNA replication (-0.0906/Gy, FDR = 0.0019), ribosomal 60S and broader ribosome subunit modules (both strongly negative with low FDR), and DNA-repair NHEJ (-0.0821/Gy, FDR = 0.0042).

This pattern is consistent with a stress-adaptation regime in which synthetic and maintenance programs are down-prioritized while stress-defense and damage-signaling programs are upregulated. Such coupling of inflammatory/stress activation with proteostasis decline is frequently discussed in vulnerability frameworks for degenerative biology [4,5], though this remains an interpretation rather than proof of pathology in this tissue context.

### Low-effect/non-informative pathways in this dataset
A conservative low-effect class (n = 57) captured pathways with near-zero GSVA slopes and no convergent significance across GSVA, FGSEA, and CAMERA. Examples included specific OXPHOS-related and prion-labeled modules with very small effect sizes and high FDR values. These were explicitly retained to avoid publication bias toward positive signals and to clarify that absence of signal here means low informativeness under current design, not biological impossibility.

### PD-focused findings as a subset of systems response
In the PD-focused phase, Familial PD genes showed a statistically robust positive dose trend (GSVA logFC/Gy = 0.04548; FDR = 0.01056), whereas broad PD-all modules remained non-significant with small effects. Thus, neurodegeneration relevance was selective and modular rather than diffuse across all PD-labeled biology.

The clean PD risk score showed a positive slope (0.04095/Gy) but did not meet conventional significance (p = 0.0743; SE = 0.02168; t = 1.8886). This supports a cautious interpretation: directional molecular tendency without strong inferential certainty for a global PD-risk score in this cohort size.

Within ND subgroup pathways, ALS-labeled modules showed the clearest GSVA significance in the subgroup panel (logFC = 0.01739; FDR = 0.04056), while AD, HD, PD-all, and prion subgroup terms were weaker or non-significant in this design.

### Cross-method concordance and donor-stability evidence
Cross-method effect concordance was high. Spearman correlation between GSVA slope and FGSEA NES was 0.908 across all tested modules and 0.904 within neurodegeneration-relevant subsets. This level of agreement indicates that major directional conclusions are not artifacts of one method.

LODO analyses further supported reproducibility for key pathways. Modules such as cyto-ribosome assembly, mtDNA replication, ISR death factors, and NHEJ retained stable slope direction under donor removal, with sign consistency reaching 1.0 in representative cases and modest slope dispersion. This provides internal robustness against dominance by any single donor profile.

### Robust versus exploratory conclusions
Robust conclusions are those supported by effect size magnitude, FDR control, cross-method coherence, and donor-stability checks. Exploratory conclusions are those with directional trends but incomplete support across these criteria. This distinction was applied throughout and is central to interpretation integrity.

## Discussion
This analysis demonstrates that radiation exposure in PBMC drives a reproducible systems-level transcriptomic shift characterized by two coordinated axes: upregulation of stress/inflammatory and cell-fate pathways, and downregulation of proteostasis/translation and selected maintenance programs. The key contribution is therefore not a binary disease claim but a coherent molecular architecture of stress vulnerability.

The PD-focused signal was real but selective. Familial PD module enrichment occurred within a broader systems response dominated by non-PD pathway remodeling. This reinforces the view that disease-labeled pathways can behave as mechanistic components of generalized stress adaptation rather than as direct evidence of disease-state equivalence in peripheral blood.

Method triangulation and reproducibility checks materially strengthen the findings. High GSVA-FGSEA concordance and stable LODO behavior for top modules suggest that the main directional biology is robust to method and donor perturbation. In small-cohort systems datasets, this form of consistency may be more informative than any single p-value threshold.

Biologically, the coexistence of inflammatory/stress amplification with proteostasis and mitochondrial-maintenance suppression is compatible with vulnerability models discussed in neurodegeneration literature [4,5]. However, this should remain an interpretation tier, not a clinical conclusion. PBMC signatures can indicate systemic stress states relevant to neurodegeneration frameworks, but they cannot establish neuronal pathology, diagnosis, or individual prognosis.

Overall, the evidence supports a conservative but meaningful claim: radiation induces a dose-structured peripheral molecular state that intersects neurodegeneration-relevant biology in a selective, reproducible, and systemically coherent manner.

## Limitations
1. Donor sample size was small (n = 5), limiting precision for weaker effects and subgroup interpretations.
2. PBMC is a peripheral compartment and cannot substitute for direct neural tissue inference.
3. Pathway overlap can induce dependence among tests and can blur semantic boundaries between modules.
4. Ortholog-sensitive interpretation was performed within available workspace evidence; comprehensive external remapping pipelines were not fully available in current runtime evidence.
5. Clinical endpoints were not available in current workspace evidence, preventing translational risk calibration.

## Conclusions
Radiation produces robust dose-responsive transcriptomic remodeling in human PBMC. The dominant pattern is systems-level: stress/inflammatory and cell-fate modules are enriched, while maintenance/proteostasis programs are impoverished. Neurodegeneration-relevant signals are present as selective modules embedded in this broader architecture, not as a global Parkinson-like state.

These findings support a reproducible molecular vulnerability framework for peripheral stress biology. They do not support diagnostic or clinical causality claims from PBMC alone. Future work should prioritize external replication, multi-timepoint designs, and cross-tissue integration to test persistence and translational relevance.

## Data and code availability statement
All analyses in this article were derived from evidence available in the current workspace, including processed tables, reports, and scripts. Core evidence includes PD-focused summary tables, global module-scan outputs, and validated reproducibility tables. Analysis scripts implementing preprocessing, donor-aware modeling, pathway triangulation, global scanning, and reproducibility checks are available in the same project repository. External clinical outcome data were not available in current workspace evidence.

## Suggested figure legends
1. Study design and analytic framework. Overview of cohort structure (5 donors, 5 doses, 25 samples), preprocessing workflow, donor-aware modeling, and pathway triangulation (GSVA, FGSEA, CAMERA), followed by reproducibility checks.
2. Global pathway effect landscape. Volcano-style display of GSVA dose slopes versus statistical significance across 236 tested modules, highlighting enriched, impoverished, and low-effect classes.
3. Category-level systems remodeling. Mean pathway effect by biological category, illustrating coordinated upregulation of stress/inflammatory domains and downregulation of maintenance/proteostasis domains.
4. Top neurodegeneration-relevant pathway effects with uncertainty. Forest plot of GSVA dose effects with 95% confidence intervals for top neurodegeneration-relevant modules; pathways not crossing zero indicate stronger directional evidence.
5. Cross-method concordance of pathway direction. Scatter of GSVA slopes versus FGSEA NES across modules, with reported Spearman correlation for all modules and neurodegeneration-relevant subset.
6. Donor-level robustness by leave-one-donor-out. Heatmap/scatter showing pathway slope stability after iteratively excluding one donor, with sign-consistency metrics and slope dispersion summaries.

## Suggested supplementary materials list
1. Full table of module-level statistics for all tested pathways (GSVA, FGSEA, CAMERA, donor summaries).
2. Operational definitions and thresholds for enriched, impoverished, and low-effect/non-informative classes.
3. PD-focused pathway evidence table with module-level donor consistency.
4. Clean PD risk-score model output with full coefficients and standard errors.
5. Neurodegeneration subgroup table (ALS, AD, HD, PD-all, prion terms) with multi-method support columns.
6. LODO full slope matrix and stability summary per pathway.
7. Sensitivity notes for pathway naming/grouping and category assignment.
8. Script manifest with execution order and expected outputs.
9. Additional diagnostic plots for p-value distributions and effect-size distributions.
10. Reproducibility metadata (software versions, seeds, environment notes where available in workspace evidence).

## Quality-control checklist
1. Internal consistency: Yes. Numeric claims are aligned with workspace-derived evidence.
2. Statistics correctly interpreted: Yes. Effect size, uncertainty, FDR, and trend interpretations are separated and conservatively phrased.
3. No clinical overclaims: Yes. No diagnostic or individual risk claims were made from PBMC evidence.
4. Reproducibility statements included: Yes. Cross-method concordance and LODO stability were explicitly reported.
5. Missing information explicitly flagged: Yes. Clinical endpoints and certain external mapping details were marked as not available in current workspace evidence.
