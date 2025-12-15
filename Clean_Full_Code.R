# ==============================================================================
# PROJECT: Radiation Analysis (GLDS-157 / Agilent Microarray)
# AUTHOR: Arthur Pelozo
# STATUS: CORRECTED (Added Log2 Transformation)
# ==============================================================================

# 1. SETUP AND LIBRARIES -------------------------------------------------------

setwd("C:/Users/arthu/Downloads/projeto_radiacao_PD")

# Load libraries
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(sva)
  library(biomaRt)
  library(fgsea)
  library(msigdbr)
  library(GSVA)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(cowplot)
  library(data.table)
})

# 2. METADATA CREATION ---------------------------------------------------------

# List files and create targets
files <- list.files("data", pattern = "\\.txt$", full.names = TRUE)
targets <- data.frame(FileName = files, BaseName = basename(files), stringsAsFactors = FALSE)

# Filter for experimental files
targets <- targets[grepl("_[A-Za-z]_[0-9\\.]+Gy", targets$BaseName), ]

# Extract Metadata
targets$Dose <- str_extract(targets$BaseName, "[0-9\\.]+Gy")
targets$DoseNum <- as.numeric(str_extract(targets$Dose, "[0-9\\.]+"))
targets$Donor <- str_extract(targets$BaseName, "_([A-Za-z])_")
targets$Donor <- gsub("_", "", targets$Donor)
targets$SampleID <- gsub("\\.txt$", "", targets$BaseName)

# Verify Design
print(table(targets$Donor, targets$DoseNum))

# 3. DATA IMPORT AND NORMALIZATION (THE FIX) -----------------------------------

message("Reading Agilent files...")
RG <- read.maimages(targets$FileName, source = "agilent", green.only = TRUE)

# A. Background Correction (NormExp + Offset)
# Offset=50 is standard to ensure we don't log(negative) or log(0) later
RG.bc <- backgroundCorrect(RG, method = "normexp", offset = 50)

# B. *** THE CRITICAL FIX: LOG2 TRANSFORMATION ***
# We log-transform BEFORE normalization to stabilize variance.
message("Applying Log2 transformation...")
E_log <- log2(RG.bc$E)

# C. Quantile Normalization
# Now we normalize the log-values, not the raw intensities.
expr_matrix <- normalizeBetweenArrays(E_log, method = "quantile")
colnames(expr_matrix) <- targets$SampleID

# D. Verification Plot (Density)
# We plot the density to confirm values are now in the 5-16 range, not 0-60,000
plotDensities(expr_matrix, main = "Log2 Expression Densities (After Fix)", legend = FALSE)

# 4. QC: PCA -------------------------------------------------------------------

expr_no_na <- expr_matrix[complete.cases(expr_matrix), ]
pca <- prcomp(t(expr_no_na), scale. = TRUE)

df_pca <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Dose = as.factor(targets$DoseNum),
  Donor = as.factor(targets$Donor)
)

ggplot(df_pca, aes(PC1, PC2, color = Dose, shape = Donor)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA — Log2 Normalized Data")

# 5. DIFFERENTIAL EXPRESSION (LIMMA) -------------------------------------------

targets$Donor <- factor(targets$Donor)
targets$DoseNum <- as.numeric(targets$DoseNum)

# Design Matrix
design <- model.matrix(~ DoseNum + Donor, data = targets)

# Fit Model
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)

# Extract Results (Dose Trend)
res_trend <- topTable(fit, coef = "DoseNum", number = Inf, sort.by = "P")

# SANITY CHECK: Look at the new logFC and AveExpr
# AveExpr should be around 7-10. logFC should be small (e.g., 0.5, 1.0, -0.8).
# If logFC is >100, something is still wrong.
print(head(res_trend))

# Save
if (!dir.exists("results")) dir.create("results")
write.csv(res_trend, file = "results/OSD157_trend_log2_probe.csv", row.names = TRUE)

# 6. ANNOTATION AND GENE COLLAPSING --------------------------------------------

probe2gene <- RG$genes$GeneName
valid <- !is.na(probe2gene) & probe2gene != "" & 
         probe2gene != "DarkCorner" & probe2gene != "GE_BrightCorner"

expr_valid <- expr_matrix[valid, ]
gene_valid <- probe2gene[valid]

# Collapse to Gene Level
expr_by_gene <- limma::avereps(expr_valid, ID = gene_valid)

# Re-run Limma on Gene Level
fit_gene <- lmFit(expr_by_gene, design)
fit_gene <- eBayes(fit_gene)
gene_ranking <- fit_gene$t[, "DoseNum"]

# Filter for HGNC
valid_hgnc <- keys(org.Hs.eg.db, keytype = "SYMBOL")
ranking_valid <- gene_ranking[names(gene_ranking) %in% valid_hgnc]
ranking_valid <- sort(ranking_valid, decreasing = TRUE)

# 7. GLOBAL GSEA (RE-CHECK) ----------------------------------------------------

m_df <- msigdbr(species = "Homo sapiens")
hallmark_list <- split(m_df[m_df$gs_collection == "H", ]$gene_symbol, 
                       m_df[m_df$gs_collection == "H", ]$gs_name)

kegg_df <- m_df[m_df$gs_subcollection %in% c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"), ]
kegg_list <- split(kegg_df$gene_symbol, kegg_df$gs_name)

message("Running GSEA on Log2 data...")
gsea_hallmark <- fgseaMultilevel(pathways = hallmark_list, stats = ranking_valid, eps = 0)
gsea_kegg <- fgseaMultilevel(pathways = kegg_list, stats = ranking_valid, eps = 0)

# Sort and Save
gsea_hallmark <- as.data.table(gsea_hallmark)[order(padj)]
gsea_kegg <- as.data.table(gsea_kegg)[order(padj)]

fwrite(gsea_hallmark, "results/GSEA_hallmark_log2.tsv", sep = "\t")
fwrite(gsea_kegg, "results/GSEA_kegg_log2.tsv", sep = "\t")

# Check Parkinson's result in the corrected data
head(gsea_hallmark, 10)
head(gsea_kegg, 10)
gsea_kegg[gsea_kegg$pathway == "KEGG_PARKINSONS_DISEASE"]


# ==============================================================================
# PART 2: TARGETED PATHWAY ANALYSIS & MITOPHAGY DISSECTION
# ==============================================================================

# 1. VERIFICATION OF INPUTS ----------------------------------------------------

# We assume 'expr_by_gene' (expression matrix) and 'ranking_valid' (stats)
# exist from the previous pipeline.
stopifnot(exists("expr_by_gene") && is.matrix(expr_by_gene))
stopifnot(exists("ranking_valid") && is.numeric(ranking_valid))

# Check overlap between expression data and ranking vector
overlap <- sum(names(ranking_valid) %in% rownames(expr_by_gene))
message(paste("Genes matching between expression and ranking:", overlap))

# 2. DEFINING SPECIFIC PD SUB-PATHWAYS -----------------------------------------

# Instead of looking at "Parkinson's" as a whole, let's break it down into
# three distinct biological modules to see where the signal comes from.

# A. Mitochondrial Complex I (Energy production, often damaged by radiation)

PD_mito_complexI <- c(
  "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS7","NDUFS8",
  "NDUFV1","NDUFV2","NDUFA2","NDUFA8","NDUFA9",
  "NDUFB9","NDUFB10","NDUFB8"
)

# B. Mitophagy (PINK1/Parkin) - The quality control mechanism

PD_mitophagy_PINK1_PRKN <- c(
  "PINK1","PRKN","PARK7","UCHL1","SNCAIP","SNCA",
  "FBXO7","ATP13A2"
)

# C. Dopaminergic Synapse - The neuron type that dies in PD

PD_dopaminergic_synapse <- c(
  "TH","SLC6A3","DDC","DRD2","SLC18A2","SLC18A1"
)

# Combine into a list for FGSEA/GSVA
custom_pd <- list(
  PD_mito_complexI = PD_mito_complexI,
  PD_mitophagy_PINK1_PRKN = PD_mitophagy_PINK1_PRKN,
  PD_dop_synapse = PD_dopaminergic_synapse
)

# Check how many genes we actually have in our data
coverage <- lapply(custom_pd, function(gs) intersect(gs, names(ranking_valid)))
print(sapply(coverage, length))

# 3. RUNNING TARGETED GSEA (fgsea) ---------------------------------------------

library(fgsea)

# Running fgsea on these 3 specific lists
# We expect Mitophagy to show signal based on general radiation literature
pd_fgsea <- fgseaMultilevel(
  pathways = custom_pd,
  stats = sort(ranking_valid, decreasing = TRUE),
  eps = 0
)

print(pd_fgsea)

# 4. SAMPLE-LEVEL ENRICHMENT (GSVA) --------------------------------------------

# We want to see how these pathways behave per sample to correlate with Dose
library(GSVA)

gsva_param <- gsvaParam(
  expr = as.matrix(expr_by_gene),
  geneSets = custom_pd,
  kcdf = "Gaussian"
)

gsva_scores <- gsva(gsva_param)

# 5. LINEAR MODELING OF PATHWAY SCORES -----------------------------------------

# Does radiation dose linearly predict the activity of these pathways?
fit_gsva <- lmFit(gsva_scores, design) # Using the design from Part 1
fit_gsva <- eBayes(fit_gsva)

# Extract Dose results
res_gsva <- topTable(fit_gsva, coef = "DoseNum", number = Inf)
print(res_gsva)

# Quick visualization of the strongest hit (Mitophagy)
plot(targets$DoseNum, gsva_scores["PD_mitophagy_PINK1_PRKN", ],
     main = "Mitophagy Score vs Radiation Dose",
     xlab = "Dose (Gy)", ylab = "GSVA Score")
abline(lm(gsva_scores["PD_mitophagy_PINK1_PRKN", ] ~ targets$DoseNum), col="red")

# 6. INTERACTION ANALYSIS (DOSE x DONOR) ---------------------------------------

# Check if different donors respond differently to radiation in these pathways
design_int <- model.matrix(~ DoseNum * Donor, data = targets)
fit_int <- lmFit(gsva_scores, design_int)
fit_int <- eBayes(fit_int)

# Look at interaction terms only
res_int <- topTable(fit_int, coef = grep("DoseNum:Donor", colnames(design_int), value = TRUE))
print(res_int)

# 7. DEEP DIVE: MITOPHAGY SUB-MODULES ------------------------------------------

# Since Mitophagy was the hit, let's dissect it further into functional steps.

PD_mitophagy_core <- c("PINK1", "PRKN", "PARK7", "FBXO7", "ATP13A2")
PD_mitophagy_regulators <- c("USP30", "USP15")
PD_mitophagy_receptors_kinases <- c("OPTN", "CALCOCO2", "TBK1")
PD_mitophagy_fusion_fission <- c("MFN1", "MFN2", "DNM1L")
PD_mitophagy_autophagy_markers <- c("SQSTM1", "MAP1LC3B")

custom_pd_mitophagy <- list(
  PD_mitophagy_core = PD_mitophagy_core,
  PD_mitophagy_regulators = PD_mitophagy_regulators,
  PD_mitophagy_receptors_kinases = PD_mitophagy_receptors_kinases,
  PD_mitophagy_fusion_fission = PD_mitophagy_fusion_fission,
  PD_mitophagy_autophagy_markers = PD_mitophagy_autophagy_markers
)

# Run GSVA on these fine-grained modules
gsva_param_mitophagy <- gsvaParam(
  expr = as.matrix(expr_by_gene),
  geneSets = custom_pd_mitophagy,
  kcdf = "Gaussian"
)
gsva_scores_mitophagy <- gsva(gsva_param_mitophagy)

# Test dose response on these sub-modules
fit_mitophagy <- lmFit(gsva_scores_mitophagy, design)
fit_mitophagy <- eBayes(fit_mitophagy)
res_mitophagy <- topTable(fit_mitophagy, coef = "DoseNum", number = Inf)
print(res_mitophagy)

# 8. VISUALIZATION OF MITOPHAGY MODULES ----------------------------------------

library(tidyverse)

# Prepare data for plotting
gsva_long <- gsva_scores_mitophagy %>%
  as.data.frame() %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(cols = -Pathway, names_to = "SampleID", values_to = "GSVA_score") %>%
  left_join(targets %>% select(SampleID, DoseNum, Donor), by = "SampleID")

# Plot 1: Dose response per module
ggplot(gsva_long, aes(x = DoseNum, y = GSVA_score)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Pathway, scales = "free_y") +
  labs(title = "Radiation dose-dependent modulation of mitophagy modules",
       x = "Radiation dose (Gy)", y = "GSVA Score") +
  theme_bw()

# Plot 2: Effect sizes
effect_df <- res_mitophagy %>% as.data.frame() %>% rownames_to_column("Pathway")

ggplot(effect_df, aes(x = logFC, y = -log10(adj.P.Val), label = Pathway)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(title = "Significance of Mitophagy Sub-modules",
       x = "LogFC (Effect Size)", y = "-log10 Adj P-value") +
  theme_bw()

# 9. RISK SCORE & MEDIATION ANALYSIS -------------------------------------------

# Hypothesis: Does radiation increase PD risk *through* Mitophagy?
# To test this, we need a "Clean" PD Risk score that DOES NOT contain the
# mitophagy genes (otherwise the correlation is mathematical, not biological).

# A. Define Lists
pd_risk_genes <- c("SNCA", "LRRK2", "PINK1", "PRKN", "PARK7", "DJ1", "ATP13A2", 
                   "FBXO7", "UCHL1", "NDUFS1", "NDUFS2", "NDUFV1", "NDUFA8", 
                   "VDAC1", "VDAC2", "GBA", "CTSD", "LAMP1", "LAMP2", 
                   "TH", "SLC6A3", "SLC18A2", "DDC")

mitophagy_core_genes <- c("PINK1", "PRKN", "PARK7", "FBXO7", "ATP13A2", "SNCA", "UCHL1", "SNCAIP")

# B. Create the "Clean" list (Set Difference)
pd_genes_clean <- setdiff(pd_risk_genes, mitophagy_core_genes)
pd_genes_clean <- intersect(pd_genes_clean, rownames(expr_by_gene))

message(paste("Calculating Clean PD Risk Score using", length(pd_genes_clean), "genes."))

# C. Calculate Scores
# Risk Score (Dependent Variable)
pd_risk_score_clean <- colMeans(expr_by_gene[pd_genes_clean, ], na.rm = TRUE)
targets$PD_Risk_Score_Clean <- pd_risk_score_clean[targets$SampleID]

# Mitophagy Score (Mediator)
targets$Mitophagy_Core_GSVA <- gsva_scores_mitophagy["PD_mitophagy_core", ]

# D. Mediation Models

# Model A: Total Effect (Radiation -> Risk)
# Does radiation increase the general vulnerability of the cell?
message("--- Model A: Radiation -> PD Risk ---")
fit_a <- lm(PD_Risk_Score_Clean ~ DoseNum + Donor, data = targets)
print(summary(fit_a))

# Model B: Mediator Association (Mitophagy -> Risk)
# Is Mitophagy activity associated with the risk profile (independent of dose)?
message("--- Model B: Mitophagy -> PD Risk ---")
fit_b <- lm(PD_Risk_Score_Clean ~ Mitophagy_Core_GSVA + Donor, data = targets)
print(summary(fit_b))

# Model C: Mediation Test (Radiation + Mitophagy -> Risk)
# If we control for Mitophagy, does the Radiation effect disappear?
# If "DoseNum" becomes less significant but "Mitophagy" stays significant, 
# it suggests Mitophagy mediates the risk.
message("--- Model C: Mediation (Rad + Mito -> Risk) ---")
fit_c <- lm(PD_Risk_Score_Clean ~ DoseNum + Mitophagy_Core_GSVA + Donor, data = targets)
print(summary(fit_c))

# 10. FINAL PLOT: THE MEDIATION LINK -------------------------------------------

ggplot(targets, aes(x = Mitophagy_Core_GSVA, y = PD_Risk_Score_Clean)) +
  geom_point(aes(color = as.factor(DoseNum)), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  theme_bw() +
  scale_color_viridis_d(name = "Dose (Gy)") +
  labs(
    title = "Mitophagy Core Predicts Independent PD Molecular Risk",
    subtitle = "Analysis using 'Clean' Score (non-overlapping genes)",
    x = "Mitophagy Core Activity (GSVA)",
    y = "PD Risk Score (Non-Mitophagy Genes)"
  )

