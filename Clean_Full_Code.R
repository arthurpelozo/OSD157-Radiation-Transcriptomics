# ==============================================================================
# PROJECT: Radiation Analysis (GLDS-157 / Agilent Microarray)
# AUTHOR: Arthur Pelozo
# STATUS: CORRECTED (Log2 Transformation + Direct Plotting)
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
RG.bc <- backgroundCorrect(RG, method = "normexp", offset = 50)

# B. *** THE CRITICAL FIX: LOG2 TRANSFORMATION ***
message("Applying Log2 transformation...")
E_log <- log2(RG.bc$E)

# C. Quantile Normalization
expr_matrix <- normalizeBetweenArrays(E_log, method = "quantile")
colnames(expr_matrix) <- targets$SampleID

# D. Verification Plot (Density)
# Direct plot to device
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

p_pca <- ggplot(df_pca, aes(PC1, PC2, color = Dose, shape = Donor)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA — Log2 Normalized Data")

print(p_pca) # Generating the plot

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

# SANITY CHECK
print(head(res_trend))

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

# Sort
gsea_hallmark <- as.data.table(gsea_hallmark)[order(padj)]
gsea_kegg <- as.data.table(gsea_kegg)[order(padj)]

# Check Parkinson's result
print(gsea_kegg[gsea_kegg$pathway == "KEGG_PARKINSONS_DISEASE"])


# ==============================================================================
# PART 2: TARGETED PATHWAY ANALYSIS & MITOPHAGY DISSECTION
# ==============================================================================

# 1. VERIFICATION OF INPUTS ----------------------------------------------------

stopifnot(exists("expr_by_gene") && is.matrix(expr_by_gene))
stopifnot(exists("ranking_valid") && is.numeric(ranking_valid))

# 2. DEFINING SPECIFIC PD SUB-PATHWAYS -----------------------------------------

PD_mito_complexI <- c(
  "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS7","NDUFS8",
  "NDUFV1","NDUFV2","NDUFA2","NDUFA8","NDUFA9",
  "NDUFB9","NDUFB10","NDUFB8"
)

PD_mitophagy_PINK1_PRKN <- c(
  "PINK1","PRKN","PARK7","UCHL1","SNCAIP","SNCA",
  "FBXO7","ATP13A2"
)

PD_dopaminergic_synapse <- c(
  "TH","SLC6A3","DDC","DRD2","SLC18A2","SLC18A1"
)

custom_pd <- list(
  PD_mito_complexI = PD_mito_complexI,
  PD_mitophagy_PINK1_PRKN = PD_mitophagy_PINK1_PRKN,
  PD_dop_synapse = PD_dopaminergic_synapse
)

# 3. RUNNING TARGETED GSEA (fgsea) ---------------------------------------------

library(fgsea)

pd_fgsea <- fgseaMultilevel(
  pathways = custom_pd,
  stats = sort(ranking_valid, decreasing = TRUE),
  eps = 0
)

print(pd_fgsea)

# 4. SAMPLE-LEVEL ENRICHMENT (GSVA) --------------------------------------------

library(GSVA)

gsva_param <- gsvaParam(
  expr = as.matrix(expr_by_gene),
  geneSets = custom_pd,
  kcdf = "Gaussian"
)

gsva_scores <- gsva(gsva_param)

# 5. LINEAR MODELING OF PATHWAY SCORES -----------------------------------------

fit_gsva <- lmFit(gsva_scores, design)
fit_gsva <- eBayes(fit_gsva)

res_gsva <- topTable(fit_gsva, coef = "DoseNum", number = Inf)
print(res_gsva)

# Quick visualization of the strongest hit (Mitophagy) - Direct Plot
plot(targets$DoseNum, gsva_scores["PD_mitophagy_PINK1_PRKN", ],
     main = "Mitophagy Score vs Radiation Dose",
     xlab = "Dose (Gy)", ylab = "GSVA Score", pch=19, col="blue")
abline(lm(gsva_scores["PD_mitophagy_PINK1_PRKN", ] ~ targets$DoseNum), col="red", lwd=2)

# 6. INTERACTION ANALYSIS (DOSE x DONOR) ---------------------------------------

design_int <- model.matrix(~ DoseNum * Donor, data = targets)
fit_int <- lmFit(gsva_scores, design_int)
fit_int <- eBayes(fit_int)
res_int <- topTable(fit_int, coef = grep("DoseNum:Donor", colnames(design_int), value = TRUE))
print(res_int)

# 7. DEEP DIVE: MITOPHAGY SUB-MODULES ------------------------------------------

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

gsva_param_mitophagy <- gsvaParam(
  expr = as.matrix(expr_by_gene),
  geneSets = custom_pd_mitophagy,
  kcdf = "Gaussian"
)
gsva_scores_mitophagy <- gsva(gsva_param_mitophagy)

fit_mitophagy <- lmFit(gsva_scores_mitophagy, design)
fit_mitophagy <- eBayes(fit_mitophagy)
res_mitophagy <- topTable(fit_mitophagy, coef = "DoseNum", number = Inf)
print(res_mitophagy)

# 8. VISUALIZATION OF MITOPHAGY MODULES ----------------------------------------

library(tidyverse)

gsva_long <- gsva_scores_mitophagy %>%
  as.data.frame() %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(cols = -Pathway, names_to = "SampleID", values_to = "GSVA_score") %>%
  left_join(targets %>% select(SampleID, DoseNum, Donor), by = "SampleID")

# Plot 1: Dose response per module
p_dose_mod <- ggplot(gsva_long, aes(x = DoseNum, y = GSVA_score)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Pathway, scales = "free_y") +
  labs(title = "Radiation dose-dependent modulation of mitophagy modules",
       x = "Radiation dose (Gy)", y = "GSVA Score") +
  theme_bw()

print(p_dose_mod) # Generating the plot

# Plot 2: Effect sizes
effect_df <- res_mitophagy %>% as.data.frame() %>% rownames_to_column("Pathway")

p_volcano_mod <- ggplot(effect_df, aes(x = logFC, y = -log10(adj.P.Val), label = Pathway)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(title = "Significance of Mitophagy Sub-modules",
       x = "LogFC (Effect Size)", y = "-log10 Adj P-value") +
  theme_bw()

print(p_volcano_mod) # Generating the plot

# 9. RISK SCORE & MEDIATION ANALYSIS -------------------------------------------

pd_risk_genes <- c("SNCA", "LRRK2", "PINK1", "PRKN", "PARK7", "DJ1", "ATP13A2", 
                   "FBXO7", "UCHL1", "NDUFS1", "NDUFS2", "NDUFV1", "NDUFA8", 
                   "VDAC1", "VDAC2", "GBA", "CTSD", "LAMP1", "LAMP2", 
                   "TH", "SLC6A3", "SLC18A2", "DDC")

mitophagy_core_genes <- c("PINK1", "PRKN", "PARK7", "FBXO7", "ATP13A2", "SNCA", "UCHL1", "SNCAIP")

pd_genes_clean <- setdiff(pd_risk_genes, mitophagy_core_genes)
pd_genes_clean <- intersect(pd_genes_clean, rownames(expr_by_gene))

# Calculate Scores
pd_risk_score_clean <- colMeans(expr_by_gene[pd_genes_clean, ], na.rm = TRUE)
targets$PD_Risk_Score_Clean <- pd_risk_score_clean[targets$SampleID]
targets$Mitophagy_Core_GSVA <- gsva_scores_mitophagy["PD_mitophagy_core", ]

# Mediation Models
message("--- Model A: Radiation -> PD Risk ---")
fit_a <- lm(PD_Risk_Score_Clean ~ DoseNum + Donor, data = targets)
print(summary(fit_a))

message("--- Model B: Mitophagy -> PD Risk ---")
fit_b <- lm(PD_Risk_Score_Clean ~ Mitophagy_Core_GSVA + Donor, data = targets)
print(summary(fit_b))

message("--- Model C: Mediation (Rad + Mito -> Risk) ---")
fit_c <- lm(PD_Risk_Score_Clean ~ DoseNum + Mitophagy_Core_GSVA + Donor, data = targets)
print(summary(fit_c))

# 10. FINAL PLOT: THE MEDIATION LINK -------------------------------------------

p_mediation <- ggplot(targets, aes(x = Mitophagy_Core_GSVA, y = PD_Risk_Score_Clean)) +
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

print(p_mediation) # Generating the plot

# ==============================================================================
# PART 3: STRATEGIC VISUALIZATIONS (Direct Plotting)
# ==============================================================================

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# FIGURE 1: THE PUSH-PULL HEATMAP
# ------------------------------------------------------------------------------

# A. Select key genes
genes_push_pull <- c("SQSTM1", "MAP1LC3B",  # Execution
                     "PINK1", "PRKN", "PARK7", "FBXO7", # Sensors
                     "USP30", "USP15")      # Regulators

# B. Prepare Data
valid_genes <- intersect(genes_push_pull, rownames(expr_by_gene))
expr_subset <- expr_by_gene[valid_genes, ]

# C. Order columns by Radiation Dose
ordem_amostras <- order(targets$DoseNum)
expr_subset <- expr_subset[, ordem_amostras]

# D. Create Annotations
annot_col <- data.frame(Dose = paste0(targets$DoseNum[ordem_amostras], "Gy"))
rownames(annot_col) <- colnames(expr_subset)

# E. Generate Heatmap (Directly to device/RMarkdown)
pheatmap(expr_subset,
         scale = "row",
         annotation_col = annot_col,
         cluster_cols = FALSE,        
         cluster_rows = TRUE,         
         show_colnames = FALSE,
         main = "Push-Pull Mechanism: Regulators (Blue) vs Sensors (Red)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# ------------------------------------------------------------------------------
# FIGURE 2: SPECIFICITY BARPLOT
# ------------------------------------------------------------------------------

# A. Create Summary Dataframe
df_spec <- data.frame(
  Pathway = c("Mitophagy (Cleaning)", "Synapse (Neuron Health)", "Complex I (Energy)"),
  NES = c(1.72, 0.92, -0.94), 
  Significance = c("Significant (p<0.05)", "Not Sig.", "Not Sig.")
)

# B. Plot
p_spec <- ggplot(df_spec, aes(x = reorder(Pathway, NES), y = NES, fill = Significance)) +
  geom_col(width = 0.6, color="black") +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(1.5, -1.5), linetype="dashed", color="gray") + 
  coord_flip() +
  scale_fill_manual(values = c("gray80", "firebrick")) +
  theme_bw() +
  labs(title = "Specificity of Radiation Response",
       y = "Normalized Enrichment Score (NES)", x = "") +
  theme(legend.position = "bottom")

print(p_spec) # Generating the plot

# ------------------------------------------------------------------------------
# FIGURE 3: CLINICAL RISK ESCALATION
# ------------------------------------------------------------------------------

p_risk <- ggplot(targets, aes(x = factor(DoseNum), y = PD_Risk_Score_Clean)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, aes(color = Donor)) +
  geom_smooth(method = "lm", aes(group=1), se=FALSE, color="red", linetype="dashed") +
  theme_classic() +
  labs(title = "The Transcriptomic Scar",
       subtitle = "Linear escalation of genetic vulnerability",
       x = "Radiation Dose (Gy)",
       y = "PD Risk Score") 

print(p_risk) # Generating the plot
