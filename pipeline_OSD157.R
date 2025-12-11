# ================================================================
# PIPELINE OSD-157 
# ================================================================
# Autor: Arthur de Melo Pelozo 
# ================================================================

###########################################
# 0. Load packages
###########################################

library(limma)
library(data.table)
library(fgsea)
library(tidyverse)

###########################################
# 1. Import sample metadata and raw files
###########################################

targets <- read.targets("data", ext = "txt", sep = "\t")
RG <- read.maimages(targets$FileName, source = "agilent", green.only = TRUE)

###########################################
# 2. Preprocessing of Agilent microarrays
###########################################

# Background correction + normalization
RG_bg   <- backgroundCorrect(RG, method = "normexp", offset = 50)
RG_norm <- normalizeBetweenArrays(RG_bg, method = "quantile")

expr_matrix <- RG_norm$E
genes_info  <- RG_norm$genes

###########################################
# 3. Remove control probes
###########################################

probe2gene <- genes_info$GeneName

valid <- !is.na(probe2gene) &
  probe2gene != "" &
  probe2gene != "DarkCorner" &
  probe2gene != "GE_BrightCorner"

expr_valid <- expr_matrix[valid, ]
gene_valid <- probe2gene[valid]

###########################################
# 4. Collapse probes → gene symbols
###########################################

expr_by_gene <- limma::avereps(expr_valid, ID = gene_valid)

###########################################
# 5. Build design matrix (dose + donor as batch)
###########################################

# Extract numeric dose (e.g., 0Gy, 0.5Gy → 0, 0.5)
dose_num <- as.numeric(gsub("Gy", "", targets$Dose))

donor <- factor(targets$Donor)

design <- model.matrix(~ dose_num + donor)
colnames(design) <- c(
  "Intercept",
  "DoseNum",
  paste0("Donor", levels(donor)[-1])
)

###########################################
# 6. Differential expression using limma
###########################################

fit <- lmFit(expr_by_gene, design)
fit <- eBayes(fit)

###########################################
# 7. Generate gene ranking for GSEA
###########################################

gene_ranking <- fit$t[, "DoseNum"]
names(gene_ranking) <- rownames(fit$t)

gene_ranking <- gene_ranking[
  !is.na(gene_ranking) & is.finite(gene_ranking)
]

gene_ranking <- sort(gene_ranking, decreasing = TRUE)

###########################################
# 8. Load and prepare MSigDB gene sets
###########################################

gmt_raw <- fread("data/msigdb_human_gene_sets.gmt",
                 sep = "\t", header = FALSE, fill = TRUE)

# Convert GMT to long format: gs_name | gene_symbol
gmt_long <- gmt_raw %>%
  pivot_longer(cols = -c(V1, V2), names_to = "n", values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  select(gs_name = V1, gene_symbol = gene)

# Build Hallmark list
hallmark_list <- gmt_long %>%
  filter(grepl("^HALLMARK_", gs_name)) %>%
  group_split(gs_name) %>%
  setNames(unique(gmt_long$gs_name[grepl("^HALLMARK_", gmt_long$gs_name)]))

# Build KEGG list
kegg_list <- gmt_long %>%
  filter(grepl("^KEGG_", gs_name)) %>%
  group_split(gs_name) %>%
  setNames(unique(gmt_long$gs_name[grepl("^KEGG_", gmt_long$gs_name)]))

# Parkinson gene set only
pd_kegg <- kegg_list["KEGG_PARKINSONS_DISEASE"]

###########################################
# 9. GSEA (fgseaMultilevel)
###########################################

gsea_hallmark <- fgseaMultilevel(
  pathways = hallmark_list,
  stats = gene_ranking,
  eps = 0
)

gsea_kegg <- fgseaMultilevel(
  pathways = kegg_list,
  stats = gene_ranking,
  eps = 0
)

gsea_pd <- fgseaMultilevel(
  pathways = pd_kegg,
  stats = gene_ranking,
  eps = 0
)

###########################################
# 10. Save results
###########################################

gsea_hallmark <- as.data.table(gsea_hallmark)[order(padj)]
gsea_kegg     <- as.data.table(gsea_kegg)[order(padj)]
gsea_pd       <- as.data.table(gsea_pd)

dir.create("results", showWarnings = FALSE)

fwrite(gsea_hallmark, "results/GSEA_Hallmark.tsv", sep = "\t")
fwrite(gsea_kegg,     "results/GSEA_KEGG.tsv", sep = "\t")
fwrite(gsea_pd,       "results/GSEA_KEGG_PD.tsv", sep = "\t")
