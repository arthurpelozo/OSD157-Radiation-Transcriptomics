Complemetary_Analysis_PD_Pathways_and_Mitophagy.R
================
arthu
2025-12-13

``` r
################################################################
#PART2 - Specific PD Pathway
################################################################
```

``` r
# 1. Ensure main expression object exists
stopifnot(exists("expr_by_gene"))
stopifnot(is.matrix(expr_by_gene))

cat("expr_by_gene dimensions:\n")
```

    ## expr_by_gene dimensions:

``` r
print(dim(expr_by_gene))
```

    ## [1] 30587    25

``` r
# 2. Ensure ranking exists
stopifnot(exists("ranking_valid"))
stopifnot(is.numeric(ranking_valid))

cat("ranking_valid length:\n")
```

    ## ranking_valid length:

``` r
print(length(ranking_valid))
```

    ## [1] 14054

``` r
# 3. Ensure rownames(gene expression) match names(ranking)
cat("Check overlap between expr and ranking names:\n")
```

    ## Check overlap between expr and ranking names:

``` r
overlap <- sum(names(ranking_valid) %in% rownames(expr_by_gene))
print(overlap)
```

    ## [1] 14054

``` r
# Show a few gene names
cat("First 20 gene names from expression:\n")
```

    ## First 20 gene names from expression:

``` r
print(head(rownames(expr_by_gene), 20))
```

    ##  [1] "APOBEC3B"     "AA085955"     "ATP11B"       "AK092846"     "DNAJA1"      
    ##  [6] "THC2741789"   "EHMT2"        "RPL23"        "T12590"       "A_24_P704878"
    ## [11] "RPS13"        "AK021474"     "HDDC3"        "PRNP"         "AK091028"    
    ## [16] "LOC150759"    "A_24_P799245" "KIAA0101"     "AK026647"     "MEGF11"

``` r
cat("First 20 gene names from ranking:\n")
```

    ## First 20 gene names from ranking:

``` r
print(head(names(ranking_valid), 20))
```

    ##  [1] "VWCE"    "PLK2"    "NUDT8"   "CLEC4D"  "MDM2"    "MYO1A"   "APBB3"   "PCNA"   
    ##  [9] "SESN2"   "REV3L"   "PHLDA3"  "CD70"    "TFRC"    "CDKN1A"  "TMEM30A" "FDXR"   
    ## [17] "SLC29A1" "LILRA3"  "ACTA2"   "TGM2"

``` r
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

cat("Custom PD pathways created.\nNames:\n")
```

    ## Custom PD pathways created.
    ## Names:

``` r
print(names(custom_pd))
```

    ## [1] "PD_mito_complexI"        "PD_mitophagy_PINK1_PRKN" "PD_dop_synapse"

``` r
coverage <- lapply(custom_pd, function(gs){
  intersect(gs, names(ranking_valid))
})

coverage_sizes <- sapply(coverage, length)

cat("Coverage counts:\n")
```

    ## Coverage counts:

``` r
print(coverage_sizes)
```

    ##        PD_mito_complexI PD_mitophagy_PINK1_PRKN          PD_dop_synapse 
    ##                      14                       7                       6

``` r
cat("\nGenes present in dataset (per pathway):\n")
```

    ## 
    ## Genes present in dataset (per pathway):

``` r
print(coverage)
```

    ## $PD_mito_complexI
    ##  [1] "NDUFS1"  "NDUFS2"  "NDUFS3"  "NDUFS4"  "NDUFS7"  "NDUFS8"  "NDUFV1"  "NDUFV2" 
    ##  [9] "NDUFA2"  "NDUFA8"  "NDUFA9"  "NDUFB9"  "NDUFB10" "NDUFB8" 
    ## 
    ## $PD_mitophagy_PINK1_PRKN
    ## [1] "PINK1"   "PARK7"   "UCHL1"   "SNCAIP"  "SNCA"    "FBXO7"   "ATP13A2"
    ## 
    ## $PD_dop_synapse
    ## [1] "TH"      "SLC6A3"  "DDC"     "DRD2"    "SLC18A2" "SLC18A1"

``` r
library(fgsea)

# Sort stats decreasing
rank_stats <- sort(ranking_valid, decreasing = TRUE)

pd_fgsea <- fgseaMultilevel(
  pathways = custom_pd,
  stats = rank_stats,
  eps = 0
)
```

    ## Warning in serialize(data, node$con): 'package:stats' pode não estar disponível ao
    ## carregar

``` r
pd_fgsea
```

    ##                    pathway       pval       padj    log2err         ES        NES
    ##                     <char>      <num>      <num>      <num>      <num>      <num>
    ## 1:          PD_dop_synapse 0.55932203 0.57986111 0.07627972  0.4241150  0.9280476
    ## 2:        PD_mito_complexI 0.57986111 0.57986111 0.06450312 -0.3261686 -0.8916441
    ## 3: PD_mitophagy_PINK1_PRKN 0.02377973 0.07133918 0.35248786  0.7092467  1.6376964
    ##     size  leadingEdge
    ##    <int>       <list>
    ## 1:     6 DRD2, TH....
    ## 2:    14 NDUFB9, ....
    ## 3:     7 SNCA, SN....

``` r
library(GSVA)

gsva_param <- gsvaParam(
  expr = as.matrix(expr_by_gene),
  geneSets = custom_pd,
  kcdf = "Gaussian"
)
gsva_scores <- gsva(gsva_param)
```

    ## ℹ GSVA version 2.4.1

    ## ℹ Searching for genes/features with constant values

    ## ℹ Calculating GSVA ranks

    ## ℹ GSVA dense (classical) algorithm

    ## ℹ Row-wise ECDF estimation with Gaussian kernels

    ## ℹ Calculating row ECDFs

    ## ℹ Calculating column ranks

    ## ℹ GSVA dense (classical) algorithm

    ## ℹ Calculating GSVA scores

    ## ✔ Calculations finished

``` r
class(gsva_param)
```

    ## [1] "gsvaParam"
    ## attr(,"package")
    ## [1] "GSVA"

``` r
dim(gsva_scores)
```

    ## [1]  3 25

``` r
gsva_scores
```

    ##                         GSM1080614_a_0.5Gy GSM1080615_a_0Gy GSM1080616_a_2Gy
    ## PD_mito_complexI                0.03251121       -0.1164269       -0.4836440
    ## PD_mitophagy_PINK1_PRKN         0.28858091       -0.1352450        0.3425630
    ## PD_dop_synapse                  0.20445938       -0.4443570        0.1966127
    ##                         GSM1080617_a_5Gy GSM1080618_a_8Gy GSM1080619_b_0.5Gy
    ## PD_mito_complexI              -0.5777041       -0.1421918         -0.5909091
    ## PD_mitophagy_PINK1_PRKN        0.6285379        0.4278432         -0.1230595
    ## PD_dop_synapse                -0.1428807       -0.2868121          0.1664806
    ##                         GSM1080620_b_0Gy GSM1080621_b_2Gy GSM1080622_b_5Gy
    ## PD_mito_complexI             -0.36397760      -0.55879602      -0.71333462
    ## PD_mitophagy_PINK1_PRKN      -0.18573095       0.06401617       0.07358654
    ## PD_dop_synapse                0.02061929      -0.11380911       0.08748867
    ##                         GSM1080623_b_8Gy GSM1080649_X_0.5Gy GSM1080650_X_0Gy
    ## PD_mito_complexI              -0.6450336          0.6453185       0.42153880
    ## PD_mitophagy_PINK1_PRKN        0.3169056          0.2590538      -0.47084047
    ## PD_dop_synapse                -0.1552925         -0.5551653      -0.02009372
    ##                         GSM1080651_X_2Gy GSM1080652_X_5Gy GSM1080653_X_8Gy
    ## PD_mito_complexI               0.4800916        0.3454355       0.59695940
    ## PD_mitophagy_PINK1_PRKN        0.2276651        0.3645095       0.56867336
    ## PD_dop_synapse                 0.3662591        0.3745481      -0.02425882
    ##                         GSM1080654_Y_0.5Gy GSM1080655_Y_0Gy GSM1080656_Y_2Gy
    ## PD_mito_complexI                 0.1436003       0.21242404       -0.1920508
    ## PD_mitophagy_PINK1_PRKN         -0.2664986      -0.51369708       -0.4313660
    ## PD_dop_synapse                  -0.2325824      -0.01837929       -0.1486883
    ##                         GSM1080657_Y_5Gy GSM1080658_Y_8Gy GSM1080659_Z_0.5Gy
    ## PD_mito_complexI              0.18885945      0.206383849         0.25163649
    ## PD_mitophagy_PINK1_PRKN      -0.37715012      0.157212212        -0.35312045
    ## PD_dop_synapse                0.02914179     -0.005872151         0.08868252
    ##                         GSM1080660_Z_0Gy GSM1080661_Z_2Gy GSM1080662_Z_5Gy
    ## PD_mito_complexI               0.5739719       0.03334292      -0.07545874
    ## PD_mitophagy_PINK1_PRKN       -0.1090582      -0.33953022      -0.28717451
    ## PD_dop_synapse                -0.3805383       0.16775059       0.10323711
    ##                         GSM1080663_Z_8Gy
    ## PD_mito_complexI              0.07742912
    ## PD_mitophagy_PINK1_PRKN       0.01260736
    ## PD_dop_synapse                0.59972994
    ## attr(,"geneSets")
    ## attr(,"geneSets")$PD_mito_complexI
    ##  [1] "NDUFS1"  "NDUFS2"  "NDUFS3"  "NDUFS4"  "NDUFS7"  "NDUFS8"  "NDUFV1"  "NDUFV2" 
    ##  [9] "NDUFA2"  "NDUFA8"  "NDUFA9"  "NDUFB9"  "NDUFB10" "NDUFB8" 
    ## 
    ## attr(,"geneSets")$PD_mitophagy_PINK1_PRKN
    ## [1] "PINK1"   "PARK7"   "UCHL1"   "SNCAIP"  "SNCA"    "FBXO7"   "ATP13A2"
    ## 
    ## attr(,"geneSets")$PD_dop_synapse
    ## [1] "TH"      "SLC6A3"  "DDC"     "DRD2"    "SLC18A2" "SLC18A1"

``` r
# 1. Fit Linear Model on GSVA Scores
fit_gsva <- lmFit(gsva_scores, design_gsva)
fit_gsva <- eBayes(fit_gsva)

res_gsva <- topTable(
  fit_gsva,
  coef = "DoseNum",
  number = Inf
)
print(res_gsva)
```

    ##                               logFC      AveExpr         t      P.Value    adj.P.Val
    ## PD_mitophagy_PINK1_PRKN  0.05785343  0.005571355  4.101675 0.0001317979 0.0003953938
    ## PD_dop_synapse           0.01961659 -0.004948791  1.230200 0.2236753015 0.2334725426
    ## PD_mito_complexI        -0.01698494 -0.010000975 -1.204243 0.2334725426 0.2334725426
    ##                                  B
    ## PD_mitophagy_PINK1_PRKN  0.7953288
    ## PD_dop_synapse          -5.7814432
    ## PD_mito_complexI        -5.8121015

``` r
# 2. Safety Check: Ensure dimensions match
stopifnot(
  ncol(gsva_scores) == nrow(design_gsva),
  all(colnames(gsva_scores) == rownames(design_gsva))
)

# 3. Visualization: Plot specific pathway vs Dose
plot(targets$DoseNum, gsva_scores["PD_mitophagy_PINK1_PRKN", ])
abline(lm(gsva_scores["PD_mitophagy_PINK1_PRKN", ] ~ targets$DoseNum))
```

![](Complemetary_Analysis_PD_Pathways_and_Mitophagy_files/figure-gfm/BLOCK%206:%20Final%20Verifications-1.png)<!-- -->

``` r
# 4. Interaction Analysis (Dose * Donor)
design_int <- model.matrix(~ DoseNum * Donor, data = targets)
print(design_int)
```

    ##    (Intercept) DoseNum Donorb DonorX DonorY DonorZ DoseNum:Donorb DoseNum:DonorX
    ## 28           1     0.5      0      0      0      0            0.0            0.0
    ## 29           1     0.0      0      0      0      0            0.0            0.0
    ## 30           1     2.0      0      0      0      0            0.0            0.0
    ## 31           1     5.0      0      0      0      0            0.0            0.0
    ## 32           1     8.0      0      0      0      0            0.0            0.0
    ## 33           1     0.5      1      0      0      0            0.5            0.0
    ## 34           1     0.0      1      0      0      0            0.0            0.0
    ## 35           1     2.0      1      0      0      0            2.0            0.0
    ## 36           1     5.0      1      0      0      0            5.0            0.0
    ## 37           1     8.0      1      0      0      0            8.0            0.0
    ## 63           1     0.5      0      1      0      0            0.0            0.5
    ## 64           1     0.0      0      1      0      0            0.0            0.0
    ## 65           1     2.0      0      1      0      0            0.0            2.0
    ## 66           1     5.0      0      1      0      0            0.0            5.0
    ## 67           1     8.0      0      1      0      0            0.0            8.0
    ## 68           1     0.5      0      0      1      0            0.0            0.0
    ## 69           1     0.0      0      0      1      0            0.0            0.0
    ## 70           1     2.0      0      0      1      0            0.0            0.0
    ## 71           1     5.0      0      0      1      0            0.0            0.0
    ## 72           1     8.0      0      0      1      0            0.0            0.0
    ## 73           1     0.5      0      0      0      1            0.0            0.0
    ## 74           1     0.0      0      0      0      1            0.0            0.0
    ## 75           1     2.0      0      0      0      1            0.0            0.0
    ## 76           1     5.0      0      0      0      1            0.0            0.0
    ## 77           1     8.0      0      0      0      1            0.0            0.0
    ##    DoseNum:DonorY DoseNum:DonorZ
    ## 28            0.0            0.0
    ## 29            0.0            0.0
    ## 30            0.0            0.0
    ## 31            0.0            0.0
    ## 32            0.0            0.0
    ## 33            0.0            0.0
    ## 34            0.0            0.0
    ## 35            0.0            0.0
    ## 36            0.0            0.0
    ## 37            0.0            0.0
    ## 63            0.0            0.0
    ## 64            0.0            0.0
    ## 65            0.0            0.0
    ## 66            0.0            0.0
    ## 67            0.0            0.0
    ## 68            0.5            0.0
    ## 69            0.0            0.0
    ## 70            2.0            0.0
    ## 71            5.0            0.0
    ## 72            8.0            0.0
    ## 73            0.0            0.5
    ## 74            0.0            0.0
    ## 75            0.0            2.0
    ## 76            0.0            5.0
    ## 77            0.0            8.0
    ## attr(,"assign")
    ##  [1] 0 1 2 2 2 2 3 3 3 3
    ## attr(,"contrasts")
    ## attr(,"contrasts")$Donor
    ## [1] "contr.treatment"

``` r
fit_int <- lmFit(gsva_scores, design_int)
fit_int <- eBayes(fit_int)

# Extract results specifically for the interaction terms
topTable(
  fit_int, 
  coef = grep("DoseNum:Donor", colnames(design_int), value = TRUE)
)
```

    ##                         DoseNum.Donorb DoseNum.DonorX DoseNum.DonorY DoseNum.DonorZ
    ## PD_dop_synapse             0.001790487     0.06068192    0.042014541     0.10943828
    ## PD_mito_complexI          -0.005545480     0.02333416    0.033970563    -0.02761591
    ## PD_mitophagy_PINK1_PRKN   -0.001434265     0.03167959    0.005106875    -0.03179120
    ##                              AveExpr         F   P.Value adj.P.Val
    ## PD_dop_synapse          -0.004948791 1.8420810 0.1176578 0.3529734
    ## PD_mito_complexI        -0.010000975 0.5287809 0.7145944 0.7694718
    ## PD_mitophagy_PINK1_PRKN  0.005571355 0.4541119 0.7694718 0.7694718

``` r
################################################################
#PART2 - Mithopaghy
################################################################
```

``` r
PD_mitophagy_core <- c(
  "PINK1", "PRKN", "PARK7", "FBXO7", "ATP13A2"
)

PD_mitophagy_regulators <- c(
  "USP30", "USP15"
)

PD_mitophagy_receptors_kinases <- c(
  "OPTN", "CALCOCO2", "TBK1"
)

PD_mitophagy_fusion_fission <- c(
  "MFN1", "MFN2", "DNM1L"
)

PD_mitophagy_autophagy_markers <- c(
  "SQSTM1", "MAP1LC3B"
)
```

``` r
custom_pd_mitophagy <- list(
  PD_mitophagy_core = PD_mitophagy_core,
  PD_mitophagy_regulators = PD_mitophagy_regulators,
  PD_mitophagy_receptors_kinases = PD_mitophagy_receptors_kinases,
  PD_mitophagy_fusion_fission = PD_mitophagy_fusion_fission,
  PD_mitophagy_autophagy_markers = PD_mitophagy_autophagy_markers
)

# Check gene coverage
coverage <- lapply(custom_pd_mitophagy, function(gs) {
  intersect(gs, rownames(expr_by_gene))
})

coverage_sizes <- sapply(coverage, length)

print(coverage_sizes)
```

    ##              PD_mitophagy_core        PD_mitophagy_regulators 
    ##                              4                              2 
    ## PD_mitophagy_receptors_kinases    PD_mitophagy_fusion_fission 
    ##                              3                              3 
    ## PD_mitophagy_autophagy_markers 
    ##                              2

``` r
print(coverage)
```

    ## $PD_mitophagy_core
    ## [1] "PINK1"   "PARK7"   "FBXO7"   "ATP13A2"
    ## 
    ## $PD_mitophagy_regulators
    ## [1] "USP30" "USP15"
    ## 
    ## $PD_mitophagy_receptors_kinases
    ## [1] "OPTN"     "CALCOCO2" "TBK1"    
    ## 
    ## $PD_mitophagy_fusion_fission
    ## [1] "MFN1"  "MFN2"  "DNM1L"
    ## 
    ## $PD_mitophagy_autophagy_markers
    ## [1] "SQSTM1"   "MAP1LC3B"

``` r
#GSVA

library(GSVA)

gsva_param_mitophagy <- gsvaParam(
  expr = as.matrix(expr_by_gene),
  geneSets = custom_pd_mitophagy,
  kcdf = "Gaussian"
)

gsva_scores_mitophagy <- gsva(gsva_param_mitophagy)
```

    ## ℹ GSVA version 2.4.1

    ## ℹ Searching for genes/features with constant values

    ## ℹ Calculating GSVA ranks

    ## ℹ GSVA dense (classical) algorithm

    ## ℹ Row-wise ECDF estimation with Gaussian kernels

    ## ℹ Calculating row ECDFs

    ## ℹ Calculating column ranks

    ## ℹ GSVA dense (classical) algorithm

    ## ℹ Calculating GSVA scores

    ## ✔ Calculations finished

``` r
dim(gsva_scores_mitophagy)
```

    ## [1]  5 25

``` r
#Dose Response
# Ensure sample alignment, 
stopifnot(all(colnames(gsva_scores_mitophagy) == rownames(design)))

fit_mitophagy <- lmFit(gsva_scores_mitophagy, design)
fit_mitophagy <- eBayes(fit_mitophagy)

res_mitophagy <- topTable(
  fit_mitophagy,
  coef = "DoseNum",
  number = Inf
)

print(res_mitophagy)
```

    ##                                       logFC       AveExpr         t      P.Value
    ## PD_mitophagy_autophagy_markers  0.065717448 -2.372134e-02  4.221867 0.0002527358
    ## PD_mitophagy_regulators        -0.081747795  6.523053e-02 -3.783917 0.0007986324
    ## PD_mitophagy_core               0.032954213 -1.239988e-02  2.602299 0.0149609191
    ## PD_mitophagy_receptors_kinases -0.030316623 -1.229879e-02 -2.486073 0.0195211284
    ## PD_mitophagy_fusion_fission     0.009059495 -4.653443e-06  0.515604 0.6104003933
    ##                                  adj.P.Val          B
    ## PD_mitophagy_autophagy_markers 0.001263679  0.3633005
    ## PD_mitophagy_regulators        0.001996581 -0.7315523
    ## PD_mitophagy_core              0.024401411 -3.4593229
    ## PD_mitophagy_receptors_kinases 0.024401411 -3.6993889
    ## PD_mitophagy_fusion_fission    0.610400393 -6.3788174

``` r
#Dose x Donor
design_int <- model.matrix(~ DoseNum * Donor, data = targets)
rownames(design_int) <- targets$SampleID

stopifnot(all(colnames(gsva_scores_mitophagy) == rownames(design_int)))

fit_int <- lmFit(gsva_scores_mitophagy, design_int)
fit_int <- eBayes(fit_int)

res_int <- topTable(
  fit_int,
  coef = grep("DoseNum:Donor", colnames(design_int), value = TRUE)
)

print(res_int)
```

    ##                                DoseNum.Donorb DoseNum.DonorX DoseNum.DonorY
    ## PD_mitophagy_regulators           0.016292110   -0.088473976    -0.14433318
    ## PD_mitophagy_fusion_fission       0.082781536    0.084699255     0.14265532
    ## PD_mitophagy_core                -0.011704453    0.018781404    -0.05914819
    ## PD_mitophagy_autophagy_markers    0.027607613    0.008321379     0.06495246
    ## PD_mitophagy_receptors_kinases   -0.009966587   -0.054103503    -0.04240591
    ##                                DoseNum.DonorZ       AveExpr         F    P.Value
    ## PD_mitophagy_regulators           -0.14030508  6.523053e-02 3.6664980 0.01630457
    ## PD_mitophagy_fusion_fission        0.02941902 -4.653443e-06 2.5133044 0.06471846
    ## PD_mitophagy_core                 -0.05843535 -1.239988e-02 1.6418568 0.19232274
    ## PD_mitophagy_autophagy_markers    -0.01884226 -2.372134e-02 0.8552081 0.50298433
    ## PD_mitophagy_receptors_kinases    -0.02733291 -1.229879e-02 0.6164044 0.65453410
    ##                                 adj.P.Val
    ## PD_mitophagy_regulators        0.08152287
    ## PD_mitophagy_fusion_fission    0.16179616
    ## PD_mitophagy_core              0.32053791
    ## PD_mitophagy_autophagy_markers 0.62873041
    ## PD_mitophagy_receptors_kinases 0.65453410

``` r
################################################################
#PART3 - Final Plots and Conclusions
################################################################

library(tidyverse)

# Convert GSVA scores to long format
gsva_long <- gsva_scores_mitophagy %>%
  as.data.frame() %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(
    cols = -Pathway,
    names_to = "SampleID",
    values_to = "GSVA_score"
  )

# Merge with metadata
gsva_long <- gsva_long %>%
  left_join(
    targets %>% select(SampleID, DoseNum, Donor),
    by = "SampleID"
  )

stopifnot(!any(is.na(gsva_long$DoseNum)))

ggplot(gsva_long, aes(x = DoseNum, y = GSVA_score)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 1
  ) +
  facet_wrap(~ Pathway, scales = "free_y") +
  labs(
    title = "Radiation dose-dependent modulation of mitophagy modules",
    x = "Radiation dose (Gy)",
    y = "GSVA enrichment score"
  ) +
  theme_bw(base_size = 13)
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](Complemetary_Analysis_PD_Pathways_and_Mitophagy_files/figure-gfm/Combined%20gene%20sets%20list-1.png)<!-- -->

``` r
ggplot(gsva_long, aes(x = DoseNum, y = GSVA_score, color = Donor)) +
  geom_point(size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  facet_wrap(~ Pathway, scales = "free_y") +
  labs(
    title = "Donor-specific radiation response across mitophagy modules",
    x = "Radiation dose (Gy)",
    y = "GSVA enrichment score"
  ) +
  theme_bw(base_size = 13)
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](Complemetary_Analysis_PD_Pathways_and_Mitophagy_files/figure-gfm/Combined%20gene%20sets%20list-2.png)<!-- -->

``` r
effect_df <- res_mitophagy %>%
  as.data.frame() %>%
  rownames_to_column("Pathway")

ggplot(effect_df, aes(x = reorder(Pathway, logFC), y = logFC)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Effect of radiation dose on mitophagy submodules",
    x = "Mitophagy module",
    y = "GSVA logFC per Gy"
  ) +
  theme_bw(base_size = 13)
```

![](Complemetary_Analysis_PD_Pathways_and_Mitophagy_files/figure-gfm/Combined%20gene%20sets%20list-3.png)<!-- -->

``` r
ggplot(effect_df, aes(x = logFC, y = -log10(adj.P.Val), label = Pathway)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(
    title = "Dose-dependent mitophagy modulation",
    x = "Effect size (logFC per Gy)",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 13)
```

![](Complemetary_Analysis_PD_Pathways_and_Mitophagy_files/figure-gfm/Combined%20gene%20sets%20list-4.png)<!-- -->

``` r
################################################################
#PART4 - PD Risk Analysis
################################################################


# DEFINE GENE LISTS

# 1. Full List of PD Risk / Vulnerability Genes
pd_risk_genes <- c(
  # Genetic PD / mitochondrial
  "SNCA", "LRRK2", "PINK1", "PRKN", "PARK7", "DJ1",
  "ATP13A2", "FBXO7", "UCHL1",
  
  # Mitochondrial stress / oxidative damage
  "NDUFS1", "NDUFS2", "NDUFV1", "NDUFA8",
  "VDAC1", "VDAC2",
  
  # Lysosomal / autophagy axis
  "GBA", "CTSD", "LAMP1", "LAMP2",
  
  # Dopaminergic vulnerability markers
  "TH", "SLC6A3", "SLC18A2", "DDC"
)

# 2. Define the Mitophagy Core Genes (The ones in your GSVA predictor)
# We must exclude these to avoid mathematical circularity
mitophagy_core_genes <- c("PINK1", "PRKN", "PARK7", "FBXO7", "ATP13A2", "SNCA", "UCHL1", "SNCAIP")

# 3. Create the "CLEAN" list (Risk genes MINUS Mitophagy Core)
pd_genes_clean <- setdiff(pd_risk_genes, mitophagy_core_genes)

# Check overlap against expression matrix
pd_genes_clean <- intersect(pd_genes_clean, rownames(expr_by_gene))

cat("Number of genes in Clean PD Risk Score:\n")
```

    ## Number of genes in Clean PD Risk Score:

``` r
print(length(pd_genes_clean))
```

    ## [1] 15

``` r
print(pd_genes_clean) # Verify that PINK1/PRKN are NOT here
```

    ##  [1] "LRRK2"   "NDUFS1"  "NDUFS2"  "NDUFV1"  "NDUFA8"  "VDAC1"   "VDAC2"   "GBA"    
    ##  [9] "CTSD"    "LAMP1"   "LAMP2"   "TH"      "SLC6A3"  "SLC18A2" "DDC"

``` r
#CALCULATE CLEAN SCORES


# Calculate score using ONLY non-mitophagy genes
pd_risk_score_clean <- colMeans(expr_by_gene[pd_genes_clean, ], na.rm = TRUE)

# Add to targets
targets$PD_Risk_Score_Clean <- pd_risk_score_clean[targets$SampleID]

# Ensure Mitophagy GSVA score is present (from previous steps)
# Assuming 'gsva_scores_mitophagy' exists from your previous run
targets$Mitophagy_Core_GSVA <- gsva_scores_mitophagy["PD_mitophagy_core", ]

# TEST RADIATION -> CLEAN PD RISK (Total Effect)

cat("\n--- MODEL A: Radiation effect on Clean PD Risk Score ---\n")
```

    ## 
    ## --- MODEL A: Radiation effect on Clean PD Risk Score ---

``` r
fit_pd_risk_clean <- lm(PD_Risk_Score_Clean ~ DoseNum + Donor, data = targets)
summary(fit_pd_risk_clean)
```

    ## 
    ## Call:
    ## lm(formula = PD_Risk_Score_Clean ~ DoseNum + Donor, data = targets)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -355.66  -69.54   12.42   84.00  281.83 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 1804.161     69.991  25.777 3.03e-16 ***
    ## DoseNum       32.457      9.454   3.433  0.00279 ** 
    ## Donorb      -539.427     89.887  -6.001 8.96e-06 ***
    ## DonorX        26.421     89.887   0.294  0.77199    
    ## DonorY       -12.667     89.887  -0.141  0.88942    
    ## DonorZ      -473.056     89.887  -5.263 4.44e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 142.1 on 19 degrees of freedom
    ## Multiple R-squared:  0.8258, Adjusted R-squared:  0.7799 
    ## F-statistic: 18.01 on 5 and 19 DF,  p-value: 1.266e-06

``` r
# TEST MITOPHAGY -> CLEAN PD RISK (The "Fire Test")

# Does mitophagy activation predict risk in OTHER pathways (lysosomal/metabolic)?

cat("\n--- MODEL B: Association between Mitophagy Core and Clean PD Risk ---\n")
```

    ## 
    ## --- MODEL B: Association between Mitophagy Core and Clean PD Risk ---

``` r
fit_assoc_clean <- lm(PD_Risk_Score_Clean ~ Mitophagy_Core_GSVA + Donor, data = targets)
summary(fit_assoc_clean)
```

    ## 
    ## Call:
    ## lm(formula = PD_Risk_Score_Clean ~ Mitophagy_Core_GSVA + Donor, 
    ##     data = targets)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -244.143  -84.217    1.908   92.688  267.199 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          1802.10      68.56  26.285  < 2e-16 ***
    ## Mitophagy_Core_GSVA   551.03     152.91   3.604 0.001893 ** 
    ## Donorb               -419.06      94.30  -4.444 0.000278 ***
    ## DonorX                -55.23      91.05  -0.607 0.551279    
    ## DonorY                305.60     124.81   2.449 0.024223 *  
    ## DonorZ               -282.45     102.83  -2.747 0.012827 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 139.4 on 19 degrees of freedom
    ## Multiple R-squared:  0.8323, Adjusted R-squared:  0.7881 
    ## F-statistic: 18.86 on 5 and 19 DF,  p-value: 8.902e-07

``` r
# MEDIATION ANALYSIS (Using Clean Score)

cat("\n--- MODEL C: Mediation Test (Does Mitophagy explain the Risk Score?) ---\n")
```

    ## 
    ## --- MODEL C: Mediation Test (Does Mitophagy explain the Risk Score?) ---

``` r
# Step 1: Rad -> Mitophagy (Already proven, fit_a in previous code)
# Step 2: Rad -> Risk (Model A above)
# Step 3: Rad + Mitophagy -> Risk (The Mediation Model)

fit_c_clean <- lm(PD_Risk_Score_Clean ~ DoseNum + Mitophagy_Core_GSVA + Donor, data = targets)
summary(fit_c_clean)
```

    ## 
    ## Call:
    ## lm(formula = PD_Risk_Score_Clean ~ DoseNum + Mitophagy_Core_GSVA + 
    ##     Donor, data = targets)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -251.208  -63.352   -0.736   70.845  225.063 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          1773.07      65.56  27.045 4.99e-16 ***
    ## DoseNum                20.29      10.29   1.971  0.06434 .  
    ## Mitophagy_Core_GSVA   369.33     169.71   2.176  0.04310 *  
    ## Donorb               -458.75      90.15  -5.089 7.66e-05 ***
    ## DonorX                -28.31      85.93  -0.329  0.74564    
    ## DonorY                200.65     127.91   1.569  0.13413    
    ## DonorZ               -345.30     100.99  -3.419  0.00306 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 129.9 on 18 degrees of freedom
    ## Multiple R-squared:  0.862,  Adjusted R-squared:  0.8161 
    ## F-statistic: 18.75 on 6 and 18 DF,  p-value: 7.636e-07

``` r
# VISUALIZATION (Clean Score)


library(ggplot2)

# Plot 1: Mitophagy predicting Clean PD Risk
plot_clean <- ggplot(targets, aes(x = Mitophagy_Core_GSVA, y = PD_Risk_Score_Clean)) +
  geom_point(aes(color = as.factor(DoseNum)), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  theme_bw() +
  scale_color_viridis_d(name = "Dose (Gy)") +
  labs(
    title = "Mitophagy Core Predicts Independent PD Molecular Risk",
    subtitle = "Analysis using 'Clean' Score (excluding PINK1/PRKN/PARK7 overlap)",
    x = "Mitophagy Core Activity (GSVA)",
    y = "PD Risk Score (Non-Mitophagy Genes)"
  )

print(plot_clean)
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](Complemetary_Analysis_PD_Pathways_and_Mitophagy_files/figure-gfm/Combined%20gene%20sets%20list-5.png)<!-- -->
