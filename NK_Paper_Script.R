# ==============================================================================
# SECTION 0. Load Required Libraries & Environment Setup
# ==============================================================================
#### 0.1 Package Dependencies --------------------------------------------------
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(gplots)
library(DESeq2)
library(dplyr)

#### 0.2 Global Options --------------------------------------------------------
options(stringsAsFactors = FALSE)
options(timeout = 600)
set.seed(1234)


# ==============================================================================
# SECTION 1. Load data
# ==============================================================================
load("./data/NK_Paper_Final_Processed.RData")


# ==============================================================================
# SECTION 2. PCA
# ==============================================================================

#### 1.1 Remove genes with zero variance across all samples (row) --------------
log2_tpm_filtered <- log2_tpm_symbol[apply(log2_tpm_symbol, 1, var) != 0, ]


pca_res <- prcomp(t(log2_tpm_filtered), center = TRUE, scale. = TRUE)


pca_df <- data.frame(PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     Condition = sample_info_sub$condition,
                     Sample = sample_info_sub$sample,
                     Cell_Line = sample_info_sub$cell_line)

var_explained <- summary(pca_res)$importance[2, 1:2] * 100  # PC1, PC2
xlab <- paste0("PC1 (", round(var_explained[1], 1), "%)")
ylab <- paste0("PC2 (", round(var_explained[2], 1), "%)")


#### 1.2 ggplot with label box -------------------------------------------------
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Cell_Line)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), size = 4, show.legend = FALSE) +
  labs(x = xlab, y = ylab,color = "Condition",       
       shape = "Cell Line") +
  scale_color_manual(values = c("Stim" = "#7B68EE", "Unstim" = "#E89E71")) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 13),   
    legend.title = element_text(size = 14)   
  )


# ==============================================================================
# SECTION 3. Heatmap
# ==============================================================================
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
col_color <- c(rep("#7B68EE", 3),rep("#E89E71", 3))
col_color <- ifelse(sample_info_sub$condition == "Stim", "#7B68EE", "#E89E71")

#### 3.1 Figure e --------------------------------------------------------------
heatmap.2(log2_tpm_filtered,
          ColSideColors = col_color,
          col = my_colors,
          density.info = "none",trace = "none",
          scale = "row",
          cexCol = 1.5,cexRow = 0.9,
          margins = c(6, 6))

legend(
  "topright",                      
  legend = c("Unstim", "Stim"),
  fill = c("#E89E71", "#7B68EE"),
  border = NA,
  bty = "n",                      
  cex = 1.2                        
)

#### 3.2 Figure f --------------------------------------------------------------
Cell_Phenotype_geneset = c('ITGAL', 'ITGB2', 'ITGAX', 'B3GAT1', 'CD2', 'CD58', 'SELL', 
                           'LILRB1','KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR3DL1',
                           'KIR2DL5A', 'KIR2DL5B', 'KIR3DL2', 'KLRC1','KLRD1','SIGLEC7',"DLA-DRA",
                           'FCGR3A', 'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DL4', 'KIR3DS1',
                           'KIR2DS5', 'KIR2DS4', 'KLRC2', 'KLRB1', 'CD226', 'CD244',
                           'KLRK1', 'NCR1', 'NCR2', 'SLAMF6', 'KLRF1',
                           'IL2RA', 'IL2RB', 'IL2RG', 'KIT', 'IL7R',
                           'CXCR1', 'CXCR3', 'CXCR4', 'CCR4', 'CCR7',
                           'IL18R1', 'CMKLR1', 'CX3CR1',
                           'FAS', 'FASLG', 'CD40LG', 'TNFSF10',
                           'LAMP1', 'TNFRSF7', 'CD27', 'CD70',"NCAM1","FCGR3B","FCGR3A","CD5",
                           "CD8A","CD8B", "CD96","GZMB","GATA3","ZEB2","NCR3","CD160","RUNX3",
                           "LOC478984","LOC486692")

col_color <- rep("#E89E71", 3)  # For RO-110, RO-111, RO-112

heatmap.2(log2_tpm_symbol[,c("RO-112","RO-110","RO-111")][Cell_Phenotype_geneset[Cell_Phenotype_geneset %in% rownames(log2_tpm_symbol)],],
          col = my_colors,
          ColSideColors = rep("#E89E71", 3) ,
          density.info = "none",trace = "none",
          cexCol = 1.5,cexRow = 1.3,
          margins = c(6, 7))


#### 3.3 Figure g --------------------------------------------------------------
Activation_related_genes <- c(
  "NCR1","PRF1","IL2RA","CD69","CXCL10","GZMA","IL10","IFNG","CD5","LOC478984","LOC486692"
)

col_color <- c(rep("#7B68EE", 3),rep("#E89E71", 3))

heatmap.2(log2_tpm_filtered[Activation_related_genes[Activation_related_genes %in% rownames(log2_tpm_filtered)],],
          col = my_colors,
          ColSideColors = col_color,
          density.info = "none",trace = "none",
          scale = "row",
          cexCol = 1.5,cexRow = 1.5,
          margins = c(6, 10))

legend(
  "topright",                      
  legend = c("Unstim", "Stim"),
  fill = c("#E89E71", "#7B68EE"),
  border = NA,
  bty = "n",                       
  cex = 1.2                        
)


# ==============================================================================
# SECTION 4. DESeq2
# ==============================================================================
#### 4.1 Make object for DESeq2 ------------------------------------------------
dds <- DESeqDataSetFromTximport(txi_final, colData = sample_info_sub, design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

#### 4.2 Run DESeq2 ------------------------------------------------------------
dds <- DESeq(dds)

#### 4.3 Save result ------------------------------------------------------------
res <- results(dds)


# ==============================================================================
# SECTION 5. Finalize DESeq2 Results (Using Ensembl Mapping)
# ==============================================================================
#### 5.1 Convert res to Data Frame ---------------------------------------------
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)

#### 5.2 Merge with Ensembl Mapping --------------------------------------------
res_annotated <- merge(res_df, mapping, by = "ensembl_gene_id", all.x = FALSE)

#### 5.3 Assign Final Gene Identifiers -----------------------------------------
res_annotated$gene <- res_annotated$external_gene_name

# Filtering: Remove rows without p-value
res_filtered <- res_annotated[!is.na(res_annotated$padj),]

#### 5.4 Duplicate by baseMean -------------------------------------------------
# We keep the one with the highest baseMean as it represents the strongest signal
res_final <- res_filtered[order(res_filtered$gene, -res_filtered$baseMean), ]
res_final <- res_final[!duplicated(res_final$gene), ]

# 
rownames(res_final) <- res_final$gene


# ==============================================================================
# SECTION 6. Volcano plot comparing these two groups (Figure d)
# ==============================================================================
#### 6.1 Identify DEGs (Criteria: padj < 0.05 & |log2FC| >= 1) -----------------
res_final$group <- ifelse(
  res_final$padj < 0.05 & res_final$log2FoldChange >= 1, "Up",
  ifelse(res_final$padj < 0.05 & res_final$log2FoldChange <= -1, "Down", "NotSig")
)

#### 6.2 Top 5 Genes by Fold Change (up-regulated) -----------------------------
top_up <- res_final[res_final$group == "Up", ]
top_up <- top_up[order(-top_up$log2FoldChange), ]
top_up <- head(top_up, 5)

#### 6.3 Top 5 Genes by Fold Change (down-regulated) ---------------------------
top_down <- res_final[res_final$group == "Down", ]
top_down <- top_down[order(top_down$log2FoldChange), ]
top_down <- head(top_down, 5)

#### 6.4 Create Labeling Dataframe for High-Interest Genes ---------------------
lab_df <- rbind(top_up, top_down)
lab_df <- lab_df[complete.cases(lab_df[, c("gene", "log2FoldChange", "padj")]), ]

#### 6.5 Volcano Plot Visualization --------------------------------------------
ggplot(res_final, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_manual(values = c("Up" = "#7B68EE", "Down" = "#E89E71", "NotSig" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = lab_df,
    aes(label = gene),
    color = "black",
    size = 4, box.padding = 0.35, point.padding = 0.25,
    max.overlaps = Inf, min.segment.length = 0, force = 2, show.legend = FALSE
  ) +
  labs(title = "Volcano Plot (Stim vs UnStim)",
       x = "log2 Fold Change", y = "-log10 adjusted p-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )

#### 6.6 Directional Validation ------------------------------------------------
# We verify that the log2FoldChange direction correctly represents Stim vs Unstim.
# Expected: Positive log2FC = Higher expression in Stimulated samples.

# Inspect top 3 up-regulated genes to ensure they follow the expected pattern
res_final %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  arrange(-log2FoldChange) %>%
  head(3)

# Manual verification with a known marker (e.g., RPL31 or IFNG)
log2_tpm_filtered["RPL31", ]

