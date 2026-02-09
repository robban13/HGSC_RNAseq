library(DESeq2)
library(EnhancedVolcano)
library(tidyverse)
library(ggplot2)
library(msigdbr)
library(decoupleR)
library(ggprism)
library(pheatmap)
library(enrichplot)
library(viridis)

######################## INPUTS ##############################################################################
counts   = read.table("salmon.merged.gene_counts_length_scaled.tsv", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
metadata = read.csv("patient_data_estimate.csv",  header=TRUE, row.names=1, check.names=FALSE)
counts = counts[, setdiff(colnames(counts), "gene_name"), drop = FALSE]

LEFT_LABEL = "Primary"
RIGHT_LABEL = "Metastasis"

######################## DESeq2 ###############################################################
metadata$condition = factor(metadata$condition, levels=c("T","M")) # To order so we have Metastasis vs Tumor 
metadata$patient =   factor(metadata$patient) 
metadata$NACT =    factor(metadata$NACT) 
metadata$HRD_test =factor(metadata$HRD_test)
str(metadata)

dds = DESeqDataSetFromMatrix(round(as.matrix(counts)),  #we need to round since we got non integer counts from salmon (should be integers but the length scaled are not!)
                             colData = metadata,
                             design = ~ NACT + patient + estimate_tumor_purity + condition)

######################## PRE-FILTERING ###########################################################
min_samples <- ceiling(0.25 * ncol(dds)) # Over 10 counts in at least 25% of  samples, (more reproducable results compareed to look at counts in n samplesi n each group since this changes depending on design)
keep <- rowSums(counts(dds) >= 10) >= min_samples # keep only genes that have at least 10 counts in at least 25% of samples included in design
dds <- dds[keep,] ; rm(keep)

######################## DESeq2 ###########################################################
# https://nbisweden.github.io/workshop-RNAseq/2111/lab_dge.html#4_Testing 

dds = DESeq(dds) 
res = results(dds, contrast=c("condition","M","T"))
resLFC = lfcShrink(dds, res = res, coef="condition_M_vs_T", type="apeglm")

######################## CREATE RANKING FOR DOWNSTREAM GSEA ####################################
resLFC$rank = resLFC$log2FoldChange * (-log10(resLFC$pvalue))  # Here we combine both logFC and p-value for ranking
resLFC = na.omit(resLFC)

res_filtered = resLFC[resLFC$padj < 0.05,] #Filter based on padj to only include significant genes for downstream analysis steps
ranks = res_filtered$rank
names(ranks) = rownames(res_filtered)

res_filtered_matrix = as.matrix(res_filtered)

######################## VOLCANO PLOT ######################################################################
EnhancedVolcano(resLFC, lab = rownames(resLFC), x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 0, title = "", subtitle = "")

source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/volcano_theme.R")
volcano_plot = volcano_theme(res = resLFC, p_cutoff = 0.05, right_label = "Metastasis", left_label = "Primary", 
                             prism_palette = "floral", right_color_index = 7, left_color_index = 1, title = "Differential expression",
                             repel_force = 0.4, repel_max_time  = 0.2, max_overlaps = 80,
                             white_spaace_above = 1.5,
                             side_label_size = 2, gene_label_size = 1.5, font_size = 8, side_label_vertical = 0)

print(volcano_plot)

#For 45x45 mm  plots, repel_force 0.1, side_label_size 1.5, font_size 6, side_label_vertical 1 
#For 90x 90 mm plots repel_force 0.4, side_label_size 2, font_size 8, side_label_vertical 0
ggsave("/Users/robertvanagas/Desktop/Fig S1A.png", plot = volcano_plot, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)
