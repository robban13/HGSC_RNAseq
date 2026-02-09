library(DESeq2)
library(EnhancedVolcano)
library(tidyverse)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(decoupleR)
library(ggprism)
library(pheatmap)
library(viridis)
######################## INPUTS ##############################################################################
counts   = read.table("salmon.merged.gene_counts_length_scaled.tsv", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
metadata = read.csv("patient_data_estimate.csv", header=TRUE, row.names=1, check.names=FALSE)

counts = counts[, setdiff(colnames(counts), "gene_name"), drop = FALSE]
metadata$patient = factor(metadata$patient)
metadata$NACT =    factor(metadata$NACT) 
metadata$HRD_test =factor(metadata$HRD_test)
str(metadata) 



## Design our analysis
sample_types_to_analyse = c("T", "M") #Since we have paired T and M samples in our dataset, we specify if we only want to look at T or M or both, all depends on our design. 

analysis_variable = "NACT" #column name in metadata to analyse for
left_side_factor  = "No"   #values inside the analysis variable
right_side_factor = "Yes"  #values inside the analysis variable
DESEQ2_DESIGN = ~ estimate_tumor_purity + patient + condition + NACT # Design parameter for DESeq2 

###### Labels for plots
LEFT_LABEL  = "Treatment naive"  
RIGHT_LABEL = "NACT"

# Colors for all plots 
PRISM_PALETTE = "floral" #select one of the ggprism palettes run this to see colors in a theme: preview_theme("floral")
RIGHT_COLOR_INDEX = 11   # color index
LEFT_COLOR_INDEX  = 4


######################## Harmonize inputs & subset (robust to NA/design) ###################
design_vars <- all.vars(DESEQ2_DESIGN) # Variables required by the DESeq2 design

metadata <- metadata %>% filter(.data[["condition"]] %in% sample_types_to_analyse)

# 1) Keep only samples whose analysis variable matches the two selected levels
keep_lvls <- c(left_side_factor, right_side_factor)
metadata <- metadata %>% filter(.data[[analysis_variable]] %in% keep_lvls)

# 2) Drop rows with NA in ANY design variable (DESeq2 requirement)
metadata <- metadata %>% filter(stats::complete.cases(pick(all_of(design_vars))))

# 3) Align counts to metadata
common_samples <- intersect(colnames(counts), rownames(metadata))
counts   <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# 4) Relevel the analysis variable to desired left/right order for contrasts
metadata[[analysis_variable]] <- factor(metadata[[analysis_variable]], levels = keep_lvls)

# Convenience strings used downstream
design_str   <- paste(deparse(DESEQ2_DESIGN), collapse = "")
contrast_var <- analysis_variable
contrast_l   <- left_side_factor
contrast_r   <- right_side_factor
coef_name    <- paste0(contrast_var, "_", contrast_r, "_vs_", contrast_l)

######################## DESeq2 ###############################################################
# Order so we have (right vs left) in results
metadata[[contrast_var]] = factor(metadata[[contrast_var]], levels = c(contrast_l, contrast_r))

dds = DESeqDataSetFromMatrix(round(as.matrix(counts)), # we need to round since we got non integer counts from salmon
                             colData = metadata,
                             design = DESEQ2_DESIGN) 

######################## PRE-FILTERING ###########################################################
min_samples <- ceiling(0.25 * ncol(dds))
keep <- rowSums(counts(dds) >= 10) >= min_samples # keep only genes that have at least 10 counts in at least 25% of samples included in design
dds <- dds[keep,] ; rm(keep)

######################## DESeq2 ###########################################################
# https://nbisweden.github.io/workshop-RNAseq/2111/lab_dge.html#4_Testing 

dds = DESeq(dds) 
res = results(dds, contrast=c(contrast_var, contrast_r, contrast_l))
resLFC = lfcShrink(dds, res = res, coef = coef_name, type = "apeglm")

######################## CREATE RANKING FOR DOWNSTREAM GSEA ####################################
resLFC$rank = resLFC$log2FoldChange * (-log10(resLFC$pvalue)) # Here we combine both logFC and p-value for ranking
resLFC = na.omit(resLFC)

resLFC_matrix = as.matrix(resLFC)
resLFC_matrix <- cbind(Gene = rownames(resLFC_matrix), resLFC_matrix)

res_filtered = resLFC
res_filtered = resLFC[resLFC$padj < 0.1 & abs(resLFC$log2FoldChange)>0,] #Filter based on padj and or Log2FC for downstream GSEA

ranks = res_filtered$rank; names(ranks) = rownames(res_filtered); ranks = na.omit(ranks)

######################## VOLCANO PLOT ######################################################################
EnhancedVolcano(resLFC, lab = rownames(resLFC), x = "log2FoldChange", y = "padj", pCutoff = 0.1, FCcutoff = 0, title = design_str, subtitle = paste0("Left: ", LEFT_LABEL, ", Right: ", RIGHT_LABEL))
EnhancedVolcano(res_filtered, lab = rownames(res_filtered), x = "log2FoldChange", y = "padj", pCutoff = 0.1, FCcutoff = 0, title = design_str, subtitle = paste0("Left: ", LEFT_LABEL, ", Right: ", RIGHT_LABEL))

source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/volcano_theme.R")
volcano_plot = volcano_theme(res = resLFC, p_cutoff = 0.1, right_label = RIGHT_LABEL, left_label = LEFT_LABEL, 
                             prism_palette = PRISM_PALETTE, right_color_index = RIGHT_COLOR_INDEX, left_color_index = LEFT_COLOR_INDEX, title = "Differential expression",
                             repel_force = 0.4, repel_max_time  = 0.2, max_overlaps = 80,
                             white_spaace_above = 1.25,
                             side_label_size = 2, gene_label_size = 1.5, font_size = 8)

print(volcano_plot)

ggsave("/Users/robertvanagas/Desktop/Fig 3A.png", plot = volcano_plot, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)



#ggsave("/Users/robertvanagas/Desktop/volcano_HRD metastasis only.png", plot = volcano_plot, width = 150, height = 150, units = "mm", bg = "white", dpi = 300)
#ggsave("~/Desktop/volcano_met_vs_tumor.svg", plot = volcano_plot, width = 150, height = 150, units = "mm", device = "svg", bg = "transparent")




######################## HEATMAP ###################
vsd = vst(dds, blind=FALSE) # variance stabilizing transformation # VST is used for other purposes than differential testing, ex clustering and heatmap visualisation of genes in samples. 
mat = assay(vsd) # convert the full gene list to a matrix
sig_genes = rownames(res_filtered) #look what genes we considered signfiicant above
mat_heatmap = assay(vsd)[sig_genes, , drop = FALSE] #keep only these genes for the heatmap. 

metadata2 = metadata
metadata2$Purity <- metadata2$estimate_tumor_purity
metadata2$Sample_site <- metadata2$condition
metadata2$Patient <- metadata2$patient

ann_want = c("Patient", "Sample_site", "HRD_test", "NACT", "Purity")
ann_keep = intersect(ann_want, colnames(metadata2))
ann_col  = metadata2[colnames(vsd), ann_keep, drop = FALSE]
# Keep numeric as numeric, others as factors
ann_col[] <- lapply(ann_col, function(x) if (is.numeric(x)) x else as.factor(x))

floral     <- ggprism::ggprism_data$fill_palettes$floral #other factors
muted_rain <- ggprism::ggprism_data$fill_palettes$muted_rainbow #for patient ids

ann_colors <- list(
  Sample_site  = c(M   = floral[7], T   = floral[1]),
  HRD_test   = c(Pos = floral[10], Neg = floral[8]),
  NACT       = c(Yes = floral[11], No  = floral[4]),
  Patient    = setNames(colorRampPalette(muted_rain)(nlevels(ann_col$Patient)), levels(ann_col$Patient)))
ann_colors$Purity <- inferno(100, begin = 0, end = 1)

colors     = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use = grDevices::colorRampPalette(colors = colors)(100)
my_breaks  = c(seq(-2, 0, length.out = ceiling(100 / 2) + 1), seq(0.05, 2, length.out = floor(100 / 2)))

heatmap_plot = pheatmap(mat_heatmap,scale = "row", annotation_col = ann_col, breaks = my_breaks, #add the breaks for better contrast, we clip gene expression to -2 and 2 
                        annotation_colors  = ann_colors, color = colors.use, clustering_method  = "ward.D2",
                        show_rownames= FALSE, main  = "",fontsize = 8, fontsize_col = 8, fontfamily = "Helvetica")

print(heatmap_plot)

ggsave("~/Desktop/Fig 3B.png", plot = heatmap_plot, width = 180, height = 180, units = "mm", bg = "white", dpi = 300)
#ggsave("~/Desktop/heatmap_met_vs_tumor.svg", plot = heatmap_plot, width = 150, height = 150, units = "mm", device = "svg", bg = "transparent")

rm(floral); rm(muted_rain); rm(ann_want); rm(ann_keep); rm(ann_colors); rm(colors); rm(my_breaks); rm(heatmap_plot);rm(colors.use); rm(ann_col)
