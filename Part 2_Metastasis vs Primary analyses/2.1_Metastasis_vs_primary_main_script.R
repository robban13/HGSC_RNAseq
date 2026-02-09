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
metadata = read.csv("patient_data_estimate.csv",  header=TRUE, row.names=1, check.names=FALSE) #patient data table with the estimate tumour purity score attached
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
                             design = ~ NACT + patient + condition)

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
ggsave("/Users/robertvanagas/Desktop/Fig_2A.png", plot = volcano_plot, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)
#ggsave("/Users/robertvanagas/Desktop/volcano_met_vs_tumor.png", plot = volcano_plot, width = 45, height = 45, units = "mm", bg = "white", dpi = 300)


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
  Sample_site  = c(M   = floral[7], T  = floral[1]),
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

ggsave("~/Desktop/Fig2B.png", plot = heatmap_plot, width = 180, height = 180, units = "mm", bg = "white", dpi = 300)


rm(floral); rm(muted_rain); rm(ann_want); rm(ann_keep); rm(ann_colors); rm(colors); rm(my_breaks); rm(heatmap_plot);rm(colors.use); rm(ann_col)



######################## PROGENy ##############################################################
# Vignette: https://saezlab.github.io/decoupleR/articles/pw_bk.html#loading-packages 

ranks_unlist = matrix(unlist(ranks), ncol=1, dimnames = list(names(ranks), "stat")) # Unlist ranks for progeny and collectri
progeny_sets = decoupleR::get_progeny(organism = "human", top = 500) # Get PROGENy sets

progeny_contrast = decoupleR::run_mlm(mat = ranks_unlist, net = progeny_sets, minsize = 10) %>% #run the MLM model 
  as_tibble() %>%
  mutate(padj = p.adjust(p_value, "BH"), direction = if_else(score > 0, "Metastasis", "Primary")) %>%
  arrange(score) %>% # Add significance labels
  mutate(stars = case_when(padj < 0.001 ~ "***", padj < 0.01  ~ "**", padj < 0.05  ~ "*", TRUE ~ ""), padj_label = if_else(padj < 0.001, "<0.001", sprintf("%.3f", padj)))


source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/barplot_theme.R")
progeny_plot = barplot_progeny(progeny_contrast, right_label = "Metastasis", left_label  = "Primary",
                               prism_palette = "floral", right_color_index = 7, left_color_index  = 1,
                               base_size = 8, padj_max_show = 0.2, show_stars = T)
print(progeny_plot)
ggsave("~/Desktop/Fig 2D.png", plot = progeny_plot, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)



######################## CLUSTERPROFILER PART 1 #######################################
library(clusterProfiler)

nice_h = function(x) str_to_title(gsub("_"," ", sub("^HALLMARK_","", x))) # create a function to remove the hallmark prefix
ranks_gsea = unlist(ranks); 
ranks_gsea = sort(ranks_gsea, decreasing = TRUE)
msig_h = msigdbr(species = "Homo sapiens", category = "H")

# Hallmark TERM2GENE / TERM2NAME
term2gene = msig_h %>% dplyr::select(term = gs_name, gene = gene_symbol) #dplyr::select only the gene symbol and what HALLMARK set it is present in.
term2name = msig_h %>% distinct(term = gs_name) %>% mutate(name = nice_h(term)) #Name labels for our HALLMARK sets


gse_h = clusterProfiler::GSEA(geneList = ranks_gsea, TERM2GENE = term2gene, TERM2NAME = term2name, minGSSize = 10, #only test gene sets with atleast 10 genes in them
                              pAdjustMethod = "BH", pvalueCutoff = 1, seed = TRUE)      #Set seed for reproducability (TRUE/FALSE)

gse_h@result$Description[gse_h@result$ID == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT" #Rename the pathway
gse_h@result$Description[gse_h@result$ID == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"] <- "ROS Pathway"

cp = gse_h@result %>%
  as_tibble() %>%
  transmute(pathway = ID, label = Description, NES, pval = pvalue, padj = p.adjust, size = setSize, direction = ifelse(NES > 0, "Metastasis", "Primary")) %>%
  arrange(NES) %>% 
  mutate(stars = case_when(padj < 0.001 ~ "***", padj < 0.01  ~ "**", padj < 0.05  ~ "*", TRUE ~ ""), padj_label = ifelse(padj < 0.001, "<0.001", sprintf("%.3f", padj)))

######################## PLOT LOGIC CLUSTERPROFILER::GSEA ############################
source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/barplot_theme.R")
gsea_plot = barplot_gsea(cp, right_label = "Metastasis", left_label  = "Primary", 
                         prism_palette = "floral", right_color_index = 7, left_color_index  = 1,
                         top_up_pathways = 100, padj_max_show = 2, top_down_pathways = 100,
                         base_size = 8, title = "Hallmark GSEA") #filter to only plot the top 10 based on NES. 
print(gsea_plot)
ggsave("~/Desktop/Fig S1C.png", plot = gsea_plot, width = 180, height = 180, units = "mm", bg = "white", dpi = 300)


source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/cnet_theme_v2.R")

cnet_plot <- cnetplot_robert_v2(gse_res = gse_h, DESeq2_results_obj = res_filtered,
                                n_terms = 4, pathway_direction  = "up", rank_by = "padj",
                                title = "GSEA: Core enriched genes for top 4 enriched hallmarks
in omental metastases after removal of overlapping adipogenesis genes")
print(cnet_plot)

ggsave("~/Desktop/Fig S1D.png", plot = cnet_plot, width = 180, height = 180, units = "mm", bg = "white", dpi = 300)




######################## CLUSTERPROFILER PART 2 AFTER REMOVING SHARED ADIPOGENESIS GENES #######################################
ranks_gsea = unlist(ranks); 
ranks_gsea = sort(ranks_gsea, decreasing = TRUE)
msig_h = msigdbr(species = "Homo sapiens", category = "H")

# Hallmark TERM2GENE / TERM2NAME
term2gene = msig_h %>% dplyr::select(term = gs_name, gene = gene_symbol) #dplyr::select only the gene symbol and what HALLMARK set it is present in.
term2name = msig_h %>% distinct(term = gs_name) %>% mutate(name = nice_h(term)) #Name labels for our HALLMARK sets

#We knoe that we have tons of adipose tissue in our samples. A few genes related to adipogenesis is also shared with other 
# pathways such as EMT and KRAS. In order to account for signal coming from adipose tissue in EMT, KRAS and other hallmark sets
# We identify what SHARED genes there are between adipogenesis and fatty acid metabolism and the remaining hallmarks, and remove them from our ranked gene list. 
# Run these three commands to remove the genes that are shared. 
excluded_terms <- c("HALLMARK_ADIPOGENESIS")
conflict_genes <- intersect(unique(term2gene$gene[term2gene$term %in% excluded_terms]),
                            unique(term2gene$gene[!(term2gene$term %in% excluded_terms)]))

ranks_gsea <- ranks_gsea[!(names(ranks_gsea) %in% conflict_genes)]


gse_h = clusterProfiler::GSEA(geneList = ranks_gsea, TERM2GENE = term2gene, TERM2NAME = term2name, minGSSize = 10, #only test gene sets with atleast 10 genes in them
                              pAdjustMethod = "BH", pvalueCutoff = 1, seed = TRUE)      #Set seed for reproducability (TRUE/FALSE)

gse_h@result$Description[gse_h@result$ID == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT" #Rename the pathway
gse_h@result$Description[gse_h@result$ID == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"] <- "ROS Pathway"

cp = gse_h@result %>%
  as_tibble() %>%
  transmute(pathway = ID, label = Description, NES, pval = pvalue, padj = p.adjust, size = setSize, direction = ifelse(NES > 0, "Metastasis", "Primary")) %>%
  arrange(NES) %>% 
  mutate(stars = case_when(padj < 0.001 ~ "***", padj < 0.01  ~ "**", padj < 0.05  ~ "*", TRUE ~ ""), padj_label = ifelse(padj < 0.001, "<0.001", sprintf("%.3f", padj)))

######################## PLOT LOGIC CLUSTERPROFILER::GSEA ############################
source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/barplot_theme.R")
gsea_plot = barplot_gsea(cp, right_label = "Metastasis", left_label  = "Primary", 
                         prism_palette = "floral", right_color_index = 7, left_color_index  = 1,
                         top_up_pathways = 100, padj_max_show = 2, top_down_pathways = 100,
                         base_size = 8, title = "Hallmark GSEA") #filter to only plot the top 10 based on NES. 
print(gsea_plot)
ggsave("~/Desktop/Fig S1E.png", plot = gsea_plot, width = 180, height = 180, units = "mm", bg = "white", dpi = 300)


source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/cnet_theme_v2.R")

cnet_plot <- cnetplot_robert_v2(gse_res = gse_h, DESeq2_results_obj = res_filtered,
                                n_terms = 5, pathway_direction  = "up", rank_by = "padj",
                                title = "GSEA: Core enriched genes for top 5 enriched hallmarks
in omental metastases after removal of overlapping adipogenesis genes")
print(cnet_plot)

ggsave("~/Desktop/Fig S1F.png", plot = cnet_plot, width = 180, height = 180, units = "mm", bg = "white", dpi = 300)


source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/barplot_theme.R")
gsea_plot = barplot_gsea(cp, right_label = "Metastasis", left_label  = "Primary", 
                         prism_palette = "floral", right_color_index = 7, left_color_index  = 1,
                         top_up_pathways = 10, padj_max_show = 2, top_down_pathways = 5, #Only 4 were down 
                         base_size = 8, title = "Hallmark GSEA") #filter to only plot the top 10 based on NES. 
print(gsea_plot)
ggsave("~/Desktop/Fig 2C.png", plot = gsea_plot, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)


#enrichplot::gseaplot2(gse_h, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
#gseaplot(gse_h, by = "all", geneSetID = "HALLMARK_ADIPOGENESIS")


######################## CollecTRI ##############################################################
# vignette: https://saezlab.github.io/decoupleR/articles/tf_bk.html 

ranks_unlist   = matrix(unlist(ranks), ncol=1, dimnames = list(names(ranks), "stat"))

collectri_sets = decoupleR::get_collectri(organism = "human") # Get the CollecTRI sets 

tf_contrast = decoupleR::run_mlm(mat = ranks_unlist, net = collectri_sets, minsize = 10) %>% # MLM (decoupleR TF vignette)
  as_tibble() %>%
  mutate(padj = p.adjust(p_value, "BH"), direction = if_else(score > 0, "Metastasis", "Primary")) %>%
  arrange(score) %>%
  rename(TF = source, contrast = condition) %>%
  mutate(stars = case_when(padj < 0.001 ~ "***", padj < 0.01  ~ "**", padj < 0.05  ~ "*", TRUE ~ ""), padj_label = if_else(padj < 0.001, "<0.001", sprintf("%.3f", padj)))


######################## Plot: TF activity) ###########################
source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/barplot_theme.R")
collectri_plot = barplot_collectri(tf_contrast, right_label = "Metastasis", left_label= "Primary", 
                                   prism_palette = "floral", right_color_index = 7, left_color_index  = 1,
                                   top_up_pathways = 15, padj_max_show = 0.2, top_down_pathways = 5,
                                   base_size = 8, title = "TF Activity")
print(collectri_plot)

ggsave("~/Desktop/Fig 2E.png", plot = collectri_plot, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)


#Visualise the expression of the TF mRNA itself from the DESEq2 results 
#Purpose: To see if the TF activity agrees with the TF mRNA expression 

#Select top n TFs to include in volcano plot, based on activity score 
n_tfs_up   <- 13
n_tfs_down <- 2  

tf_up <- tf_contrast %>% filter(score > 0) %>% slice_max(order_by = abs(score), n = n_tfs_up)
tf_down <- tf_contrast %>% filter(score < 0) %>% slice_max(order_by = abs(score), n = n_tfs_down)
tf_top <- bind_rows(tf_up, tf_down)

tfs_of_interest <- unique(tf_top$TF)
tfs_in_res <- intersect(tfs_of_interest, rownames(resLFC)) #Intersect TFs with DESeq2 result rows
resLFC_tfs_only <- resLFC[tfs_in_res, , drop = FALSE] # Subset DESeq2 results to TFs only


source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/volcano_theme.R")
tf_volcano <- volcano_theme(res= resLFC_tfs_only, p_cutoff  = 0.1,
                            right_label = "Metastasis",left_label = "Primary", prism_palette = "floral",
                            right_color_index = 7, left_color_index  = 1, title  = "Differential Expression
of TF genes",
                            repel_force  = 0.4, repel_max_time  = 0.2, max_overlaps = 80,
                            white_spaace_above = 1.25, side_label_size  = 2, gene_label_size  = 1.5, font_size = 8)

print(tf_volcano) # The TFAP2B gene is not in or resLFC since it gets pre-filtered due to low expression! 

ggsave("/Users/robertvanagas/Desktop/Fig S1B.png", plot  = tf_volcano, width  = 90, height = 90, units  = "mm", bg = "white", dpi = 300)
