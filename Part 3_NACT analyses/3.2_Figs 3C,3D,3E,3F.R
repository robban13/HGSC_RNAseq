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

LEFT_LABEL = "Treatment Naive"
RIGHT_LABEL = "Post-NACT"

######################## DESeq2 ###############################################################
metadata$condition = factor(metadata$condition, levels=c("T","M")) # To order so we have Metastasis vs Tumor 
metadata$patient =   factor(metadata$patient) 
metadata$NACT =    factor(metadata$NACT) 
metadata$HRD_test =factor(metadata$HRD_test)
str(metadata)

dds = DESeqDataSetFromMatrix(round(as.matrix(counts)),  #we need to round since we got non integer counts from salmon (should be integers but the length scaled are not!)
                             colData = metadata,
                             design = ~ estimate_tumor_purity + patient + condition + NACT + condition:NACT)

######################## PRE-FILTERING ###########################################################
min_samples <- ceiling(0.25 * ncol(dds)) # Over 10 counts in at least 25% of  samples, (more reproducable results compareed to look at counts in n samplesi n each group since this changes depending on design)
keep <- rowSums(counts(dds) >= 10) >= min_samples # keep only genes that have at least 10 counts in at least 25% of samples included in design
dds <- dds[keep,] ; rm(keep)

######################## DESeq2 ###########################################################
# https://nbisweden.github.io/workshop-RNAseq/2111/lab_dge.html#4_Testing 

dds = DESeq(dds) 
resultsNames(dds)



#################### PRIMARY 

res_primary = results(dds, name = "NACT_Yes_vs_No") #Since T (primary tumour is our default level, running this will investigate genes upregulated post-NACT in PRIMARY TUMOURS ONLY) based on resultsnames
res_primary_LFC = lfcShrink(dds, res = res_primary, coef = "NACT_Yes_vs_No" , type = "apeglm")

res_primary_LFC$rank = res_primary_LFC$log2FoldChange * (-log10(res_primary_LFC$pvalue))  # Here we combine both logFC and p-value for ranking
res_primary_LFC = na.omit(res_primary_LFC)

res_primary_filtered = res_primary_LFC[res_primary_LFC$padj < 0.1 & abs(res_primary_LFC$log2FoldChange) > 0.2,] #Filter based on padj to only include significant genes for downstream analysis steps
ranks_primary = res_primary_filtered$rank
names(ranks_primary) = rownames(res_primary_filtered)


EnhancedVolcano(res_primary_LFC, lab = rownames(res_primary_LFC), x = "log2FoldChange", y = "padj", pCutoff = 0.1, FCcutoff = 0, title = "", subtitle = "")
EnhancedVolcano(res_primary_filtered, lab = rownames(res_primary_filtered), x = "log2FoldChange", y = "padj", pCutoff = 0.1, FCcutoff = 0, title = "", subtitle = "")



up_primary <- rownames(res_primary_LFC)[res_primary_LFC$padj < 0.1 & res_primary_LFC$log2FoldChange > 0.2]
down_primary <- rownames(res_primary_LFC)[res_primary_LFC$padj < 0.1 & res_primary_LFC$log2FoldChange < 0.2]
#################### Omental 



res_omental = results(dds, name = "conditionM.NACTYes")
res_omental_LFC = lfcShrink(dds, res = res_omental, coef = "conditionM.NACTYes" , type = "apeglm")

res_omental_LFC$rank = res_omental_LFC$log2FoldChange * (-log10(res_omental_LFC$pvalue))  # Here we combine both logFC and p-value for ranking
res_omental_LFC = na.omit(res_omental_LFC)

res_omental_filtered = res_omental_LFC[res_omental_LFC$padj < 0.1 & abs(res_omental_LFC$log2FoldChange) > 0.2,] #Filter based on padj to only include significant genes for downstream analysis steps
ranks_omental = res_omental_filtered$rank
names(ranks_omental) = rownames(res_omental_filtered)


EnhancedVolcano(res_omental_LFC, lab = rownames(res_omental_LFC), x = "log2FoldChange", y = "padj", pCutoff = 0.1, FCcutoff = 0, title = "", subtitle = "")
EnhancedVolcano(res_omental_filtered, lab = rownames(res_omental_filtered), x = "log2FoldChange", y = "padj", pCutoff = 0.1, FCcutoff = 0, title = "", subtitle = "")

up_omental <- rownames(res_omental_LFC)[res_omental_LFC$padj < 0.1 & res_omental_LFC$log2FoldChange > 0.2]
down_omental <- rownames(res_omental_LFC)[res_omental_LFC$padj < 0.1 & res_omental_LFC$log2FoldChange < 0.2]

##########

venn_list_4 <- list(
  "Upregulated
in Primary"    = up_primary,
  "Downregulated
in Primary"  = down_primary,
  "Upregulated
in Omental"   =  up_omental,
  "Downregulated
in Omental" =  down_omental)

library(eulerr)
fit <- euler(venn_list_4)
plot(fit, quantities = TRUE, fills = TRUE)

library(scales) # Required for alpha/lightening

# 1. Extract colors from ggprism 'floral' palette
# Index 1 is for Primary, Index 7 is for Omental
primary_base <- prism_color_pal("floral")(12)[1] # Coral/Red
omental_base <- prism_color_pal("floral")(12)[7] # Purple

# 2. Create the color vector
# Order: Up-Primary, Down-Primary, Up-Omental, Down-Omental
# We use alpha = 0.4 for downregulated to make them "tints"
venn_colors <- c(
  primary_base,                # Upregulated Primary (Full)
  alpha(primary_base, 0.6),    # Downregulated Primary (Light)
  omental_base,                # Upregulated Omental (Full)
  alpha(omental_base, 0.6)     # Downregulated Omental (Light)
)
png(filename = "/Users/robertvanagas/Desktop/Advanced Cancer biology/Report/Figures/Raw figures/Fig 3E.png", width = 90, height = 90, units = "mm", res = 300)

# 3. Plot with Helvetica and correct styling
plot(fit, 
     fills = list(fill = venn_colors, alpha = 0.6), # Overall transparency for overlaps
     labels = list(fontfamily = "Helvetica", cex = 0.5),
     quantities = list(fontfamily = "Helvetica", cex = 0.5),
     edges = list(col = "black", lwd = 1)) # White edges look very clean with ggprism colors


dev.off()

genes_up_primary  = unique(up_primary)
genes_down_primary = unique(down_primary)

background_genes = unique(rownames(res_primary_LFC)) 

con = pipe("pbcopy", "w"); writeLines(genes_up_primary, con); close(con)   # copy the gene list to clipboard
con = pipe("pbcopy", "w"); writeLines(genes_down_primary, con); close(con) 
con = pipe("pbcopy", "w"); writeLines(background_genes, con); close(con)



genes_up_omental  = unique(up_omental)
genes_down_omental = unique(down_omental)

background_genes = unique(rownames(res_omental_LFC)) 

con = pipe("pbcopy", "w"); writeLines(genes_up_omental, con); close(con)   # copy the gene list to clipboard
con = pipe("pbcopy", "w"); writeLines(genes_down_omental, con); close(con) 
con = pipe("pbcopy", "w"); writeLines(background_genes, con); close(con)




source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/volcano_theme.R")
volcano_plot_primary = volcano_theme(res = res_primary_filtered, p_cutoff = 0.1, right_label = "Post-NACT", left_label = "Treatment Naive" ,
                             prism_palette = "floral", right_color_index = 11, left_color_index = 4, title = "Primary Tumours",
                             repel_force = 0.005, repel_max_time  = 0.2, max_overlaps = 80,
                             white_spaace_above = 1.25,
                             side_label_size = 2, gene_label_size = 1.5, font_size = 8)

print(volcano_plot_primary)
ggsave("/Users/robertvanagas/Desktop/Fig 3C.png", plot = volcano_plot_primary, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)


source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/volcano_theme.R")
volcano_plot_omental = volcano_theme(res = res_omental_filtered, p_cutoff = 0.1, right_label = "Post-NACT", left_label = "Treatment Naive", 
                             prism_palette = "floral", right_color_index = 11, left_color_index = 4, title = "Omental Metastases",
                             repel_force = 0.005, repel_max_time  = 0.2, max_overlaps = 80,
                             white_spaace_above = 1.25,
                             side_label_size = 2, gene_label_size = 1.5, font_size = 8)

print(volcano_plot_omental)


ggsave("/Users/robertvanagas/Desktop/Fig 3D.png", plot = volcano_plot_omental, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)


######################## RANKS ##############################################################

ranks = ranks_primary



######################## CLUSTERPROFILER #######################################
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
  transmute(pathway = ID, label = Description, NES, pval = pvalue, padj = p.adjust, size = setSize, direction = ifelse(NES > 0, "Post-NACT", "Treatment Naive")) %>%
  arrange(NES) %>% 
  mutate(stars = case_when(padj < 0.001 ~ "***", padj < 0.01  ~ "**", padj < 0.05  ~ "*", TRUE ~ ""), padj_label = ifelse(padj < 0.001, "<0.001", sprintf("%.3f", padj)))

######################## PLOT LOGIC CLUSTERPROFILER::GSEA ############################
source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/barplot_theme.R")
gsea_plot = barplot_gsea(cp, right_label = "Post-NACT", left_label  = "Treatment Naive", 
                         prism_palette = "floral", right_color_index = 11, left_color_index  = 4,
                         top_up_pathways = 10, padj_max_show = 1, top_down_pathways = 5,
                         base_size = 10, title = "Hallmark GSEA") #filter to only plot the top 10 based on NES. 
print(gsea_plot)

ggsave("~/Desktop/Fig S2D.png", plot = gsea_plot, width = 90, height = 90, units = "mm", bg = "white", dpi = 300)


source("/Users/robertvanagas/Desktop/OSCAT_bulk_rnaseq/scripts/cnetplot_robert.R")
cnetplot_robert_up <- cnetplot_robert(gse_res = gse_h, DESeq2_results_obj = res_primary_filtered, 
                                      pathway_direction = "up", rank_by = "padj", n_terms = 5,
                                      # skip_pathways = c("HALLMARK_ADIPOGENESIS"), # "HALLMARK_FATTY_ACID_METABOLISM")
                                      label_nonoverlap_genes = TRUE, label_overlap_genes = TRUE, title = "Core Enriched Genes",
                                      subtitle = "", alpha = 0.01,   # DE significance for rings (padj threshold)
                                      abs_lfc_thr = 1.5,             # additional |log2FC| threshold for the "high" ring)
                                      ring_color_base = "black",     # ring colour when padj < alpha (any |LFC|)
                                      ring_color_high = "green")     # ring colour when padj < alpha & |LFC| >= abs_lfc_thr

plot(cnetplot_robert_up)

gseaplot(gse_h, by = "all", geneSetID = "HALLMARK_APOPTOSIS")
gseaplot(gse_h, by = "all", geneSetID = "HALLMARK_P53_PATHWAY")
gseaplot(gse_h, by = "all", geneSetID = "HALLMARK_HYPOXIA")
gseaplot(gse_h, by = "all", geneSetID = "HALLMARK_INFLAMMATORY_RESPONSE")
gseaplot(gse_h, by = "all", geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

enrichplot::gseaplot2(gse_h, geneSetID = "HALLMARK_APOPTOSIS")



#Create a cnet plot with all genes being included, not only the core enriched genes defined by gsea!
library(ggraph) 

# Assume 'gse_h' and 'res_primary_LFC' exist.
# Prepare Colors (Squished LogFC)
logFC_vector <- res_primary_LFC$log2FoldChange
names(logFC_vector) <- rownames(res_primary_LFC)
logFC_squished <- ifelse(logFC_vector > 2, 2, ifelse(logFC_vector < -2, -2, logFC_vector))

gse_h_all <- gse_h
all_gene_sets <- gse_h@geneSets
my_ranked_genes <- names(ranks_gsea)

gse_h_all@result$core_enrichment <- sapply(gse_h_all@result$ID, function(id) {
  pathway_genes <- all_gene_sets[[id]]
  genes_in_my_data <- intersect(pathway_genes, my_ranked_genes)
  paste(genes_in_my_data, collapse = "/")
})

# Select Top n Names
top_n_names <- gse_h_all@result %>%
  as_tibble() %>%
  arrange(desc(abs(NES))) %>%
  head(10) %>%
  pull(Description)

set.seed(123) 

library(ggplot2)
library(RColorBrewer)


# --- 1. CREATE CNET PLOT ---
# We generate the plot with node_label = "none" so we can manually add custom layers
cnet_plot_primary_after_NACT <- cnetplot(gse_h_all, foldChange = logFC_squished, showCategory = top_n_names,  node_label = "none", 
                                         categorySize = "pvalue", layout = "fr",  niter = 1000000, color_category = "black") 

# --- 2. APPLY CUSTOM STYLING ---
cnet_plot_primary_after_NACT <- cnet_plot_primary_after_NACT +

  scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", 
                        midpoint = 0, limit = c(-2, 2), oob = scales::squish, name = expression(Log[2]*FC)) +
  
  # Gene Labels
  geom_node_text(aes(label = name, filter = !name %in% top_n_names), 
                 repel = TRUE,
                 size = 2.5,
                 color = "black",
                 family = "Helvetica",
                 bg.color = "white",
                 bg.r = 0.1,
                 segment.color = "black",
                 segment.size = 0.4,
                 min.segment.length = 1,   # Always draw the line if set to 0
                 force = 1,
                 box.padding = 0.4,
                 point.padding = 0.2,
                 max.overlaps = Inf) +
  
  # Hallmark Labels
  geom_node_text(aes(label = name, filter = name %in% top_n_names), 
                 repel = TRUE,
                 size = 3.5,
                 family = "Helvetica",
                 fontface = "bold",
                 bg.color = "white",
                 bg.r = 0,
                 box.padding = 0.5) +
  
  ggtitle("GSEA: Core enrichment genes for enriched Hallmark genesets
          post-NACT in Primary Tumours") +
  theme_void() + 
  theme(
    legend.position = "right", 
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Helvetica", size = 12),
    text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica", size = 8),
    legend.text = element_text(family = "Helvetica", size = 7))

print(cnet_plot_primary_after_NACT)

ggsave("~/Desktop/Fig 3F.png", plot = cnet_plot_primary_after_NACT, width = 180, height = 180, units = "mm", bg = "white", dpi = 600)


