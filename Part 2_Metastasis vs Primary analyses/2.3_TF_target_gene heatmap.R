library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(tibble)
library(ggprism) # Needed for the floral palette

########################################################################
## 1. PREPARE DATA & CALCULATE CLUSTERING (The Missing Piece)
########################################################################

# 1. Define TFs and Matrix
tfs_to_plot <- tf_use_strong
vsd_mat     <- assay(vsd)

# 2. Identify ALL unique genes involved (TF drivers + their targets)
# We need this to calculate the clustering correctly. We cannot cluster 
# on the final matrix because it will have duplicated rows.
targets_in_scope <- tf_target_DE %>%
  filter(TF %in% tfs_to_plot) %>%
  pull(target) %>%
  unique()

all_unique_genes <- unique(c(tfs_to_plot, targets_in_scope))

# 3. Filter to genes actually present in our VST matrix
all_unique_genes <- intersect(all_unique_genes, rownames(vsd_mat))

# 4. Create a subset matrix for clustering calculation
mat_for_clustering <- vsd_mat[all_unique_genes, , drop = FALSE]

# 5. Scale and Cluster (Generate hc_cols)
# This mimics the logic from your "correction" script
mat_cluster_scaled <- t(scale(t(mat_for_clustering)))
mat_cluster_scaled[is.na(mat_cluster_scaled)] <- 0

dist_mat_cols <- dist(t(mat_cluster_scaled), method = "euclidean")
hc_cols       <- hclust(dist_mat_cols, method = "ward.D2")

########################################################################
## 2. PREPARE TARGET DATA FOR PLOTTING
########################################################################

target_data_subset <- tf_target_DE %>%
  dplyr::select(TF, target, log2FoldChange, padj, mor)

# Filter target data to only include valid genes found in our matrix
target_data_valid <- target_data_subset %>%
  filter(target %in% rownames(vsd_mat))

########################################################################
## 3. BUILD THE VISUAL MATRIX (The "Fat" Rows)
########################################################################

matrix_list     <- list()
annot_list      <- list()
display_labels  <- c() 
n_repeats       <- 5 

for(tf in tfs_to_plot) {
  
  # --- A. The TF Driver (FAT VERSION) ---
  if(tf %in% rownames(vsd_mat)) {
    # 1. Duplicate row
    tf_vec   <- vsd_mat[tf, , drop=FALSE]
    tf_block <- tf_vec[rep(1, n_repeats), , drop=FALSE]
    
    # 2. Create IDs & Annotation
    row_ids  <- paste0(tf, "_Driver_", 1:n_repeats)
    rownames(tf_block) <- row_ids
    
    matrix_list[[paste0(tf, "_Driver")]] <- tf_block
    annot_list[[paste0(tf, "_Driver")]]  <- data.frame(rowname = row_ids, TF_Module = tf, Group = "TF Driver")
    
    # 3. Label only the middle row
    labels_vec <- rep("", n_repeats)
    labels_vec[ceiling(n_repeats / 2)] <- tf
    names(labels_vec) <- row_ids
    display_labels <- c(display_labels, labels_vec)
  }
  
  # --- B. Activated Targets ---
  targets_up <- target_data_valid %>% filter(TF == tf, log2FoldChange > 0) %>% pull(target) %>% unique()
  
  if(length(targets_up) > 0) {
    mat_up <- vsd_mat[targets_up, , drop=FALSE]
    new_ids <- paste0(rownames(mat_up), "_", tf, "_Up")
    rownames(mat_up) <- new_ids
    
    matrix_list[[paste0(tf, "_Up")]] <- mat_up
    annot_list[[paste0(tf, "_Up")]]  <- data.frame(rowname = new_ids, TF_Module = tf, Group = "Activated Targets" )
    
    display_labels <- c(display_labels, setNames(rep("", length(new_ids)), new_ids))
  }
  
  # --- C. Repressed Targets ---
  targets_down <- target_data_valid %>% filter(TF == tf, log2FoldChange < 0) %>% pull(target) %>% unique()
  
  if(length(targets_down) > 0) {
    mat_down <- vsd_mat[targets_down, , drop=FALSE]
    new_ids  <- paste0(rownames(mat_down), "_", tf, "_Down")
    rownames(mat_down) <- new_ids
    
    matrix_list[[paste0(tf, "_Down")]] <- mat_down
    annot_list[[paste0(tf, "_Down")]]  <- data.frame(rowname = new_ids, TF_Module = tf, Group = "Repressed Targets")
    
    display_labels <- c(display_labels, setNames(rep("", length(new_ids)), new_ids))
  }
}

# Consolidate Matrix & Annotation
final_mat  <- do.call(rbind, matrix_list)
full_annot <- do.call(rbind, annot_list) 

# Create Final Annotation (removing Group so it has no legend)
final_annot <- full_annot %>% 
  remove_rownames() %>% 
  column_to_rownames("rowname") %>% 
  dplyr::select(-Group)

########################################################################
## 4. DEFINE GAPS & COLORS
########################################################################

# Calculate gaps using the full annotation
group_vec   <- paste0(full_annot$TF_Module, full_annot$Group)
gap_indices <- which(group_vec[-1] != group_vec[-length(group_vec)])

# Re-establish palettes (in case they were removed)
floral     <- ggprism::ggprism_data$fill_palettes$floral
muted_rain <- ggprism::ggprism_data$fill_palettes$muted_rainbow

# Define Colors
ann_colors <- list() # Reset to be safe
ann_colors$TF_Module <- setNames(viridis(length(tfs_to_plot)), tfs_to_plot)

# Scale Matrix
final_mat_scaled <- t(scale(t(final_mat)))
final_mat_scaled[is.na(final_mat_scaled)] <- 0 

########################################################################
## 5. UPDATE METADATA (for plot labels)
########################################################################

# Re-create ann_col based on current metadata to ensure it exists
ann_col_plot <- metadata %>%
  dplyr::select(condition, HRD_test, NACT, patient, estimate_tumor_purity) %>%
  dplyr::rename(
    `Sample site` = condition,
    HRD = HRD_test,
    Patient = patient,
    Purity = estimate_tumor_purity
  ) %>%
  mutate(`Sample site` = recode(`Sample site`, "T" = "Primary", "M" = "Metastasis"))

# Update Annotation Colors
ann_colors$`Sample site` = c(Metastasis = floral[7], Primary = floral[1])
ann_colors$HRD           = c(Pos = floral[10], Neg = floral[8])
ann_colors$NACT          = c(Yes = floral[11], No = floral[4])
ann_colors$Patient       = setNames(colorRampPalette(muted_rain)(nlevels(factor(ann_col_plot$Patient))), levels(factor(ann_col_plot$Patient)))
ann_colors$Purity        = inferno(100, begin = 0, end = 1)

# Define breaks and colors
colors     = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use = grDevices::colorRampPalette(colors = colors)(100)
my_breaks  = c(seq(-2, 0, length.out = ceiling(100 / 2) + 1), seq(0.05, 2, length.out = floor(100 / 2)))

########################################################################
# 6. PLOT
########################################################################

p_genes <- pheatmap(
  final_mat_scaled, 
  cluster_rows      = FALSE,      
  cluster_cols      = hc_cols,     
  annotation_col    = ann_col_plot, 
  annotation_row    = final_annot,
  annotation_colors = ann_colors, 
  color             = colors.use, 
  breaks            = my_breaks,
  gaps_row          = gap_indices, 
  labels_row        = display_labels[rownames(final_mat)], 
  show_rownames     = TRUE, 
  fontsize          = 12,
  fontsize_row      = 12,     
  main              = "",
  border_color      = NA)

print(p_genes)
# OBS!!! Clustering is skewed due to row duplication! This means that we can not draw conclusions from this clustering!!
# It is still included due to visual reasons, but manually remove the clustering from final plot in the report

ggsave("~/Desktop/FigS2C.png", plot = p_genes, width = 200, height = 200, units = "mm", bg = "white", dpi = 300)
