library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggprism)
library(ggnewscale)
library(ggrepel)

######################## DATA PREPARATION ########################

counts   <- read.table("salmon.merged.gene_counts_length_scaled.tsv", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
metadata <- read.csv("patient_data_estimate.csv",  header=TRUE, row.names=1, check.names=FALSE)
counts   <- counts[, setdiff(colnames(counts), "gene_name"), drop = FALSE]

metadata$condition <- factor(metadata$condition, levels=c("T","M")) # Sample site is called named condition. T = Primary tumor, M = Omental metastasis
metadata$patient   <- factor(metadata$patient)

dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), colData = metadata, design = ~ 1)
min_samples <- ceiling(0.25 * ncol(dds))
keep <- rowSums(counts(dds) >= 10) >= min_samples
dds <- dds[keep,]
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup = c("condition","NACT","estimate_tumor_purity", "HRD_test"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$sample   <- rownames(pcaData)
pcaData$HRD_test <- factor(pcaData$HRD_test, levels = c("Neg", "Pos"))

######################## PCA PLOTTING ########################

pcaplot_purity <- ggplot(pcaData, aes(PC1, PC2)) +
  
  geom_point(aes(shape = condition, fill = estimate_tumor_purity, colour = NACT), size = 3.8, stroke = 1.8) +
  scale_shape_manual(
    values = c(T = 21, M = 24), 
    labels = c(T = "Primary", M = "Metastasis"),
    name   = "Sample site") +
  scale_fill_viridis_c(option = "inferno",  end    = 1,  name   = "Tumour purity") +
  scale_colour_manual( values = c(No = "black", Yes = "red"), labels = c(No = "No", Yes = "NACT"), name   = "NACT") +
  guides(
    colour = guide_legend(override.aes = list(fill = NA, shape = 21, size = 4.5, stroke = 2.2), order = 1),
    shape  = guide_legend(override.aes = list(fill = NA, colour = "grey40"), order = 2),
    fill   = guide_colourbar(title.position = "top", order = 3)) +
  new_scale_color() + 
  geom_text_repel(
    aes(x = PC1, y = PC2, label = sample, colour = HRD_test), 
    inherit.aes = FALSE, show.legend = FALSE, size = 3, fontface = "bold", family = "Helvetica",
    box.padding = 0.6, point.padding = 0.15, 
    min.segment.length = grid::unit(0.05, "lines"), max.overlaps = Inf) +
  geom_point(aes(x = PC1, y = PC2, colour = HRD_test), alpha = 0, size = 0, inherit.aes = FALSE) +
  scale_colour_manual(name   = "HRD Status", values = c("Pos" = "blue", "Neg" = "black"), guide  = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 15), order = 4)) +
  labs(title = "PCA", x = paste0("PC1 (", percentVar[1], "%)"), y = paste0("PC2 (", percentVar[2], "%)")) +
  theme_prism(base_size = 12, base_family = "Helvetica") +
  theme(
    plot.margin       = margin(18, 12, 12, 12),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 2),
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    axis.line         = element_blank(),
    axis.ticks        = element_line(linewidth = 0.8),
    axis.ticks.length = unit(3, "pt"),
    legend.title      = element_text(face = "bold"),
    legend.text       = element_text(face = "bold"))

print(pcaplot_purity)

ggsave("/Users/robertvanagas/Desktop/Fig 1C.png", plot = pcaplot_purity, width = 200, height = 175, units  = "mm", bg  = "white", dpi= 300)


library(variancePartition)
library(limma)
library(edgeR)
library(ggplot2)
library(BiocParallel)
param = MulticoreParam(workers = 4)
register(param, default = TRUE)


## 1. Prepare expression matrix and metadata
# Expression matrix: genes x samples
expr_filt = assay(vsd)
meta = as.data.frame(colData(vsd))

# Covariates to include
vars_to_test = c("condition", "NACT", "estimate_tumor_purity", "patient", "RIN", "HRD_test", "Age_at_surg")
meta_subset = meta[, vars_to_test, drop = FALSE]

# Set variable types

# Categorical
meta_subset$patient   = factor(meta_subset$patient)
meta_subset$condition = factor(meta_subset$condition)
meta_subset$NACT      = factor(meta_subset$NACT)
meta_subset$HRD_test  = factor(meta_subset$HRD_test)

# Numeric
meta_subset$estimate_tumor_purity = as.numeric(meta_subset$estimate_tumor_purity)
meta_subset$RIN = as.numeric(meta_subset$RIN)
meta_subset$Age_at_surg = as.numeric(meta_subset$Age_at_surg)

## 2. Define variancePartition formula

# All categorical variables are modelled as random effects: (1|...)
# Numeric variables remain as fixed effects (see documentation)

form = ~ (1 | patient) + (1 | condition) + (1 | NACT) + (1 | HRD_test) +
  estimate_tumor_purity + RIN + Age_at_surg

## 4. Fit variancePartition model
varPart = fitExtractVarPartModel(expr_filt, form, meta_subset, BPPARAM = param) #will take a few minutes


## 5. PLOT LOGIC

vp <- sortCols(varPart)  # reorder columns for nicer plotting

vp_renamed <- vp
colnames(vp_renamed) <- c(
  "patient"               = "Patient
ID",
  "Residuals"             = "Residuals",
  "estimate_tumor_purity" = "Tumour
purity",
  "Age_at_surg"           = "Age at
surgery",
  "RIN"                   = "RIN",
  "condition"             = "Sample
site",
  "NACT"                  = "NACT",
  "HRD_test"              = "HRD
status"
)[colnames(vp_renamed)]

vp_violin <- plotVarPart(vp_renamed)
cov_levels <- levels(vp_violin$data$variable)  
floral_pal <- ggprism::ggprism_data$fill_palettes$floral
pal_use    <- setNames(floral_pal[seq_along(cov_levels)], cov_levels)
pal_use["Residuals"] <- "grey70" 


vp_violin <- vp_violin +
  ggprism::theme_prism(base_size = 12, base_family = "Helvetica") +
  scale_fill_manual(values = pal_use) +
  theme(legend.position = "none") +
  labs(title = "variancePartition")

vp_violin

ggsave("/Users/robertvanagas/Desktop/Fig 1B.png", plot = vp_violin, width = 180, height = 90, units = "mm", bg = "white", dpi = 300)

