library(tidyverse)
library(ggplot2)
library(ggprism)
library(rstatix)
library(immunedeconv)

############### 1. Inputs & Deconvolution
metadata <- read.csv('patient_data.csv', header=TRUE, row.names=1, check.names=FALSE) %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, patient, condition)

tpm = read.table("salmon.merged.gene_tpm.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
tpm = tpm[, setdiff(colnames(tpm), "gene_name"), drop = FALSE]

res_estimate = deconvolute(tpm, method = "estimate")

# Rename map for ESTIMATE specific results
name_map <- c("estimate score" = "Estimate Score", "immune score" = "Immune Score", 
              "stroma score" = "Stroma Score", "tumor purity" = "Tumour Purity")

############### 2. Tidy results & Identify Pairing
td_raw <- res_estimate %>%
  pivot_longer(-cell_type, names_to = "sample", values_to = "score") %>%
  left_join(metadata, by = "sample") %>%
  filter(condition %in% c("T", "M")) %>%
  mutate(
    condition = factor(condition, levels = c("T", "M"), labels = c("Primary", "Metastasis")),
    cell_type = recode(cell_type, !!!name_map)
  ) %>%

filter(cell_type %in% c("Stroma Score", "Tumour Purity")) #plot only these in final plot

# Identify orphans vs pairs
valid_patients <- td_raw %>%
  group_by(cell_type, patient) %>%
  summarise(is_paired = n_distinct(condition) == 2, .groups = "drop")

td_raw <- td_raw %>%
  left_join(valid_patients, by = c("cell_type", "patient")) %>%
  mutate(pairing_status = ifelse(is_paired, "Paired", "Orphan"))

############### 3. Calculate Stats & Variable Bracket Coordinates
stats_df <- td_raw %>%
  filter(pairing_status == "Paired") %>%
  group_by(cell_type) %>%
  summarise(
    p = wilcox.test(score[condition == "Primary"], 
                    score[condition == "Metastasis"], 
                    paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  )

# Calculate variable heights per facet
bracket_df <- td_raw %>%
  group_by(cell_type) %>%
  summarise(
    # Get max of each group to define tip endings
    max_p = max(score[condition == "Primary"], na.rm = TRUE),
    max_m = max(score[condition == "Metastasis"], na.rm = TRUE),
    # Global max per facet to define horizontal bar height
    facet_max = max(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(stats_df, by = "cell_type") %>%
  mutate(
    bracket_y = facet_max + (facet_max * 0.15),
    # Tips reach down to 5% above the highest point of that specific group
    tip_primary_y = max_p + (facet_max * 0.05),
    tip_met_y = max_m + (facet_max * 0.05),
    label_y = bracket_y + (facet_max * 0.02),
    label = paste0("p = ", formatC(p, format = "g", digits = 1))
  )

# Data for connecting lines (Pairs only)
lines_data <- td_raw %>%
  filter(pairing_status == "Paired") %>%
  select(patient, cell_type, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score)

############### 4. Plotting
floral_colors <- ggprism::prism_colour_pal(palette = "floral")(12)
col_primary <- floral_colors[1] 
col_met     <- floral_colors[7] 

estimate_plot <- ggplot(td_raw, aes(x = condition, y = score)) +
  # 1. Boxplots
  geom_boxplot(aes(fill = condition), width = 0.5, outlier.shape = NA, alpha = 0.6) +
  
  # 2. Connecting Lines
  geom_line(data = td_raw %>% filter(pairing_status == "Paired"),
            aes(group = patient), color = "grey50", alpha = 0.4, linewidth = 0.3) +
  
  # 3. Raw Points
  geom_point(aes(fill = condition, shape = pairing_status),
             color = "black", stroke = 0.2, size = 2.5, alpha = 0.9) +
  
  # 4. Significance Brackets (Variable tips reaching toward data)
  # Horizontal bar
  geom_segment(data = bracket_df, aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y)) +
  # Left tip (reaches down to Primary data)
  geom_segment(data = bracket_df, aes(x = 1, xend = 1, y = bracket_y, yend = tip_primary_y)) +
  # Right tip (reaches down to Metastasis data)
  geom_segment(data = bracket_df, aes(x = 2, xend = 2, y = bracket_y, yend = tip_met_y)) +
  # P-value text
  geom_text(data = bracket_df, aes(x = 1.5, y = label_y, label = label), 
            size = 3, vjust = 0) +
  
  # Faceting & Scales (Crucial for Tumour Purity vs Scores)
  facet_wrap(~cell_type, scales = "free_y") +
  scale_fill_manual(values = c("Primary" = col_primary, "Metastasis" = col_met)) +
  scale_shape_manual(values = c("Paired" = 21, "Orphan" = 22)) +
  
  # Aesthetics
  guides(shape = "none") +
  theme_prism(base_size = 12) +
  labs(title = "ESTIMATE", y = "Score / Purity", x = "") +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

print(estimate_plot)

ggsave("~/Desktop/Fig 1A.pdf", plot = estimate_plot, width = 180, height = 90, units = "mm", bg = "white", dpi = 300)


######## Attach tumor purity & stroma scores from estimate to metadata and save

# res_estimate: rows = c("stroma score","immune score","estimate score","tumor purity"), cols = sample IDs
scores_tbl = res_estimate %>%
  filter(cell_type %in% c("tumor purity", "stroma score")) %>%
  mutate(cell_type = recode(cell_type, "tumor purity" = "estimate_tumor_purity", "stroma score" = "estimate_stroma_score")) %>%
  pivot_longer(-cell_type, names_to = "sample", values_to = "value") %>%
  mutate(value = suppressWarnings(as.numeric(value))) %>%
  pivot_wider(names_from = cell_type, values_from = value)

# Join onto metadata (rownames are sample IDs)
metadata_updated = metadata %>%
  rownames_to_column("sample") %>%
  left_join(scores_tbl, by = "sample") %>%
  mutate(patient = factor(patient), condition = factor(condition, levels = c("T","M"))) %>%
  column_to_rownames("sample")

write.csv(metadata_updated, "patient_data_estimate.csv")

