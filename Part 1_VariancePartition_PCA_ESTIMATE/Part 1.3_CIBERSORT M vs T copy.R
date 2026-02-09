library(tidyverse)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(rstatix)

fix_ids <- function(x) gsub("^X", "", x) #Sample names have an X before any zero otherwise it breaks cibersort, remove this X downstream.

# ==============================================================================
# 0) USER SETTINGS
# ==============================================================================
P_VAL_DISPLAY_THRESHOLD <- 0.2  
SELECTED_CELL_TYPES <- c("B cells naive", "Macrophages M2", "T cells CD8") #Subset to only plot these ones for report. 

# ==============================================================================
# 1) Inputs & Data Prep
# ==============================================================================
metadata <- read.csv("patient_data_estimate.csv", header = TRUE, row.names = 1, check.names = FALSE) %>%
  mutate(sample = fix_ids(rownames(.))) %>%
  dplyr::select(sample, patient, condition) 

cib <- read.csv("CIBERSORTx_Job5_Results_absolute_scaled.csv", header = TRUE, row.names = 1, check.names = FALSE) %>%
  rownames_to_column("sample") %>%
  mutate(sample = fix_ids(sample))

cell_cols <- setdiff(colnames(cib), c("sample","P-value","Correlation","RMSE","Absolute score (sig.score)"))

# --- RAW DATA (Renaming Tumor -> Primary) ---
td_raw <- cib %>%
  pivot_longer(all_of(cell_cols), names_to = "cell_type", values_to = "score") %>%
  mutate(score = suppressWarnings(as.numeric(score))) %>%
  left_join(metadata, by = "sample") %>%
  filter(!is.na(score)) %>%
  filter(condition %in% c("T", "M")) %>%
  mutate(condition = factor(condition, levels = c("T","M"), labels = c("Primary","Metastasis")))

if (!is.null(SELECTED_CELL_TYPES) && length(SELECTED_CELL_TYPES) > 0) {
  td_raw <- td_raw %>% filter(cell_type %in% SELECTED_CELL_TYPES)
}

# --- IDENTIFY ORPHANS VS PAIRS ---
valid_patients <- td_raw %>%
  group_by(patient) %>%
  summarise(
    has_primary = any(condition == "Primary"),
    has_met     = any(condition == "Metastasis")
  )

td_raw <- td_raw %>%
  left_join(valid_patients, by = "patient") %>%
  mutate(pairing_status = ifelse(has_primary & has_met, "Paired", "Orphan"))

# --- AVERAGED DATA (For Stats & Connecting Lines) ---
td_averaged <- td_raw %>%
  group_by(cell_type, patient, condition, pairing_status) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

# ==============================================================================
# 2) Plotting Function
# ==============================================================================

plot_paired_comparison <- function(data_raw, data_avg) {
  
  # --- STEP A: DEFINE VALID CELL TYPES ---
  abund_threshold <- 0.01
  
  valid_cells <- data_raw %>%
    group_by(cell_type) %>%
    summarise(mean_score = mean(score, na.rm = TRUE)) %>%
    filter(mean_score > abund_threshold) %>%
    pull(cell_type)
  
  clean_raw <- data_raw %>% filter(cell_type %in% valid_cells)
  clean_avg <- data_avg %>% filter(cell_type %in% valid_cells)
  
  if(nrow(clean_raw) == 0) return(NULL)
  
  # --- STEP B: PREPARE GEOMETRIES ---
  
  # 1. DOTS (Use RAW Data)
  dots_data <- clean_raw %>%
    mutate(
      x_base = as.numeric(factor(cell_type, levels = levels(factor(clean_raw$cell_type)))),
      x_dodge = ifelse(condition == "Primary", -0.2, 0.2), 
      x_final = x_base + x_dodge 
    )
  
  # 2. LINES (Use AVERAGED Data)
  lines_data <- clean_avg %>%
    filter(pairing_status == "Paired") %>%
    select(patient, cell_type, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    mutate(
      x_base = as.numeric(factor(cell_type, levels = levels(factor(clean_raw$cell_type)))),
      x_primary_final = x_base - 0.2, 
      x_met_final     = x_base + 0.2
    )
  
  # 3. STATS (Use AVERAGED Data -> Valid Paired Test)
  stats_input <- clean_avg %>% filter(pairing_status == "Paired")
  
  box_maxes <- clean_raw %>%
    group_by(cell_type, condition) %>%
    summarise(max_val = max(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = max_val)
  
  p_values_df <- stats_input %>%
    group_by(cell_type) %>%
    summarise(
      p = tryCatch({
        pair_wide <- cur_data() %>% pivot_wider(names_from = condition, values_from = score)
        if(nrow(pair_wide) >= 3) {
          wilcox.test(pair_wide$Primary, pair_wide$Metastasis, paired = TRUE, exact = FALSE)$p.value
        } else { NA_real_ }
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    adjust_pvalue(method = "BH") %>% 
    add_significance("p")
  
  stats_df <- p_values_df %>%
    left_join(box_maxes, by = "cell_type") %>%
    filter(p <= P_VAL_DISPLAY_THRESHOLD) %>% 
    mutate(
      bracket_y = pmax(Primary, Metastasis, na.rm = TRUE) + (max(clean_raw$score, na.rm = TRUE) * 0.10),
      label = case_when(
        p < 0.001 ~ "p < 0.001",
        TRUE ~ paste0("p = ", sprintf("%.3f", p))
      ),
      x_center = as.numeric(factor(cell_type, levels = levels(factor(clean_raw$cell_type)))),
      x_primary = x_center - 0.2, 
      x_met     = x_center + 0.2
    )
  
  # --- COLORS ---
  floral_colors <- ggprism::prism_colour_pal(palette = "floral")(12)
  col_primary <- floral_colors[1] 
  col_met     <- floral_colors[7] 
  
  # --- PLOT ---
  p <- ggplot() +
    
    # 1. Boxplots (Background)
    geom_boxplot(data = clean_raw, aes(x = cell_type, y = score, fill = condition),
                 width = 0.7, position = position_dodge(width = 0.8), 
                 outlier.shape = NA, alpha = 0.6) +
    
    # 2. LINES (Averaged connections)
    geom_segment(data = lines_data,
                 aes(x = x_primary_final, xend = x_met_final, y = Primary, yend = Metastasis),
                 color = "grey50", alpha = 0.4, linewidth = 0.3) +
    
    # 3. POINTS (Raw Data)
    geom_point(data = dots_data,
               aes(x = x_final, y = score, fill = condition, shape = pairing_status),
               color = "black", stroke = 0.2, 
               size = 2, alpha = 0.9) +
    
    # 4. Stats Brackets
    geom_segment(data = stats_df, aes(x = x_primary, xend = x_met, y = bracket_y, yend = bracket_y)) +
    geom_segment(data = stats_df, aes(x = x_primary, xend = x_primary, y = bracket_y, yend = Primary + (max(clean_raw$score)*0.05))) +
    geom_segment(data = stats_df, aes(x = x_met, xend = x_met, y = bracket_y, yend = Metastasis + (max(clean_raw$score)*0.05))) +
    geom_text(data = stats_df, aes(x = x_center, y = bracket_y + (max(clean_raw$score)*0.02), label = label), vjust = 0, size = 3.5) +
    
    scale_fill_manual(values = c("Primary" = col_primary, "Metastasis" = col_met)) +
    scale_shape_manual(values = c("Paired" = 21, "Orphan" = 22)) +
    
    # >>> HIDE SHAPE LEGEND <<<
    guides(shape = "none") +
    
    theme_prism(base_size = 16) +
    labs(title = "CIBERSORTx",
         y = "Absolute Score", x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "top",
          plot.subtitle = element_text(face = "italic", size = 10, colour = "grey40"))
  
  return(p)
}

# ==============================================================================
# 4) Generate Plot
# ==============================================================================

cibersort_paired_plot <- plot_paired_comparison(td_raw, td_averaged)
print(cibersort_paired_plot)

ggsave("~/Desktop/Fig 2F.pdf", plot = cibersort_paired_plot, width = 180, height = 180, units = "mm", bg = "white", dpi = 300)
