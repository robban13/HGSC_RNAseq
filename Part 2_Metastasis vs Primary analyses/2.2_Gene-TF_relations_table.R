library(dplyr)
library(tibble)
library(clusterProfiler)
library(msigdbr)
library(stringr)

short_label <- function(x, max_chars = 80) {
  x <- gsub("^GOMF_", "", x)
  x <- gsub("_", " ", x)
  x <- stringr::str_to_title(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "…"), x)
} #Helper function to tidy GO labels

msig_gomf     <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF")
term2gene_gomf <- msig_gomf %>% dplyr::select(term = gs_name, gene = gene_symbol)
term2name_gomf <- msig_gomf %>% dplyr::distinct(term = gs_name) %>% dplyr::mutate(name = short_label(term), name = make.unique(name))

## 1. parameters

key_tfs  <- tfs_in_res               # all TFs you had in DE results / Collectri
progeny_sources_keep <- c("TGFb")    # Progeny sets that we want to keep based on significance
hallmark_labels_keep <- c("EMT", "Kras Signaling Up", "Allograft Rejection", "Fatty Acid Metabolism", "E2f Targets", "Adipogenesis")

tf_padj_thr      <- 0.1   # TF transcript padj cutoff (NOTE that even if we set padj = 1 to include all TFs whose activity is significant, we will still filter out Tfs whose activity and expression disagree!)
# see tf_use_strong <- tf_info %>% section!!
tf_lfc_thr       <- 0     # TF transcript |log2FC| cutoff
target_padj_thr  <- 0.05  # gene padj cutoff
target_lfc_thr   <- 0   # gene |log2FC| cutoff from DE results
tf_activity_padj_thr <- 1 # TF activity padj cutoff 

## 2. DE table -----------------------------------------------------------------

deg_df <- as.data.frame(resLFC) %>% rownames_to_column("gene")

## 3. TF mRNA DE info ----------------------------------------------------------

tf_de_table <- deg_df %>%
  filter(gene %in% key_tfs) %>%
  mutate(
    is_TF_DE      = padj < tf_padj_thr & abs(log2FoldChange) >= tf_lfc_thr,
    tf_mrna_group = if_else(log2FoldChange > 0, "Metastasis", "Tumor"),
    tf_mrna_sign  = sign(log2FoldChange))

## 3b. TF activity info (decoupleR CollecTRI) ----------------------------------

tf_activity_df <- tf_contrast %>%
  filter(TF %in% key_tfs) %>%
  mutate(
    tf_group  = direction,                   # "Metastasis" or "Tumor"
    tf_active = padj < tf_activity_padj_thr ) # check so the TF activity passes the threshold for significance 


tf_info <- tf_de_table %>%
  select(TF = gene, is_TF_DE, tf_mrna_group, tf_mrna_sign) %>%
  inner_join(tf_activity_df, by = "TF") %>%
  mutate(
    agree_mrna_activity = tf_mrna_group == tf_group,
    tf_dir_sign         = if_else(tf_group == "Metastasis", 1L, -1L))

# Compact sanity table: TF transcript vs activity
tf_diag <- tf_de_table %>%
  select(TF = gene, tf_mrna_log2FC = log2FoldChange, tf_mrna_padj = padj,
         is_TF_DE, tf_mrna_group) %>%
  inner_join(tf_activity_df %>%
               select(TF, tf_activity_score = score, tf_activity_padj = padj,  tf_group, tf_active), by = "TF") %>%
  mutate(
    agree_mrna_activity     = tf_mrna_group == tf_group,
    mrna_activity_agreement = if_else(agree_mrna_activity, "Agree", "Disagree"),
    tf_DE_flag              = if_else(is_TF_DE, "Pass", "Fail"),
    agree_and_pass          = if_else(agree_mrna_activity & is_TF_DE, "Yes", "No"))

cat("\nTF transcript vs activity (all TFs):\n")
print(
  tf_diag %>%
    arrange(TF) %>%
    select(TF, tf_mrna_log2FC, tf_mrna_padj, tf_mrna_group, tf_activity_score, tf_activity_padj, tf_group, tf_DE_flag, tf_active, mrna_activity_agreement, agree_and_pass))

## pick final TFs used downstream 
tf_use_strong <- tf_info %>%
  filter(is_TF_DE, tf_active, agree_mrna_activity) %>%
  pull(TF)

tf_use_strong

## 4. TF → target evidence (CollecTRI) ----------------------------------------

tf_target_DE <- collectri_sets %>%
  filter(source %in% tf_use_strong) %>%
  rename(TF = source, target = target) %>%
  inner_join(deg_df, by = c("target" = "gene")) %>%
  left_join(tf_info %>% select(TF, tf_group, tf_dir_sign), by = "TF") %>%
  filter(!is.na(padj), padj < target_padj_thr, abs(log2FoldChange) >= target_lfc_thr) %>%
  mutate(concordant = sign(mor * tf_dir_sign * log2FoldChange) > 0)

tf_gene_scores <- tf_target_DE %>%
  group_by(target) %>%
  summarise(
    n_TFs_any        = n_distinct(TF),
    n_TFs_concordant = n_distinct(TF[concordant]),
    TFs_concordant   = paste(sort(unique(TF[concordant])), collapse = ", "),
    n_TFs_discordant = n_distinct(TF[!concordant]),
    TFs_discordant   = paste(sort(unique(TF[!concordant])), collapse = ", "),
    min_padj         = min(padj),
    .groups = "drop") %>%
  mutate(
    TFs_concordant = na_if(TFs_concordant, ""),
    TFs_discordant = na_if(TFs_discordant, ""))

## 5. PROGENy evidence --------------------------------------------------------

progeny_res_sub <- progeny_contrast %>%
  filter(source %in% progeny_sources_keep) %>%
  transmute(
    Progeny       = source,
    progeny_score = score,
    progeny_padj  = padj,
    progeny_dir   = direction,
    path_sign     = sign(score))

progeny_gene_df <- progeny_sets %>%
  filter(source %in% progeny_sources_keep) %>%
  rename(Progeny = source) %>%
  inner_join(deg_df, by = c("target" = "gene")) %>%
  left_join(progeny_res_sub, by = "Progeny") %>%
  mutate(
    progeny_is_DE    = padj < target_padj_thr & abs(log2FoldChange) >= target_lfc_thr,
    progeny_conc_flag = progeny_is_DE & !is.na(path_sign) &
      sign(weight * log2FoldChange * path_sign) > 0)

progeny_per_gene <- progeny_gene_df %>%
  filter(progeny_is_DE) %>%
  group_by(target) %>%
  summarise(
    progeny_selected     = paste(sort(unique(Progeny)), collapse = "; "),
    n_progeny_any        = n_distinct(Progeny),
    n_progeny_concordant = n_distinct(Progeny[progeny_conc_flag]),
    n_progeny_discordant = n_distinct(Progeny[progeny_is_DE & !progeny_conc_flag]),
    progeny_concordant   = paste(sort(unique(Progeny[progeny_conc_flag])), collapse = ", "),
    progeny_discordant   = paste(sort(unique(Progeny[progeny_is_DE & !progeny_conc_flag])), collapse = ", "),
    .groups = "drop" ) %>%
  mutate(
    progeny_concordant = na_if(progeny_concordant, ""),
    progeny_discordant = na_if(progeny_discordant, ""))

## 6. Hallmark membership -----------------------------------------------------

hallmarks_keep <- cp %>%
  filter(label %in% hallmark_labels_keep) %>%
  pull(pathway) %>%
  unique()

term2gene_sub <- term2gene %>% filter(term %in% hallmarks_keep)

hallmarks_per_gene <- term2gene_sub %>%
  left_join(cp %>% select(pathway, label), by = c("term" = "pathway")) %>%
  transmute(target = gene, hallmark_label = label) %>%
  group_by(target) %>%
  summarise(hallmarks_selected = paste(sort(unique(hallmark_label)), collapse = "; "), .groups = "drop")

## 7. GO:MF membership --------------------------------------------------------

go_mf_per_gene <- term2gene_gomf %>%
  left_join(term2name_gomf, by = "term") %>%
  group_by(gene) %>%
  summarise(go_mf_names = paste(sort(unique(name)), collapse = "; "), .groups = "drop")

## 8. Gene-level DE direction -------------------------------------------------

dir_df <- deg_df %>%
  transmute(
    target      = gene,
    gene_group  = if_else(log2FoldChange > 0, "Metastasis", "Tumor"),
    gene_log2FC = log2FoldChange)

## 9. Final integration -------------------------------------------------------

full_gene_summary <- tf_gene_scores %>%
  full_join(progeny_per_gene,   by = "target") %>%
  left_join(dir_df,             by = "target") %>%
  left_join(hallmarks_per_gene, by = "target") %>%
  left_join(go_mf_per_gene,     by = c("target" = "gene")) %>%
  mutate(
    n_TFs_any            = ifelse(is.na(n_TFs_any),            0L, n_TFs_any),
    n_TFs_concordant     = ifelse(is.na(n_TFs_concordant),     0L, n_TFs_concordant),
    n_TFs_discordant     = ifelse(is.na(n_TFs_discordant),     0L, n_TFs_discordant),
    n_progeny_any        = ifelse(is.na(n_progeny_any),        0L, n_progeny_any),
    n_progeny_concordant = ifelse(is.na(n_progeny_concordant), 0L, n_progeny_concordant),
    n_progeny_discordant = ifelse(is.na(n_progeny_discordant), 0L, n_progeny_discordant)) %>%
  select(
    target,
    hallmarks_selected,
    go_mf_names,
    gene_group, gene_log2FC, min_padj,
    n_TFs_any, n_TFs_concordant, TFs_concordant,
    n_TFs_discordant, TFs_discordant,
    n_progeny_any, n_progeny_concordant, progeny_concordant,
    n_progeny_discordant, progeny_discordant) %>%
  arrange(gene_group, desc(n_TFs_concordant), desc(n_progeny_concordant), min_padj)

view(full_gene_summary)
