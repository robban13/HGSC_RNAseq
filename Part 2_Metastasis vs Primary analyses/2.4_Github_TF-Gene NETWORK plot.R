library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(grid) 
library(scales) 

# Config section

show_all_targets <- F #Set to TRUE to show ALL target genes. Set to FALSE to only show the specific list below.

# Custom List: If show_all_targets is FALSE, only these genes will be labeled.
custom_gene_list <- c("CD36", "FABP4", "ADH1B", "LPL" ,#From FABP4 and CD34 publications
                      "IL8", "CXCR1", "CXCR2",        #IL8 and the receptor parts
                      "IL6", "IL6R" ,                 #IL6 and receptors
                      "CCL2", "CCR2", "CCR4",         #MCP-1 = CCL2 and receptors CCR2 CCR4
                      "TIMP1", "CD63",
                      "PLIN1", "ADIPOQ", "MT2A", "FASN",
                      "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "STAT3",
                      "SERPINE1", "MMP2", "VCAM1", "FN1",
                      "SGK1", "SCD", "GPX4", "ACSL1"
)  #


# Data preparation for network
# (Assuming target_data_subset, tf_use_strong, and resLFC exist)

network_edges <- target_data_subset %>%
  filter(TF %in% tf_use_strong) %>%
  rename(from = TF, to = target) %>%
  mutate(interaction_type = ifelse(log2FoldChange > 0, "Activated", "Repressed")) %>%
  dplyr::select(from, to, interaction_type, log2FoldChange)

tf_nodes     <- unique(network_edges$from)
target_nodes <- unique(network_edges$to)
all_nodes    <- unique(c(tf_nodes, target_nodes))

res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)

# Create Node Data with Smart Labeling Logic
nodes_df <- data.frame(name = all_nodes) %>%
  left_join(res_df %>% dplyr::select(gene, log2FoldChange), by = c("name" = "gene")) %>%
  mutate(
    node_type = ifelse(name %in% tf_nodes, "TF Driver", "Target Gene"),
    size      = ifelse(node_type == "TF Driver", 7, 2.5), 
    
    # Logic for TF Labels (Always Show)
    label_tf   = ifelse(node_type == "TF Driver", name, NA),
    
    # Logic for Target Gene Labels (Dependent on config)
    label_gene = case_when(
      node_type == "Target Gene" & show_all_targets == TRUE ~ name, # Show all if Toggle is ON
      node_type == "Target Gene" & name %in% custom_gene_list ~ name, # Show only match if Toggle is OFF
      TRUE ~ NA_character_ # Otherwise hide
    )
  )

# CREATE GRAPH
graph_obj <- tbl_graph(nodes = nodes_df, edges = network_edges) %>%
  mutate(degree = centrality_degree(mode = 'all'))

# Define Weights (Keep 5x weight for leaf nodes to keep structure clean)
E(graph_obj)$weight <- sapply(1:ecount(graph_obj), function(i) {
  ends <- ends(graph_obj, i, names = FALSE) 
  deg1 <- degree(graph_obj)[ends[1]]
  deg2 <- degree(graph_obj)[ends[2]]
  if(deg1 == 1 | deg2 == 1) return(5) else return(1)
})

set.seed(123)
coords <- layout_with_fr(graph_obj, weights = E(graph_obj)$weight, niter = 50000)
graph_obj <- graph_obj %>% mutate(x = coords[,1], y = coords[,2])


# PLOT 
p_balanced <- ggraph(graph_obj, layout = "manual", x = x, y = y) + 
  
  # A. Edges
  geom_edge_link(aes(filter = interaction_type == "Activated"),
                 arrow = arrow(length = unit(2.5, 'mm'), type = "closed"), 
                 end_cap = circle(2, 'mm'),  color = "grey60",  alpha = 0.6, width = 0.4) +
  geom_edge_link(aes(filter = interaction_type == "Repressed"),
                 arrow = arrow(length = unit(2.5, 'mm'), angle = 90, ends = "last"), 
                 end_cap = circle(2, 'mm'),  color = "grey60",  alpha = 0.6, width = 0.4) +
  
  # B. Nodes
  geom_node_point(aes(fill = log2FoldChange, size = size, shape = node_type), color = "black", stroke = 0.6) +
  
  # C.1 Target Labels
  geom_node_text(aes(label = label_gene), 
                 repel = TRUE,              # Activates ggrepel logic
                 size = 3,              
                 color = "black", 
                 bg.color = "white", 
                 bg.r = 0.1,                # Slight background radius to clear space around text
                 
                 # --- LINE SETTINGS ---
                 segment.color = "black",   # Color of the connector line
                 segment.size = 0.4,        # Thickness of the line
                 min.segment.length = 0,    # ALWAYS draw the line even if close
                 
                 # --- SPACING SETTINGS ---
                 force = 10,                 # Strength of push away from points
                 box.padding = 0.8,         # Distance between label and point
                 point.padding = 0,       # Distance between point and start of line
                 max.overlaps = Inf) +  
  
  # C.2 TF Labels
  geom_node_text(aes(label = label_tf), 
                 repel = TRUE, 
                 fontface = "bold", 
                 size = 3, 
                 bg.color = "white", 
                 bg.r = 0,
                 box.padding = 0.5) +
  
  # D. Scales and customisation
  scale_edge_color_manual(values = c("Activated" = "#D73027", "Repressed" = "#4575B4")) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limit = c(-2, 2), oob = scales::squish, name = "Log2FC") +
  scale_shape_manual(values = c("TF Driver" = 23, "Target Gene" = 21), guide = "none") + 
  scale_size_identity() +
  theme_graph(base_family = "Helvetica") +
  theme(legend.position = "right") +
  labs(title = "", subtitle = "")

print(p_balanced)

ggsave("~/Desktop/Fig2G.png", plot = p_balanced, width = 180, height = 180, units = "mm", bg = "white", dpi = 600)
