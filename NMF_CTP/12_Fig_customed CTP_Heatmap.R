# ============================================================
# Plot CTP discovery results in the style of the reference study
#
# Required outputs from the completed CTP discovery pipeline:
#   output/CTP_analysis/rds/merged_A_matrices.rds
#   output/CTP_analysis/rds/factor_pearson_correlation_matrix.rds
#   output/CTP_analysis/rds/Final_igraph_object.rds
#   output/CTP_analysis/rds/ctp_aggweights_matrix.rds
#   output/CTP_analysis/tables/nodes_annotated.csv
#
# Outputs:
#   Fig_S1B_real_factor_correlation_heatmap.pdf
#   Fig_S1C_network_by_dataset.pdf
#   Fig_S1C_network_by_community.pdf
#   Fig_S1D_community_factor_counts.pdf
#   Fig_S1E_CTP_marker_weight_heatmap.pdf
# ============================================================

setwd("D:/HCC_ICB/cuda_nmf_results")

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
required_packages <- c(
  "tidyverse",
  "igraph",
  "tidygraph",
  "ggraph",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer",
  "grid"
)

to_install <- setdiff(required_packages, rownames(installed.packages()))

if (length(to_install) > 0) {
  install.packages(to_install)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

set.seed(123)
# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
analysis_dir <- file.path("output", "CTP_analysis")
rds_dir      <- file.path(analysis_dir, "rds")
table_dir    <- file.path(analysis_dir, "tables")
plot_dir     <- file.path(analysis_dir, "reference_style_plots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  file.path(rds_dir, "merged_A_matrices.rds"),
  file.path(rds_dir, "factor_pearson_correlation_matrix.rds"),
  file.path(rds_dir, "Final_igraph_object.rds"),
  file.path(rds_dir, "ctp_aggweights_matrix.rds"),
  file.path(table_dir, "nodes_annotated.csv")
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}

MP          <- readRDS(file.path(rds_dir, "merged_A_matrices.rds"))
cor_matrix  <- readRDS(file.path(rds_dir, "factor_pearson_correlation_matrix.rds"))
graph       <- readRDS(file.path(rds_dir, "Final_igraph_object.rds"))
agg_weights <- readRDS(file.path(rds_dir, "ctp_aggweights_matrix.rds"))
nodes_csv   <- read_csv(
  file.path(table_dir, "nodes_annotated.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Consistent dataset colors
# ------------------------------------------------------------
datasets <- sort(unique(nodes_csv$dataset))

base_dataset_colors <- c(
  "#D73027",
  "#4575B4",
  "#1A9850",
  "#984EA3",
  "#FF7F00",
  "#A65628",
  "#F781BF",
  "#999999"
)

dataset_colors <- setNames(
  rep(base_dataset_colors, length.out = length(datasets)),
  datasets
)

dataset_shapes <- setNames(
  rep(c(16, 17, 15, 18, 8, 3, 7, 4), length.out = length(datasets)),
  datasets
)

message("Datasets and colors:")
print(dataset_colors)


# ============================================================
# Figure S1B-like: all-factor Pearson correlation heatmap
# ============================================================

# Use one common dendrogram for rows and columns.
hc <- hclust(
  as.dist(1 - cor_matrix),
  method = "average"
)

factor_dataset <- str_remove(
  colnames(cor_matrix),
  "_k[0-9]+_P[0-9]+$"
)

factor_k <- as.integer(
  str_match(
    colnames(cor_matrix),
    "_k([0-9]+)_P[0-9]+$"
  )[, 2]
)

factor_annotation <- data.frame(
  Dataset = factor_dataset,
  K = factor_k,
  row.names = colnames(cor_matrix)
)

ha_factor <- HeatmapAnnotation(
  Dataset = factor_annotation$Dataset,
  col = list(Dataset = dataset_colors),
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(2.5, "mm")
)

ra_factor <- rowAnnotation(
  Dataset = factor_annotation$Dataset,
  col = list(Dataset = dataset_colors),
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(2.5, "mm")
)

cor_col_fun <- colorRamp2(
  c(-0.3, 0, 0.3, 0.5, 0.75, 1),
  c(
    "#313695",
    "#FFFFFF",
    "#FEE8C8",
    "#FDBB84",
    "#E34A33",
    "#7F0000"
  )
)

ht_cor <- Heatmap(
  cor_matrix,
  name = "Pearson r",
  col = cor_col_fun,
  cluster_rows = as.dendrogram(hc),
  cluster_columns = as.dendrogram(hc),
  top_annotation = ha_factor,
  left_annotation = ra_factor,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  border = FALSE,
  use_raster = TRUE,
  raster_quality = 4,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(35, "mm")
  )
)

pdf(
  file.path(plot_dir, "Fig_S1B_real_factor_correlation_heatmap.pdf"),
  width = 6.2,
  height = 6
)

draw(
  ht_cor,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()


# ============================================================
# Prepare graph attributes
# ============================================================

# Some reference code expects label_kp. Create it when absent.
if (!"label_kp" %in% vertex_attr_names(graph)) {
  V(graph)$label_kp <- paste0(
    V(graph)$dataset,
    "\nk",
    V(graph)$k,
    "-P",
    V(graph)$pattern
  )
}

# Ensure CTP_strict is available as a graph attribute.
if (!"CTP_strict" %in% vertex_attr_names(graph)) {
  strict_map <- setNames(nodes_csv$CTP_strict, nodes_csv$name)
  V(graph)$CTP_strict <- unname(strict_map[V(graph)$name])
}

# Store one deterministic FR layout and reuse it.
set.seed(123)
layout_full <- create_layout(
  as_tbl_graph(graph),
  layout = "fr",
  weights = E(graph)$cor,
  niter = 2000
)


# ============================================================
# Figure S1C-like panel 1: network colored by dataset
# ============================================================

p_network_dataset <- ggraph(layout_full) +
  geom_edge_link(
    aes(
      edge_alpha = cor,
      edge_width = cor
    ),
    colour = "grey70",
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(
      colour = dataset,
      shape = dataset,
      size = centrality_degree
    ),
    alpha = 0.9,
    stroke = 0.25
  ) +
  scale_colour_manual(values = dataset_colors) +
  scale_shape_manual(values = dataset_shapes) +
  scale_size_continuous(range = c(1.2, 5.5)) +
  scale_edge_alpha(range = c(0.15, 0.75)) +
  scale_edge_width(range = c(0.15, 1.25)) +
  coord_equal() +
  theme_void(base_size = 10) +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    colour = "Dataset",
    shape = "Dataset",
    size = "Degree"
  )

ggsave(
  file.path(plot_dir, "Fig_S1C_network_by_dataset.pdf"),
  p_network_dataset,
  width = 6,
  height = 5.5,
  device = cairo_pdf
)

# ============================================================
# Figure S1C-like panel 2: candidate CTP network by community
# ============================================================

graph_ctp <- as_tbl_graph(graph) %>%
  activate(nodes) %>%
  filter(CTP_strict == "yes")

if (gorder(graph_ctp) == 0) {
  stop("No CTP_strict == 'yes' nodes found in the graph.")
}

ctp_communities <- sort(
  unique(
    as.integer(
      as_tibble(graph_ctp, active = "nodes")$comm_infomap
    )
  )
)

n_communities <- length(ctp_communities)

community_palette <- if (n_communities <= 12) {
  brewer.pal(max(3, n_communities), "Paired")[seq_len(n_communities)]
} else {
  colorRampPalette(brewer.pal(12, "Paired"))(n_communities)
}

community_colors <- setNames(
  community_palette,
  as.character(ctp_communities)
)

set.seed(123)
layout_ctp <- create_layout(
  graph_ctp,
  layout = "fr",
  weights = edge_attr(as.igraph(graph_ctp), "cor"),
  niter = 2000
)

p_network_ctp <- ggraph(layout_ctp) +
  geom_edge_link(
    aes(
      edge_alpha = cor,
      edge_width = cor
    ),
    colour = "grey72",
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(
      colour = factor(comm_infomap),
      shape = dataset,
      size = centrality_degree
    ),
    alpha = 0.92,
    stroke = 0.25
  ) +
  scale_colour_manual(values = community_colors) +
  scale_shape_manual(values = dataset_shapes) +
  scale_size_continuous(range = c(1.5, 6)) +
  scale_edge_alpha(range = c(0.2, 0.85)) +
  scale_edge_width(range = c(0.2, 1.4)) +
  coord_equal() +
  theme_void(base_size = 10) +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    colour = "Infomap\ncommunity",
    shape = "Dataset",
    size = "Degree"
  )

ggsave(
  file.path(plot_dir, "c"),
  p_network_ctp,
  width = 6.5,
  height = 5.8,
  device = cairo_pdf
)


# ============================================================
# Figure S1D-like: factor count per community by dataset
# ============================================================

community_counts <- nodes_csv %>%
  mutate(
    community_label = paste0("C", comm_infomap)
  ) %>%
  count(
    comm_infomap,
    community_label,
    dataset,
    CTP_strict,
    name = "n_factors"
  ) %>%
  group_by(
    comm_infomap,
    community_label
  ) %>%
  mutate(
    total_factors = sum(n_factors),
    is_candidate_ctp = any(CTP_strict == "yes")
  ) %>%
  ungroup()

community_levels <- community_counts %>%
  distinct(
    comm_infomap,
    community_label,
    total_factors,
    is_candidate_ctp
  ) %>%
  arrange(
    desc(is_candidate_ctp),
    desc(total_factors),
    comm_infomap
  ) %>%
  pull(community_label)

community_counts <- community_counts %>%
  mutate(
    community_label = factor(
      community_label,
      levels = community_levels
    )
  )

community_totals <- community_counts %>%
  group_by(
    community_label
  ) %>%
  summarise(
    total_factors = first(total_factors),
    is_candidate_ctp = first(is_candidate_ctp),
    .groups = "drop"
  )

p_community_counts <- ggplot(
  community_counts,
  aes(
    x = community_label,
    y = n_factors,
    fill = dataset
  )
) +
  geom_col(
    width = 0.78,
    colour = "black",
    linewidth = 0.15
  ) +
  geom_text(
    data = filter(community_totals, is_candidate_ctp),
    aes(
      x = community_label,
      y = total_factors + max(total_factors) * 0.035,
      label = "*"
    ),
    inherit.aes = FALSE,
    size = 5,
    vjust = 0
  ) +
  scale_fill_manual(values = dataset_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.13))
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "right"
  ) +
  labs(
    x = "Infomap community",
    y = "Number of NMF factors",
    fill = "Dataset",
    caption = "* Multi-dataset candidate CTP"
  )

ggsave(
  file.path(plot_dir, "Fig_S1D_community_factor_counts.pdf"),
  p_community_counts,
  width = max(7, length(community_levels) * 0.34),
  height = 3.5,
  device = cairo_pdf
)


# ============================================================
# Figure S1E-like: top marker-gene aggregate-weight heatmap
# ============================================================

agg_weights <- as.matrix(agg_weights)
storage.mode(agg_weights) <- "double"

# Use the same candidate-community ordering as in the count plot.
candidate_order <- nodes_csv %>%
  filter(CTP_strict == "yes") %>%
  distinct(community_id, comm_infomap) %>%
  mutate(
    plot_order = match(
      paste0("C", comm_infomap),
      community_levels
    )
  ) %>%
  arrange(plot_order) %>%
  pull(community_id)

candidate_order <- intersect(
  candidate_order,
  colnames(agg_weights)
)

if (length(candidate_order) == 0) {
  stop(
    "No candidate community IDs match columns of agg_weights. ",
    "Check colnames(agg_weights) and nodes_csv$community_id."
  )
}

agg_plot <- agg_weights[, candidate_order, drop = FALSE]

n_top_per_ctp <- 100

top_gene_list <- vector("list", length(candidate_order))
names(top_gene_list) <- candidate_order

for (community in candidate_order) {
  ranked_genes <- rownames(agg_plot)[
    order(
      agg_plot[, community],
      decreasing = TRUE,
      na.last = NA
    )
  ]
  
  top_gene_list[[community]] <- head(
    unique(ranked_genes),
    n_top_per_ctp
  )
}

# Preserve repeated genes across different CTP blocks, as in the reference.
heatmap_gene_vector <- unlist(
  top_gene_list,
  use.names = FALSE
)

row_split_vector <- rep(
  names(top_gene_list),
  lengths(top_gene_list)
)

marker_mat <- agg_plot[
  heatmap_gene_vector,
  candidate_order,
  drop = FALSE
]

# Label a small, objective set of genes:
# top five genes from each community.
genes_to_label <- unique(
  unlist(
    lapply(top_gene_list, head, 5),
    use.names = FALSE
  )
)

label_positions <- which(
  rownames(marker_mat) %in% genes_to_label
)

label_names <- rownames(marker_mat)[label_positions]

# Repeated labels are shown in dark red.
duplicated_label_names <- names(
  table(label_names)[table(label_names) > 1]
)

label_colors <- ifelse(
  label_names %in% duplicated_label_names,
  "#B2182B",
  "black"
)

left_mark <- rowAnnotation(
  Gene = anno_mark(
    at = label_positions,
    labels = label_names,
    side = "left",
    labels_gp = gpar(
      fontsize = 5.5,
      col = label_colors
    ),
    link_width = unit(3, "mm"),
    padding = 0.5
  )
)

row_split_factor <- factor(
  row_split_vector,
  levels = candidate_order
)

positive_values <- marker_mat[
  is.finite(marker_mat) &
    marker_mat > 0
]

if (length(positive_values) == 0) {
  stop("No positive finite values found in marker_mat.")
}

weight_breaks <- unique(
  as.numeric(
    quantile(
      positive_values,
      probs = c(0, 0.5, 0.9, 0.99, 1),
      na.rm = TRUE
    )
  )
)

# Guard against duplicated quantile cutoffs.
if (length(weight_breaks) < 3) {
  weight_breaks <- seq(
    min(positive_values),
    max(positive_values),
    length.out = 5
  )
}

weight_colors <- colorRampPalette(
  c(
    "#FFFFFF",
    "#FFF7BC",
    "#FEC44F",
    "#EC7014",
    "#8C2D04"
  )
)(length(weight_breaks))

weight_col_fun <- colorRamp2(
  weight_breaks,
  weight_colors
)

ht_markers <- Heatmap(
  marker_mat,
  name = "Aggregate\nweight",
  col = weight_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 7),
  row_split = row_split_factor,
  row_title = NULL,
  row_gap = unit(0.6, "mm"),
  left_annotation = left_mark,
  border = TRUE,
  use_raster = TRUE,
  raster_quality = 4,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 7)
  )
)

marker_height <- max(
  6,
  0.024 * nrow(marker_mat)
)

pdf(
  file.path(plot_dir, "Fig_S1E_CTP_marker_weight_heatmap.pdf"),
  width = max(5, 0.35 * ncol(marker_mat) + 3),
  height = marker_height
)

draw(
  ht_markers,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()


# ============================================================
# 7. Optional: export figure metadata
# ============================================================

write_csv(
  community_counts,
  file.path(plot_dir, "Fig_S1D_community_factor_counts_data.csv")
)

marker_gene_table <- imap_dfr(
  top_gene_list,
  function(genes, community) {
    tibble(
      community_id = community,
      rank = seq_along(genes),
      gene = genes,
      aggregate_weight = agg_plot[genes, community]
    )
  }
)

write_csv(
  marker_gene_table,
  file.path(plot_dir, "Fig_S1E_top_marker_genes.csv")
)

message(
  "\nFinished. Figures were written to:\n",
  normalizePath(plot_dir, mustWork = FALSE)
)
