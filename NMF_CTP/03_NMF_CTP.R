# ============================================================
# CTP discovery pipeline
#
# Input:
#   gene × factor Ascaled matrix generated from Python
#
# Main steps:
#   1. Read merged Ascaled matrix
#   2. Parse factor information
#   3. Calculate factor-factor Pearson correlations
#   4. Construct correlation network
#   5. Infomap community detection
#   6. Identify cross-dataset communities
#   7. Aggregate gene weights for each community
#   8. Export top genes and diagnostic plots
# ============================================================
setwd("D:\\HCC_ICB\\cuda_nmf_results")
getwd()
[1] "D:/HCC_ICB/cuda_nmf_results"
# ============================================================
# 0. Install and load packages
# ============================================================

required_packages <- c(
  "tidyverse",
  "data.table",
  "Hmisc",
  "igraph",
  "tidygraph",
  "ggraph",
  "ggrepel",
  "patchwork",
  "pheatmap",
  "RColorBrewer"
)

packages_to_install <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Hmisc)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
})


# ============================================================
# 1. Paths and parameters
# ============================================================
getwd()
input_file <- file.path(
  "data",
  "all_groups_merged_Ascaled_K10_K25.tsv"
)

output_dir <- file.path(
  "output",
  "CTP_analysis"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
rds_dir <- file.path(output_dir, "rds")
topgene_dir <- file.path(output_dir, "top_genes")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(topgene_dir, recursive = TRUE, showWarnings = FALSE)


# Correlation threshold used by the original workflow
cor_threshold <- 0.50

# Infomap trials
infomap_trials <- 200

# Number of top genes retained per community
n_top_genes <- 100

# Rank offset used by the original comm_agg_weights function
rank_buffer <- 9

# The original script used FALSE.
# FALSE = each factor contributes equally.
# TRUE  = each dataset contributes approximately equally within a community.
weight_dataset_equal <- FALSE

# Minimum number of different datasets required for a candidate CTP
min_datasets_per_ctp <- 2

# Optional stricter filter:
# require at least this many factors in each represented dataset.
min_factors_per_dataset <- 1

set.seed(123)


# ============================================================
# 2. Read merged gene × factor matrix
# ============================================================

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    normalizePath(input_file, mustWork = FALSE)
  )
}

message("Reading matrix: ", input_file)

MP_df <- fread(
  input_file,
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(MP_df) < 2) {
  stop("The input file contains fewer than two columns.")
}

# The first column exported by pandas contains gene names
gene_column <- colnames(MP_df)[1]

message("Gene-name column detected as: ", gene_column)

gene_names <- as.character(MP_df[[gene_column]])

MP_df[[gene_column]] <- NULL

MP <- as.matrix(MP_df)
storage.mode(MP) <- "double"

rownames(MP) <- gene_names

colnames(MP) <- str_replace(
  colnames(MP),
  "_Pattern([0-9]+)$",
  "_P\\1"
)

head(colnames(MP), 20)

rm(MP_df)
gc()


# ============================================================
# 3. Basic input checks
# ============================================================

cat("\n")
cat("Merged Ascaled matrix dimensions\n")
cat("--------------------------------\n")
cat("Genes:   ", nrow(MP), "\n")
cat("Factors: ", ncol(MP), "\n\n")

expected_factors <- 3 * sum(10:25)

cat("Expected factors for 3 datasets and K=10:25: ",
    expected_factors, "\n")

if (ncol(MP) != expected_factors) {
  warning(
    "The factor number is ", ncol(MP),
    ", rather than the expected ", expected_factors,
    ". The pipeline can still run, but check whether all K values were included."
  )
}

if (anyDuplicated(rownames(MP)) > 0) {
  duplicated_genes <- unique(
    rownames(MP)[duplicated(rownames(MP))]
  )
  
  stop(
    "Duplicated gene names detected. Examples: ",
    paste(head(duplicated_genes, 10), collapse = ", ")
  )
}

if (anyDuplicated(colnames(MP)) > 0) {
  duplicated_factors <- unique(
    colnames(MP)[duplicated(colnames(MP))]
  )
  
  stop(
    "Duplicated factor names detected. Examples: ",
    paste(head(duplicated_factors, 10), collapse = ", ")
  )
}

if (anyNA(MP)) {
  stop(
    "NA values detected in the merged Ascaled matrix. ",
    "Please resolve missing values before calculating correlations."
  )
}

if (any(!is.finite(MP))) {
  stop("Infinite or non-finite values detected in the matrix.")
}

if (any(MP < 0)) {
  warning(
    "Negative values were detected. ",
    "Ascaled/NMF gene loadings are normally non-negative."
  )
}

zero_variance_factors <- apply(
  MP,
  2,
  function(x) sd(x) == 0
)

if (any(zero_variance_factors)) {
  warning(
    sum(zero_variance_factors),
    " zero-variance factors will be removed."
  )
  
  MP <- MP[, !zero_variance_factors, drop = FALSE]
}

saveRDS(
  MP,
  file.path(rds_dir, "merged_A_matrices.rds")
)


# ============================================================
# 4. Parse factor names
# ============================================================

# Expected factor-name format:
#   DATASET_k10_P1
#
# Dataset names may contain letters, numbers and underscores.
# Parsing is performed from the right side of the string.

factor_info <- tibble(
  factor = colnames(MP)
) %>%
  mutate(
    dataset = str_remove(
      factor,
      "_k[0-9]+_P[0-9]+$"
    ),
    k = as.integer(
      str_match(factor, "_k([0-9]+)_P[0-9]+$")[, 2]
    ),
    pattern = as.integer(
      str_match(factor, "_k[0-9]+_P([0-9]+)$")[, 2]
    )
  )

invalid_factor_names <- factor_info %>%
  filter(
    is.na(k) |
      is.na(pattern) |
      dataset == factor
  )

if (nrow(invalid_factor_names) > 0) {
  print(head(invalid_factor_names, 20))
  
  stop(
    "Some factor names do not follow DATASET_kK_PPATTERN format. ",
    "Please rename them before continuing."
  )
}

cat("\nDatasets detected:\n")
print(
  factor_info %>%
    count(dataset, name = "number_of_factors")
)

cat("\nK values detected:\n")
print(
  factor_info %>%
    count(dataset, k) %>%
    arrange(dataset, k)
)

write_csv(
  factor_info,
  file.path(table_dir, "factor_information.csv")
)


# ============================================================
# 5. Factor-level normalization check
# ============================================================

factor_qc <- tibble(
  factor = colnames(MP),
  minimum = apply(MP, 2, min),
  maximum = apply(MP, 2, max),
  mean = colMeans(MP),
  sd = apply(MP, 2, sd),
  sum = colSums(MP),
  nonzero_genes = colSums(MP > 0)
) %>%
  left_join(factor_info, by = "factor")

write_csv(
  factor_qc,
  file.path(table_dir, "factor_QC.csv")
)

p_factor_sum <- ggplot(
  factor_qc,
  aes(
    x = dataset,
    y = sum,
    fill = dataset
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.35,
    size = 0.8
  ) +
  theme_classic(base_size = 12) +
  guides(fill = "none") +
  labs(
    title = "Distribution of total Ascaled loading per factor",
    x = NULL,
    y = "Column sum"
  )

ggsave(
  file.path(plot_dir, "factor_loading_sum_by_dataset.pdf"),
  p_factor_sum,
  width = 7,
  height = 5
)


# ============================================================
# 6. Pearson correlation between factors
# ============================================================

message("Calculating pairwise Pearson correlations...")

# Rows are genes and columns are factors.
# Therefore this calculates factor-factor correlations.
cor_res <- Hmisc::rcorr(
  MP,
  type = "pearson"
)

cor_matrix <- cor_res$r
p_matrix <- cor_res$P
n_matrix <- cor_res$n

saveRDS(
  cor_matrix,
  file.path(rds_dir, "factor_pearson_correlation_matrix.rds")
)

fwrite(
  as.data.frame(cor_matrix, check.names = FALSE) %>%
    rownames_to_column("factor"),
  file.path(table_dir, "factor_pearson_correlation_matrix.tsv"),
  sep = "\t"
)


# ============================================================
# 7. Convert correlation matrix into an edge table
# ============================================================

upper_index <- which(
  upper.tri(cor_matrix),
  arr.ind = TRUE
)

edges_all <- tibble(
  pattern1 = rownames(cor_matrix)[upper_index[, 1]],
  pattern2 = colnames(cor_matrix)[upper_index[, 2]],
  cor = cor_matrix[upper_index],
  pval = p_matrix[upper_index],
  n_genes = n_matrix[upper_index]
) %>%
  left_join(
    factor_info %>%
      rename(
        pattern1 = factor,
        dataset1 = dataset,
        k1 = k,
        factor_number1 = pattern
      ),
    by = "pattern1"
  ) %>%
  left_join(
    factor_info %>%
      rename(
        pattern2 = factor,
        dataset2 = dataset,
        k2 = k,
        factor_number2 = pattern
      ),
    by = "pattern2"
  ) %>%
  mutate(
    cross_dataset = dataset1 != dataset2
  )

write_csv(
  edges_all,
  file.path(table_dir, "all_factor_pair_correlations.csv")
)

edges <- edges_all %>%
  filter(
    is.finite(cor),
    cor >= cor_threshold
  )

cat("\nCorrelation network\n")
cat("-------------------\n")
cat("Correlation threshold: ", cor_threshold, "\n")
cat("Retained edges:        ", nrow(edges), "\n")
cat(
  "Cross-dataset edges:  ",
  sum(edges$cross_dataset),
  "\n\n"
)

if (nrow(edges) == 0) {
  stop(
    "No edges remain at correlation threshold ",
    cor_threshold,
    ". Try lowering the threshold to 0.50 or 0.55."
  )
}

write_csv(
  edges,
  file.path(table_dir, "edges_cor_threshold.csv")
)


# ============================================================
# 8. Correlation diagnostics
# ============================================================

p_cor_distribution <- ggplot(
  edges_all,
  aes(
    x = cor,
    fill = cross_dataset
  )
) +
  geom_histogram(
    bins = 100,
    alpha = 0.65,
    position = "identity"
  ) +
  geom_vline(
    xintercept = cor_threshold,
    linetype = 2,
    linewidth = 0.8
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Distribution of factor-factor Pearson correlations",
    subtitle = paste0(
      "Dashed line: correlation threshold = ",
      cor_threshold
    ),
    x = "Pearson correlation",
    y = "Number of factor pairs",
    fill = "Cross-dataset"
  )

ggsave(
  file.path(plot_dir, "factor_correlation_distribution.pdf"),
  p_cor_distribution,
  width = 8,
  height = 5
)


# Heatmap for factors having at least one retained edge
connected_factors <- union(
  edges$pattern1,
  edges$pattern2
)

cor_connected <- cor_matrix[
  connected_factors,
  connected_factors,
  drop = FALSE
]

factor_annotation <- factor_info %>%
  filter(factor %in% connected_factors) %>%
  select(factor, dataset, k) %>%
  column_to_rownames("factor")

pdf(
  file.path(plot_dir, "connected_factor_correlation_heatmap.pdf"),
  width = 14,
  height = 13
)

pheatmap(
  cor_connected,
  annotation_row = factor_annotation,
  annotation_col = factor_annotation,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_method = "average",
  border_color = NA,
  breaks = seq(-1, 1, length.out = 101),
  main = paste0(
    "Connected NMF factors, Pearson r ≥ ",
    cor_threshold
  )
)

dev.off()


# ============================================================
# 9. Build igraph object
# ============================================================

vertices <- factor_info %>%
  transmute(
    name = factor,
    dataset = dataset,
    k = k,
    pattern = pattern
  )

graph <- graph_from_data_frame(
  d = edges %>%
    transmute(
      from = pattern1,
      to = pattern2,
      cor = cor,
      pval = pval,
      cross_dataset = cross_dataset
    ),
  directed = FALSE,
  vertices = vertices
)

# Remove isolated nodes because they do not belong to a correlation community
graph <- delete_vertices(
  graph,
  which(degree(graph) == 0)
)

cat("Graph nodes: ", vcount(graph), "\n")
cat("Graph edges: ", ecount(graph), "\n")


# ============================================================
# 10. Infomap community detection
# ============================================================

set.seed(123)

infomap_result <- cluster_infomap(
  graph,
  e.weights = E(graph)$cor,
  nb.trials = infomap_trials
)

V(graph)$comm_infomap <- membership(infomap_result)

cat(
  "Infomap communities: ",
  length(unique(V(graph)$comm_infomap)),
  "\n"
)


# ============================================================
# 11. Node centralities and CTP annotation
# ============================================================

V(graph)$centrality_degree <- degree(
  graph,
  mode = "all",
  normalized = FALSE
)

V(graph)$centrality_strength <- strength(
  graph,
  mode = "all",
  weights = E(graph)$cor
)

V(graph)$centrality_betweenness <- betweenness(
  graph,
  directed = FALSE,
  weights = 1 / pmax(E(graph)$cor, 1e-8),
  normalized = TRUE
)

V(graph)$centrality_closeness <- closeness(
  graph,
  mode = "all",
  weights = 1 / pmax(E(graph)$cor, 1e-8),
  normalized = TRUE
)

nodes <- as_data_frame(
  graph,
  what = "vertices"
) %>%
  as_tibble() %>%
  group_by(comm_infomap) %>%
  mutate(
    community_size = n(),
    community_n_datasets = n_distinct(dataset),
    CTP = if_else(
      community_n_datasets >= min_datasets_per_ctp,
      "yes",
      "no"
    )
  ) %>%
  ungroup() %>%
  mutate(
    community_id = paste0("CTP_", comm_infomap)
  )


# Check minimum contribution per represented dataset
community_dataset_counts <- nodes %>%
  count(
    comm_infomap,
    community_id,
    dataset,
    name = "n_factors"
  )

community_min_dataset_count <- community_dataset_counts %>%
  group_by(
    comm_infomap,
    community_id
  ) %>%
  summarise(
    minimum_factors_in_dataset = min(n_factors),
    .groups = "drop"
  )

nodes <- nodes %>%
  left_join(
    community_min_dataset_count,
    by = c("comm_infomap", "community_id")
  ) %>%
  mutate(
    CTP_strict = if_else(
      CTP == "yes" &
        minimum_factors_in_dataset >= min_factors_per_dataset,
      "yes",
      "no"
    )
  )

write_csv(
  nodes,
  file.path(table_dir, "nodes_annotated.csv")
)


# ============================================================
# 12. Community summary
# ============================================================

community_summary <- nodes %>%
  count(
    comm_infomap,
    community_id,
    dataset,
    name = "n_factors"
  ) %>%
  group_by(
    comm_infomap,
    community_id
  ) %>%
  mutate(
    total_factors = sum(n_factors),
    fraction = n_factors / total_factors
  ) %>%
  summarise(
    total_factors = first(total_factors),
    n_datasets = n_distinct(dataset),
    largest_dataset_fraction = max(fraction),
    smallest_dataset_count = min(n_factors),
    dataset_composition = paste0(
      dataset,
      ":",
      n_factors,
      collapse = "; "
    ),
    CTP = if_else(
      n_datasets >= min_datasets_per_ctp,
      "yes",
      "no"
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(CTP),
    desc(n_datasets),
    desc(total_factors)
  )

write_csv(
  community_summary,
  file.path(table_dir, "community_summary.csv")
)

print(community_summary, n = Inf)


# ============================================================
# 13. Annotate edges with community information
# ============================================================

edge_table_graph <- as_data_frame(
  graph,
  what = "edges"
) %>%
  as_tibble() %>%
  rename(
    pattern1 = from,
    pattern2 = to
  ) %>%
  left_join(
    nodes %>%
      select(
        name,
        community_id,
        comm_infomap,
        CTP,
        centrality_degree,
        centrality_strength,
        centrality_betweenness,
        centrality_closeness
      ) %>%
      rename_with(
        ~ paste0(.x, "1"),
        -name
      ),
    by = c("pattern1" = "name")
  ) %>%
  left_join(
    nodes %>%
      select(
        name,
        community_id,
        comm_infomap,
        CTP,
        centrality_degree,
        centrality_strength,
        centrality_betweenness,
        centrality_closeness
      ) %>%
      rename_with(
        ~ paste0(.x, "2"),
        -name
      ),
    by = c("pattern2" = "name")
  ) %>%
  left_join(
    factor_info %>%
      rename(
        pattern1 = factor,
        dataset1 = dataset,
        k1 = k,
        pattern_number1 = pattern
      ),
    by = "pattern1"
  ) %>%
  left_join(
    factor_info %>%
      rename(
        pattern2 = factor,
        dataset2 = dataset,
        k2 = k,
        pattern_number2 = pattern
      ),
    by = "pattern2"
  ) %>%
  mutate(
    cross_dataset = dataset1 != dataset2,
    same_community = comm_infomap1 == comm_infomap2
  )

write_csv(
  edge_table_graph,
  file.path(table_dir, "edges_annotated_all.csv")
)

write_csv(
  edge_table_graph %>%
    filter(cross_dataset),
  file.path(table_dir, "edges_annotated_cross_dataset.csv")
)

saveRDS(
  graph,
  file.path(rds_dir, "Final_igraph_object.rds")
)

saveRDS(
  infomap_result,
  file.path(rds_dir, "infomap_result.rds")
)


# ============================================================
# 14. Network plot
# ============================================================

graph_tbl <- as_tbl_graph(graph)

set.seed(123)

p_network <- ggraph(
  graph_tbl,
  layout = "fr"
) +
  geom_edge_link(
    aes(
      alpha = cor,
      width = cor
    ),
    colour = "grey70",
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(
      colour = dataset,
      shape = factor(comm_infomap),
      size = centrality_degree
    ),
    alpha = 0.85
  ) +
  scale_size_continuous(
    range = c(1.5, 7)
  ) +
  guides(
    shape = "none"
  ) +
  theme_void() +
  labs(
    title = "NMF factor correlation network",
    subtitle = paste0(
      "Edges: Pearson r ≥ ",
      cor_threshold,
      "; communities: Infomap"
    ),
    colour = "Dataset",
    size = "Degree"
  )

ggsave(
  file.path(plot_dir, "factor_network_by_dataset.pdf"),
  p_network,
  width = 12,
  height = 10
)


p_network_community <- ggraph(
  graph_tbl,
  layout = "fr"
) +
  geom_edge_link(
    aes(
      alpha = cor,
      width = cor
    ),
    colour = "grey75",
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(
      colour = factor(comm_infomap),
      size = centrality_degree
    ),
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(1.5, 7)
  ) +
  guides(
    colour = guide_legend(
      title = "Infomap community",
      override.aes = list(size = 4)
    )
  ) +
  theme_void() +
  labs(
    title = "Infomap communities of NMF factors",
    subtitle = paste0(
      "Pearson correlation threshold = ",
      cor_threshold
    ),
    size = "Degree"
  )

ggsave(
  file.path(plot_dir, "factor_network_by_community.pdf"),
  p_network_community,
  width = 12,
  height = 10
)


# ============================================================
# 15. Community composition plots
# ============================================================

community_dataset_plot_df <- nodes %>%
  count(
    community_id,
    dataset,
    name = "n_factors"
  ) %>%
  group_by(community_id) %>%
  mutate(
    total_factors = sum(n_factors),
    fraction = n_factors / total_factors
  ) %>%
  ungroup()

community_order <- community_summary %>%
  arrange(
    desc(n_datasets),
    desc(total_factors)
  ) %>%
  pull(community_id)

community_dataset_plot_df$community_id <- factor(
  community_dataset_plot_df$community_id,
  levels = community_order
)

p_community_composition <- ggplot(
  community_dataset_plot_df,
  aes(
    x = community_id,
    y = n_factors,
    fill = dataset
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(
    title = "Dataset composition of Infomap communities",
    x = "Community",
    y = "Number of factors",
    fill = "Dataset"
  )

ggsave(
  file.path(plot_dir, "community_dataset_composition_counts.pdf"),
  p_community_composition,
  width = 9,
  height = max(5, length(community_order) * 0.3)
)


p_community_fraction <- ggplot(
  community_dataset_plot_df,
  aes(
    x = community_id,
    y = fraction,
    fill = dataset
  )
) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(
    labels = scales::percent
  ) +
  theme_classic(base_size = 11) +
  labs(
    title = "Relative dataset composition of Infomap communities",
    x = "Community",
    y = "Fraction of factors",
    fill = "Dataset"
  )

ggsave(
  file.path(plot_dir, "community_dataset_composition_fraction.pdf"),
  p_community_fraction,
  width = 9,
  height = max(5, length(community_order) * 0.3)
)


# ============================================================
# 16. Community aggregation function
# ============================================================

comm_agg_weights <- function(
    df,
    merged_A_matrices,
    weight_dataset_equal = FALSE,
    use_random_patterns = FALSE,
    rank_buffer = 9,
    ntopgenes = 100,
    random_seed = 123) {
  
  required_columns <- c(
    "name",
    "community_id",
    "dataset"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in df: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  communities <- unique(df$community_id)
  
  agg_weights <- matrix(
    0,
    nrow = nrow(merged_A_matrices),
    ncol = length(communities),
    dimnames = list(
      rownames(merged_A_matrices),
      communities
    )
  )
  
  set.seed(random_seed)
  
  for (community in communities) {
    
    pattern_names <- unique(
      df$name[df$community_id == community]
    )
    
    missing_patterns <- setdiff(
      pattern_names,
      colnames(merged_A_matrices)
    )
    
    if (length(missing_patterns) > 0) {
      stop(
        "Patterns missing from merged matrix: ",
        paste(head(missing_patterns, 10), collapse = ", ")
      )
    }
    
    if (use_random_patterns) {
      pattern_names <- sample(
        unique(df$name),
        length(pattern_names),
        replace = FALSE
      )
    }
    
    dataset_vector <- df$dataset[
      match(pattern_names, df$name)
    ]
    
    small_A <- merged_A_matrices[
      ,
      pattern_names,
      drop = FALSE
    ]
    
    # Rank 1 = highest loading gene
    rank_A <- apply(
      small_A,
      2,
      function(x) {
        rank(
          -x,
          ties.method = "average"
        )
      }
    )
    
    if (is.null(dim(rank_A))) {
      rank_A <- matrix(
        rank_A,
        ncol = 1,
        dimnames = list(
          rownames(small_A),
          colnames(small_A)
        )
      )
    }
    
    rownames(rank_A) <- rownames(small_A)
    colnames(rank_A) <- colnames(small_A)
    
    weights_A <- 1 / (
      rank_A + rank_buffer
    )
    
    if (weight_dataset_equal) {
      
      number_of_datasets <- length(
        unique(dataset_vector)
      )
      
      factor_counts <- table(dataset_vector)
      
      dataset_weights <- 1 /
        as.numeric(factor_counts[dataset_vector]) /
        number_of_datasets
      
      # Multiply each pattern column by its dataset-specific weight
      weights_A <- sweep(
        weights_A,
        2,
        dataset_weights,
        FUN = "*"
      )
    }
    
    # Preserve the original workflow:
    # community weight = mean across member factors
    agg_weights[, community] <- rowMeans(
      weights_A
    )
  }
  
  agg_ranks <- apply(
    agg_weights,
    2,
    function(x) {
      rank(
        -x,
        ties.method = "average"
      )
    }
  )
  
  if (is.null(dim(agg_ranks))) {
    agg_ranks <- matrix(
      agg_ranks,
      ncol = 1,
      dimnames = dimnames(agg_weights)
    )
  }
  
  rownames(agg_ranks) <- rownames(agg_weights)
  colnames(agg_ranks) <- colnames(agg_weights)
  
  number_to_select <- min(
    ntopgenes,
    nrow(agg_weights)
  )
  
  top_genes <- lapply(
    colnames(agg_weights),
    function(community) {
      
      gene_order <- order(
        agg_weights[, community],
        decreasing = TRUE
      )
      
      selected_genes <- rownames(
        agg_weights
      )[gene_order[seq_len(number_to_select)]]
      
      tibble(
        gene = selected_genes,
        aggWeight = agg_weights[
          selected_genes,
          community
        ],
        rank = seq_len(number_to_select)
      )
    }
  )
  
  names(top_genes) <- colnames(agg_weights)
  
  list(
    aggWeights = agg_weights,
    topGenes = top_genes,
    ranks = agg_ranks
  )
}


# ============================================================
# 17. Select multi-dataset communities
# ============================================================

ctp_nodes <- nodes %>%
  filter(
    CTP_strict == "yes"
  )

if (nrow(ctp_nodes) == 0) {
  stop(
    "No multi-dataset communities were detected. ",
    "Check the correlation threshold and community composition."
  )
}

candidate_ctps <- sort(
  unique(ctp_nodes$community_id)
)

cat("\nCandidate multi-dataset CTPs:\n")
print(candidate_ctps)

cat(
  "\nNumber of candidate CTP communities: ",
  length(candidate_ctps),
  "\n"
)


# ============================================================
# 18. Aggregate marker-gene weights
# ============================================================

message("Aggregating gene weights within CTP communities...")

agg <- comm_agg_weights(
  df = ctp_nodes,
  merged_A_matrices = MP,
  weight_dataset_equal = weight_dataset_equal,
  use_random_patterns = FALSE,
  rank_buffer = rank_buffer,
  ntopgenes = n_top_genes,
  random_seed = 123
)

agg_weights <- agg$aggWeights
agg_ranks <- agg$ranks
top_genes <- agg$topGenes

saveRDS(
  agg,
  file.path(rds_dir, "ctp_aggregation_complete.rds")
)

saveRDS(
  agg_weights,
  file.path(rds_dir, "ctp_aggweights_matrix.rds")
)

saveRDS(
  top_genes,
  file.path(rds_dir, "ctp_topGenes.rds")
)


# ============================================================
# 19. Export aggregate-weight matrix
# ============================================================

agg_weights_df <- as.data.frame(
  agg_weights,
  check.names = FALSE
) %>%
  rownames_to_column("gene")

write_csv(
  agg_weights_df,
  file.path(table_dir, "ctp_aggweights.csv")
)

fwrite(
  agg_weights_df,
  file.path(table_dir, "ctp_aggweights.tsv"),
  sep = "\t"
)

agg_ranks_df <- as.data.frame(
  agg_ranks,
  check.names = FALSE
) %>%
  rownames_to_column("gene")

fwrite(
  agg_ranks_df,
  file.path(table_dir, "ctp_gene_ranks.tsv"),
  sep = "\t"
)


# ============================================================
# 20. Export top genes
# ============================================================

top_genes_long <- bind_rows(
  top_genes,
  .id = "community_id"
)

write_csv(
  top_genes_long,
  file.path(table_dir, "ctp_top100_genes_long.csv")
)

for (community in names(top_genes)) {
  
  safe_name <- str_replace_all(
    community,
    "[^A-Za-z0-9_-]",
    "_"
  )
  
  write_csv(
    top_genes[[community]],
    file.path(
      topgene_dir,
      paste0(safe_name, "_top", n_top_genes, "_genes.csv")
    )
  )
}


# ============================================================
# 21. Top-gene heatmap
# ============================================================

# Use the union of the top 30 genes from every community
n_heatmap_genes_per_ctp <- 30

heatmap_genes <- top_genes %>%
  lapply(
    function(x) {
      head(
        x$gene,
        n_heatmap_genes_per_ctp
      )
    }
  ) %>%
  unlist() %>%
  unique()

heatmap_matrix <- agg_weights[
  heatmap_genes,
  ,
  drop = FALSE
]

# Row-wise z-score only for visualization

heatmap_matrix_z <- t(
  scale(
    t(heatmap_matrix)
  )
)

heatmap_matrix_z[
  !is.finite(heatmap_matrix_z)
] <- 0

pdf(
  file.path(plot_dir, "ctp_top_gene_aggregate_weight_heatmap.pdf"),
  width = 10,
  height = max(10, length(heatmap_genes) * 0.08)
)

pheatmap(
  heatmap_matrix_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  border_color = NA,
  main = "Aggregated marker-gene weights of candidate CTPs"
)

dev.off()


# ============================================================
# 22. CTP member-factor table
# ============================================================

ctp_member_factors <- ctp_nodes %>%
  arrange(
    community_id,
    dataset,
    k,
    pattern
  ) %>%
  select(
    community_id,
    comm_infomap,
    name,
    dataset,
    k,
    pattern,
    centrality_degree,
    centrality_strength,
    centrality_betweenness,
    centrality_closeness
  )

write_csv(
  ctp_member_factors,
  file.path(table_dir, "ctp_member_factors.csv")
)


# ============================================================
# 23. Representative factor per community
# ============================================================

representative_factors <- ctp_nodes %>%
  group_by(community_id) %>%
  arrange(
    desc(centrality_strength),
    desc(centrality_degree)
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    community_id,
    representative_factor = name,
    dataset,
    k,
    pattern,
    centrality_degree,
    centrality_strength,
    centrality_betweenness,
    centrality_closeness
  )

write_csv(
  representative_factors,
  file.path(table_dir, "ctp_representative_factors.csv")
)


# ============================================================
# 24. Generate a blank annotation template
# ============================================================

annotation_template <- community_summary %>%
  filter(
    community_id %in% candidate_ctps
  ) %>%
  select(
    community_id,
    total_factors,
    n_datasets,
    largest_dataset_fraction,
    dataset_composition
  ) %>%
  mutate(
    proposed_name = NA_character_,
    biological_annotation = NA_character_,
    evidence = NA_character_,
    retain = TRUE
  )

write_csv(
  annotation_template,
  file.path(table_dir, "ctp_manual_annotation_template.csv")
)


# ============================================================
# 25. Save analysis parameters
# ============================================================

analysis_parameters <- list(
  input_file = input_file,
  matrix_dimensions = dim(MP),
  expected_factors = expected_factors,
  correlation_threshold = cor_threshold,
  infomap_trials = infomap_trials,
  n_top_genes = n_top_genes,
  rank_buffer = rank_buffer,
  weight_dataset_equal = weight_dataset_equal,
  min_datasets_per_ctp = min_datasets_per_ctp,
  min_factors_per_dataset = min_factors_per_dataset,
  random_seed = 123,
  candidate_ctps = candidate_ctps,
  run_time = Sys.time()
)

saveRDS(
  analysis_parameters,
  file.path(rds_dir, "analysis_parameters.rds")
)


# ============================================================
# 26. Session information
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)


# ============================================================
# 27. Final summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CTP analysis completed\n")
cat("============================================================\n")
cat("Genes in merged matrix:       ", nrow(MP), "\n")
cat("Factors in merged matrix:     ", ncol(MP), "\n")
cat("Correlation threshold:        ", cor_threshold, "\n")
cat("Network nodes:                ", vcount(graph), "\n")
cat("Network edges:                ", ecount(graph), "\n")
cat(
  "Infomap communities:         ",
  length(unique(nodes$comm_infomap)),
  "\n"
)
cat(
  "Candidate multi-dataset CTPs:",
  length(candidate_ctps),
  "\n"
)
cat("Top genes per CTP:            ", n_top_genes, "\n")
cat("Output directory:             ", output_dir, "\n")
cat("============================================================\n")


#####================================================================================



