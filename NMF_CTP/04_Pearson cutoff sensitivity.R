# ============================================================
# Pearson cutoff sensitivity analysis
# thresholds: 0.40, 0.45, 0.50, 0.55
# ============================================================

library(dplyr)
library(igraph)
library(purrr)
library(tidyr)
library(ggplot2)
library(tibble)

thresholds <- c(
  0.40,
  0.45,
  0.50,
  0.55
)

# Infomap重复次数，与主流程保持一致
infomap_trials <- 200

# ============================================================
# 1. Prepare vertex table
# ============================================================

# igraph要求顶点名称列叫name，并且位于第一列
vertex_table <- factor_info %>%
  transmute(
    name = as.character(factor),
    dataset = as.character(dataset),
    k = as.integer(k),
    pattern = as.integer(pattern)
  ) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )

threshold_results <- vector(
  mode = "list",
  length = length(thresholds)
)

names(threshold_results) <- as.character(thresholds)


# ============================================================
# 2. Run network analysis at each threshold
# ============================================================

for (i in seq_along(thresholds)) {
  
  thr <- thresholds[i]
  
  cat("\n")
  cat("==================================================\n")
  cat("Pearson threshold =", thr, "\n")
  cat("==================================================\n")
  
  # ----------------------------------------------------------
  # 2.1 Filter edges
  # ----------------------------------------------------------
  
  edge_use <- edges_all %>%
    filter(
      is.finite(cor),
      cor >= thr
    ) %>%
    transmute(
      from = as.character(pattern1),
      to = as.character(pattern2),
      weight = as.numeric(cor),
      cross_dataset = as.logical(cross_dataset)
    ) %>%
    as.data.frame(
      stringsAsFactors = FALSE
    )
  
  cat("Retained edges:", nrow(edge_use), "\n")
  
  # 如果该阈值没有任何边
  if (nrow(edge_use) == 0) {
    
    threshold_results[[i]] <- tibble(
      threshold = thr,
      n_edges = 0L,
      n_nodes_total = nrow(vertex_table),
      n_nodes_connected = 0L,
      isolated_nodes = nrow(vertex_table),
      isolated_fraction = 1,
      network_density = NA_real_,
      largest_component = 0L,
      largest_component_fraction = NA_real_,
      n_communities = 0L,
      largest_community = 0L,
      largest_community_fraction = NA_real_,
      cross_dataset_edges = 0L,
      cross_dataset_fraction = NA_real_,
      n_candidate_ctps_2datasets = 0L,
      n_candidate_ctps_3datasets = 0L
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 2.2 Construct graph including all factors
  # ----------------------------------------------------------
  
  g_all <- igraph::graph_from_data_frame(
    d = edge_use,
    directed = FALSE,
    vertices = vertex_table
  )
  
  degree_all <- igraph::degree(
    g_all,
    mode = "all"
  )
  
  isolated_nodes <- sum(
    degree_all == 0
  )
  
  isolated_fraction <-
    isolated_nodes /
    igraph::vcount(g_all)
  
  # 删除孤立节点后再进行community detection
  g <- igraph::delete_vertices(
    g_all,
    igraph::V(g_all)[degree_all == 0]
  )
  
  cat("Total nodes:", igraph::vcount(g_all), "\n")
  cat("Connected nodes:", igraph::vcount(g), "\n")
  cat("Isolated nodes:", isolated_nodes, "\n")
  cat("Graph edges:", igraph::ecount(g), "\n")
  
  # ----------------------------------------------------------
  # 2.3 Extract numeric edge weights safely
  # ----------------------------------------------------------
  
  edge_weights <- igraph::edge_attr(
    g,
    "weight"
  )
  
  edge_weights <- as.numeric(
    unlist(
      edge_weights,
      use.names = FALSE
    )
  )
  
  if (length(edge_weights) != igraph::ecount(g)) {
    stop(
      "Edge weight length does not match edge number at threshold ",
      thr
    )
  }
  
  if (
    anyNA(edge_weights) ||
    any(!is.finite(edge_weights))
  ) {
    stop(
      "NA or non-finite edge weights detected at threshold ",
      thr
    )
  }
  
  cat(
    "Edge-weight range:",
    paste(
      round(
        range(edge_weights),
        4
      ),
      collapse = "–"
    ),
    "\n"
  )
  
  # ----------------------------------------------------------
  # 2.4 Connected components
  # ----------------------------------------------------------
  
  component_result <- igraph::components(g)
  
  largest_component <- max(
    component_result$csize
  )
  
  largest_component_fraction <-
    largest_component /
    igraph::vcount(g)
  
  # ----------------------------------------------------------
  # 2.5 Infomap community detection
  # ----------------------------------------------------------
  
  set.seed(123)
  
  infomap_result <- igraph::cluster_infomap(
    graph = g,
    e.weights = edge_weights,
    nb.trials = infomap_trials
  )
  
  membership_vector <- igraph::membership(
    infomap_result
  )
  
  community_sizes <- as.numeric(
    table(membership_vector)
  )
  
  # ----------------------------------------------------------
  # 2.6 Community composition across datasets
  # ----------------------------------------------------------
  
  node_result <- tibble(
    name = igraph::V(g)$name,
    dataset = as.character(
      igraph::vertex_attr(
        g,
        "dataset"
      )
    ),
    community = as.integer(
      membership_vector[
        igraph::V(g)$name
      ]
    )
  )
  
  community_result <- node_result %>%
    group_by(community) %>%
    summarise(
      n_factors = n(),
      n_datasets = n_distinct(dataset),
      dataset_composition = paste0(
        dataset,
        collapse = ";"
      ),
      .groups = "drop"
    )
  
  n_candidate_ctps_2datasets <- sum(
    community_result$n_datasets >= 2
  )
  
  n_candidate_ctps_3datasets <- sum(
    community_result$n_datasets >= 3
  )
  
  # ----------------------------------------------------------
  # 2.7 Cross-dataset edge statistics
  # ----------------------------------------------------------
  
  cross_dataset_edges <- sum(
    edge_use$cross_dataset,
    na.rm = TRUE
  )
  
  cross_dataset_fraction <- mean(
    edge_use$cross_dataset,
    na.rm = TRUE
  )
  
  # ----------------------------------------------------------
  # 2.8 Save threshold-level metrics
  # ----------------------------------------------------------
  
  threshold_results[[i]] <- tibble(
    threshold = thr,
    
    n_edges = igraph::ecount(g),
    
    n_nodes_total = igraph::vcount(g_all),
    
    n_nodes_connected = igraph::vcount(g),
    
    isolated_nodes = isolated_nodes,
    
    isolated_fraction = isolated_fraction,
    
    network_density = igraph::edge_density(
      g,
      loops = FALSE
    ),
    
    largest_component = largest_component,
    
    largest_component_fraction =
      largest_component_fraction,
    
    n_communities = length(
      community_sizes
    ),
    
    largest_community = max(
      community_sizes
    ),
    
    largest_community_fraction =
      max(community_sizes) /
      sum(community_sizes),
    
    cross_dataset_edges =
      cross_dataset_edges,
    
    cross_dataset_fraction =
      cross_dataset_fraction,
    
    n_candidate_ctps_2datasets =
      n_candidate_ctps_2datasets,
    
    n_candidate_ctps_3datasets =
      n_candidate_ctps_3datasets
  )
  
  cat(
    "Infomap communities:",
    length(community_sizes),
    "\n"
  )
  
  cat(
    "Largest community:",
    max(community_sizes),
    "\n"
  )
  
  cat(
    "Cross-dataset edges:",
    cross_dataset_edges,
    "\n"
  )
  
  cat(
    "CTPs supported by >=2 datasets:",
    n_candidate_ctps_2datasets,
    "\n"
  )
  
  cat(
    "CTPs supported by all 3 datasets:",
    n_candidate_ctps_3datasets,
    "\n"
  )
}


# ============================================================
# 3. Combine results
# ============================================================

threshold_summary <- bind_rows(
  threshold_results
)

print(
  threshold_summary,
  width = Inf
)


# ============================================================
# 4. Save result table
# ============================================================

write.csv(
  threshold_summary,
  file.path(
    table_dir,
    "pearson_cutoff_sensitivity.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 5. Select metrics for plotting
# ============================================================

metrics_to_plot <- c(
  "n_edges",
  "isolated_fraction",
  "largest_component_fraction",
  "largest_community_fraction",
  "cross_dataset_edges",
  "n_candidate_ctps_2datasets"
)

threshold_plot_df <- threshold_summary %>%
  select(
    threshold,
    all_of(metrics_to_plot)
  ) %>%
  pivot_longer(
    cols = -threshold,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      
      n_edges =
        "Retained edges",
      
      isolated_fraction =
        "Isolated-node fraction",
      
      largest_component_fraction =
        "Largest-component fraction",
      
      largest_community_fraction =
        "Largest-community fraction",
      
      cross_dataset_edges =
        "Cross-dataset edges",
      
      n_candidate_ctps_2datasets =
        "Multi-dataset CTPs"
    ),
    
    threshold = factor(
      threshold,
      levels = thresholds,
      labels = sprintf(
        "%.2f",
        thresholds
      )
    )
  )


# ============================================================
# 6. Nature-style sensitivity plot
# ============================================================

p_cutoff <- ggplot(
  threshold_plot_df,
  aes(
    x = threshold,
    y = value,
    group = 1
  )
) +
  geom_line(
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_point(
    size = 2.4,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.6
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    x = "Pearson correlation cutoff",
    y = NULL,
    title = "Sensitivity analysis of factor-correlation cutoff"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 9,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 9.5,
      colour = "black"
    ),
    
    axis.text = element_text(
      size = 8,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.35,
      colour = "black"
    ),
    
    panel.grid = element_blank(),
    
    panel.spacing = unit(
      0.9,
      "lines"
    ),
    
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

p_cutoff


# ============================================================
# 7. Save plot
# ============================================================

ggsave(
  filename = file.path(
    plot_dir,
    "pearson_cutoff_sensitivity.pdf"
  ),
  plot = p_cutoff,
  width = 7.2,
  height = 7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    plot_dir,
    "pearson_cutoff_sensitivity.tiff"
  ),
  plot = p_cutoff,
  width = 7.2,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
