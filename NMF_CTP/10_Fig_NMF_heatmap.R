# ============================================================
# Nature-style Fig. S1B
# Consensus Pearson-correlation heatmap of CUDA-NMF factors
#
# Input:
#   output/CTP_analysis/rds/merged_A_matrices.rds
#
# Optional input:
#   output/CTP_analysis/rds/rand_merged_A_matrices.rds
#
# Output:
#   output/CTP_analysis/reference_style_plots/
#       Fig_S1B_Nature_real.pdf
#       Fig_S1B_Nature_real.png
#       Fig_S1B_Nature_real_vs_random.pdf   (if random matrix exists)
#       Fig_S1B_factor_order.csv
# ============================================================

setwd("D:/HCC_ICB/cuda_nmf_results")

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
required_packages <- c(
  "ComplexHeatmap",
  "circlize",
  "tidyverse",
  "RColorBrewer",
  "grid"
)

to_install <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(to_install) > 0) {
  install.packages(to_install)
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(tidyverse)
  library(RColorBrewer)
  library(grid)
})

set.seed(123)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
rds_dir <- file.path(
  "output",
  "CTP_analysis",
  "rds"
)

plot_dir <- file.path(
  "output",
  "CTP_analysis",
  "reference_style_plots"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

real_file <- file.path(
  rds_dir,
  "merged_A_matrices.rds"
)

random_file <- file.path(
  rds_dir,
  "rand_merged_A_matrices.rds"
)

if (!file.exists(real_file)) {
  stop(
    "Cannot find real factor matrix:\n",
    normalizePath(real_file, mustWork = FALSE)
  )
}

# ------------------------------------------------------------
# 2. Read and validate factor matrix
# ------------------------------------------------------------
MP <- readRDS(real_file)
MP <- as.matrix(MP)
storage.mode(MP) <- "double"

if (is.null(rownames(MP)) || is.null(colnames(MP))) {
  stop("MP must have gene row names and factor column names.")
}

if (anyNA(MP) || any(!is.finite(MP))) {
  stop("MP contains NA or non-finite values.")
}

zero_sd <- apply(
  MP,
  2,
  function(x) sd(x) == 0
)

if (any(zero_sd)) {
  warning(
    sum(zero_sd),
    " zero-variance factors were removed."
  )
  MP <- MP[, !zero_sd, drop = FALSE]
}

message(
  "Real matrix: ",
  nrow(MP),
  " genes × ",
  ncol(MP),
  " factors"
)

# ------------------------------------------------------------
# 3. Parse factor metadata
#
# Expected format:
#   DATASET_k10_P1
#
# Dataset names may contain underscores.
# ------------------------------------------------------------
parse_factor_metadata <- function(factor_names) {
  out <- tibble(
    factor = factor_names,
    dataset = str_remove(
      factor_names,
      "_k[0-9]+_P[0-9]+$"
    ),
    k = as.integer(
      str_match(
        factor_names,
        "_k([0-9]+)_P[0-9]+$"
      )[, 2]
    ),
    pattern = as.integer(
      str_match(
        factor_names,
        "_k[0-9]+_P([0-9]+)$"
      )[, 2]
    )
  )
  
  invalid <- out %>%
    filter(
      is.na(k) |
        is.na(pattern) |
        dataset == factor
    )
  
  if (nrow(invalid) > 0) {
    print(invalid)
    stop(
      "Some factor names do not follow ",
      "DATASET_kK_PPATTERN."
    )
  }
  
  out
}

factor_meta <- parse_factor_metadata(
  colnames(MP)
)

datasets <- unique(factor_meta$dataset)

# ------------------------------------------------------------
# 4. Color design
#
# Dataset colors are saturated but restrained.
# K is represented by a continuous grey scale.
# Correlation colors deliberately emphasize r >= 0.50.
# ------------------------------------------------------------
dataset_palette_base <- c(
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#009E73", # green
  "#CC79A7", # purple-pink
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#999999", # grey
  "#000000"  # black
)

dataset_cols <- setNames(
  rep(
    dataset_palette_base,
    length.out = length(datasets)
  ),
  datasets
)

k_min <- min(factor_meta$k)
k_max <- max(factor_meta$k)

k_col_fun <- colorRamp2(
  c(k_min, mean(c(k_min, k_max)), k_max),
  c("#F2F2F2", "#A6A6A6", "#262626")
)

# White is used for weak correlations.
# The biologically relevant network threshold r = 0.50
# is deliberately placed at the beginning of the strong red range.
cor_col_fun <- colorRamp2(
  c(-0.25, 0, 0.25, 0.50, 0.70, 1.00),
  c(
    "#3B4CC0",
    "#F7F7F7",
    "#F7F7F7",
    "#FDD0A2",
    "#EF6548",
    "#7F0000"
  )
)

# ------------------------------------------------------------
# 5. Correlation and clustering
# ------------------------------------------------------------
cor_real <- cor(
  MP,
  method = "pearson",
  use = "pairwise.complete.obs"
)

if (anyNA(cor_real)) {
  stop("The real correlation matrix contains NA values.")
}

# Correlation distance gives one identical order for rows and columns.
# Negative correlations are treated as distant.
distance_real <- as.dist(
  1 - cor_real
)


hc_real <- hclust(
  distance_real,
  method = "average"
)

factor_order <- colnames(cor_real)[hc_real$order]

factor_meta_ordered <- factor_meta %>%
  mutate(
    cluster_order = match(
      factor,
      factor_order
    )
  ) %>%
  arrange(cluster_order)

write_csv(
  factor_meta_ordered,
  file.path(
    plot_dir,
    "Fig_S1B_factor_order.csv"
  )
)

# ------------------------------------------------------------
# 6. Annotation objects
# ------------------------------------------------------------
annotation_df <- factor_meta %>%
  select(
    factor,
    dataset,
    k
  ) %>%
  column_to_rownames("factor")

top_anno <- HeatmapAnnotation(
  Dataset = annotation_df$dataset,
  Rank = annotation_df$k,
  col = list(
    Dataset = dataset_cols,
    Rank = k_col_fun
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(
    fontsize = 7,
    fontface = "bold"
  ),
  simple_anno_size = unit(2.2, "mm"),
  gap = unit(0.5, "mm"),
  border = FALSE,
  show_annotation_name = TRUE,
  show_legend = FALSE
)


# ------------------------------------------------------------
# 7. Heatmap constructor
# ------------------------------------------------------------
make_nature_heatmap <- function(
    cor_mat,
    hc,
    title = NULL,
    show_dataset_legend = TRUE
) {
  Heatmap(
    cor_mat,
    name = "Pearson r",
    col = cor_col_fun,
    
    cluster_rows = as.dendrogram(hc),
    cluster_columns = as.dendrogram(hc),
    
    top_annotation = top_anno,
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    row_title = NULL,
    column_title = title,
    column_title_gp = gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    
    border = FALSE,
    rect_gp = gpar(
      col = NA
    ),
    
    use_raster = TRUE,
    raster_quality = 5,
    raster_device = "png",
    
    heatmap_legend_param = list(
      title = "Pearson r",
      title_gp = gpar(
        fontsize = 8,
        fontface = "bold"
      ),
      labels_gp = gpar(
        fontsize = 7
      ),
      at = c(
        -0.25,
        0,
        0.25,
        0.50,
        0.70,
        1.00
      ),
      labels = c(
        "-0.25",
        "0",
        "0.25",
        "0.50",
        "0.70",
        "1.00"
      ),
      legend_height = unit(
        31,
        "mm"
      ),
      grid_width = unit(
        3,
        "mm"
      )
    )
  )
}

ht_real <- make_nature_heatmap(
  cor_mat = cor_real,
  hc = hc_real,
  title = "Real-data CUDA-NMF factors"
)

# ------------------------------------------------------------
# 8. Dataset legend
# ------------------------------------------------------------
dataset_legend <- Legend(
  title = "Dataset",
  labels = names(dataset_cols),
  legend_gp = gpar(
    fill = unname(dataset_cols),
    col = NA
  ),
  title_gp = gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = gpar(
    fontsize = 7
  ),
  grid_width = unit(
    3,
    "mm"
  ),
  grid_height = unit(
    3,
    "mm"
  ),
  gap = unit(
    1.2,
    "mm"
  )
)

rank_legend <- Legend(
  title = "NMF rank",
  col_fun = k_col_fun,
  at = unique(
    round(
      c(
        k_min,
        mean(c(k_min, k_max)),
        k_max
      )
    )
  ),
  title_gp = gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = gpar(
    fontsize = 7
  ),
  legend_height = unit(
    23,
    "mm"
  ),
  grid_width = unit(
    3,
    "mm"
  )
)

# ------------------------------------------------------------
# 9. Export real-data heatmap
# ------------------------------------------------------------
pdf(
  file.path(
    plot_dir,
    "Fig_S1B_Nature_real.pdf"
  ),
  width = 6,
  height = 5,
  useDingbats = FALSE
)

draw(
  ht_real,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(
    dataset_legend,
    rank_legend
  ),
  merge_legends = FALSE,
  padding = unit(
    c(3, 3, 3, 3),
    "mm"
  )
)

dev.off()

png(
  file.path(
    plot_dir,
    "Fig_S1B_Nature_real.png"
  ),
  width = 2200,
  height = 2000,
  res = 400,
  bg = "white"
)

draw(
  ht_real,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(
    dataset_legend,
    rank_legend
  ),
  merge_legends = FALSE,
  padding = unit(
    c(3, 3, 3, 3),
    "mm"
  )
)

dev.off()

# ------------------------------------------------------------
# 10. Optional real-versus-random comparison
#
# The random matrix must be generated by:
#   gene-wise expression permutation
#   -> rerun NMF across all ranks and cohorts
#   -> merge the resulting factors
#
# Do not replace it with a shuffled correlation matrix.
# ------------------------------------------------------------
if (file.exists(random_file)) {
  rand_MP <- readRDS(random_file)
  rand_MP <- as.matrix(rand_MP)
  storage.mode(rand_MP) <- "double"
  
  if (anyNA(rand_MP) || any(!is.finite(rand_MP))) {
    stop("rand_MP contains NA or non-finite values.")
  }
  
  rand_zero_sd <- apply(
    rand_MP,
    2,
    function(x) sd(x) == 0
  )
  
  if (any(rand_zero_sd)) {
    rand_MP <- rand_MP[
      ,
      !rand_zero_sd,
      drop = FALSE
    ]
  }
  
  random_meta <- parse_factor_metadata(
    colnames(rand_MP)
  )
  
  if (!setequal(
    random_meta$dataset,
    factor_meta$dataset
  )) {
    warning(
      "Real and random matrices contain different dataset names."
    )
  }
  
  cor_random <- cor(
    rand_MP,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  
  distance_random <- as.dist(
    pmax(
      0,
      1 - cor_random
    )
  )
  
  hc_random <- hclust(
    distance_random,
    method = "average"
  )
  
  # Build annotations specific to the random factors.
  random_annotation_df <- random_meta %>%
    select(
      factor,
      dataset,
      k
    ) %>%
    column_to_rownames("factor")
  
  top_anno_random <- HeatmapAnnotation(
    Dataset = random_annotation_df$dataset,
    Rank = random_annotation_df$k,
    col = list(
      Dataset = dataset_cols,
      Rank = k_col_fun
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(
      fontsize = 7,
      fontface = "bold"
    ),
    simple_anno_size = unit(
      2.2,
      "mm"
    ),
    gap = unit(
      0.5,
      "mm"
    ),
    border = FALSE,
    show_annotation_name = TRUE,
    show_legend = FALSE
  )
  
  left_anno_random <- rowAnnotation(
    Dataset = random_annotation_df$dataset,
    col = list(
      Dataset = dataset_cols
    ),
    simple_anno_size = unit(
      2.2,
      "mm"
    ),
    border = FALSE,
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  
  ht_random <- Heatmap(
    cor_random,
    name = "Pearson r",
    col = cor_col_fun,
    
    cluster_rows = as.dendrogram(
      hc_random
    ),
    cluster_columns = as.dendrogram(
      hc_random
    ),
    
    top_annotation = top_anno_random,
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    column_title = "Gene-wise permuted controls",
    column_title_gp = gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    
    border = FALSE,
    rect_gp = gpar(
      col = NA
    ),
    
    use_raster = TRUE,
    raster_quality = 5,
    raster_device = "png",
    
    show_heatmap_legend = FALSE
  )
  
  pdf(
    file.path(
      plot_dir,
      "Fig_S1B_Nature_real_vs_random.pdf"
    ),
    width = 9.2,
    height = 4.8,
    useDingbats = FALSE
  )
  
  draw(
    ht_real + ht_random,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(
      dataset_legend,
      rank_legend
    ),
    merge_legends = FALSE,
    padding = unit(
      c(3, 3, 3, 3),
      "mm"
    )
  )
  
  dev.off()
  
  message(
    "Random control detected: real-versus-random panel was created."
  )
} else {
  message(
    "No random factor matrix found. ",
    "Only the real-data panel was created."
  )
}

message(
  "\nFinished.\nOutput directory:\n",
  normalizePath(
    plot_dir,
    mustWork = FALSE
  )
)

