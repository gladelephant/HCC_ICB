# ============================================================
# Eight selected CTP heatmaps
#
# Panel A:
#   CTP marker aggregate-weight heatmap
#   rows    = top marker genes
#   columns = 8 selected CTPs
#
# Panel B:
#   Bulk sample expression heatmap in the style of the reference Fig. 1B
#   rows    = top marker genes, split by CTP
#   columns = tumor samples
#
# IMPORTANT:
# Replace selected_ctps with the exact community IDs or final names
# used in your own analysis.
# ============================================================

setwd("D:/HCC_ICB/cuda_nmf_results")

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
required_packages <- c(
  "tidyverse",
  "data.table",
  "ComplexHeatmap",
  "circlize",
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
  library(tidyverse)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(123)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
analysis_dir <- file.path(
  "output",
  "CTP_analysis"
)

rds_dir <- file.path(
  analysis_dir,
  "rds"
)

table_dir <- file.path(
  analysis_dir,
  "tables"
)

plot_dir <- file.path(
  analysis_dir,
  "selected_8_CTP_heatmaps"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

agg_file <- file.path(
  rds_dir,
  "ctp_aggweights_matrix.rds"
)

topgene_file <- file.path(
  rds_dir,
  "ctp_topGenes.rds"
)

if (!file.exists(agg_file)) {
  stop("Missing file: ", agg_file)
}

if (!file.exists(topgene_file)) {
  stop("Missing file: ", topgene_file)
}

agg_weights <- readRDS(agg_file)
top_genes <- readRDS(topgene_file)

agg_weights <- as.matrix(agg_weights)
storage.mode(agg_weights) <- "double"

# ------------------------------------------------------------
# 2. Select the eight biologically valuable CTPs
#
# Replace these examples with your exact 8 CTP IDs.
#
# First inspect:
#   colnames(agg_weights)
#   names(top_genes)
# ------------------------------------------------------------
print(colnames(agg_weights))
print(names(top_genes))

selected_ctps <- c(
  "CTP_10",
  "CTP_25",
  "CTP_5",
  "CTP_14",
  "CTP_11",
  "CTP_3",
  "CTP_8",
  "CTP_1"
)
missing_ctps <- setdiff(
  selected_ctps,
  colnames(agg_weights)
)

if (length(missing_ctps) > 0) {
  stop(
    "These selected CTPs are absent from agg_weights:\n",
    paste(missing_ctps, collapse = ", "),
    "\n\nAvailable CTPs are:\n",
    paste(colnames(agg_weights), collapse = ", ")
  )
}

selected_ctps <- intersect(
  selected_ctps,
  colnames(agg_weights)
)

# Optional publication labels.
# Replace the right-hand labels after biological annotation.
ctp_labels <- setNames(
  selected_ctps,
  selected_ctps
)

# Example:
# ctp_labels <- c(
#   CTP_3  = "Stem-like",
#   CTP_7  = "EMT",
#   CTP_9  = "IFN response",
#   CTP_12 = "Myeloid niche",
#   CTP_15 = "T-cell inflamed",
#   CTP_18 = "Angiogenic",
#   CTP_21 = "Proliferative",
#   CTP_24 = "Metabolic"
# )

# ------------------------------------------------------------
# 3. Parameters
# ------------------------------------------------------------
n_marker_weight <- 50
n_marker_expression <- 30

# For repeated genes:
# TRUE  = preserve the same gene in multiple CTP row blocks,
#         closest to the reference study.
# FALSE = show every gene only once.
preserve_repeated_genes <- TRUE

# ------------------------------------------------------------
# 4. Construct top-gene lists from aggregate weights
# ------------------------------------------------------------
get_ranked_genes <- function(
    ctp,
    n_genes
) {
  x <- agg_weights[, ctp]
  
  x <- x[
    is.finite(x)
  ]
  
  genes <- names(
    sort(
      x,
      decreasing = TRUE
    )
  )
  
  head(
    unique(genes),
    n_genes
  )
}

marker_lists_weight <- setNames(
  lapply(
    selected_ctps,
    get_ranked_genes,
    n_genes = n_marker_weight
  ),
  selected_ctps
)

marker_lists_expression <- setNames(
  lapply(
    selected_ctps,
    get_ranked_genes,
    n_genes = n_marker_expression
  ),
  selected_ctps
)

# ============================================================
# PANEL A
# Eight-CTP aggregate marker-weight heatmap
# ============================================================

if (preserve_repeated_genes) {
  marker_vector_weight <- unlist(
    marker_lists_weight,
    use.names = FALSE
  )
  
  marker_group_weight <- rep(
    names(marker_lists_weight),
    lengths(marker_lists_weight)
  )
} else {
  marker_vector_weight <- unique(
    unlist(
      marker_lists_weight,
      use.names = FALSE
    )
  )
  
  marker_group_weight <- vapply(
    marker_vector_weight,
    function(gene) {
      hits <- names(marker_lists_weight)[
        vapply(
          marker_lists_weight,
          function(x) gene %in% x,
          logical(1)
        )
      ]
      
      hits[1]
    },
    character(1)
  )
}

weight_mat <- agg_weights[
  marker_vector_weight,
  selected_ctps,
  drop = FALSE
]

colnames(weight_mat) <- unname(
  ctp_labels[colnames(weight_mat)]
)

row_split_weight <- factor(
  unname(
    ctp_labels[marker_group_weight]
  ),
  levels = unname(
    ctp_labels[selected_ctps]
  )
)

positive_weights <- weight_mat[
  is.finite(weight_mat) &
    weight_mat > 0
]

if (length(positive_weights) == 0) {
  stop("No positive aggregate weights were found.")
}

weight_breaks <- unique(
  as.numeric(
    quantile(
      positive_weights,
      probs = c(0, 0.50, 0.85, 0.97, 1),
      na.rm = TRUE
    )
  )
)

if (length(weight_breaks) < 3) {
  weight_breaks <- seq(
    min(positive_weights),
    max(positive_weights),
    length.out = 5
  )
}

weight_col_fun <- colorRamp2(
  weight_breaks,
  colorRampPalette(
    c(
      "#FFFFFF",
      "#FFF7BC",
      "#FEC44F",
      "#EC7014",
      "#8C2D04"
    )
  )(length(weight_breaks))
)

# Objectively label top 3 genes from every selected CTP.
label_genes <- unique(
  unlist(
    lapply(
      marker_lists_weight,
      head,
      10
    ),
    use.names = FALSE
  )
)

label_at <- which(
  rownames(weight_mat) %in% label_genes
)

left_mark <- rowAnnotation(
  Gene = anno_mark(
    at = label_at,
    labels = rownames(weight_mat)[label_at],
    side = "left",
    labels_gp = gpar(
      fontsize = 6
    ),
    link_width = unit(
      3,
      "mm"
    ),
    padding = 0.5
  )
)

h_weight <- Heatmap(
  weight_mat,
  name = "Aggregate\nweight",
  col = weight_col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_split = row_split_weight,
  row_gap = unit(
    0.8,
    "mm"
  ),
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  
  row_title_side = "right",
  row_title_rot = 0,
  row_title_gp = gpar(
    fontsize = 7,
    fontface = "bold"
  ),
  
  left_annotation = left_mark,
  
  border = FALSE,
  use_raster = TRUE,
  raster_quality = 5,
  
  heatmap_legend_param = list(
    title_gp = gpar(
      fontsize = 8,
      fontface = "bold"
    ),
    labels_gp = gpar(
      fontsize = 7
    )
  )
)

pdf(
  file.path(
    plot_dir,
    "Eight_CTP_marker_aggregate_weight_heatmap.pdf"
  ),
  width = 6,
  height = max(
    7,
    nrow(weight_mat) * 0.025
  ),
  useDingbats = FALSE
)

draw(
  h_weight,
  heatmap_legend_side = "right",
  padding = unit(
    c(3, 3, 3, 3),
    "mm"
  )
)

dev.off()

# ============================================================
# PANEL B
# Eight-CTP bulk-expression heatmap
#
# Expected files from your prior bulk workflow:
#   HCC_ICB_altas_452cases.csv
#       rows    = genes
#       columns = samples
#
#   Total_metadata_452例.csv
#       must contain a sample ID column and cohort/group columns.
#
# Adjust the four settings immediately below when necessary.
# ============================================================

expression_file <- "HCC_ICB_altas_452cases.csv"
metadata_file <- "Total_metadata_452例.csv"

sample_id_column <- "sampleID"
cohort_column <- "Group"

# TRUE when expression_file contains TPM values.
# FALSE when it is already log2(TPM + 1).
expression_is_raw_tpm <- TRUE

if (
  file.exists(expression_file) &&
  file.exists(metadata_file)
) {
  expression_df <- fread(
    expression_file,
    data.table = FALSE,
    check.names = FALSE
  )
  
  gene_column <- colnames(expression_df)[1]
  gene_names <- as.character(
    expression_df[[gene_column]]
  )
  
  expression_df[[gene_column]] <- NULL
  
  expr <- as.matrix(
    expression_df
  )
  
  storage.mode(expr) <- "double"
  rownames(expr) <- gene_names
  
  if (expression_is_raw_tpm) {
    expr <- log2(
      expr + 1
    )
  }
  
  meta <- read.csv(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!sample_id_column %in% colnames(meta)) {
    stop(
      "Metadata does not contain sample ID column: ",
      sample_id_column
    )
  }
  
  if (!cohort_column %in% colnames(meta)) {
    stop(
      "Metadata does not contain cohort column: ",
      cohort_column
    )
  }
  
  common_samples <- intersect(
    meta[[sample_id_column]],
    colnames(expr)
  )
  
  if (length(common_samples) < 2) {
    stop(
      "Too few shared sample IDs between expression and metadata."
    )
  }
  
  meta <- meta[
    match(
      common_samples,
      meta[[sample_id_column]]
    ),
    ,
    drop = FALSE
  ]
  
  expr <- expr[
    ,
    common_samples,
    drop = FALSE
  ]
  
  expression_genes <- intersect(
    unique(
      unlist(
        marker_lists_expression,
        use.names = FALSE
      )
    ),
    rownames(expr)
  )
  
  if (length(expression_genes) == 0) {
    stop(
      "None of the selected CTP marker genes were found in expression data."
    )
  }
  
  # Preserve repeated genes in different CTP row blocks.
  expression_rows <- unlist(
    lapply(
      marker_lists_expression,
      function(x) {
        x[
          x %in% rownames(expr)
        ]
      }
    ),
    use.names = FALSE
  )
  
  expression_row_group <- rep(
    names(marker_lists_expression),
    vapply(
      marker_lists_expression,
      function(x) {
        sum(
          x %in% rownames(expr)
        )
      },
      integer(1)
    )
  )
  
  expr_sub <- expr[
    expression_rows,
    ,
    drop = FALSE
  ]
  
  # Cohort-specific gene z-scaling.
  # This reduces platform/cohort scale effects while preserving
  # relative expression patterns within each cohort.
  cohort_values <- as.character(
    meta[[cohort_column]]
  )
  
  cohort_levels <- unique(
    cohort_values
  )
  
  scaled_parts <- lapply(
    cohort_levels,
    function(cohort) {
      sample_ids <- meta[
        cohort_values == cohort,
        sample_id_column
      ]
      
      m <- expr_sub[
        ,
        sample_ids,
        drop = FALSE
      ]
      
      z <- t(
        scale(
          t(m)
        )
      )
      
      z[
        !is.finite(z)
      ] <- 0
      
      z
    }
  )
  
  expr_z <- do.call(
    cbind,
    scaled_parts
  )
  
  # Restore metadata sample order.
  expr_z <- expr_z[
    ,
    meta[[sample_id_column]],
    drop = FALSE
  ]
  
  # Limit extreme values for visualization only.
  expr_z[expr_z > 2.5] <- 2.5
  expr_z[expr_z < -2.5] <- -2.5
  
  # Sample clustering from selected marker-gene expression.
  sample_cor <- cor(
    expr_z,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  
  sample_hc <- hclust(
    as.dist(
      1 - sample_cor
    ),
    method = "average"
  )
  
  sample_order <- colnames(
    expr_z
  )[sample_hc$order]
  
  expr_z <- expr_z[
    ,
    sample_order,
    drop = FALSE
  ]
  
  meta_plot <- meta[
    match(
      sample_order,
      meta[[sample_id_column]]
    ),
    ,
    drop = FALSE
  ]
  
  cohort_palette <- setNames(
    rep(
      c(
        "#0072B2",
        "#D55E00",
        "#009E73",
        "#CC79A7",
        "#E69F00",
        "#56B4E9",
        "#999999"
      ),
      length.out = length(
        unique(
          meta_plot[[cohort_column]]
        )
      )
    ),
    unique(
      meta_plot[[cohort_column]]
    )
  )
  
  top_anno_expression <- HeatmapAnnotation(
    Cohort = meta_plot[[cohort_column]],
    col = list(
      Cohort = cohort_palette
    ),
    simple_anno_size = unit(
      2.5,
      "mm"
    ),
    annotation_name_gp = gpar(
      fontsize = 8,
      fontface = "bold"
    ),
    show_annotation_name = TRUE
  )
  
  expression_col_fun <- colorRamp2(
    c(
      -2.5,
      0,
      2.5
    ),
    c(
      "#2166AC",
      "#F7F7F7",
      "#B2182B"
    )
  )
  
  expression_split <- factor(
    unname(
      ctp_labels[expression_row_group]
    ),
    levels = unname(
      ctp_labels[selected_ctps]
    )
  )
  
  h_expression <- Heatmap(
    expr_z,
    name = "Scaled\nexpression",
    col = expression_col_fun,
    
    top_annotation = top_anno_expression,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    row_split = expression_split,
    row_gap = unit(
      0.8,
      "mm"
    ),
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    row_title_side = "right",
    row_title_rot = 0,
    row_title_gp = gpar(
      fontsize = 7,
      fontface = "bold"
    ),
    
    border = FALSE,
    use_raster = TRUE,
    raster_quality = 5,
    
    heatmap_legend_param = list(
      at = c(
        -2,
        0,
        2
      ),
      title_gp = gpar(
        fontsize = 8,
        fontface = "bold"
      ),
      labels_gp = gpar(
        fontsize = 7
      )
    )
  )
  
  pdf(
    file.path(
      plot_dir,
      "Eight_CTP_bulk_marker_expression_heatmap.pdf"
    ),
    width = 8,
    height = max(
      7,
      nrow(expr_z) * 0.025
    ),
    useDingbats = FALSE
  )
  
  draw(
    h_expression,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(
      c(3, 3, 3, 3),
      "mm"
    )
  )
  
  dev.off()
  
  write_csv(
    tibble(
      sampleID = sample_order,
      heatmap_order = seq_along(
        sample_order
      )
    ) %>%
      left_join(
        meta,
        by = setNames(
          sample_id_column,
          "sampleID"
        )
      ),
    file.path(
      plot_dir,
      "Eight_CTP_expression_heatmap_sample_order.csv"
    )
  )
  
  message(
    "Panel B was generated successfully."
  )
} else {
  message(
    "Expression or metadata file was not found.\n",
    "Panel A was generated, but Panel B was skipped.\n",
    "Expected expression file: ",
    expression_file,
    "\nExpected metadata file: ",
    metadata_file
  )
}

# ------------------------------------------------------------
# 5. Export selected CTP marker table
# ------------------------------------------------------------
selected_marker_table <- imap_dfr(
  marker_lists_weight,
  function(genes, ctp) {
    tibble(
      CTP_ID = ctp,
      CTP_label = unname(
        ctp_labels[ctp]
      ),
      rank = seq_along(
        genes
      ),
      gene = genes,
      aggregate_weight = agg_weights[
        genes,
        ctp
      ]
    )
  }
)

write_csv(
  selected_marker_table,
  file.path(
    plot_dir,
    "Eight_CTP_selected_marker_genes.csv"
  )
)

message(
  "\nFinished.\nOutput directory:\n",
  normalizePath(
    plot_dir,
    mustWork = FALSE
  )
)
