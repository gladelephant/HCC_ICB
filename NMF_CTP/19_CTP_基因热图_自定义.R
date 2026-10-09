
#1. 加载R包 #####================================================================================

library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("G:/08_HCC_靶免治疗/CTP_V2")


#2. 设置路径与代表基因 #####================================================================================

rds_dir <- "output/CTP_analysis/rds"
out_dir <- "output/CTP_analysis/selected_8_CTP_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CTP_marker_list <- list(
  CTP_1  = c("COL11A1", "PDPN", "COMP", "FGF7", "SFRP4"),
  CTP_3  = c("IGF2", "IGF2BP1", "TRIM71", "MYCN", "PEG3"),
  CTP_5  = c("AXIN2", "NOTUM", "NKD1", "LGR5", "SP5"),
  CTP_8  = c("AKR1B10", "NQO1", "SRXN1", "KRT23", "S100P"),
  CTP_10 = c("CD8B", "PDCD1", "CXCL13", "TIGIT", "GZMB"),
  CTP_11 = c("STAB2", "FCN2", "OIT3", "MARCO", "RSPO3"),
  CTP_14 = c("CYP2A7", "THRSP", "GCGR", "SLC25A47", "MOGAT2"),
  CTP_25 = c("GAGE2E", "PAGE2B", "XAGE1A", "SSX1", "IL13RA2")
)

selected_ctps <- names(CTP_marker_list)
marker_genes <- unlist(CTP_marker_list, use.names = FALSE)

row_groups <- factor(
  rep(selected_ctps, lengths(CTP_marker_list)),
  levels = selected_ctps
)


#3. Panel A：CTP aggregate weight热图 #####================================================================================

agg_weights <- as.matrix(
  readRDS(file.path(rds_dir, "ctp_aggweights_matrix.rds"))
)

weight_mat <- agg_weights[
  marker_genes,
  selected_ctps,
  drop = FALSE
]

# weight_col <- colorRamp2(
#   quantile(weight_mat, c(0, 0.5, 0.9, 1), na.rm = TRUE),
#   c("#FFFFFF", "#FFF7BC", "#FEC44F", "#B35806")
# )
# 
# library(circlize)
# library(RColorBrewer)
# #双向spectral配色#####################################
# weight_col <- colorRamp2(
#   seq(
#     min(weight_mat, na.rm = TRUE),
#     max(weight_mat, na.rm = TRUE),
#     length.out = 7
#   ),
#   rev(brewer.pal(7, "Spectral"))
# )
# 
# #单向蓝色配色#####################################
# weight_col <- colorRamp2(
#   seq(
#     min(weight_mat, na.rm = TRUE),
#     max(weight_mat, na.rm = TRUE),
#     length.out = 5
#   ),
#   c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
# )

#Morpheus风格单向红色配色 #####================================================================================

weight_col <- colorRamp2(
  seq(
    min(weight_mat, na.rm = TRUE),
    max(weight_mat, na.rm = TRUE),
    length.out = 5
  ),
  c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D")
)

h_weight <- Heatmap(
  weight_mat,
  name = "Aggregate weight",
  col = weight_col,
  
  row_split = row_groups,
  row_gap = unit(1.2, "mm"),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 45,
  
  row_title_side = "right",
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  
  border = FALSE,
  use_raster = FALSE
)

pdf(
  file.path(out_dir, "CTP_marker_aggregate_weight.pdf"),
  width = 7,
  height = 9,
  useDingbats = FALSE
)

draw(h_weight, heatmap_legend_side = "right")
dev.off()
