
setwd("F:/01_HCC_ICB_paper/5基因得分与CTP得分相关性分析")

#1. 加载R包 #####================================================================================

library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(scales)


#2. 设置路径与参数 #####================================================================================

ctp_dir <- "G:/08_HCC_靶免治疗/CTP_V2/output/CTP_analysis/top_genes"
enrich_file <- "F:/01_HCC_ICB_paper/5基因得分与CTP得分相关性分析/enrichr_self_enrichment_items.csv"

enrich_dir <- file.path(ctp_dir, "enrichment")
out_dir <- "CTP_Nature_style"

dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

selected_file <- "CTP_selected_pathways.csv"
max_pathways <- Inf

ctp_order <- c(
  "CTP_1", "CTP_3", "CTP_5", "CTP_8",
  "CTP_10", "CTP_11", "CTP_14"
)

ctp_colors <- c(
  CTP_1  = "#B65C58",
  CTP_3  = "#637DA7",
  CTP_5  = "#468B7A",
  CTP_8  = "#9473A6",
  CTP_10 = "#C38B43",
  CTP_11 = "#4E91A2",
  CTP_14 = "#7686B0"
)


#3. 读取CTP Top50基因 #####================================================================================

files <- list.files(
  ctp_dir,
  pattern = "^CTP_[0-9]+_top100_genes\\.csv$",
  full.names = TRUE
)

gene_list <- setNames(
  lapply(files, function(f) {
    read_csv(f, show_col_types = FALSE) |>
      arrange(rank) |>
      pull(gene) |>
      as.character() |>
      unique() |>
      head(50)
  }),
  str_remove(basename(files), "_top100_genes\\.csv$")
)

gene_list <- gene_list[
  order(as.integer(str_remove(names(gene_list), "CTP_")))
]

CTP_top50_df <- as_tibble(gene_list)

write_csv(
  CTP_top50_df,
  file.path(ctp_dir, "All_CTP_top50_genes.csv")
)


#4. 导出CTP代表基因 #####================================================================================

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

CTP_marker_df <- tibble(
  CTP = rep(names(CTP_marker_list), lengths(CTP_marker_list)),
  gene = unlist(CTP_marker_list, use.names = FALSE)
)

write_csv(CTP_marker_df, "CTP_40_marker_genes.csv")


#5. 批量富集分析 #####================================================================================

TERM2GENE <- read_csv(
  enrich_file,
  show_col_types = FALSE
) |>
  transmute(
    term = gs_name,
    gene = gene_symbol
  ) |>
  filter(!is.na(term), !is.na(gene)) |>
  distinct()

cc <- compareCluster(
  geneCluster = gene_list,
  fun = "enricher",
  TERM2GENE = TERM2GENE,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

enrich_all <- as.data.frame(cc)

write_csv(
  enrich_all,
  file.path(enrich_dir, "CTP_enrichment_50genes.csv")
)


#6. 整理人工筛选的富集通路 #####================================================================================

df <- read_csv(
  selected_file,
  show_col_types = FALSE
) |>
  mutate(
    Cluster = factor(Cluster, levels = ctp_order),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio)),
    GeneRatio_num = GeneCount / as.numeric(sub(".*/", "", GeneRatio))
  ) |>
  filter(!is.na(Cluster)) |>
  arrange(Cluster, p.adjust)

if (is.finite(max_pathways)) {
  df <- df |>
    group_by(Cluster) |>
    slice_head(n = as.integer(max_pathways)) |>
    ungroup()
}

df <- df |>
  ungroup() |>
  mutate(
    row_id = row_number(),
    group_id = as.integer(Cluster),
    y = -row_id - (group_id - 1) * 0.9
  )

group_centers <- df |>
  group_by(Cluster) |>
  summarise(
    y_min = min(y) - 0.35,
    y_max = max(y) + 0.35,
    y_center = mean(range(y)),
    .groups = "drop"
  )

x_max <- max(
  ceiling(max(df$GeneRatio_num) * 10) / 10,
  0.05
)

bar_x <- x_max * 1.12
label_x <- x_max * 1.18
x_total_max <- x_max * 1.45


#7. 绘制CTP富集Lollipop图 #####================================================================================

p <- ggplot(df) +
  
  geom_segment(
    aes(
      x = 0, xend = GeneRatio_num,
      y = y, yend = y,
      color = Cluster
    ),
    linewidth = 0.42,
    alpha = 0.48
  ) +
  
  geom_point(
    aes(
      x = GeneRatio_num,
      y = y,
      size = GeneCount,
      fill = Cluster
    ),
    shape = 21,
    color = "white",
    stroke = 0.18
  ) +
  
  geom_segment(
    data = group_centers,
    aes(
      x = bar_x, xend = bar_x,
      y = y_min, yend = y_max,
      color = Cluster
    ),
    inherit.aes = FALSE,
    linewidth = 2.8
  ) +
  
  geom_text(
    data = group_centers,
    aes(
      x = label_x,
      y = y_center,
      label = Cluster
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.8,
    fontface = "bold",
    color = "#303030"
  ) +
  
  scale_color_manual(values = ctp_colors, guide = "none") +
  scale_fill_manual(values = ctp_colors, guide = "none") +
  
  scale_size_area(
    max_size = 5,
    breaks = c(5, 10, 15, 20),
    limits = c(0, max(20, df$GeneCount)),
    name = "Gene count"
  ) +
  
  scale_x_continuous(
    limits = c(0, x_total_max),
    breaks = seq(0, x_max, by = 0.1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = 0)
  ) +
  
  scale_y_continuous(
    breaks = df$y,
    labels = df$Description,
    expand = expansion(add = 0.7)
  ) +
  
  guides(
    size = guide_legend(
      title.position = "left",
      direction = "horizontal",
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "black",
        color = "black",
        stroke = 0,
        alpha = 1
      )
    )
  ) +
  
  labs(x = "Gene ratio", y = NULL) +
  
  theme_classic(base_size = 9) +
  
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text.y = element_text(
      size = 7.5,
      color = "#303030",
      margin = margin(r = 5)
    ),
    
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key = element_blank(),
    legend.key.width = unit(10, "mm"),
    legend.key.height = unit(7, "mm"),
    
    plot.margin = margin(10, 12, 10, 10)
  )

print(p)


#8. 保存图片 #####================================================================================

ggsave(
  file.path(out_dir, "CTP_Nature_lollipop_final.pdf"),
  p,
  width = 6.7,
  height = 6,
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(out_dir, "CTP_Nature_lollipop_final.png"),
  p,
  width = 6.7,
  height = 6,
  dpi = 600,
  bg = "white"
)
