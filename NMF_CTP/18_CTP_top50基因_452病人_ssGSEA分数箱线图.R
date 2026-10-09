
#1. 加载R包 #####================================================================================

library(tidyverse)
library(GSVA)
library(ggplot2)


#2. 设置路径与参数 #####================================================================================

setwd("G:/08_HCC_靶免治疗/CTP_V2")

expr_file <- "F:/01_HCC_ICB_paper/Supplemental_tables/HCC_ICB_altas_452cases_log2TPM.csv"

ctp_gene_file <- "G:/08_HCC_靶免治疗/CTP_V2/output/CTP_analysis/rds/ctp_topGenes.rds"

metadata_file <- "G:/08_HCC_靶免治疗/HCC免疫治疗bulk集合/TPM_清洗后/HCC_ICB_total/Total_metadata_452例.csv"

out_dir <- "output/CTP_analysis/CTP_Top50_ssGSEA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

selected_ctps <- c(
  "CTP_1", "CTP_3", "CTP_5", "CTP_8",
  "CTP_10", "CTP_11", "CTP_14", "CTP_25"
)

top_n <- 50

response_colors <- c(
  R  = "#438C79",
  NR = "#C76B70"
)


#3. 读取452例bulk表达矩阵 #####================================================================================

expr_df <- read_csv(
  expr_file,
  show_col_types = FALSE
)

expr <- expr_df |>
  column_to_rownames(names(expr_df)[1]) |>
  as.matrix()

storage.mode(expr) <- "numeric"

# 删除缺失基因名
expr <- expr[
  !is.na(rownames(expr)) & rownames(expr) != "",
  ,
  drop = FALSE
]

# 重复基因取平均表达量
gene_count <- table(rownames(expr))

expr <- rowsum(
  expr,
  group = rownames(expr),
  reorder = FALSE
)

expr <- sweep(
  expr,
  1,
  as.numeric(gene_count[rownames(expr)]),
  "/"
)


#4. 提取8个CTP的Top50基因 #####================================================================================

top_genes <- readRDS(ctp_gene_file)

ctp_gene_sets <- lapply(selected_ctps, function(ctp) {
  
  gene_df <- top_genes[[ctp]] |>
    as_tibble() |>
    filter(
      !is.na(gene),
      !is.na(aggWeight)
    ) |>
    arrange(desc(aggWeight)) |>
    distinct(gene, .keep_all = TRUE) |>
    slice_head(n = top_n)
  
  intersect(
    gene_df$gene,
    rownames(expr)
  )
  
})

names(ctp_gene_sets) <- selected_ctps

gene_summary <- tibble(
  CTP = selected_ctps,
  TopN = top_n,
  Matched_Genes = lengths(ctp_gene_sets)
)

print(gene_summary, n = Inf)

write_csv(
  gene_summary,
  file.path(out_dir, "CTP_Top50_gene_summary.csv")
)

# 保存实际参与ssGSEA计算的基因集合
gene_set_df <- enframe(
  ctp_gene_sets,
  name = "CTP",
  value = "gene"
) |>
  unnest(gene)

write_csv(
  gene_set_df,
  file.path(out_dir, "CTP_Top50_gene_sets.csv")
)


#5. 计算Top50基因ssGSEA评分 #####================================================================================

ssgsea_param <- GSVA::ssgseaParam(
  exprData = expr,
  geneSets = ctp_gene_sets,
  minSize = 5,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea_scores <- GSVA::gsva(
  ssgsea_param,
  verbose = TRUE
)

saveRDS(
  ssgsea_scores,
  file.path(out_dir, "CTP_Top50_ssGSEA_scores.rds")
)


#6. 整理452例患者的ssGSEA评分 #####================================================================================

score_df <- t(ssgsea_scores) |>
  as.data.frame() |>
  rownames_to_column("sample")

meta <- read_csv(
  metadata_file,
  show_col_types = FALSE
) |>
  transmute(
    sample = as.character(sample),
    Response = case_when(
      ORR == "Responder" ~ "R",
      ORR == "Non_responder" ~ "NR",
      TRUE ~ NA_character_
    )
  )

df <- score_df |>
  inner_join(meta, by = "sample") |>
  filter(!is.na(Response)) |>
  distinct(sample, .keep_all = TRUE)

plot_df <- df |>
  select(
    sample,
    Response,
    all_of(selected_ctps)
  ) |>
  pivot_longer(
    cols = all_of(selected_ctps),
    names_to = "CTP",
    values_to = "Score"
  ) |>
  mutate(
    CTP = factor(
      CTP,
      levels = selected_ctps
    ),
    Response = factor(
      Response,
      levels = c("R", "NR")
    )
  ) |>
  filter(is.finite(Score))



#7. 8个CTP整体Wilcoxon检验 #####================================================================================

stat_df <- plot_df |>
  group_by(CTP) |>
  summarise(
    n_R = sum(Response == "R"),
    n_NR = sum(Response == "NR"),
    
    median_R = median(Score[Response == "R"]),
    median_NR = median(Score[Response == "NR"]),
    
    pvalue = wilcox.test(
      Score ~ Response,
      exact = FALSE
    )$p.value,
    
    y_max = max(Score),
    y_range = diff(range(Score)),
    
    .groups = "drop"
  ) |>
  mutate(
    Difference = median_NR - median_R,
    
    Direction = case_when(
      Difference > 0 ~ "NR > R",
      Difference < 0 ~ "R > NR",
      TRUE ~ "Equal"
    ),
    
    # 防止个别CTP分数范围过小时括号重叠
    y_span = if_else(
      y_range > 0,
      y_range,
      pmax(abs(y_max) * 0.1, 0.1)
    ),
    
    # 括号横线位置
    bracket_y = y_max + 0.10 * y_span,
    
    # 括号两端短竖线长度
    tip_length = 0.025 * y_span,
    
    # P值位于括号上方
    y_position = y_max + 0.25 * y_span,
    
    label = paste0(
      "p = ",
      format.pval(
        pvalue,
        digits = 2,
        eps = 0.001
      )
    )
  )

print(stat_df, n = Inf)


#8. 绘制8个CTP整体箱线图（4列×2行，带P值括号） #####================================================================================

set.seed(123)

p <- ggplot(
  plot_df,
  aes(
    x = Response,
    y = Score,
    fill = Response
  )
) +
  
  # 箱线图
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    linewidth = 0.45,
    color = "#404040",
    alpha = 0.65
  ) +
  
  # 患者散点
  geom_jitter(
    aes(color = Response),
    width = 0.17,
    height = 0,
    size = 0.65,
    alpha = 0.28
  ) +
  
  # 显著性括号：水平线
  geom_segment(
    data = stat_df,
    aes(
      x = 1,
      xend = 2,
      y = bracket_y,
      yend = bracket_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "#454545"
  ) +
  
  # 显著性括号：左侧短竖线
  geom_segment(
    data = stat_df,
    aes(
      x = 1,
      xend = 1,
      y = bracket_y,
      yend = bracket_y - tip_length
    ),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "#454545"
  ) +
  
  # 显著性括号：右侧短竖线
  geom_segment(
    data = stat_df,
    aes(
      x = 2,
      xend = 2,
      y = bracket_y,
      yend = bracket_y - tip_length
    ),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "#454545"
  ) +
  
  # P值居中显示在括号上方
  geom_text(
    data = stat_df,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.2,
    color = "#303030"
  ) +
  
  # 8个CTP：上4张，下4张
  facet_wrap(
    ~ CTP,
    ncol = 4,
    nrow = 2,
    scales = "free_y"
  ) +
  
  # R与NR配色
  scale_fill_manual(
    values = response_colors
  ) +
  
  scale_color_manual(
    values = response_colors
  ) +
  
  # 为P值和括号预留空间
  scale_y_continuous(
    expand = expansion(
      mult = c(0.06, 0.10)
    )
  ) +
  
  labs(
    x = NULL,
    y = "CTP Top50 ssGSEA score"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 8)
    ),
    
    axis.line = element_line(
      linewidth = 0.4,
      color = "#444444"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.35
    ),
    
    panel.spacing = unit(9, "mm"),
    
    legend.position = "none",
    
    plot.margin = margin(
      10, 12, 10, 10
    )
  )

print(p)


#9. 保存图片与统计结果 #####================================================================================

# PDF矢量图
ggsave(
  filename = file.path(
    out_dir,
    "CTP_Top50_ssGSEA_overall_Wilcox.pdf"
  ),
  plot = p,
  width = 6.2,
  height = 4.3,
  device = cairo_pdf,
  bg = "white"
)

# PNG高清图
ggsave(
  filename = file.path(
    out_dir,
    "CTP_Top50_ssGSEA_overall_Wilcox.png"
  ),
  plot = p,
  width = 12,
  height = 6,
  dpi = 600,
  bg = "white"
)

# Wilcoxon统计结果
write_csv(
  stat_df,
  file.path(
    out_dir,
    "CTP_Top50_ssGSEA_Wilcox_results.csv"
  )
)

# ssGSEA评分矩阵
write_csv(
  score_df,
  file.path(
    out_dir,
    "CTP_Top50_ssGSEA_scores.csv"
  )
)

# 绘图数据
write_csv(
  plot_df,
  file.path(
    out_dir,
    "CTP_Top50_ssGSEA_plot_data.csv"
  )
)
