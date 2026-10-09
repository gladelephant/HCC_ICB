
#1. 加载R包 #####================================================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

#2. 设置路径 #####================================================================================

ctp_score_file <- "G:\\08_HCC_靶免治疗\\CTP_V2\\output\\CTP_analysis\\rds\\CTP_scores.rds"
  
metadata_file <- "G:/08_HCC_靶免治疗/HCC免疫治疗bulk集合/TPM_清洗后/HCC_ICB_total/Total_metadata_452例.csv"

out_dir <- "output/CTP_analysis/selected_8_CTP_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#3. 读取CTP评分与疗效信息 #####================================================================================

scores <- read_rds(
  ctp_score_file
)

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

df <- scores |>
  mutate(sample = as.character(sample)) |>
  inner_join(meta, by = "sample") |>
  filter(!is.na(Response)) |>
  distinct(sample, .keep_all = TRUE)


#2. 提取8个CTP评分 #####================================================================================

selected_ctps <- c(
  "CTP_1", "CTP_3", "CTP_5", "CTP_8",
  "CTP_10", "CTP_11", "CTP_14", "CTP_25"
)

score_df <- df |>
  select(all_of(selected_ctps)) |>
  mutate(across(everything(), as.numeric))


#3. 计算Spearman相关性 #####================================================================================

cor_mat <- cor(
  score_df,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# 计算两两相关的P值
pairs <- combn(selected_ctps, 2, simplify = FALSE)

cor_results <- purrr::map_dfr(pairs, function(pair) {
  
  x <- score_df[[pair[1]]]
  y <- score_df[[pair[2]]]
  
  test <- cor.test(
    x, y,
    method = "spearman",
    exact = FALSE
  )
  
  tibble(
    CTP_1 = pair[1],
    CTP_2 = pair[2],
    rho = unname(test$estimate),
    pvalue = test$p.value
  )
}) |>
  mutate(FDR = p.adjust(pvalue, method = "BH"))

write_csv(
  cor_results,
  file.path(out_dir, "CTP_score_correlations.csv")
)


#4. 绘制相关性热图 #####================================================================================

cor_col <- colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "#FFFFFF", "#B2182B")
)

# 只显示下三角（包含对角线）
cor_plot <- cor_mat
cor_plot[upper.tri(cor_plot)] <- NA

h_cor <- Heatmap(
  cor_plot,
  name = "Spearman rho",
  col = cor_col,
  
  na_col = "white",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_side = "left",
  column_names_rot = 45,
  
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  
  rect_gp = gpar(col = "white", lwd = 1),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    if (i >= j) {
      grid.text(
        sprintf("%.2f", cor_mat[i, j]),
        x, y,
        gp = gpar(
          fontsize = 8,
          col = if (abs(cor_mat[i, j]) > 0.65)
            "white" else "#303030"
        )
      )
    }
  },
  
  width = unit(75, "mm"),
  height = unit(75, "mm"),
  
  heatmap_legend_param = list(
    at = c(-1, -0.5, 0, 0.5, 1)
  )
)


#5. 保存相关性热图 #####================================================================================

pdf(
  file.path(out_dir, "PanelC_CTP_score_correlation.pdf"),
  width = 5.5,
  height = 5,
  useDingbats = FALSE
)

draw(h_cor, heatmap_legend_side = "right")
dev.off()
