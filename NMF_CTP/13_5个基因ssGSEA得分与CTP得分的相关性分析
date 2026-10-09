
library(tidyverse)
library(GSVA)
library(meta)
library(ggsci)

# 1. 数据读取 ============================================================
expr_bulk <- read.csv(
  "F:/01_HCC_ICB_paper/Supplemental_tables/HCC_ICB_altas_452cases_log2TPM.csv",
  row.names = 1, check.names = FALSE
) |> as.matrix()

meta_bulk <- read.csv(
  "F:/01_HCC_ICB_paper/Supplemental_tables/Total_metadata_452例.csv"
)

ctp_scores <- read_rds(
  "G:/08_HCC_靶免治疗/CTP_V2/output/CTP_analysis/rds/CTP_scores.rds"
)

gene5 <- c("REG3A", "GLUD1", "AIFM2", "CTNNA2", "CACFD1")
gene4 <- setdiff(gene5, "REG3A")

# 2. 样本匹配与ssGSEA =====================================================
samples <- Reduce(intersect, list(
  colnames(expr_bulk), meta_bulk$sample, ctp_scores$sample
))

expr_bulk <- expr_bulk[, samples, drop = FALSE]

scores <- gsva(
  ssgseaParam(
    exprData = expr_bulk,
    geneSets = list(Gene5 = gene5, Gene4 = gene4),
    normalize = TRUE
  ),
  verbose = FALSE
) |>
  t() |>
  as.data.frame() |>
  rownames_to_column("sample")

df <- scores |>
  left_join(meta_bulk, by = "sample") |>
  left_join(ctp_scores, by = "sample")

ctp_cols <- grep("^CTP_[0-9]+$", names(df), value = TRUE)

# 3. Gene5与所有CTP的Spearman相关性 =======================================
cor_all <- map_dfr(ctp_cols, \(ctp) {
  test <- cor.test(df$Gene5, df[[ctp]],
                   method = "spearman", exact = FALSE)
  
  tibble(
    CTP = ctp,
    rho = unname(test$estimate),
    p = test$p.value
  )
}) |>
  mutate(FDR = p.adjust(p, "BH")) |>
  arrange(desc(rho))

print(cor_all, n = Inf)

# 4. 八个CTP叠加散点图 ====================================================
plot_df <- df |>
  select(Gene5, all_of(ctp_cols)) |>
  pivot_longer(-Gene5, names_to = "CTP", values_to = "Score")

cols <- c(
  CTP_5  = "#E64B35",
  CTP_8  = "#4DBBD5",
  CTP_14 = "#B09C85",
  CTP_3  = "#00A087",
  CTP_10 = "#3C5488",
  CTP_25 = "#8491B4",
  CTP_11 = "#F39B7F",
  CTP_1  = "#7E6148"
)

p1 <- ggplot(plot_df, aes(Gene5, Score, color = CTP)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  scale_color_manual(
    values = cols,
    breaks = cor_all$CTP,
    labels = sprintf("%s (r = %.3f)", cor_all$CTP, cor_all$rho)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 12)
  ) +
  labs(x = "Five-gene ssGSEA score", y = "CTP score")

print(p1)
ggsave("Eight_CTP_regression.pdf", p1, width = 4.5, height = 3)

# 5. Gene5与CTP_5的分队列Spearman相关性 ==================================
cor_dataset <- df |>
  group_by(Dataset) |>
  group_modify(~ {
    dat <- .x |> drop_na(Gene5, CTP_5)
    test <- cor.test(dat$Gene5, dat$CTP_5,
                     method = "spearman", exact = FALSE)
    
    tibble(
      n = nrow(dat),
      rho = unname(test$estimate),
      p = test$p.value
    )
  }) |>
  ungroup()

print(cor_dataset)

# 6. Gene5与CTP_5的分队列散点图 ==========================================
legend_labels <- setNames(
  sprintf("%s (r = %.3f)", cor_dataset$Dataset, cor_dataset$rho),
  cor_dataset$Dataset
)

p2 <- ggplot(df, aes(Gene5, CTP_5, color = Dataset)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  scale_color_npg(labels = legend_labels) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 12)
  ) +
  labs(x = "Five-gene ssGSEA score", y = "CTP_5 score")

print(p2)
ggsave("Gene5_CTP5_scatter.pdf", p2, width = 5, height = 3)

# 7. Spearman相关系数随机效应Meta分析 ====================================
fit <- metacor(
  cor = rho,
  n = n,
  studlab = Dataset,
  data = cor_dataset,
  sm = "ZCOR",
  method.tau = "REML",
  method.random.ci = "HK",
  common = FALSE,
  random = TRUE
)

summary(fit)

meta_result <- tibble(
  rho = tanh(fit$TE.random),
  CI_low = tanh(fit$lower.random),
  CI_high = tanh(fit$upper.random),
  p = fit$pval.random,
  I2 = 100 * fit$I2,
  tau2 = fit$tau2
)

print(meta_result)

# 8. Nature风格森林图 =====================================================
forest_df <- bind_rows(
  tibble(
    Dataset = as.character(fit$studlab),
    n = fit$n,
    rho = tanh(fit$TE),
    ci_low = tanh(fit$lower),
    ci_high = tanh(fit$upper),
    Type = "Cohort"
  ),
  tibble(
    Dataset = "Pooled (REML)",
    n = sum(cor_dataset$n),
    rho = meta_result$rho,
    ci_low = meta_result$CI_low,
    ci_high = meta_result$CI_high,
    Type = "Pooled"
  )
) |>
  mutate(
    Dataset = factor(Dataset, levels = rev(unique(Dataset))),
    Type = factor(Type, levels = c("Cohort", "Pooled"))
  )

p3 <- ggplot(forest_df,
             aes(x = rho, y = Dataset, color = Type)) +
  geom_vline(
    xintercept = 0, linetype = 2,
    color = "grey75", linewidth = 0.4
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = Dataset),
    linewidth = 0.8
  ) +
  geom_point(aes(size = n), shape = 16) +
  scale_color_manual(values = c(
    Cohort = "#3C5488",
    Pooled = "#E64B35"
  )) +
  scale_size(range = c(2.5, 5)) +
  scale_x_continuous(
    breaks = seq(-0.5, 1, 0.25)
  ) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11),
    plot.subtitle = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  ) +
  labs(
    subtitle = sprintf(
      "pooled r = %.2f, p = %.2g, I² = %.1f%%",
      meta_result$rho,
      meta_result$p,
      meta_result$I2
    ),
    x = "Spearman correlation (95% CI)",
    y = NULL
  )

print(p3)

ggsave(
  "forest_plot_Gene5_CTP5_score_association.pdf",
  p3, width = 4.6, height = 2.3, bg = "white"
)
