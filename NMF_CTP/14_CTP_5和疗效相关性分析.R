
library(tidyverse)
library(broom)
library(metafor)

# 1. 数据读取与样本匹配 ====================================================
expr <- read.csv(
  "F:/01_HCC_ICB_paper/Supplemental_tables/HCC_ICB_altas_452cases_log2TPM.csv",
  row.names = 1, check.names = FALSE
) |> as.matrix()

storage.mode(expr) <- "double"

meta <- read.csv(
  "F:/01_HCC_ICB_paper/Supplemental_tables/Total_metadata_452例.csv"
)

top_genes <- readRDS(
  "G:/08_HCC_靶免治疗/CTP_V2/output/CTP_analysis/rds/ctp_topGenes.rds"
)

samples <- intersect(colnames(expr), meta$sample)
expr <- expr[, samples, drop = FALSE]
meta <- meta[match(samples, meta$sample), ]

# 2. CTP加权评分（weighted = FALSE） ======================================
CTP.marker.gene.list <- lapply(top_genes, \(x) {
  data.frame(aggWeight = x$aggWeight, row.names = x$gene)
})

CTPscorev2 <- function(marker_list, mat) {
  scores <- lapply(marker_list, \(x) {
    w <- setNames(x$aggWeight, rownames(x))
    w <- w[names(w) %in% rownames(mat)]
    colSums(sweep(mat[names(w), , drop = FALSE], 1, w, "*"))
  })
  do.call(rbind, scores)
}

ctp.score.raw <- CTPscorev2(CTP.marker.gene.list, expr)

ctp.score.df <- t(ctp.score.raw) |>
  as.data.frame() |>
  rownames_to_column("sample")

ctp.meta.raw <- inner_join(meta, ctp.score.df, by = "sample")

saveRDS(ctp.score.raw, "HCC_CTP_raw_score.rds")
write.csv(ctp.meta.raw, "HCC_CTP_raw_score_metadata.csv",
          row.names = FALSE)

# 3. 队列内Z-score标准化 ===================================================
ctp.names <- rownames(ctp.score.raw)

ctp.meta.z <- ctp.meta.raw |>
  group_by(Group) |>
  mutate(across(all_of(ctp.names), ~ as.numeric(scale(.x)))) |>
  ungroup()

ctp.long <- ctp.meta.z |>
  pivot_longer(
    cols = all_of(ctp.names),
    names_to = "CTP",
    values_to = "score"
  )

# 4. 应答变量 ==============================================================
ctp.analysis <- ctp.long |>
  filter(ORR %in% c("Responder", "Non_responder")) |>
  mutate(
    NR = as.integer(ORR == "Non_responder"),
    ORR = factor(ORR, levels = c("Responder", "Non_responder"))
  )

# 5. 各队列Logistic回归 ====================================================
logistic.result <- ctp.analysis |>
  group_by(Group, CTP) |>
  group_modify(~ {
    dat <- filter(.x, is.finite(score))
    n_NR <- sum(dat$NR == 1)
    n_R <- sum(dat$NR == 0)
    
    if (nrow(dat) < 6 || n_NR < 2 || n_R < 2) {
      return(tibble(
        n = nrow(dat), n_NR = n_NR, n_R = n_R,
        beta = NA_real_, SE = NA_real_, OR = NA_real_,
        CI_low = NA_real_, CI_high = NA_real_, pvalue = NA_real_
      ))
    }
    
    fit <- glm(NR ~ score, family = binomial(), data = dat)
    x <- tidy(fit) |> filter(term == "score")
    
    tibble(
      n = nrow(dat), n_NR = n_NR, n_R = n_R,
      beta = x$estimate,
      SE = x$std.error,
      OR = exp(x$estimate),
      CI_low = exp(x$estimate - 1.96 * x$std.error),
      CI_high = exp(x$estimate + 1.96 * x$std.error),
      pvalue = x$p.value
    )
  }) |>
  ungroup()

write.csv(logistic.result, "CTP_ORR_logistic_by_cohort.csv",
          row.names = FALSE)

# 6. REML随机效应Meta分析 ==================================================
meta.result <- logistic.result |>
  filter(is.finite(beta), is.finite(SE), SE > 0) |>
  group_by(CTP) |>
  group_modify(~ {
    if (nrow(.x) < 2) {
      return(tibble(
        k = nrow(.x), beta_meta = NA_real_, SE_meta = NA_real_,
        OR_meta = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        p_meta = NA_real_, tau2 = NA_real_, I2 = NA_real_
      ))
    }
    
    fit <- rma.uni(
      yi = beta, sei = SE, data = .x,
      method = "REML", test = "z"
    )
    
    tibble(
      k = fit$k,
      beta_meta = as.numeric(fit$b),
      SE_meta = fit$se,
      OR_meta = exp(as.numeric(fit$b)),
      CI_low = exp(fit$ci.lb),
      CI_high = exp(fit$ci.ub),
      p_meta = fit$pval,
      tau2 = fit$tau2,
      I2 = fit$I2
    )
  }) |>
  ungroup() |>
  mutate(FDR = p.adjust(p_meta, method = "BH"))

# 7. 方向一致性 ===========================================================
direction.result <- logistic.result |>
  filter(is.finite(beta)) |>
  group_by(CTP) |>
  summarise(
    n_dataset = n(),
    n_resistance = sum(beta > 0),
    n_response = sum(beta < 0),
    direction_consistency = pmax(n_resistance, n_response) / n_dataset,
    .groups = "drop"
  )

meta.result <- meta.result |>
  left_join(direction.result, by = "CTP") |>
  arrange(beta_meta)

write.csv(meta.result, "CTP_ORR_meta_results.csv",
          row.names = FALSE)

print(meta.result, n = Inf)

# 8. Figure A：所有CTP的Meta效应 ==========================================
plot.meta <- meta.result |>
  filter(is.finite(OR_meta), is.finite(FDR)) |>
  mutate(
    log2OR = log2(OR_meta),
    significance = case_when(
      FDR < 0.05 & OR_meta > 1 ~ "Resistance",
      FDR < 0.05 & OR_meta < 1 ~ "Response",
      TRUE ~ "Not significant"
    ),
    CTP = reorder(CTP, log2OR)
  )

pA <- ggplot(plot.meta, aes(x = log2OR, y = CTP)) +
  geom_vline(xintercept = 0, linetype = 2,
             color = "grey55", linewidth = 0.45) +
  geom_segment(
    aes(x = log2(CI_low), xend = log2(CI_high),
        yend = CTP, color = significance),
    linewidth = 0.65
  ) +
  geom_point(aes(color = significance), size = 2.4) +
  scale_color_manual(values = c(
    "Resistance" = "#B2182B",
    "Response" = "#2166AC",
    "Not significant" = "grey70"
  )) +
  labs(
    x = expression(log[2]*"(OR for non-response per 1-SD increase)"),
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(pA)
ggsave("CTP_Meta_effect_landscape.pdf", pA, width = 3.1, height = 2.4)

# 9. Figure B：CTP_5与CTP_10森林图 =========================================
key.ctps <- c("CTP_5", "CTP_10")

forest.cohort <- logistic.result |>
  filter(CTP %in% key.ctps, is.finite(OR)) |>
  transmute(
    CTP, Dataset = as.character(Group),
    OR, CI_low, CI_high, Type = "Cohort"
  )

forest.meta <- meta.result |>
  filter(CTP %in% key.ctps) |>
  transmute(
    CTP, Dataset = "Meta-analysis",
    OR = OR_meta, CI_low, CI_high, Type = "Meta"
  )

dataset.order <- c(unique(as.character(logistic.result$Group)),
                   "Meta-analysis")

forest.plot.df <- bind_rows(forest.cohort, forest.meta) |>
  mutate(Dataset = factor(Dataset, levels = rev(dataset.order)))

pB <- ggplot(forest.plot.df, aes(x = OR, y = Dataset)) +
  geom_vline(xintercept = 1, linetype = 2,
             color = "grey55", linewidth = 0.45) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, yend = Dataset),
    linewidth = 0.65
  ) +
  geom_point(aes(shape = Type), size = 2.4) +
  scale_shape_manual(values = c(Cohort = 16, Meta = 18)) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  facet_wrap(~ CTP, nrow = 1) +
  labs(
    x = "Odds ratio for ICB non-response per 1-SD increase",
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(pB)
ggsave("CTP5_CTP10_forest.pdf", pB, width = 6, height = 3.5)

# 10. Figure C：CTP_5全集小提琴图 + Meta P值 ===============================
distribution.df <- ctp.analysis |>
  filter(CTP == "CTP_5", is.finite(score)) |>
  mutate(
    ORR = factor(
      ORR,
      levels = c("Responder", "Non_responder"),
      labels = c("Responder", "Non-responder")
    )
  )

p_value <- meta.result |>
  filter(CTP == "CTP_5") |>
  pull(p_meta)

y_max <- max(distribution.df$score)
h <- diff(range(distribution.df$score)) * 0.09
y_pos <- y_max + 2 * h

pC <- ggplot(distribution.df, aes(x = ORR, y = score, fill = ORR)) +
  geom_violin(width = 0.85, trim = FALSE,
              alpha = 0.35, linewidth = 0.35) +
  geom_boxplot(width = 0.22, outlier.shape = NA,
               alpha = 0.8, linewidth = 0.45) +
  geom_jitter(width = 0.08, size = 0.7, alpha = 0.35) +
  annotate("segment", x = 1, xend = 2,
           y = y_pos, yend = y_pos, linewidth = 0.4) +
  annotate("segment", x = c(1, 2), xend = c(1, 2),
           y = y_pos, yend = y_pos - h / 2, linewidth = 0.4) +
  annotate("text", x = 1.5, y = y_pos + h,
           label = sprintf("P = %.2g", p_value), size = 4) +
  scale_fill_manual(values = c(
    "Responder" = "#2166AC",
    "Non-responder" = "#B2182B"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "CTP_5 activity (Z-score)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(linewidth = 0.4)
  )

print(pC)
ggsave("CTP5_violin_all_patients_P.pdf", pC,
       width = 3.5, height = 4)

# 11. Figure D：跨队列效应一致性热图 =======================================
heat.df <- logistic.result |>
  filter(is.finite(beta)) |>
  left_join(
    meta.result |> select(CTP, beta_meta),
    by = "CTP"
  ) |>
  mutate(CTP = reorder(CTP, beta_meta))

max.beta <- max(abs(heat.df$beta), na.rm = TRUE)

pD <- ggplot(heat.df, aes(x = Group, y = CTP, fill = beta)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-max.beta, max.beta),
    name = "log(OR)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(pD)
ggsave("CTP_cohort_effect_heatmap.pdf", pD,
       width = 2, height = 2.5)
