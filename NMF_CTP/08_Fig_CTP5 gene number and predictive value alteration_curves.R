library(dplyr)
library(ggplot2)

top_n_values <- c(
  50, 100, 200, 300, 500,
  800, 1000, 1500, 2000,
  3000, 5000
)

ctp5_topn_result <- data.frame()

for (top_n in top_n_values) {
  
  # 避免 top_n 超过总基因数
  n_use <- min(
    top_n,
    nrow(weight.use)
  )
  
  # CTP_5 按权重排序，取前 n_use 个基因
  genes <- names(
    sort(
      weight.use[, "CTP_5"],
      decreasing = TRUE
    )
  )[1:n_use]
  
  # 提取对应权重
  weights <- weight.use[
    genes,
    "CTP_5"
  ]
  
  # 权重归一化
  weights <- weights / sum(weights)
  
  # 计算每个样本的 CTP_5 score
  score <- colSums(
    expr.rank[
      genes,
      ,
      drop = FALSE
    ] * weights
  )
  
  # 合并 metadata
  score.df <- data.frame(
    sample = names(score),
    score = score
  ) %>%
    inner_join(
      meta,
      by = "sample"
    ) %>%
    group_by(Group) %>%
    mutate(
      score = as.numeric(
        scale(score)
      )
    ) %>%
    ungroup()
  
  # 两组中位数
  median_r <- median(
    score.df$score[
      score.df$ORR == "Responder"
    ],
    na.rm = TRUE
  )
  
  median_nr <- median(
    score.df$score[
      score.df$ORR == "Non_responder"
    ],
    na.rm = TRUE
  )
  
  # Wilcoxon 检验
  p <- wilcox.test(
    score ~ ORR,
    data = score.df,
    exact = FALSE
  )$p.value
  
  ctp5_topn_result <- rbind(
    ctp5_topn_result,
    data.frame(
      TopN = n_use,
      median_Responder = median_r,
      median_Non_responder = median_nr,
      difference = median_r - median_nr,
      pvalue = p,
      neg_log10_p = -log10(p)
    )
  )
}



p1 <- ggplot(
  ctp5_topn_result,
  aes(
    x = TopN,
    y = neg_log10_p
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2
  ) +
  scale_x_log10(
    breaks = top_n_values
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    x = "Number of top CTP_5 genes",
    y = expression(-log[10](P)),
    title = "CTP_5 signature-size sensitivity",
    subtitle = "Dashed line: P = 0.05"
  )

p1
