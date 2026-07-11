library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
# Project aggregated CTP signatures to bulk RNA-seq ============================

# 将前面得到的每个CTP聚合基因权重（agg_weights）
# 投影回452例bulk RNA-seq样本，
# 计算每个样本的CTP活性分数，
# 并比较Responder与Non-responder之间的差异


meta <- read.csv(
  "D:/HCC_ICB/Total_metadata_452例.csv",
  check.names = FALSE
)
meta <- meta[, -c(1, 5)]
expr <- read.csv(
  "D:/HCC_ICB/HCC_ICB_altas_452cases_log2TPM.csv",
  row.names = 1,
  check.names = FALSE
)

common_genes <- intersect(
  rownames(expr),
  rownames(agg_weights)
)

length(common_genes)

expr.use <- as.matrix(
  expr[common_genes, ]
)

weight.use <- as.matrix(
  agg_weights[common_genes, ]
)

#把每个样本的基因表达转成样本内百分位排名
#不再直接使用不同样本之间的绝对表达值，而是看某个CTP基因在每个样本内部是否相对高表达
#对测序深度和总体表达尺度更稳健；
#对不同cohort之间的技术差异不那么敏感；
#与ssGSEA、AUCell、rank-based scoring的思路接近。

expr.rank <- apply(
  expr.use,
  2,
  function(x) {
    rank(
      x,
      ties.method = "average"
    ) / length(x)
  }
)
dim(expr.rank)

#计算每个样本的CTP score
ctp.score <- matrix(
  NA,
  nrow = ncol(expr.rank),
  ncol = ncol(weight.use)
)

rownames(ctp.score) <- colnames(expr.rank)
colnames(ctp.score) <- colnames(weight.use)

for (ctp in colnames(weight.use)) {
  
  genes <- rownames(weight.use)
  
  # 取对应权重
  weights <- weight.use[
    genes,
    ctp
  ]
  
  # 权重归一化
  weights <- weights / sum(weights)
  
  # 每个样本的加权 CTP score
  ctp.score[, ctp] <- colSums(
    expr.rank[genes, ] * weights
  )
}



ctp.score.df <- as.data.frame(
  ctp.score
)

ctp.score.df$sample <- rownames(
  ctp.score.df
)

ctp.meta <- meta %>%
  inner_join(
    ctp.score.df,
    by = "sample"
  )


ctp.names <- colnames(
  agg_weights
)
#在每个Group内部进行Z-score标准化
#对每一个CTP，在每个Group内部单独做均值为0、标准差为1的标准化，避免不同Group之间的批次效应影响CTP score的比较
ctp.meta.z <- ctp.meta %>%
  group_by(Group) %>%
  mutate(
    across(
      all_of(ctp.names),
      ~ as.numeric(scale(.x))
    )
  ) %>%
  ungroup()

ctp.long <- ctp.meta.z %>%
  pivot_longer(
    cols = all_of(ctp.names),
    names_to = "CTP",
    values_to = "score"
  )

orr.result <- ctp.long %>%
  group_by(CTP) %>%
  summarise(
    median_Responder = median(
      score[ORR == "Responder"],
      na.rm = TRUE
    ),
    
    median_Non_responder = median(
      score[ORR == "Non_responder"],
      na.rm = TRUE
    ),
    
    difference =
      median_Responder -
      median_Non_responder,
    
    pvalue = wilcox.test(
      score ~ ORR,
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(
      pvalue,
      method = "BH"
    )
  ) %>%
  arrange(FDR)


orr.result


ggplot(
  ctp.long,
  aes(
    x = ORR,
    y = score,
    fill = ORR
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 0.6,
    alpha = 0.4
  ) +
  facet_wrap(
    ~ CTP,
    scales = "free_y",
    ncol = 4
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  ) +
  labs(
    x = NULL,
    y = "CTP score",
    title = "CTP activity in Responder and Non-responder"
  )


top_genes[["CTP_10"]]
top_genes[["CTP_25"]]
top_genes[["CTP_5"]]
top_genes[["CTP_14"]]
