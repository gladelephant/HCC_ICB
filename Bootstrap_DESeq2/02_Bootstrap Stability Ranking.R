# ============================================================
# Direction-specific Bootstrap Stability Ranking
# NR-up: 耐药相关候选基因
# R-up : 敏感相关候选基因
# ============================================================

# 如果重新打开RStudio，可先读取：
# res_list <- readRDS("JCO_82_bootstrap_results.rds")

n_success <- sum(
  !vapply(res_list, is.null, logical(1))
)

stopifnot(n_success > 0)

genes_all <- rownames(count82)

# 1. 分别统计正向和负向入选次数
n_NR <- setNames(integer(length(genes_all)), genes_all)
n_R  <- setNames(integer(length(genes_all)), genes_all)

sum_NR <- setNames(numeric(length(genes_all)), genes_all)
sum_R  <- setNames(numeric(length(genes_all)), genes_all)

for (beta in res_list) {
  
  if (is.null(beta)) next
  
  # NR高表达：耐药相关
  beta_NR <- beta[beta > 0]
  
  # R高表达：敏感相关
  beta_R <- beta[beta < 0]
  
  if (length(beta_NR) > 0) {
    g <- names(beta_NR)
    
    n_NR[g] <- n_NR[g] + 1L
    sum_NR[g] <- sum_NR[g] + beta_NR
  }
  
  if (length(beta_R) > 0) {
    g <- names(beta_R)
    
    n_R[g] <- n_R[g] + 1L
    sum_R[g] <- sum_R[g] + beta_R
  }
}

# 2. 生成方向性总表
direction_stability <- data.frame(
  gene = genes_all,
  
  NR_n = as.integer(n_NR),
  NR_frequency = as.numeric(n_NR / n_success),
  NR_mean_log2FC = as.numeric(
    ifelse(n_NR > 0, sum_NR / pmax(n_NR, 1), NA_real_)
  ),
  
  R_n = as.integer(n_R),
  R_frequency = as.numeric(n_R / n_success),
  R_mean_log2FC = as.numeric(
    ifelse(n_R > 0, sum_R / pmax(n_R, 1), NA_real_)
  ),
  
  stringsAsFactors = FALSE
)

# 3. 耐药相关基因：NR高表达
resistance_genes <- subset(
  direction_stability,
  NR_n > 0
)

resistance_genes <- resistance_genes[
  order(-resistance_genes$NR_frequency,
        resistance_genes$gene),
]

rownames(resistance_genes) <- NULL
resistance_genes$rank <- seq_len(nrow(resistance_genes))

resistance_genes <- resistance_genes[
  , c("rank", "gene",
      "NR_n", "NR_frequency", "NR_mean_log2FC",
      "R_n", "R_frequency")
]

# 4. 敏感相关基因：R高表达
sensitivity_genes <- subset(
  direction_stability,
  R_n > 0
)

sensitivity_genes <- sensitivity_genes[
  order(-sensitivity_genes$R_frequency,
        sensitivity_genes$gene),
]

rownames(sensitivity_genes) <- NULL
sensitivity_genes$rank <- seq_len(nrow(sensitivity_genes))

sensitivity_genes <- sensitivity_genes[
  , c("rank", "gene",
      "R_n", "R_frequency", "R_mean_log2FC",
      "NR_n", "NR_frequency")
]

# 5. 保存
write.csv(
  resistance_genes,
  "JCO_82_Resistance_Genes_Rank.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_genes,
  "JCO_82_Sensitivity_Genes_Rank.csv",
  row.names = FALSE
)

write.csv(
  direction_stability,
  "JCO_82_Direction_Stability_All.csv",
  row.names = FALSE
)

# 6. 查看结果
head(resistance_genes, 20)
head(sensitivity_genes, 20)
