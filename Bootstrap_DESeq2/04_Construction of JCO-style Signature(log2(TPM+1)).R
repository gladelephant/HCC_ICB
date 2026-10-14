expr82 <- log2(tpm82 + 1)
expr35 <- log2(tpm35 + 1)
#####================================================================================

library(pROC)
library(ggplot2)
library(dplyr)

set.seed(123)

dir.create("JCO_Resistance_TPM", showWarnings = FALSE)

# ============================================================
# 1. 共同基因、样本匹配、耐药基因排名
# ============================================================

common_genes <- intersect(rownames(expr82), rownames(expr35))
expr82 <- expr82[common_genes, , drop = FALSE]
expr35 <- expr35[common_genes, , drop = FALSE]

meta82 <- meta82[match(colnames(expr82), meta82$sample), ]
meta35 <- meta35[match(colnames(expr35), meta35$sample), ]

y82 <- as.integer(meta82$ORR == "Non-responder")
y35 <- as.integer(meta35$ORR == "Non-responder")

rank_df <- resistance_genes %>%
  filter(
    gene %in% common_genes,
    NR_frequency > R_frequency,
    is.finite(NR_mean_log2FC),
    NR_mean_log2FC > 0
  ) %>%
  arrange(desc(NR_frequency), gene)

max_genes <- min(200, nrow(rank_df))
if (max_genes < 3) stop("可用耐药基因不足3个")

ranked_genes <- rank_df$gene[1:max_genes]
weights <- rank_df$NR_mean_log2FC[1:max_genes]
gene_numbers <- 2:max_genes

# ============================================================
# 2. Top 2-200累积加权评分
# ============================================================

compute_all_scores <- function(expr, genes, weights) {
  
  X <- expr[genes, , drop = FALSE]
  WX <- sweep(X, 1, weights, "*")
  
  cumulative <- apply(WX, 2, cumsum)
  
  cumulative <- matrix(
    cumulative,
    nrow = length(genes),
    ncol = ncol(expr)
  )
  
  scores <- sweep(
    cumulative, 1, cumsum(weights), "/"
  )
  
  scores <- t(scores)
  
  rownames(scores) <- colnames(expr)
  colnames(scores) <- paste0("Top", seq_along(genes))
  
  scores
}

score82_all <- compute_all_scores(
  expr82, ranked_genes, weights
)[, paste0("Top", gene_numbers), drop = FALSE]

score35_all <- compute_all_scores(
  expr35, ranked_genes, weights
)[, paste0("Top", gene_numbers), drop = FALSE]

# ============================================================
# 3. Bootstrap AUROC
# ============================================================

bootstrap_auc <- function(scores, y, B = 2000) {
  
  idx0 <- which(y == 0)
  idx1 <- which(y == 1)
  
  boot_auc <- matrix(
    NA_real_,
    nrow = B,
    ncol = ncol(scores)
  )
  
  for (b in seq_len(B)) {
    
    ix <- c(
      sample(idx0, length(idx0), replace = TRUE),
      sample(idx1, length(idx1), replace = TRUE)
    )
    
    for (j in seq_len(ncol(scores))) {
      
      boot_auc[b, j] <- as.numeric(
        auc(
          roc(
            y[ix],
            scores[ix, j],
            levels = c(0, 1),
            direction = "<",
            quiet = TRUE
          )
        )
      )
    }
  }
  
  q <- apply(
    boot_auc, 2, quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  
  list(
    summary = data.frame(
      low = q[1, ],
      median = q[2, ],
      high = q[3, ],
      n_genes = gene_numbers
    ),
    bootstrap = boot_auc
  )
}

boot82 <- bootstrap_auc(score82_all, y82)
boot35 <- bootstrap_auc(score35_all, y35)

bootstrap_conf_df <- bind_rows(
  mutate(boot82$summary, data = "Train"),
  mutate(boot35$summary, data = "Validation")
) %>%
  select(data, low, median, high, n_genes)

write.csv(
  bootstrap_conf_df,
  "JCO_Resistance_TPM/bootstrap_conf_df.csv",
  row.names = FALSE
)

saveRDS(
  list(train = boot82, validation = boot35),
  "JCO_Resistance_TPM/bootstrap_AUROC.rds"
)

# ============================================================
# 4. JCO式AUROC曲线
# ============================================================

p1 <- ggplot(
  bootstrap_conf_df,
  aes(n_genes, median, color = data)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 0.7) +
  scale_color_manual(
    values = c(
      Train = "#2471A3",
      Validation = "#C0392B"
    )
  ) +
  labs(
    x = "Number of genes",
    y = "Median AUROC across bootstraps",
    color = "Cohort"
  ) +
  theme_classic(base_size = 13)

print(p1)

ggsave(
  "JCO_Resistance_TPM/JCO_AUROC_Curve.pdf",
  p1, width = 9, height = 5.5
)

# ============================================================
# 5. 确定Signature规模
# ============================================================

valid <- bootstrap_conf_df %>%
  filter(data == "Validation") %>%
  arrange(n_genes)

# A. 最大Bootstrap中位AUROC
max_auc <- max(valid$median)

max_k <- min(
  valid$n_genes[valid$median == max_auc]
)

# B. 99%性能平台期：最少基因达到最高AUROC的99%
plateau_k <- min(
  valid$n_genes[
    valid$median >= 0.99 * max_auc
  ]
)

# C. 平滑曲线的几何Elbow
smooth_y <- predict(
  smooth.spline(
    valid$n_genes,
    valid$median,
    spar = 0.65
  ),
  x = valid$n_genes
)$y

x <- valid$n_genes
xn <- (x - min(x)) / (max(x) - min(x))
yn <- (smooth_y - min(smooth_y)) /
  (max(smooth_y) - min(smooth_y))

baseline <- yn[1] + (yn[length(yn)] - yn[1]) * xn

elbow_k <- x[which.max(yn - baseline)]

# 三种候选规模
candidate_sizes <- data.frame(
  Method = c(
    "Maximum AUROC",
    "99% plateau",
    "Geometric elbow"
  ),
  n_genes = c(max_k, plateau_k, elbow_k)
)

candidate_sizes$AUROC <- valid$median[
  match(candidate_sizes$n_genes, valid$n_genes)
]

candidate_sizes$Delta_from_max <- max_auc -
  candidate_sizes$AUROC

print(candidate_sizes)

write.csv(
  candidate_sizes,
  "JCO_Resistance_TPM/Signature_Size_Candidates.csv",
  row.names = FALSE
)

# ============================================================
# 6. 标注候选规模
# ============================================================

p2 <- ggplot(valid, aes(n_genes, median)) +
  geom_line(color = "#2471A3", linewidth = 0.9) +
  geom_line(
    aes(y = smooth_y),
    color = "#E67E22",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0.99 * max_auc,
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = plateau_k,
    color = "#C0392B",
    linetype = "dashed"
  ) +
  geom_point(
    data = candidate_sizes,
    aes(x = n_genes, y = AUROC),
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = "Resistance Signature Size Selection",
    subtitle = paste0(
      "Maximum: ", max_k,
      " | 99% plateau: ", plateau_k,
      " | Elbow: ", elbow_k
    ),
    x = "Number of genes",
    y = "Bootstrap median AUROC"
  ) +
  theme_classic(base_size = 13)

print(p2)

ggsave(
  "JCO_Resistance_TPM/Signature_Size_Selection.pdf",
  p2, width = 9, height = 5.5
)

# ============================================================
# 7. Bootstrap基因规模选择稳定性
# ============================================================

selected_k_boot <- apply(
  boot35$bootstrap,
  1,
  function(a) {
    min(gene_numbers[a >= 0.99 * max(a)])
  }
)

selection_frequency <- as.data.frame(
  table(selected_k_boot)
)

colnames(selection_frequency) <- c(
  "n_genes", "n_boot"
)

selection_frequency$n_genes <- as.integer(
  as.character(selection_frequency$n_genes)
)

selection_frequency$frequency <-
  selection_frequency$n_boot / length(selected_k_boot)

selection_frequency <- selection_frequency %>%
  arrange(desc(frequency))

write.csv(
  selection_frequency,
  "JCO_Resistance_TPM/Signature_Size_Stability.csv",
  row.names = FALSE
)

p3 <- ggplot(
  selection_frequency,
  aes(n_genes, frequency)
) +
  geom_col(fill = "#2471A3") +
  labs(
    x = "Selected number of genes",
    y = "Bootstrap selection frequency"
  ) +
  theme_classic(base_size = 13)

print(p3)

ggsave(
  "JCO_Resistance_TPM/Signature_Size_Stability.pdf",
  p3, width = 9, height = 5.5
)

# ============================================================
# 8. 保存候选模型
# ============================================================

# 暂时采用99%平台期规则
selected_k <- plateau_k

final_signature <- rank_df[1:selected_k, ]

final_model <- list(
  genes = final_signature$gene,
  weights = final_signature$NR_mean_log2FC,
  n_genes = selected_k,
  expression = "log2(TPM+1)",
  direction = "Non-responder high",
  selection_method = "99% bootstrap median AUROC plateau"
)

write.csv(
  final_signature,
  "JCO_Resistance_TPM/Final_Signature.csv",
  row.names = FALSE
)

saveRDS(
  final_model,
  "JCO_Resistance_TPM/Final_Model.rds"
)

print(candidate_sizes)
print(head(selection_frequency, 15))
print(final_signature)


