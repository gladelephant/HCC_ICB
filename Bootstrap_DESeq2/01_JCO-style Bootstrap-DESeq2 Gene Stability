# ============================================================
# HCC-ICB | JCO-style Bootstrap-DESeq2 Gene Stability
# Training cohort: 82 patients
# Bootstrap: 500 iterations
# Parallel workers: 50
# ============================================================

# 0. CPU线程设置
# 建议在全新R会话中、加载R包之前执行
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

library(DESeq2)
library(ashr)
library(future)
library(future.apply)


# 1. 参数
B <- 500
min_fdr <- 0.1
n_workers <- 50

meta82 <- read.csv("RNA_82sample_metadata.csv", check.names = FALSE)
count82 <- read.csv("Inhouse_82例_RNA_count_matrix_filt.csv",
                    row.names = 1, check.names = FALSE)
colnames(meta82)[1:3] <- c("sample", "ORR", "imaging")
meta82$ORR <- factor(meta82$ORR,
                     levels = c("Responder", "Non-responder"))

# 2. 数据检查与匹配
stopifnot(
  !anyDuplicated(meta82$sample),
  !anyDuplicated(colnames(count82)),
  !anyDuplicated(rownames(count82)),
  all(colnames(count82) %in% meta82$sample)
)

meta82 <- meta82[
  match(colnames(count82), meta82$sample),
  ,
  drop = FALSE
]

rownames(meta82) <- colnames(count82)

meta82$ORR <- factor(
  meta82$ORR,
  levels = c("Responder", "Non-responder")
)

stopifnot(
  identical(rownames(meta82), colnames(count82)),
  !anyNA(meta82$ORR),
  all(table(meta82$ORR) >= 2),
  is.numeric(as.matrix(count82)),
  all(is.finite(as.matrix(count82))),
  all(as.matrix(count82) >= 0)
)

# 3. 单次Bootstrap函数
bootstrap_fun <- function(b) {
  
  tryCatch({
    
    # 按ORR分层、有放回抽样
    ix <- unlist(
      lapply(
        split(seq_len(nrow(meta82)), meta82$ORR),
        function(i) sample(i, length(i), replace = TRUE)
      ),
      use.names = FALSE
    )
    
    cb <- count82[, ix, drop = FALSE]
    mb <- meta82[ix, , drop = FALSE]
    
    # 重复抽中的患者赋予唯一列名
    colnames(cb) <- paste0("boot_", seq_along(ix))
    rownames(mb) <- colnames(cb)
    
    # DESeq2
    dds <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(cb)),
      colData = mb,
      design = ~ ORR
    )
    
    dds <- DESeq(dds, quiet = TRUE)
    
    # Non-responder vs Responder
    res <- results(
      dds,
      contrast = c(
        "ORR",
        "Non-responder",
        "Responder"
      )
    )
    
    # ashr效应量收缩
    shr <- lfcShrink(
      dds,
      contrast = c(
        "ORR",
        "Non-responder",
        "Responder"
      ),
      type = "ashr",
      quiet = TRUE
    )
    
    # 筛选显著差异基因
    selected <- rownames(res)[
      !is.na(res$padj) &
        res$padj < min_fdr
    ]
    
    # 提取收缩log2FC，并保留基因名
    beta <- as.numeric(
      shr[selected, "log2FoldChange"]
    )
    
    names(beta) <- selected
    
    beta <- beta[is.finite(beta)]
    
    return(beta)
    
  }, error = function(e) {
    
    message(
      "Bootstrap ", b, " failed: ",
      conditionMessage(e)
    )
    
    return(NULL)
  })
}

# 4. 50进程并行运行
plan(multisession, workers = n_workers)

set.seed(123)
t0 <- Sys.time()

res_list <- future_lapply(
  seq_len(B),
  bootstrap_fun,
  future.seed = TRUE,
  future.chunk.size = 1
)

t1 <- Sys.time()

plan(sequential)

# 5. 第一时间保存Bootstrap原始结果
saveRDS(
  res_list,
  "JCO_82_bootstrap_results.rds"
)

# 6. 检查成功次数
n_success <- sum(
  !vapply(res_list, is.null, logical(1))
)

if (n_success == 0) {
  stop("所有Bootstrap均失败，请检查错误信息。")
}

if (n_success < B) {
  warning(
    "仅 ", n_success, "/", B,
    " 次Bootstrap成功。建议检查失败原因。"
  )
}

# 7. 汇总基因稳定性与效应量
selection <- setNames(
  integer(nrow(count82)),
  rownames(count82)
)

effect_sum <- setNames(
  numeric(nrow(count82)),
  rownames(count82)
)

effect_n <- setNames(
  integer(nrow(count82)),
  rownames(count82)
)

for (beta in res_list) {
  
  if (is.null(beta)) next
  
  genes <- names(beta)
  
  selection[genes] <- selection[genes] + 1L
  
  effect_sum[genes] <- effect_sum[genes] + beta
  
  effect_n[genes] <- effect_n[genes] + 1L
}

# 8. 生成稳定性结果
stability <- data.frame(
  gene = names(selection),
  selection_n = as.integer(selection),
  frequency = as.numeric(selection / n_success),
  mean_log2FC = as.numeric(
    ifelse(
      effect_n > 0,
      effect_sum / pmax(effect_n, 1),
      NA_real_
    )
  ),
  stringsAsFactors = FALSE
)

# 按Bootstrap入选频率降序排列
# 并列时按基因名称排序，确保可重复
stability <- stability[
  order(
    -stability$frequency,
    stability$gene
  ),
  ,
  drop = FALSE
]

rownames(stability) <- NULL

stability$rank <- seq_len(nrow(stability))

stability <- stability[
  ,
  c(
    "rank",
    "gene",
    "selection_n",
    "frequency",
    "mean_log2FC"
  )
]

# 9. 提取频率>=50%的稳定基因
stable_genes <- subset(
  stability,
  frequency >= 0.50 &
    is.finite(mean_log2FC)
)

# 10. 保存结果
write.csv(
  stability,
  "JCO_82_gene_stability.csv",
  row.names = FALSE
)

write.csv(
  stable_genes,
  "JCO_82_stable_genes_F50.csv",
  row.names = FALSE
)

# 11. 输出运行结果
cat(
  "Bootstrap成功次数：",
  n_success, "/", B, "\n"
)

cat(
  "至少入选一次的基因数：",
  sum(selection > 0), "\n"
)

cat(
  "入选频率>=50%的基因数：",
  nrow(stable_genes), "\n"
)

cat(
  "总运行时间（分钟）：",
  as.numeric(
    difftime(t1, t0, units = "mins")
  ),
  "\n"
)

head(stability, 20)
head(stable_genes, 20)

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


expr82 <- log2(tpm82 + 1)
expr35 <- log2(tpm35 + 1)
