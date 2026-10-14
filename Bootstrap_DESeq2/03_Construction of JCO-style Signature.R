
# ============================================================
# HCC-ICB | Resistance-only JCO-style Signature
# Training: 82 patients
# Internal validation: 35 patients
# Bootstrap: 500 completed iterations
#
# NOTE:
# This is a VST-based adaptation.
# It is NOT an exact replication of JCO's log2(TPM+1).
# ============================================================

library(DESeq2)
library(pROC)
library(ggplot2)

set.seed(123)

dir.create("JCO_Resistance_Validation",
           showWarnings = FALSE)

outdir <- "JCO_Resistance_Validation"

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

count82 <- read.csv(
  "Inhouse_82例_RNA_count_matrix_filt.csv",
  row.names = 1,
  check.names = FALSE
)

count35 <- read.csv(
  "Inhouse_35例_RNA_count_matrix_filt.csv",
  row.names = 1,
  check.names = FALSE
)

meta82 <- read.csv(
  "RNA_82sample_metadata.csv",
  check.names = FALSE
)

meta35 <- read.csv(
  "RNA_35sample_metadata.csv",
  check.names = FALSE
)

colnames(meta82)[1:3] <- c(
  "sample", "ORR", "imaging"
)

colnames(meta35)[1:3] <- c(
  "sample", "ORR", "imaging"
)

res_list <- readRDS(
  "JCO_82_bootstrap_results.rds"
)

# ------------------------------------------------------------
# 2. Data checks
# ------------------------------------------------------------

prepare_cohort <- function(count, meta) {
  
  count <- as.matrix(count)
  
  stopifnot(
    !anyDuplicated(rownames(count)),
    !anyDuplicated(colnames(count)),
    !anyDuplicated(meta$sample),
    all(colnames(count) %in% meta$sample)
  )
  
  storage.mode(count) <- "numeric"
  
  stopifnot(
    all(is.finite(count)),
    all(count >= 0),
    all(abs(count - round(count)) < 1e-8)
  )
  
  meta <- meta[
    match(colnames(count), meta$sample),
    ,
    drop = FALSE
  ]
  
  rownames(meta) <- colnames(count)
  
  meta$ORR <- factor(
    meta$ORR,
    levels = c("Responder", "Non-responder")
  )
  
  stopifnot(
    identical(rownames(meta), colnames(count)),
    !anyNA(meta$ORR),
    all(table(meta$ORR) > 0)
  )
  
  list(count = count, meta = meta)
}

d82 <- prepare_cohort(count82, meta82)
d35 <- prepare_cohort(count35, meta35)

count82 <- d82$count
count35 <- d35$count

meta82 <- d82$meta
meta35 <- d35$meta

cat("Training ORR:\n")
print(table(meta82$ORR))

cat("Validation ORR:\n")
print(table(meta35$ORR))

# ------------------------------------------------------------
# 3. Common gene universe
# ------------------------------------------------------------

common_genes <- intersect(
  rownames(count82),
  rownames(count35)
)

if (length(common_genes) < 1000) {
  stop("Common genes too few. Check gene IDs.")
}

count82 <- count82[common_genes, , drop = FALSE]
count35 <- count35[common_genes, , drop = FALSE]

cat("Common genes:", length(common_genes), "\n")

# ------------------------------------------------------------
# 4. Direction-specific Bootstrap ranking
#    Only genes with positive shrunken log2FC
#    Non-responder > Responder
# ------------------------------------------------------------

n_success <- sum(
  !vapply(res_list, is.null, logical(1))
)

stopifnot(n_success > 0)

n_NR <- setNames(
  integer(length(common_genes)),
  common_genes
)

n_R <- n_NR

sum_NR <- setNames(
  numeric(length(common_genes)),
  common_genes
)

for (beta in res_list) {
  
  if (is.null(beta)) next
  
  stopifnot(
    is.numeric(beta),
    !is.null(names(beta))
  )
  
  beta <- beta[
    names(beta) %in% common_genes &
      is.finite(beta)
  ]
  
  pos <- beta[beta > 0]
  neg <- beta[beta < 0]
  
  if (length(pos) > 0) {
    g <- names(pos)
    n_NR[g] <- n_NR[g] + 1L
    sum_NR[g] <- sum_NR[g] + pos
  }
  
  if (length(neg) > 0) {
    g <- names(neg)
    n_R[g] <- n_R[g] + 1L
  }
}

resistance_rank <- data.frame(
  gene = common_genes,
  NR_n = as.integer(n_NR),
  NR_frequency = as.numeric(n_NR / n_success),
  R_n = as.integer(n_R),
  R_frequency = as.numeric(n_R / n_success),
  mean_log2FC = as.numeric(
    ifelse(
      n_NR > 0,
      sum_NR / pmax(n_NR, 1),
      NA_real_
    )
  ),
  stringsAsFactors = FALSE
)

# Primary rule: NR-up frequency
# Exclude genes whose opposite-direction
# frequency is >= NR-up frequency.

resistance_rank <- subset(
  resistance_rank,
  NR_n > 0 &
    NR_frequency > R_frequency &
    is.finite(mean_log2FC) &
    mean_log2FC > 0
)

resistance_rank <- resistance_rank[
  order(
    -resistance_rank$NR_frequency,
    resistance_rank$gene
  ),
  ,
  drop = FALSE
]

rownames(resistance_rank) <- NULL
resistance_rank$rank <- seq_len(
  nrow(resistance_rank)
)

resistance_rank <- resistance_rank[
  ,
  c(
    "rank", "gene",
    "NR_n", "NR_frequency",
    "R_n", "R_frequency",
    "mean_log2FC"
  )
]

write.csv(
  resistance_rank,
  file.path(outdir, "Resistance_Gene_Rank.csv"),
  row.names = FALSE
)

cat("Resistance candidate genes:",
    nrow(resistance_rank), "\n")

print(head(resistance_rank, 20))

# ------------------------------------------------------------
# 5. VST expression
#
# Each cohort is normalized independently
# WITHOUT using ORR labels.
#
# VST is a practical exploratory substitute,
# not identical to log2(TPM+1).
# ------------------------------------------------------------

get_vst <- function(count, meta) {
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(count),
    colData = meta,
    design = ~ 1
  )
  
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  
  dds <- estimateSizeFactors(dds)
  
  vst_obj <- varianceStabilizingTransformation(
    dds,
    blind = TRUE
  )
  
  assay(vst_obj)
}

expr82 <- get_vst(count82, meta82)
expr35 <- get_vst(count35, meta35)

# ------------------------------------------------------------
# 6. Keep genes available in both expression matrices
# ------------------------------------------------------------

available <- intersect(
  rownames(expr82),
  rownames(expr35)
)

resistance_rank <- resistance_rank[
  resistance_rank$gene %in% available,
  ,
  drop = FALSE
]

rownames(resistance_rank) <- NULL
resistance_rank$rank <- seq_len(
  nrow(resistance_rank)
)

# ------------------------------------------------------------
# 7. Build cumulative weighted signatures
#
# Score = sum(weight * expression) / sum(weight)
#
# All weights are positive because this is
# resistance-only.
# Higher score = higher predicted non-response.
# ------------------------------------------------------------

max_genes <- min(
  200,
  nrow(resistance_rank)
)

if (max_genes < 2) {
  stop("Fewer than two resistance genes available.")
}

gene_numbers <- 2:max_genes

ranked_genes <- resistance_rank$gene[
  seq_len(max_genes)
]

weights <- resistance_rank$mean_log2FC[
  seq_len(max_genes)
]

stopifnot(all(weights > 0))

weighted_scores <- function(expr, genes, weights) {
  
  X <- expr[genes, , drop = FALSE]
  
  as.numeric(
    colSums(
      sweep(X, 1, weights, "*")
    ) / sum(weights)
  )
}

# ------------------------------------------------------------
# 8. Evaluate AUROC for each signature size
# ------------------------------------------------------------

y35 <- as.integer(
  meta35$ORR == "Non-responder"
)

results_list <- vector(
  "list",
  length(gene_numbers)
)

score35_all <- matrix(
  NA_real_,
  nrow = ncol(expr35),
  ncol = length(gene_numbers),
  dimnames = list(
    colnames(expr35),
    paste0("Top", gene_numbers)
  )
)

for (i in seq_along(gene_numbers)) {
  
  k <- gene_numbers[i]
  
  genes_k <- ranked_genes[seq_len(k)]
  weights_k <- weights[seq_len(k)]
  
  score35 <- weighted_scores(
    expr35,
    genes_k,
    weights_k
  )
  
  score35_all[, i] <- score35
  
  roc_obj <- roc(
    response = y35,
    predictor = score35,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  results_list[[i]] <- data.frame(
    n_genes = k,
    AUROC = as.numeric(auc(roc_obj))
  )
}

auc_curve <- do.call(
  rbind,
  results_list
)

# ------------------------------------------------------------
# 9. Select signature size
#
# Rule:
# Smallest model within 0.01 AUROC
# of the observed maximum.
#
# This is an explicit parsimony rule,
# NOT a literal implementation of JCO's
# visually determined elbow method.
# ------------------------------------------------------------

best_auc <- max(auc_curve$AUROC)

auc_tolerance <- 0.01

candidate_sizes <- auc_curve$n_genes[
  auc_curve$AUROC >=
    best_auc - auc_tolerance
]

selected_k <- min(candidate_sizes)

selected_auc <- auc_curve$AUROC[
  auc_curve$n_genes == selected_k
]

cat("\nMaximum validation AUROC:",
    round(best_auc, 4), "\n")

cat("Selected gene number:",
    selected_k, "\n")

cat("Selected validation AUROC:",
    round(selected_auc, 4), "\n")

# ------------------------------------------------------------
# 10. Plot AUROC versus gene number
# ------------------------------------------------------------

p <- ggplot(
  auc_curve,
  aes(x = n_genes, y = AUROC)
) +
  geom_line(
    color = "#2471A3",
    linewidth = 0.8
  ) +
  geom_point(
    color = "#2471A3",
    size = 1
  ) +
  geom_vline(
    xintercept = selected_k,
    linetype = "dashed",
    color = "#C0392B"
  ) +
  geom_point(
    data = subset(
      auc_curve,
      n_genes == selected_k
    ),
    color = "#C0392B",
    size = 3
  ) +
  labs(
    title = "Resistance Signature: Internal Validation",
    subtitle = paste0(
      "Selected Top ", selected_k,
      " genes | AUROC = ",
      round(selected_auc, 3)
    ),
    x = "Number of resistance genes",
    y = "AUROC"
  ) +
  theme_classic(base_size = 13)

print(p)

ggsave(
  file.path(outdir, "AUROC_vs_Gene_Number.pdf"),
  p,
  width = 8,
  height = 5
)

write.csv(
  auc_curve,
  file.path(outdir, "AUROC_by_Gene_Number.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Freeze final signature
# ------------------------------------------------------------

final_signature <- resistance_rank[
  seq_len(selected_k),
  ,
  drop = FALSE
]

write.csv(
  final_signature,
  file.path(outdir, "Final_Resistance_Signature.csv"),
  row.names = FALSE
)

final_genes <- final_signature$gene
final_weights <- final_signature$mean_log2FC

# ------------------------------------------------------------
# 12. Apply fixed signature to 82 and 35
# ------------------------------------------------------------

score82 <- weighted_scores(
  expr82,
  final_genes,
  final_weights
)

score35 <- weighted_scores(
  expr35,
  final_genes,
  final_weights
)

evaluate_roc <- function(meta, score) {
  
  y <- as.integer(
    meta$ORR == "Non-responder"
  )
  
  roc_obj <- roc(
    y,
    score,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  ci <- ci.auc(
    roc_obj,
    method = "bootstrap",
    boot.n = 2000,
    boot.stratified = TRUE
  )
  
  list(
    roc = roc_obj,
    auc = as.numeric(auc(roc_obj)),
    ci_low = as.numeric(ci[1]),
    ci_high = as.numeric(ci[3])
  )
}

roc82 <- evaluate_roc(meta82, score82)
roc35 <- evaluate_roc(meta35, score35)

performance <- data.frame(
  Cohort = c("Training_82", "Validation_35"),
  N = c(nrow(meta82), nrow(meta35)),
  AUROC = c(roc82$auc, roc35$auc),
  CI_low = c(roc82$ci_low, roc35$ci_low),
  CI_high = c(roc82$ci_high, roc35$ci_high)
)

print(performance)

write.csv(
  performance,
  file.path(outdir, "Final_Model_Performance.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Save patient-level scores
# ------------------------------------------------------------

patient_scores <- rbind(
  data.frame(
    cohort = "Training_82",
    sample = meta82$sample,
    ORR = as.character(meta82$ORR),
    Resistance_Score = score82
  ),
  data.frame(
    cohort = "Validation_35",
    sample = meta35$sample,
    ORR = as.character(meta35$ORR),
    Resistance_Score = score35
  )
)

write.csv(
  patient_scores,
  file.path(outdir, "Patient_Resistance_Scores.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    genes = final_genes,
    weights = final_weights,
    n_genes = selected_k,
    auc_tolerance = auc_tolerance,
    expression_method = "VST",
    score_formula = "weighted mean",
    resistance_rank = resistance_rank
  ),
  file.path(outdir, "Final_Resistance_Model.rds")
)

cat("\nAnalysis completed.\n")
cat("Selected genes:", selected_k, "\n")
print(performance)



auc_curve[
  auc_curve$n_genes %in%
    c(5, 8, 10, 15, 20, 50, 100, 200),
]

final_signature[
  ,
  c(
    "rank",
    "gene",
    "NR_n",
    "NR_frequency",
    "R_n",
    "R_frequency",
    "mean_log2FC"
  )
]

library(ggplot2)

df35 <- data.frame(
  ORR = meta35$ORR,
  Resistance_Score = score35
)

ggplot(df35, aes(x = ORR, y = Resistance_Score,
                 fill = ORR)) +
  geom_boxplot(width = 0.5, alpha = 0.3,
               outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2) +
  theme_classic(base_size = 13) +
  labs(
    title = "Top 8 Resistance Signature",
    x = NULL,
    y = "Resistance Score (VST)"
  ) +
  guides(fill = "none")
