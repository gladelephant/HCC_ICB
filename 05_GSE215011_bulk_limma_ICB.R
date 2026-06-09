matrix<-fread("GSE215011_gene_description_human_samples.txt")
summary(as.numeric(as.matrix(matrix)))

matrix <- subset(
  matrix,
  gene_biotype == "protein_coding"
)

colnames(matrix)

## 重复基因保留表达counts最大的行
count2 <- matrix %>%
  mutate(total_count = rowSums(across(where(is.numeric)))) %>%
  group_by(gene_name) %>%
  slice_max(order_by = total_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()
rownames(count2) <-count2$gene_name

count2 <- count2[, c(
  "A12_949T",
  "A16_557T",
  "A17_152T",
  "A19_16T",
  "A19_171T",
  "A11_385T",
  "A18_440T",
  "A19_1T",
  "A19_121T",
  "A19_174T"
)]

##保留蛋白编码基因
library("qs")
gtf <- qread("G:\\HCC_靶免治疗\\bulk_in_house\\Homo_sapiens_gtf")
proteincoding <- gtf %>% filter(gene_biotype == "protein_coding")

# 剔除 RPS/RPL 核糖体蛋白基因
genes.use <- grep(
  pattern = "^(RPS|RPL)",
  x = rownames(count2),
  value = TRUE,
  invert = TRUE
)

# 再保留 protein coding genes
genename_filt <- genes.use[genes.use %in% proteincoding$gene_name]

# 过滤 count matrix
count2_filt <- count2[rownames(count2) %in% genename_filt, ]

keep <- rowSums(count2_filt >= 1) >= 3

count2_filt <- count2_filt[keep, ]
expr_log2 <- log2(count2_filt + 1)

colnames(expr_log2)

group <- factor(c(
  "NR",
  "NR",
  "NR",
  "NR",
  "R",
  "NR",
  "R",
  "R",
  "R",
  "R"
))

metadata<-data.frame(
  sample = colnames(expr_log2),
  group = group
)

library(limma)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_log2, design)

contrast.matrix <- makeContrasts(
  R_vs_NR = R - NR,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(
  fit2,
  coef = "R_vs_NR",
  number = Inf,
  adjust.method = "BH"
)


deg$gene<-row.names(deg)
getwd()

obj_list <- list(
  original_matrix = count2,
  protein_coding_filtered_matrix = count2_filt,
  deg_list = deg,
  metadata = metadata
)

saveRDS(
  obj_list,
  "GSE215011_HCC_ICB.rds"
)
