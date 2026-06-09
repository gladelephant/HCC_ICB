library(readxl)
df <- read_excel("GSE181946_count_HCC免疫治疗.xlsx")
View(df)

install.packages(
  "org.Hs.eg.db_3.23.1.tar.gz",
  repos = NULL,
  type = "source"
)

df$GeneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = df$Tag,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

##保留蛋白编码基因
library("qs")
gtf <- qread("G:\\HCC_靶免治疗\\bulk_in_house\\Homo_sapiens_gtf")
proteincoding <- gtf %>% filter(gene_biotype == "protein_coding")
# 剔除 RPS/RPL 核糖体蛋白基因
genes.use <- grep(
  pattern = "^(RPS|RPL)",
  x = df$GeneSymbol,
  value = TRUE,
  invert = TRUE
)

# 再保留 protein coding genes
genename_filt <- genes.use[genes.use %in% proteincoding$gene_name]

# 过滤 count matrix
df_filt <- df[
  !is.na(df$GeneSymbol) &
    df$GeneSymbol %in% genename_filt,
]

## 2. 读取 count 矩阵
count <- df_filt
count$Tag <-NULL
## 假设第一列是 gene id
## 重复基因保留表达counts最大的行
count2 <- count %>%
  mutate(total_count = rowSums(across(where(is.numeric)))) %>%
  group_by(GeneSymbol) %>%
  slice_max(order_by = total_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()
rownames(count2) <-count2$GeneSymbol
count2$total_count<-NULL
count2$GeneSymbol<-NULL

gset = fread("GSE181946_series_matrix.txt",data.table = F,sep = "\t",skip=30)

meta <- data.frame(
  GSM = c(
    "GSM5514677","GSM5514678","GSM5514679","GSM5514680",
    "GSM5514681","GSM5514682","GSM5514683","GSM5514684",
    "GSM5514685","GSM5514686","GSM5514687","GSM5514688",
    "GSM5514689","GSM5514690","GSM5514691","GSM5514692",
    "GSM5514693"
  ),
  Sample = c(
    "327","332","333","334","335",
    "337","338","343","358","360",
    "342","355","251","445",
    "480","500","605"
  ),
  Response = c(
    rep("PD", 10),
    rep("PR", 7)
  ),
  stringsAsFactors = FALSE
)


count <- count2

## 6. 过滤低表达基因
keep <- rowSums(count >= 10) >= 5

count_filt <- count[keep, ]

## 7. 设置分组因子
table(meta$Response)
meta$Response <- factor(meta$Response, levels = c("PR", "PD"))

## 8. 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = count_filt,
  colData = meta,
  design = ~ Response
)

## 9. 差异分析
dds <- DESeq(dds)

## Resistant vs Sensitive
res <- results(dds, contrast = c("Response", "PR", "PD"))

## 10. 整理结果
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df %>%
  arrange(padj) %>%
  mutate(
    change = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Stable"
    )
  )

write.csv(res_df, "DESeq2_DEG_results.csv", row.names = FALSE)
## 11. 提取显著差异基因
deg <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

write.csv(deg, "DEG_padj0.1_log2FC0.5.csv", row.names = FALSE)

## 8. Save analysis objects -------------------------------------------------

saveRDS(
  list(
    expression_matrix = count,
    filter_matrix = count_filt,
    metadata = meta,
    deg_all = res_df,
    deg_sig = deg
  ),
  file = "GSE181946_HCC免疫治疗.rds")

