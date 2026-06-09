data_dir <- "G:\\HCC_靶免治疗\\肝癌免疫治疗数据集\\GSE279750_RAW"
# 获取所有xlsx文件
files <- list.files(
  data_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE
)

# 读取所有文件
data_list <- lapply(files, function(f){
  
  df <- read_excel(f)
  
  # 第一列必须是gene_id
  colnames(df)[1] <- "gene_id"
  
  # 表达列改成文件名
  sample_name <- tools::file_path_sans_ext(basename(f))
  colnames(df)[2] <- sample_name
  
  # 去重复gene
  df <- df %>%
    distinct(gene_id, .keep_all = TRUE)
  
  return(df)
})

# 逐个inner_join（只保留交集基因）
merged_matrix <- reduce(
  data_list,
  inner_join,
  by = "gene_id"
)

# gene_id作为行名
merged_matrix <- as.data.frame(merged_matrix)

rownames(merged_matrix) <- merged_matrix$gene_id

merged_matrix$gene_id <- NULL

summary(as.numeric(as.matrix(merged_matrix)))

df1 <- data.frame(
  sample = colnames(merged_matrix)
)
write.csv(df1,"df1.csv")
#####================================================================================
expr <- as.matrix(merged_matrix)
mode(expr) <- "numeric"
expr_filt <- expr[
  rowSums(expr > 1) >= 3,
]
expr_log2 <- log2(expr_filt + 1)


group <- factor(c(
  "NR", "NR", "NR", "NR",
  "R", "R", "R", "R", "R", "R"
))

install.packages("BiocManager")
BiocManager::install("limma")

library(limma)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

##保留蛋白编码基因
library("qs")
gtf <- qread("G:\\HCC_靶免治疗\\bulk_in_house\\Homo_sapiens_gtf")
proteincoding <- gtf %>% filter(gene_biotype == "protein_coding")
# 剔除 RPS/RPL 核糖体蛋白基因
genes.use <- grep(
  pattern = "^(RPS|RPL)",
  x = rownames(expr_log2),
  value = TRUE,
  invert = TRUE
)

# 再保留 protein coding genes
genename_filt <- genes.use[genes.use %in% proteincoding$gene_name]

# 过滤 count matrix
expr_log2 <- expr_log2[rownames(expr_log2) %in% genename_filt, ]

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
  original_matrix = merged_matrix,
  protein_coding_filtered_matrix = expr_log2,
  deg_list = deg
)

saveRDS(
  obj_list,
  "GSE279750_HCC_ICB.rds"
)
