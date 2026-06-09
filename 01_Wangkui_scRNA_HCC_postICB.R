library(Seurat)
library(SeuratDisk)
library(Matrix)
library(tidyverse)
rm(list = setdiff(ls(), c("")))
# 你的总目录
base_dir <- "G:\\HCC_靶免治疗\\王葵-靶免"

# 获取所有样本文件夹
sample_dirs <- list.dirs(base_dir, recursive = FALSE)

# 读取每个样本
seu_list <- lapply(sample_dirs, function(x){
  
  # 文件夹名作为sample名
  sample_name <- basename(x)
  
  # 读取10X矩阵
  mat <- Read10X(data.dir = x)
  
  # 创建Seurat对象
  seu <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  # 添加sample信息
  seu$sample <- sample_name
  
  return(seu)
})

# 给list命名
names(seu_list) <- basename(sample_dirs)

# 合并
combined <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  add.cell.ids = names(seu_list),
  project = "Merged_scRNA"
)

# 去掉下划线前面的部分
combined@meta.data$lesion <-combined@meta.data$sample
combined@meta.data$lesion <- sub(".*_", "", combined@meta.data$lesion)
# 重命名
combined@meta.data$lesion <- recode(combined@meta.data$lesion,
                     "adj" = "ANT",
                     "N" = "Normal",
                     "T" = "Tumor")
table(combined@meta.data$lesion )

combined@meta.data$response<-combined@meta.data$orig.ident
#采用RECIST1.1标准
combined@meta.data$response <- ifelse(
  grepl("^P[13]", combined@meta.data$response),
  "Responder",
  "Non_responder"
)
#排除1例ICC病人
combined <- combined[, !grepl("^P2_", combined$orig.ident)]

table(combined@meta.data$orig.ident)

library(readxl)
library(dplyr)
library(magrittr)

df <- read_xlsx("G:\\HCC_靶免治疗\\Wangkui_dataset\\meteadata.xlsx")

combined@meta.data <- combined@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  left_join(df, by = c("orig.ident" = "sample")) %>%
  tibble::column_to_rownames("cell")

colnames(combined@meta.data)

combined@meta.data <- combined@meta.data %>%
  select(-ends_with(".x")) %>%
  rename_with(~ gsub("\\.y$", "", .x), ends_with(".y"))

colnames(combined@meta.data)[colnames(combined@meta.data) == "Response"] <- "response"
write_rds(combined,"Wangkui_HCC_post_ICB.rds")

step 1: factor to character, or else your factor will be number in adata
i <- sapply(combined@meta.data, is.factor)
combined@meta.data[i] <- lapply(combined@meta.data[i], as.character)
# step 2: change to h5ad
SaveH5Seurat(combined,filename="scobj.h5seurat", overwrite = TRUE)
Convert("scobj.h5seurat", dest = "h5ad", assay="RNA", overwrite = TRUE)
