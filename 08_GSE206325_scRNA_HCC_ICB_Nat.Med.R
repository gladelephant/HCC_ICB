rm(list = setdiff(ls(), c("")))
library(SingleCellExperiment)
library(Seurat)
library(SeuratObject)
#####================================================================================
# 假设你的对象叫 sce
sce

# 看有哪些 assay
assayNames(sce)

# 转 Seurat
seu <- as.Seurat(
  sce,
  counts = "counts",
  data = NULL
)

meta2 <- as.data.frame(colData(sce))
length((unique(meta2$cell_to_sample_ID)))
head(meta)
length((unique(meta2$cell_to_sample_ID)))
length((unique(meta$sample_ID)))

#####================================================================================
meta<-read.csv("GSE206325_sample_annots_Liver_Treated_patients.csv")

library(dplyr)

seu@meta.data <- seu@meta.data %>%
  mutate(cell_to_sample_ID = as.character(cell_to_sample_ID)) %>%
  left_join(
    meta %>%
      mutate(sample_ID = as.character(sample_ID)),
    by = c("cell_to_sample_ID" = "sample_ID")
  )
row.names(seu@meta.data)<-seu@meta.data$barcodes
table(seu@meta.data$prep)
#####================================================================================
meta3<-read.csv("HCC_immunotherapy_clinical_metadata.csv")
seu@meta.data <- seu@meta.data %>%
  mutate(patient_ID = as.character(patient_ID)) %>%
  left_join(
    meta3 %>%
      mutate(ID = as.character(ID)),
    by = c("patient_ID" = "ID"))
row.names(seu@meta.data)<-seu@meta.data$barcodes

length(unique(meta$patient_ID))
#####================================================================================
missing_patients <- setdiff(
  unique(meta$patient_ID),
  unique(seu@meta.data$patient_ID)
)
missing_patients
length(missing_patients)

write_rds(seu,"GSE206325_HCC_ICB_Nat_Med.rds")
