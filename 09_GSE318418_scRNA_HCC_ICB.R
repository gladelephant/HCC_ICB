library(Seurat)
library(stringr)
library(Matrix)

data_dir <- x

matrix_files <- list.files(
  data_dir,
  pattern = "matrix.mtx.gz$",
  recursive = TRUE,
  full.names = TRUE
)

seurat_list <- list()

for (mtx in matrix_files) {
  
  sample_dir <- dirname(mtx)
  sample_name <- basename(sample_dir)
  
  cat("\nReading:", sample_name, "\n")
  
  files <- list.files(sample_dir, full.names = TRUE)
  
  barcode_file <- files[str_detect(basename(files), "barcodes.tsv.gz$")]
  feature_file <- files[str_detect(basename(files), "features.tsv.gz$|genes.tsv.gz$")]
  
  mat <- readMM(mtx)
  
  features <- read.delim(
    gzfile(feature_file),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  barcodes <- read.delim(
    gzfile(barcode_file),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  cat("matrix:", dim(mat), "\n")
  cat("features:", nrow(features), "\n")
  cat("barcodes:", nrow(barcodes), "\n")
  
  if (nrow(mat) != nrow(features)) {
    stop("matrix 行数和 features 行数不一致: ", sample_name)
  }
  
  if (ncol(mat) != nrow(barcodes)) {
    stop("matrix 列数和 barcodes 行数不一致: ", sample_name)
  }
  
  gene_names <- features[[1]]
  cell_names <- barcodes[[1]]
  
  gene_names[is.na(gene_names) | gene_names == ""] <- paste0(
    "UnknownGene_",
    which(is.na(gene_names) | gene_names == "")
  )
  
  rownames(mat) <- make.unique(gene_names)
  colnames(mat) <- cell_names
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  obj$sample <- sample_name
  obj$GSM <- str_extract(sample_name, "GSM[0-9]+")
  obj$tissue <- ifelse(
    str_detect(sample_name, "_T$|_T_"),
    "Tumor",
    ifelse(str_detect(sample_name, "_N$|_N_"), "Normal", NA)
  )
  
  seurat_list[[sample_name]] <- obj
  
  cat("成功:", sample_name, "\n")
}

cat("\n成功读取样本数:", length(seurat_list), "\n")

combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "HCC_combined"
)


remove_samples <- c(
  "GSM9494103",
  "GSM9494105",
  "GSM9494107",
  "GSM9494109",
  "GSM9494111",
  "GSM9494113"
)

cells_remove <- grep(
  paste(remove_samples, collapse = "|"),
  colnames(combined),
  value = TRUE
)

combined_subset <- subset(
  combined,
  cells = setdiff(colnames(combined), cells_remove)
)

table(combined_subset@meta.data$orig.ident,combined_subset$tissue)

combined_subset$tissue[
  combined_subset$orig.ident == "GSM9494102_A_18N_Live"
] <- "Normal"

combined_subset$tissue[
  combined_subset$orig.ident == "GSM9494110_A_23_Live"
] <- "Tumor"

saveRDS(
  combined,
  file = file.path(data_dir, "HCC_combined_raw.rds")
)
