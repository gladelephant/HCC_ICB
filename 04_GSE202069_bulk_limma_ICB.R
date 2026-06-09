############################################################
## Differential expression analysis of GSE202069 bulk RNA-seq
## Input: log2(TPM + 1) expression matrix
## Method: limma empirical Bayes linear model
############################################################

## 1. Load packages ---------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(limma)
  library(qs)
})

## 2. Set paths -------------------------------------------------------------

data_dir <- "./GSE202069_bulk_rna_seq"
expr_file <- file.path(data_dir, "GSE202069_gene_tpm_expression.txt")
meta_file <- file.path(data_dir, "GSE202069_metadata.csv")
gtf_file  <- "G:/HCC_靶免治疗/bulk_in_house/Homo_sapiens_gtf"

out_dir <- file.path(data_dir, "DEG_limma_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 3. Read expression matrix and metadata ----------------------------------

## Expression values were provided as log2(TPM + 1) by the original authors
expr <- read.table(
  expr_file,
  row.names = 1,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

metadata <- read.csv(
  meta_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(all(c("sample", "response") %in% colnames(metadata)))

rownames(metadata) <- metadata$sample

## Keep only samples present in both expression matrix and metadata
common_samples <- intersect(colnames(expr), metadata$sample)

expr <- expr[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

stopifnot(identical(colnames(expr), rownames(metadata)))
colnames(expr)
rownames(metadata)
## 4. Gene filtering --------------------------------------------------------

gtf <- qread(gtf_file)

protein_coding_genes <- gtf %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(gene_name) %>%
  unique()

## Remove ribosomal protein genes: RPS and RPL families
non_ribo_genes <- rownames(expr)[
  !grepl("^RP[SL][0-9A-Z]", rownames(expr))
]

genes_keep <- intersect(non_ribo_genes, protein_coding_genes)

expr <- expr[genes_keep, , drop = FALSE]

## Remove genes with zero variance across samples
expr <- expr[apply(expr, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

## 5. Define groups ---------------------------------------------------------

metadata$response <- factor(
  metadata$response,
  levels = c("Non_responders", "Responders")
)

stopifnot(!any(is.na(metadata$response)))

design <- model.matrix(~ 0 + response, data = metadata)

colnames(design) <- levels(metadata$response)
rownames(design) <- rownames(metadata)

stopifnot(identical(colnames(expr), rownames(design)))

## 6. Differential expression analysis -------------------------------------

fit <- lmFit(expr, design)

contrast_matrix <- makeContrasts(
  Responders_vs_Non_responders = Responders - Non_responders,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE)

deg_all <- topTable(
  fit2,
  coef = "Responders_vs_Non_responders",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg_all <- deg_all %>%
  rownames_to_column("gene") %>%
  mutate(
    comparison = "Responders_vs_Non_responders",
    regulation = case_when(
      adj.P.Val < 0.1 & logFC >  0.5 ~ "Up_in_Responders",
      adj.P.Val < 0.1 & logFC < -0.5 ~ "Down_in_Responders",
      TRUE ~ "Not_significant"
    )
  )

## 7. Save results ----------------------------------------------------------

write.csv(
  deg_all,
  file = file.path(out_dir, "GSE202069_limma_all_genes.csv"),
  row.names = FALSE
)

deg_sig <- deg_all %>%
  filter(adj.P.Val < 0.1, abs(logFC) > 0.5)

write.csv(
  deg_sig,
  file = file.path(out_dir, "GSE202069_limma_DEGs_FDR0.05_logFC1.csv"),
  row.names = FALSE
)

## 8. Save analysis objects -------------------------------------------------

saveRDS(
  list(
    expression_matrix = expr,
    metadata = metadata,
    design = design,
    contrast_matrix = contrast_matrix,
    fit = fit2,
    deg_all = deg_all,
    deg_sig = deg_sig
  ),
  file = file.path(out_dir, "GSE202069_limma_analysis_objects.rds")
)

## 9. Record session information -------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)
