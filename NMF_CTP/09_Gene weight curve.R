ctp_name <- "CTP_5"
ctp_weights <- read.delim(
  file.path(
    "output",
    "CTP_analysis",
    "tables",
    "ctp_aggweights.tsv"
  ),
  row.names = 1,
  check.names = FALSE
)

dim(ctp_weights)
head(ctp_weights)
colnames(ctp_weights)
ctp_sorted <- sort(
  ctp_weights[, ctp_name],
  decreasing = TRUE
)

plot(
  seq_along(ctp_sorted),
  ctp_sorted,
  type = "l",
  log = "x",
  xlab = "Gene rank",
  ylab = "Aggregated weight",
  main = paste0(ctp_name, " ranked gene weights")
)+
abline(v = c(30, 50, 100, 200,500,1000), lty = 2)
