library(kneedle)
# https://github.com/etam4260/kneedle

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
rank <- seq_along(ctp_sorted)
x <- log10(rank)
y <- as.numeric(ctp_sorted)
library(kneedle)
# 自动判断 decreasing 和 concavity
kneedle_result <- kneedle(x, y)

knee_index  <- as.integer(kneedle_result[1])
elbow_index <- as.integer(kneedle_result[2])

#####================================================================================
ctp_sorted <- sort(
  ctp_weights[, ctp_name],
  decreasing = TRUE,
  na.last = NA
)

rank <- seq_along(ctp_sorted)
y <- as.numeric(ctp_sorted)

kneedle_result <- kneedle(
  x = rank,
  y = y,
  decreasing = TRUE,
  concave = FALSE,
  sensitivity = 1
)

kneedle_result

knee_x <- kneedle_result[1]
knee_y <- kneedle_result[2]

# 转换为最接近的整数排名
knee_index <- which.min(abs(rank - knee_x))

#####================================================================================
plot(
  rank,
  y,
  type = "l",
  log = "x",
  lwd = 1.8,
  lty = 1,
  xlab = "Gene rank",
  ylab = "Aggregated weight",
  main = paste0(ctp_name, " ranked gene weights")
)

abline(
  v = knee_index,
  col = "red",
  lwd = 2,
  lty = 3
)

points(
  knee_index,
  y[knee_index],
  col = "blue",
  pch = 19
)

legend(
  "topright",
  legend = paste0("Kneedle cutoff = ", knee_index),
  col = "blue",
  lwd = 2,
  pch = 19,
  bty = "n"
)

