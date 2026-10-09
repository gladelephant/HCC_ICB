
library(tidyverse)
library(data.table)
library(igraph)

setwd("D:/HCC_ICB/cuda_nmf_results")

# 1. 读取NMF Ascaled矩阵
MP_df <- fread(
  "G:\\08_HCC_靶免治疗\\CTP_V2\\all_groups_merged_Ascaled_K10_K25.tsv",
  data.table = FALSE,
  check.names = FALSE
)

MP <- as.matrix(MP_df[, -1])
storage.mode(MP) <- "double"
rownames(MP) <- MP_df[[1]]
colnames(MP) <- str_replace(
  colnames(MP), "_Pattern([0-9]+)$", "_P\\1"
)

# 2. 解析NMF factor信息
factor_info <- tibble(
  factor = colnames(MP),
  dataset = str_remove(factor, "_k[0-9]+_P[0-9]+$")
)

# 3. Pearson相关网络 + Infomap
cor_matrix <- cor(MP, method = "pearson")
idx <- which(upper.tri(cor_matrix) & cor_matrix >= 0.50,
             arr.ind = TRUE)

edges <- tibble(
  from = colnames(MP)[idx[, 1]],
  to = colnames(MP)[idx[, 2]],
  weight = cor_matrix[idx]
)

graph <- graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = factor_info |> rename(name = factor)
)

graph <- delete_vertices(graph, which(degree(graph) == 0))

set.seed(123)
community <- cluster_infomap(
  graph,
  e.weights = E(graph)$weight,
  nb.trials = 200
)

# 4. 筛选跨数据集CTP
ctp_nodes <- as_data_frame(graph, what = "vertices") |>
  as_tibble() |>
  mutate(
    community_id = paste0("CTP_", membership(community))
  ) |>
  group_by(community_id) |>
  filter(n_distinct(dataset) >= 2) |>
  ungroup()

# 5. 计算5基因的factor-level loading-rank percentile
genes <- intersect(
  c("REG3A", "GLUD1", "AIFM2", "CACFD1", "CTNNA2"),
  rownames(MP)
)

rank_matrix <- apply(
  MP, 2,
  \(x) 1 - rank(-x, ties.method = "average") / length(x)
)

factor_score <- colMeans(rank_matrix[genes, , drop = FALSE])

factor_result <- ctp_nodes |>
  mutate(NR29_score = as.numeric(factor_score[name])) |>
  select(community_id, dataset, name, NR29_score)

# 6. Figure 7：Cross-dataset factor consistency
nature_cols <- c(
  "#3C5488", "#00A087", "#E64B35",
  "#8491B4", "#F39B7F", "#91D1C2",
  "#7E6148", "#B09C85"
)

p <- ggplot(
  factor_result,
  aes(
    x = community_id,
    y = NR29_score,
    colour = dataset
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.15,
    size = 0.5,
    alpha = 0.55
  ) +
  coord_flip() +
  scale_color_manual(values = nature_cols) +
  labs(
    x = NULL,
    y = "",
    title = "Loading-rank of resistant signature"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11),
    axis.line = element_line(linewidth = 0.4),
    legend.title = element_blank(),
    legend.position = "right"
  )

print(p)

ggsave(
  "Cross_dataset_factor_consistency.pdf",
  p,
  width = 3.5,
  height = 3.5
)
