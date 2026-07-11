library(dplyr)
library(ggplot2)
#画的是每个CTP的gene loading（基因权重）排序曲线
#一个CTP是由少数几个高权重基因驱动，还是由大量中等权重基因共同构成？

plot.ctps <- c(
  "CTP_5",
  "CTP_10",
  "CTP_25",
  "CTP_14"
)

plot.df <- data.frame()

for(ctp in plot.ctps){
  
  x <- sort(
    weight.use[, ctp],
    decreasing = TRUE
  )
  
  tmp <- data.frame(
    Rank = 1:length(x),
    Loading = x,
    CTP = ctp
  )
  
  plot.df <- rbind(
    plot.df,
    tmp
  )
}

ggplot(
  plot.df,
  aes(
    Rank,
    Loading,
    color = CTP
  )
) +
  geom_line(size=1) +
  theme_classic(base_size = 14) +
  labs(
    x="Gene rank",
    y="Gene loading",
    title="Gene loading distribution"
  )



ggplot(
  plot.df,
  aes(
    Rank,
    Loading,
    color=CTP
  )
)+
  geom_line(size=1)+
  scale_x_log10()+
  theme_classic(base_size=14)+
  labs(
    x="Gene rank (log10)",
    y="Loading",
    title="Gene loading distribution"
  )
