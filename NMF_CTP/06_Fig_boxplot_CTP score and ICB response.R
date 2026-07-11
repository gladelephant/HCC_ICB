#CTP score and ICB response 
#===========================================================
# CTP score：Responder vs Non-responder
# Nature风格箱线图 + FDR标注
#===========================================================

library(ggplot2)
library(dplyr)
library(ggpubr)
library(grid)

#-----------------------------------------------------------
# 1. 设置分组顺序
#-----------------------------------------------------------

ctp.long <- ctp.long %>%
  mutate(
    ORR = factor(
      ORR,
      levels = c("Non_responder", "Responder")
    )
  )

library(stringr)

ctp.order <- ctp.long %>%
  distinct(CTP) %>%
  mutate(
    number = as.numeric(str_extract(CTP, "\\d+"))
  ) %>%
  arrange(number) %>%
  pull(CTP)

ctp.long$CTP <- factor(
  ctp.long$CTP,
  levels = ctp.order
)

stat.df$CTP <- factor(
  stat.df$CTP,
  levels = ctp.order
)

#-----------------------------------------------------------
# 2. 使用已有的 orr.result 构建FDR标注数据
#-----------------------------------------------------------

stat.df <- ctp.long %>%
  group_by(CTP) %>%
  summarise(
    ymax = max(score, na.rm = TRUE),
    ymin = min(score, na.rm = TRUE),
    yrange = ymax - ymin,
    .groups = "drop"
  ) %>%
  left_join(
    orr.result %>%
      select(CTP, FDR),
    by = "CTP"
  ) %>%
  mutate(
    # 必须与ctp.long中ORR的实际名称一致
    group1 = "Non_responder",
    group2 = "Responder",
    
    # 防止某个CTP的取值范围为0
    yrange = ifelse(
      yrange == 0,
      abs(ymax) * 0.1 + 0.1,
      yrange
    ),
    
    # 横线位置
    y.position = ymax + yrange * 0.13,
    
    # FDR显示格式
    label = case_when(
      FDR < 0.001 ~ paste0(
        " ",
        formatC(
          FDR,
          format = "e",
          digits = 2
        )
      ),
      TRUE ~ paste0(
        " ",
        formatC(
          FDR,
          format = "f",
          digits = 3
        )
      )
    )
  )

# 查看标注数据
print(stat.df)

#-----------------------------------------------------------
# 3. 绘图
#-----------------------------------------------------------

p.ctp.orr <- ggplot(
  ctp.long,
  aes(
    x = ORR,
    y = score,
    fill = ORR
  )
) +
  
  geom_boxplot(
    width = 0.60,
    linewidth = 0.5,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  
  geom_jitter(
    aes(color = ORR),
    width = 0.15,
    height = 0,
    size = 0.2,
    alpha = 0.45
  ) +
  
  facet_wrap(
    ~ CTP,
    scales = "free_y",
    ncol = 4
  ) +
  
  # 添加横线、两端短竖线和FDR
  stat_pvalue_manual(
    stat.df,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.025,
    bracket.size = 0.45,
    size = 3.1,
    vjust = -0.15,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      "Non_responder" = "#4C78A8",
      "Responder" = "#E68613"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Non_responder" = "#4C78A8",
      "Responder" = "#E68613"
    )
  ) +
  
  # 给顶部FDR标注留出空间
  scale_y_continuous(
    expand = expansion(
      mult = c(0.06, 0.22)
    )
  ) +
  
  labs(
    x = NULL,
    y = "CTP score"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "none",
    
    # 去除CTP标题黑色边框
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 11,
      colour = "black"
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 12,
      colour = "black"
    ),
    
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1
    ),
    
    axis.line = element_line(
      linewidth = 0.6,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    panel.spacing = unit(
      1,
      "lines"
    ),
    
    plot.margin = margin(
      t = 12,
      r = 8,
      b = 8,
      l = 8
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  )

p.ctp.orr

ggsave(
  filename = "CTP_ORR_boxplot_FDR.pdf",
  plot = p.ctp.orr,
  width = 6,
  height = 6,
  device = cairo_pdf
)

