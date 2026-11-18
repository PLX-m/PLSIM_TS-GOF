library(ggplot2)
library(patchwork)

df <- data.frame(
  c = rep(c(0, 0.25, 0.5, 0.75, 1), times = 6),
  value = c(
    # n=800, d=s=10：Z1, Z2, CCT
    0.048, 0.649, 0.833, 0.858, 0.871,
    0.048, 0.628, 0.837, 0.856, 0.858,
    0.054, 0.732, 0.894, 0.911, 0.909,
    
    # n=800, d=s=20：Z1, Z2, CCT
    0.039, 0.404, 0.605, 0.652, 0.700,
    0.035, 0.387, 0.616, 0.663, 0.682,
    0.035, 0.488, 0.714, 0.747, 0.772
  ),
  Method = rep(c("Z1", "Z2", "Tc"), each = 5, times = 2),
  scenario = rep(c("d10", "d20"), each = 15)
)

# ==============================
# 字体大小参数
# ==============================
base_size        <- 16  # theme_bw 的基础字体大小
title_size       <- 16  # 子图标题字体大小
axis_title_size  <- 16  # 横纵坐标标题字体大小
axis_text_size   <- 14  # 坐标轴刻度字体大小
legend_text_size <- 14  # 图例文字字体大小
legend_title_size <- 16 # 图例标题字体大小


# ==============================
# 颜色
# ==============================
my_cols <- c(
  "Z1" = "#3B82B8",
  "Z2" = "#E69F00",
  "Tc" = "#D55E5E"
)

# ==============================
# 线型
# ==============================
my_linetypes <- c(
  "Z1" = "solid",
  "Z2" = "solid",
  "Tc" = "solid"
)

# ==============================
# 点的形状
# ==============================
my_shapes <- c(
  "Z1" = 16,   # 圆点
  "Z2" = 17,   # 三角形
  "Tc" = 15    # 方块
)

# 图例标签
my_labels <- c(
  expression(Z[n]^{(1)}),
  expression(Z[n]^{(2)}),
  expression(T[C])
)


# ==============================
# d=10 图
# ==============================
p1 <- ggplot(
  subset(df, scenario == "d10"),
  aes(
    x = c, y = value,
    color = Method,
    linetype = Method,
    shape = Method
   )
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = my_cols,
    breaks = c("Z1", "Z2", "Tc"),
    labels = my_labels
  ) + 
  scale_linetype_manual(
    values = my_linetypes,
    breaks = c("Z1", "Z2", "Tc"),
    labels = my_labels
  ) +
  scale_shape_manual(
    values = my_shapes,
    breaks = c("Z1", "Z2", "Tc"),
    labels = my_labels
  ) +
  # 5% 显著性水平
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "c",   
    y = "Probability of rejection",
    title = "p=q=10"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    plot.title   = element_text(size = title_size),        # 子图标题
    axis.title   = element_text(size = axis_title_size),   # 横纵坐标标题
    axis.text    = element_text(size = axis_text_size),    # 坐标轴刻度
    legend.text  = element_text(size = legend_text_size),  # 图例文字
    legend.title = element_text(size = legend_title_size)  # 图例标题
    # panel.grid = element_blank()   # 去掉网格线
    
  )

# ==============================
# d=20 图
# ==============================
p2 <- ggplot(
  subset(df, scenario == "d20"),
  aes(
    x = c, y = value,
    color = Method,
    linetype = Method,
    shape = Method
    )
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = my_cols,
    breaks = c("Z1", "Z2", "Tc"),
    labels = my_labels
  ) + 
  scale_linetype_manual(
    values = my_linetypes,
    breaks = c("Z1", "Z2", "Tc"),
    labels = my_labels
  ) +
  scale_shape_manual(
    values = my_shapes,
    breaks = c("Z1", "Z2", "Tc"),
    labels = my_labels
  ) +
  # 5% 显著性水平
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "c",
    y = "Probability of rejection",
    title = "p=q=20"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    plot.title   = element_text(size = title_size),        # 子图标题
    axis.title   = element_text(size = axis_title_size),   # 横纵坐标标题
    axis.text    = element_text(size = axis_text_size),    # 坐标轴刻度
    legend.text  = element_text(size = legend_text_size),  # 图例文字
    legend.title = element_text(size = legend_title_size)  # 图例标题
    # panel.grid = element_blank()   # 去掉网格线
    
  )



# ==============================
# 拼接（横向）
# ==============================
p1 + p2

p_final <- (p1 + p2) + plot_layout(guides = "collect")
 
ggsave(
  filename = "plot/power_plot_05exp(-t_val^2)_shapes.pdf",
  plot = p_final,
  width = 12,
  height = 5,
  units = "in",
  device = cairo_pdf   # 字体更清晰
)


