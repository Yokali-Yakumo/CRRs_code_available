
<!--
Helper: GO-term dotplots for up/down-regulated modules from GO result matrices (GO/results/*.tsv).
--> 
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 读取数据
df7_raw <- read.table("GO/results/d7_d0_res_up_module_mat.tsv", check.names = FALSE, stringsAsFactors = FALSE,header = T,sep = "\t")
df14_raw <- read.csv("GO/results/d14-d7_res_up_module_mat.tsv", check.names = FALSE, stringsAsFactors = FALSE,header = T,sep = "\t")

# 统一前两列列名
colnames(df7_raw)[1:2] <- c("id", "module_label")
colnames(df14_raw)[1:2] <- c("id", "module_label")

# 增加时间维度标识
df7_raw$Timepoint <- "d7"
df14_raw$Timepoint <- "d14"

# 2. 【核心步骤】合并为一个大表并转为长格式
# 合并后，所有的计算都基于这个 single_df
single_df <- bind_rows(df7_raw, df14_raw) %>%
  dplyr::select(module_label, Timepoint, matches("^cluster[1-7]$")) %>%
  pivot_longer(cols = matches("^cluster[1-7]$"), 
               names_to = "Cluster", 
               values_to = "Value") %>%
  dplyr::filter(!is.na(Value)) # 只保留有值的数据

# 3. 计算全局排序权重
# 逻辑：对于每个 Module，找到它在所有时间点中，最早出现的 Cluster 编号
module_order <- single_df %>%
  mutate(cluster_num = as.numeric(gsub("cluster", "", Cluster))) %>%
  group_by(module_label) %>%
  summarise(
    first_idx = min(cluster_num),      # 阶梯的主排序依据：最早出现的列
    max_score = max(Value, na.rm = TRUE) # 次要依据：数值大小
  ) %>%
  arrange(first_idx, desc(max_score)) %>%
  pull(module_label)

# 4. 应用排序并处理绘图属性
offset <- 0.08 

plot_data <- single_df %>%
  mutate(
    # 转换为因子
    module_label = factor(module_label, levels = rev(module_order)),
    Cluster = factor(Cluster, levels = paste0("cluster", 1:7)),
    # 转换为数值以便计算偏移
    y_num = as.numeric(module_label),
    x_num = as.numeric(Cluster)
  ) %>%
  mutate(
    # 逻辑：
    # d14 在左上：X减小，Y增大
    # d7  在右下：X增大，Y减小
    x_final = ifelse(Timepoint == "d14", x_num + offset, x_num - offset),
    y_final = ifelse(Timepoint == "d14", y_num - offset, y_num + offset)
  )

# 4. 绘图
# ... 前面计算坐标偏移和排序的代码保持不变 ...

# 4. 绘图部分
p <- ggplot() +
  # 第二层：绘制 d14 (左上，形状 21)
  geom_point(
    data = dplyr::filter(plot_data, Timepoint == "d14"),
    aes(x = x_final, y = y_final, size = Value, color = Value),
    shape = 21,     # 带边框的圆
    fill = "grey95", # 边框颜色固定为黑色（或深红色 "#a50f15"）
    stroke = 1.5    # 这里可以自由控制空心圆的描边宽度
  ) +
  # 第一层：绘制 d7 (右下，实心点)
  geom_point(
    data = dplyr::filter(plot_data, Timepoint == "d7"),
    aes(x = x_final, y = y_final, size = Value, color = Value),
    shape = 16,     # 纯实心圆
    alpha = 1     # 稍微加一点透明度，防止完全遮挡
  ) +
  

  
  # 统一颜色映射：将 color 和 fill 映射到同一个红色渐变色阶
  scale_color_gradientn(colors = c("#fb9a99", "#e31a1c", "#800026"), name = "Score") +
  #scale_fill_gradientn(colors = c("#fee0d2", "#fc9272", "#de2d26"), name = "Score") +
  
  # 比例尺
  scale_size_continuous(range = c(5,8)) +
  
  # 坐标轴还原
  scale_x_continuous(breaks = 1:7, labels = paste0("cluster", 1:7)) +
  scale_y_continuous(breaks = 1:length(module_order), labels = rev(module_order)) +
  
  # 样式
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,color="black"),
    axis.text.y = element_text(color="black"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    legend.box = "vertical"
  ) +    labs(title = "",x = "", y = "")

ggsave("outplots/fig3/fig3_d.go.dotplots.pdf",p,width=9,height=10,dpi=300)



library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 读取数据
df7_raw <- read.csv("GO/results/d7_d0_res_down_module_mat.tsv", check.names = FALSE, stringsAsFactors = FALSE,header = T,sep = "\t")
df14_raw <- data.frame(matrix(nrow=1,ncol=ncol(df7_raw)))
colnames(df14_raw) <- colnames(df7_raw)

# 统一前两列列名
colnames(df7_raw)[1:2] <- c("id", "module_label")
colnames(df14_raw)[1:2] <- c("id", "module_label")

# 增加时间维度标识
df7_raw$Timepoint <- "d7"
df14_raw$Timepoint <- "d14"

# 2. 【核心步骤】合并为一个大表并转为长格式
# 合并后，所有的计算都基于这个 single_df
single_df <- bind_rows(df7_raw, df14_raw) %>%
  dplyr::select(module_label, Timepoint, matches("^cluster[1-7]$")) %>%
  pivot_longer(cols = matches("^cluster[1-7]$"), 
               names_to = "Cluster", 
               values_to = "Value") %>%
  dplyr::filter(!is.na(Value)) # 只保留有值的数据

# 3. 计算全局排序权重
# 逻辑：对于每个 Module，找到它在所有时间点中，最早出现的 Cluster 编号
module_order <- single_df %>%
  mutate(cluster_num = as.numeric(gsub("cluster", "", Cluster))) %>%
  group_by(module_label) %>%
  summarise(
    first_idx = min(cluster_num),      # 阶梯的主排序依据：最早出现的列
    max_score = max(Value, na.rm = TRUE) # 次要依据：数值大小
  ) %>%
  arrange(first_idx, desc(max_score)) %>%
  pull(module_label)

# 4. 应用排序并处理绘图属性
offset <- 0.08 

plot_data <- single_df %>%
  mutate(
    # 转换为因子
    module_label = factor(module_label, levels = rev(module_order)),
    Cluster = factor(Cluster, levels = paste0("cluster", 1:7)),
    # 转换为数值以便计算偏移
    y_num = as.numeric(module_label),
    x_num = as.numeric(Cluster)
  ) %>%
  mutate(
    # 逻辑：
    # d14 在左上：X减小，Y增大
    # d7  在右下：X增大，Y减小
    x_final = ifelse(Timepoint == "d14", x_num + offset, x_num - offset),
    y_final = ifelse(Timepoint == "d14", y_num - offset, y_num + offset)
  )

# 4. 绘图
# ... 前面计算坐标偏移和排序的代码保持不变 ...

# 4. 绘图部分
p <- ggplot() +
  # 第二层：绘制 d14 (左上，形状 21)
  geom_point(
    data = dplyr::filter(plot_data, Timepoint == "d14"),
    aes(x = x_final, y = y_final, size = Value, color = Value),
    shape = 21,     # 带边框的圆
    fill = "grey95", # 边框颜色固定为黑色（或深红色 "#a50f15"）
    stroke = 1.5    # 这里可以自由控制空心圆的描边宽度
  ) +
  # 第一层：绘制 d7 (右下，实心点)
  geom_point(
    data = dplyr::filter(plot_data, Timepoint == "d7"),
    aes(x = x_final, y = y_final, size = Value, color = Value),
    shape = 16,     # 纯实心圆
    alpha = 1     # 稍微加一点透明度，防止完全遮挡
  ) +
  

  
  # 统一颜色映射：将 color 和 fill 映射到同一个红色渐变色阶
  scale_color_gradientn(colors = c("#4eb3d3", "#0868ac", "#084081"), name = "Score") +
  #scale_fill_gradientn(colors = c("#fee0d2", "#fc9272", "#de2d26"), name = "Score") +
  
  # 比例尺
  scale_size_continuous(range = c(5,8)) +
  
  # 坐标轴还原
  scale_x_continuous(breaks = 1:7, labels = paste0("cluster", 1:7)) +
  scale_y_continuous(breaks = 1:length(module_order), labels = rev(module_order)) +
  
  # 样式
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,color="black"),
    axis.text.y = element_text(color="black"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    legend.box = "vertical"
  ) +    labs(title = "",x = "", y = "")

ggsave("outplots/fig3/fig3_d.go.down.dotplots.pdf",p,width=9,height=10,dpi=300)
