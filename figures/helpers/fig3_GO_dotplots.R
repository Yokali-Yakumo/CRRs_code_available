library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Read data
df7_raw <- read.table("GO/results/d7_d0_res_up_module_mat.tsv", check.names = FALSE, stringsAsFactors = FALSE,header = T,sep = "\t")
df14_raw <- read.csv("GO/results/d14-d7_res_up_module_mat.tsv", check.names = FALSE, stringsAsFactors = FALSE,header = T,sep = "\t")

# Standardise the names of the first two columns
colnames(df7_raw)[1:2] <- c("id", "module_label")
colnames(df14_raw)[1:2] <- c("id", "module_label")

# Add the time-point label
df7_raw$Timepoint <- "d7"
df14_raw$Timepoint <- "d14"

# 2. Key step: merge into one table and reshape to long format
# After merging, all calculations are based on this single_df
single_df <- bind_rows(df7_raw, df14_raw) %>%
  dplyr::select(module_label, Timepoint, matches("^cluster[1-7]$")) %>%
  pivot_longer(cols = matches("^cluster[1-7]$"), 
               names_to = "Cluster", 
               values_to = "Value") %>%
  dplyr::filter(!is.na(Value)) # keep only entries that have a value

# 3. Compute the global ordering keys
# Logic: for each module, find the earliest cluster index at which it appears across all time points
module_order <- single_df %>%
  mutate(cluster_num = as.numeric(gsub("cluster", "", Cluster))) %>%
  group_by(module_label) %>%
  summarise(
    first_idx = min(cluster_num),      # primary ordering key: earliest column of the staircase
    max_score = max(Value, na.rm = TRUE) # secondary key: magnitude of the score
  ) %>%
  arrange(first_idx, desc(max_score)) %>%
  pull(module_label)

# 4. Apply the ordering and set up plotting attributes
offset <- 0.08 

plot_data <- single_df %>%
  mutate(
    # convert to factors
    module_label = factor(module_label, levels = rev(module_order)),
    Cluster = factor(Cluster, levels = paste0("cluster", 1:7)),
    # convert to numeric so the offsets can be computed
    y_num = as.numeric(module_label),
    x_num = as.numeric(Cluster)
  ) %>%
  mutate(
    # Logic:
    # d14 sits upper left: X decreases, Y increases
    # d7  sits lower right: X increases, Y decreases
    x_final = ifelse(Timepoint == "d14", x_num + offset, x_num - offset),
    y_final = ifelse(Timepoint == "d14", y_num - offset, y_num + offset)
  )

# 4. Plot
# ... the preceding code that computes the coordinate offsets and ordering is unchanged ...

# 4. Plotting
p <- ggplot() +
  # Second layer: plot d14 (upper left, shape 21)
  geom_point(
    data = dplyr::filter(plot_data, Timepoint == "d14"),
    aes(x = x_final, y = y_final, size = Value, color = Value),
    shape = 21,     # circle with a border
    fill = "grey95", # border colour fixed to black (or dark red "#a50f15")
    stroke = 1.5    # stroke width of the hollow circle can be set freely here
  ) +
  # First layer: plot d7 (lower right, filled point)
  geom_point(
    data = dplyr::filter(plot_data, Timepoint == "d7"),
    aes(x = x_final, y = y_final, size = Value, color = Value),
    shape = 16,     # solid filled circle
    alpha = 1     # a little transparency prevents complete occlusion
  ) +
  

  
  # Unified colour mapping: map color and fill to the same gradient scale
  scale_color_gradientn(colors = c("#fb9a99", "#e31a1c", "#800026"), name = "Score") +
  #scale_fill_gradientn(colors = c("#fee0d2", "#fc9272", "#de2d26"), name = "Score") +
  
  # Size scale
  scale_size_continuous(range = c(5,8)) +
  
  # Restore the axes
  scale_x_continuous(breaks = 1:7, labels = paste0("cluster", 1:7)) +
  scale_y_continuous(breaks = 1:length(module_order), labels = rev(module_order)) +
  
  # Styling
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

# 1. Read data
df7_raw <- read.csv("GO/results/d7_d0_res_down_module_mat.tsv", check.names = FALSE, stringsAsFactors = FALSE,header = T,sep = "\t")
df14_raw <- data.frame(matrix(nrow=1,ncol=ncol(df7_raw)))
colnames(df14_raw) <- colnames(df7_raw)

# Standardise the names of the first two columns
colnames(df7_raw)[1:2] <- c("id", "module_label")
colnames(df14_raw)[1:2] <- c("id", "module_label")

# Add the time-point label
df7_raw$Timepoint <- "d7"
df14_raw$Timepoint <- "d14"

# 2. Key step: merge into one table and reshape to long format
# After merging, all calculations are based on this single_df
single_df <- bind_rows(df7_raw, df14_raw) %>%
  dplyr::select(module_label, Timepoint, matches("^cluster[1-7]$")) %>%
  pivot_longer(cols = matches("^cluster[1-7]$"), 
               names_to = "Cluster", 
               values_to = "Value") %>%
  dplyr::filter(!is.na(Value)) # keep only entries that have a value

# 3. Compute the global ordering keys
# Logic: for each module, find the earliest cluster index at which it appears across all time points
module_order <- single_df %>%
  mutate(cluster_num = as.numeric(gsub("cluster", "", Cluster))) %>%
  group_by(module_label) %>%
  summarise(
    first_idx = min(cluster_num),      # primary ordering key: earliest column of the staircase
    max_score = max(Value, na.rm = TRUE) # secondary key: magnitude of the score
  ) %>%
  arrange(first_idx, desc(max_score)) %>%
  pull(module_label)

# 4. Apply the ordering and set up plotting attributes
offset <- 0.08 

plot_data <- single_df %>%
  mutate(
    # convert to factors
    module_label = factor(module_label, levels = rev(module_order)),
    Cluster = factor(Cluster, levels = paste0("cluster", 1:7)),
    # convert to numeric so the offsets can be computed
    y_num = as.numeric(module_label),
    x_num = as.numeric(Cluster)
  ) %>%
  mutate(
    # Logic:
    # d14 sits upper left: X decreases, Y increases
    # d7  sits lower right: X increases, Y decreases
    x_final = ifelse(Timepoint == "d14", x_num + offset, x_num - offset),
    y_final = ifelse(Timepoint == "d14", y_num - offset, y_num + offset)
  )

# 4. Plot
# ... the preceding code that computes the coordinate offsets and ordering is unchanged ...

# 4. Plotting
p <- ggplot() +
  # Second layer: plot d14 (upper left, shape 21)
  geom_point(
    data = dplyr::filter(plot_data, Timepoint == "d14"),
    aes(x = x_final, y = y_final, size = Value, color = Value),
    shape = 21,     # circle with a border
    fill = "grey95", # border colour fixed to black (or dark red "#a50f15")
    stroke = 1.5    # stroke width of the hollow circle can be set freely here
  ) +
  # First layer: plot d7 (lower right, filled point)
  geom_point(
    data = dplyr::filter(plot_data, Timepoint == "d7"),
    aes(x = x_final, y = y_final, size = Value, color = Value),
    shape = 16,     # solid filled circle
    alpha = 1     # a little transparency prevents complete occlusion
  ) +
  

  
  # Unified colour mapping: map color and fill to the same gradient scale
  scale_color_gradientn(colors = c("#4eb3d3", "#0868ac", "#084081"), name = "Score") +
  #scale_fill_gradientn(colors = c("#fee0d2", "#fc9272", "#de2d26"), name = "Score") +
  
  # Size scale
  scale_size_continuous(range = c(5,8)) +
  
  # Restore the axes
  scale_x_continuous(breaks = 1:7, labels = paste0("cluster", 1:7)) +
  scale_y_continuous(breaks = 1:length(module_order), labels = rev(module_order)) +
  
  # Styling
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,color="black"),
    axis.text.y = element_text(color="black"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    legend.box = "vertical"
  ) +    labs(title = "",x = "", y = "")

ggsave("outplots/fig3/fig3_d.go.down.dotplots.pdf",p,width=9,height=10,dpi=300)
