
<!--
Helper: chromHMM state-composition summary for a set of genomic windows grouped by cluster/time (returns a plottable data.frame).
--> 
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

chromHMM_composition <- function(
    win,
    chromHMM.annot,
    time,
    cluster,
    annotation_col = 4,
    use_bp = TRUE
) {

  ## chromHMM.annot中定义的全部可能注释
  all.annotations <- unique(
    trimws(as.character(chromHMM.annot$Annotation))
  )

  ## 提取注释列
  if (is.numeric(annotation_col)) {
    annotation <- win[[annotation_col]]
  } else {
    if (!annotation_col %in% colnames(win)) {
      stop("输入数据中不存在注释列：", annotation_col)
    }
    annotation <- win[[annotation_col]]
  }

  annotation <- trimws(as.character(annotation))

  ## 检查NA
  if (any(is.na(annotation) | annotation == "")) {
    stop(
      time, " ", cluster,
      "中存在NA或空白ChromHMM注释。"
    )
  }

  ## 检查是否出现chromHMM.annot之外的注释
  unknown.annotations <- setdiff(
    unique(annotation),
    all.annotations
  )

  if (length(unknown.annotations) > 0) {
    stop(
      time, " ", cluster,
      "中存在未定义注释：",
      paste(unknown.annotations, collapse = ", ")
    )
  }

  dat <- data.frame(
    start = as.numeric(win[[2]]),
    end = as.numeric(win[[3]]),
    Annotation = annotation,
    stringsAsFactors = FALSE
  )

  ## use_bp=TRUE：按覆盖碱基数统计
  ## use_bp=FALSE：每个窗口权重均为1，相当于table(V4)
  if (use_bp) {
    dat$Weight <- dat$end - dat$start
  } else {
    dat$Weight <- 1
  }

  if (any(dat$Weight <= 0, na.rm = TRUE)) {
    stop(
      time, " ", cluster,
      "中存在终点小于或等于起点的窗口。"
    )
  }

  result <- dat %>%
    group_by(Annotation) %>%
    summarise(
      Count = n(),
      Covered_bp = sum(end - start),
      Weight = sum(Weight),
      .groups = "drop"
    ) %>%

    ## 补齐当前数据中没有出现的注释
    complete(
      Annotation = all.annotations,
      fill = list(
        Count = 0,
        Covered_bp = 0,
        Weight = 0
      )
    ) %>%

    mutate(
      Proportion = Weight / sum(Weight),
      Time = time,
      Cluster = cluster
    ) %>%

    select(
      Time,
      Cluster,
      Annotation,
      Count,
      Covered_bp,
      Proportion
    )

  return(result)
}



make_chromHMM_cluster_panel <- function(
    data,
    cluster.name,
    show.y.axis = FALSE
) {

  p <- data %>%
    filter(
      Cluster == cluster.name
    ) %>%

    ggplot(
      aes(
        x = Time,
        y = Proportion,
        fill = Annotation
      )
    ) +

    geom_col(
      width = 0.68,
      color = "white",
      linewidth = 0.25,
      position = position_stack(
        reverse = TRUE
      )
    ) +

    ## 每个单独图仍使用facet strip展示cluster名称
    facet_wrap(
      ~ Cluster,
      nrow = 1
    ) +

    scale_fill_manual(
      values = chromHMM.colors,
      breaks = annotation.order,
      drop = FALSE,
      name = "ChromHMM annotation"
    ) +

    scale_y_continuous(
      labels = percent_format(
        accuracy = 1
      ),
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      limits = c(0, 1),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +

    labs(
      x = NULL,
      y = if (show.y.axis) {
        "ChromHMM annotation proportion"
      } else {
        NULL
      }
    ) +

    guides(
      fill = guide_legend(
        ncol = 1,
        byrow = TRUE
      )
    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(
      strip.background = element_rect(
        fill = "grey95",
        color = "black",
        linewidth = 0.5
      ),

      strip.text = element_text(
        face = "bold",
        size = 10,
        color = "black"
      ),

      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = "black",
        size = 8.5
      ),

      axis.text.y = element_text(
        color = "black",
        size = 9
      ),

      axis.title.y = element_text(
        size = 10,
        margin = margin(
          r = 7
        )
      ),

      axis.ticks.x = element_blank(),

      panel.grid = element_blank(),

      legend.position = "right",

      legend.title = element_text(
        face = "bold",
        size = 9
      ),

      legend.text = element_text(
        size = 8.5
      ),

      legend.key.width = unit(
        0.45,
        "cm"
      ),

      legend.key.height = unit(
        0.38,
        "cm"
      ),

      legend.spacing.y = unit(
        0.05,
        "cm"
      ),

      plot.margin = margin(
        t = 5,
        r = 4,
        b = 5,
        l = 4
      )
    )

  ## 除每行第一个面板外，不重复显示Y轴
  if (!show.y.axis) {

    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      )
  }

  return(p)
}