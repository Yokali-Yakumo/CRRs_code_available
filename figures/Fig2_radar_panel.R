library(dplyr)
library(fmsb)
library(scales)

# ============================================================
# 增强时间差异的雷达图
#
# 标准化方式：
# 对每个 mark，在所有 cluster 和所有时间点之间进行全局 min-max scaling
#
# 保证：
# 1. 同一 mark 不同时间点的差异更明显
# 2. 不同 cluster 的雷达图仍使用完全相同的尺度
# 3. 保留原代码块1的颜色、线型和填充样式
# ============================================================

plot_cluster_radars <- function(
    orig_df,
    cluster_labels,
    outdir = "outplots/radar_publication",
    mark_order = c(
      "H3K27ac", "DNase", "H3K4me1", "H3K9ac",
      "H3K36me3", "H3K4me3", "H3K27me3", "H3K9me3"
    ),
    axis_expand = 1.05,
    font_family = "Arial",
    pdf_width = 5,
    pdf_height = 5,
    show_legend = FALSE
) {

  # ----------------------------------------------------------
  # 1. 输入检查
  # ----------------------------------------------------------
  if (length(cluster_labels) != nrow(orig_df)) {
    stop(
      "cluster_labels 长度为 ", length(cluster_labels),
      "，orig_df 行数为 ", nrow(orig_df),
      "，两者不一致。"
    )
  }

  # 优先按照名称对齐 cluster_labels
  if (
    !is.null(names(cluster_labels)) &&
    !is.null(rownames(orig_df)) &&
    all(rownames(orig_df) %in% names(cluster_labels))
  ) {
    orig_df$cluster <- unname(
      cluster_labels[rownames(orig_df)]
    )
  } else {
    orig_df$cluster <- unname(cluster_labels)
  }

  # ----------------------------------------------------------
  # 2. 识别 coverage 列
  # ----------------------------------------------------------
  cov_cols <- grep(
    "_d[0-9]+_coverage$",
    colnames(orig_df),
    value = TRUE
  )

  if (length(cov_cols) == 0) {
    stop(
      "未找到形如 H3K27ac_d0_coverage 的列，",
      "请检查列名格式。"
    )
  }

  # ----------------------------------------------------------
  # 3. 提取时间点
  # ----------------------------------------------------------
  times <- unique(
    sub(".*_(d[0-9]+)_coverage$", "\\1", cov_cols)
  )

  times <- times[
    order(as.numeric(sub("^d", "", times)))
  ]

  if (length(times) > 3) {
    stop(
      "当前颜色和线型按照最多3个时间点设置，",
      "但检测到：",
      paste(times, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 4. 确定 mark 及其排列顺序
  # ----------------------------------------------------------
  available_marks <- unique(
    sub("_d[0-9]+_coverage$", "", cov_cols)
  )

  marks <- mark_order[
    mark_order %in% available_marks
  ]

  if (length(marks) < 3) {
    stop("实际识别出的 mark 少于3个，无法绘制雷达图。")
  }

  missing_marks <- setdiff(
    mark_order,
    available_marks
  )

  if (length(missing_marks) > 0) {
    message(
      "以下 mark 未在数据中找到，将被忽略：",
      paste(missing_marks, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 5. 检查所有 mark × time 组合是否完整
  # ----------------------------------------------------------
  expected_cols <- unlist(
    lapply(
      marks,
      function(m) {
        paste0(m, "_", times, "_coverage")
      }
    )
  )

  missing_cols <- setdiff(
    expected_cols,
    colnames(orig_df)
  )

  if (length(missing_cols) > 0) {
    stop(
      "以下 coverage 列缺失：\n",
      paste(missing_cols, collapse = "\n")
    )
  }

  # ----------------------------------------------------------
  # 6. 按 cluster 聚合
  # ----------------------------------------------------------
  safe_mean <- function(x) {

    if (all(is.na(x))) {
      return(NA_real_)
    }

    mean(x, na.rm = TRUE)
  }

  cluster_cov <- orig_df %>%
    dplyr::select(cluster, all_of(expected_cols)) %>%
    group_by(cluster) %>%
    summarise(
      across(
        all_of(expected_cols),
        safe_mean
      ),
      .groups = "drop"
    ) %>%
    as.data.frame()

  # ----------------------------------------------------------
  # 7. 检查缺失值
  # ----------------------------------------------------------
  coverage_matrix <- as.matrix(
    cluster_cov[, expected_cols, drop = FALSE]
  )

  if (any(!is.finite(coverage_matrix))) {

    bad_cols <- expected_cols[
      colSums(!is.finite(coverage_matrix)) > 0
    ]

    stop(
      "cluster 均值中存在 NA、NaN 或 Inf。\n",
      "涉及列：\n",
      paste(bad_cols, collapse = "\n"),
      "\n不建议将缺失值直接当作0。"
    )
  }

  # ----------------------------------------------------------
  # 8. 计算每个 mark 的全局 minimum 和 maximum
  #
  # 注意：
  # 范围跨越所有 cluster 和所有时间点
  #
  # 因此不同 cluster 的图仍然可以相互比较
  # ----------------------------------------------------------
  global_min_mark <- setNames(
    numeric(length(marks)),
    marks
  )

  global_max_mark <- setNames(
    numeric(length(marks)),
    marks
  )

  for (m in marks) {

    mark_cols <- paste0(
      m, "_", times, "_coverage"
    )

    mark_values <- unlist(
      cluster_cov[, mark_cols, drop = FALSE],
      use.names = FALSE
    )

    global_min_mark[m] <- min(
      mark_values,
      na.rm = TRUE
    )

    global_max_mark[m] <- max(
      mark_values,
      na.rm = TRUE
    )
  }

  # 每个 mark 的全局动态范围
  global_range_mark <- global_max_mark - global_min_mark

  # 防止某个 mark 所有值完全一致，导致除以0
  zero_range_marks <- names(global_range_mark)[
    !is.finite(global_range_mark) |
      global_range_mark <= .Machine$double.eps
  ]

  if (length(zero_range_marks) > 0) {

    warning(
      "以下 mark 的全局动态范围为0，",
      "缩放分母将设为1：",
      paste(zero_range_marks, collapse = ", ")
    )

    global_range_mark[zero_range_marks] <- 1
  }

  # ----------------------------------------------------------
  # 9. 创建输出目录并保存原始数据及尺度参数
  # ----------------------------------------------------------
  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    cluster_cov,
    file = file.path(
      outdir,
      "cluster_mean_raw_coverage.csv"
    ),
    row.names = FALSE
  )

  scale_table <- data.frame(
    mark = marks,
    global_minimum = unname(
      global_min_mark[marks]
    ),
    global_maximum = unname(
      global_max_mark[marks]
    ),
    global_range = unname(
      global_range_mark[marks]
    ),
    scaling_method = paste0(
      "(value - global minimum) / ",
      "(global maximum - global minimum)"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    scale_table,
    file = file.path(
      outdir,
      "radar_scaling_parameters.csv"
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 10. 绘图样式
  #
  # 完全保留代码块1的颜色和线型
  # ----------------------------------------------------------

    # 色盲友好配色：
    # Day 0 = 蓝；Day 2 = 橙；Day 7 = 紫红
    time_colors <- c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#CC79A7"   # reddish purple
    )[seq_along(times)]

    # 不使用填充色
    fill_colors <- rep(
    NA_character_,
    length(times)
    )

    # 保持线型差异；即使灰度打印也能区分
    line_types <- c(
    1,  # solid
    2,  # dashed
    3   # dotted
    )[seq_along(times)]

    # 不同时间点采用不同点形状
    point_types <- c(
    16,  # circle
    17,  # triangle
    15   # square
    )[seq_along(times)]


  time_labels <- paste0(
    "Day ",
    sub("^d", "", times)
  )

  # ----------------------------------------------------------
  # 11. 单个 cluster 的绘图函数
  # ----------------------------------------------------------
  draw_radar <- function(cluster_index) {

    cluster_name <- as.character(
      cluster_cov$cluster[cluster_index]
    )

    # 原始 time × mark 数据矩阵
    radar_raw <- matrix(
      NA_real_,
      nrow = length(times),
      ncol = length(marks),
      dimnames = list(times, marks)
    )

    for (ti in seq_along(times)) {

      for (mj in seq_along(marks)) {

        colname <- paste0(
          marks[mj],
          "_",
          times[ti],
          "_coverage"
        )

        radar_raw[ti, mj] <- cluster_cov[
          cluster_index,
          colname
        ]
      }
    }

    # --------------------------------------------------------
    # 全局 min-max scaling
    #
    # 所有 cluster 使用同一组 global_min/global_max
    # --------------------------------------------------------
    radar_scaled <- sweep(
      radar_raw,
      MARGIN = 2,
      STATS = global_min_mark[marks],
      FUN = "-"
    )

    radar_scaled <- sweep(
      radar_scaled,
      MARGIN = 2,
      STATS = global_range_mark[marks],
      FUN = "/"
    )

    # 浮点误差保护
    radar_scaled[radar_scaled < 0] <- 0
    radar_scaled[radar_scaled > 1] <- 1

    # fmsb 要求第一行是最大值，第二行是最小值
    plot_data <- as.data.frame(
      rbind(
        max = rep(
          axis_expand,
          length(marks)
        ),
        min = rep(
          0,
          length(marks)
        ),
        radar_scaled
      )
    )

    colnames(plot_data) <- marks

    # --------------------------------------------------------
    # 页面边距
    #
    # 左右边距略加大，减少 H3K4me1 等标签被截断
    # --------------------------------------------------------
    # --------------------------------------------------------
    # 页面边距
    # --------------------------------------------------------
    if (show_legend) {

    par(
        mar = c(5.2, 5.5, 4.2, 5.5),
        pty = "s",
        xpd = NA,
        family = font_family
    )

    } else {

    par(
        mar = c(4.5, 5.5, 4.2, 5.5),
        pty = "s",
        xpd = NA,
        family = font_family
    )
    }

    # --------------------------------------------------------
    # 雷达图
    #
    # 保留原代码块1的样式
    # --------------------------------------------------------
    radarchart(
    plot_data,

    # 不显示径向数值刻度，使图更简洁
    axistype = 0,

    # ------------------------------------------------------
    # 数据曲线：无填充，线条和点更突出
    # ------------------------------------------------------
    pcol = time_colors,
    pfcol = fill_colors,
    plwd = 3.8,
    plty = line_types,

    # 数据点：稍大、稍醒目
    pty = point_types,
    pcex = 1.35,

    # ------------------------------------------------------
    # 网格：减少圈数、显著弱化
    # seg = 3 即仅保留较少的内部同心多边形
    # ------------------------------------------------------
    cglcol = "#D9D9D9",
    cglwd = 0.7,
    seg = 2,
    cglty = 1,


    # 轴线采用更浅灰色
    axislabcol = "#666666",

    # ------------------------------------------------------
    # Marker 标签：全部使用黑色，不与曲线竞争
    # ------------------------------------------------------
    vlabels = marks,
    vlcex = 1.05
    )

    # 可选图例
    if (show_legend) {

      legend(
        x = "bottom",
        inset = c(0, -0.16),
        legend = time_labels,
        bty = "n",
        pch = 16,
        lty = line_types,
        lwd = 3,
        col = time_colors,
        cex = 0.9,
        pt.cex = 1.2,
        horiz = TRUE,
        xpd = NA
      )
    }
  }

  # ----------------------------------------------------------
  # 12. 输出每个 cluster 的 PDF 和 PNG
  # ----------------------------------------------------------
  for (i in seq_len(nrow(cluster_cov))) {

    cluster_name <- as.character(
      cluster_cov$cluster[i]
    )

    safe_cluster_name <- gsub(
      "[^A-Za-z0-9_.-]",
      "_",
      cluster_name
    )

    pdf_file <- file.path(
      outdir,
      paste0(
        "cluster_",
        safe_cluster_name,
        "_radar.pdf"
      )
    )

    png_file <- file.path(
      outdir,
      paste0(
        "cluster_",
        safe_cluster_name,
        "_radar.png"
      )
    )

    # 矢量PDF
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = pdf_width,
      height = pdf_height,
      family = font_family
    )

    draw_radar(i)

    dev.off()

    # 高分辨率PNG
    png(
      filename = png_file,
      width = 2000,
      height = 2000,
      res = 400,
      type = "cairo",
      family = font_family
    )

    draw_radar(i)

    dev.off()

    message("已保存：", pdf_file)
    message("已保存：", png_file)
  }

  # ----------------------------------------------------------
  # 13. 多页PDF
  # ----------------------------------------------------------
  all_pdf <- file.path(
    outdir,
    "all_clusters_radar.pdf"
  )

  grDevices::cairo_pdf(
    filename = all_pdf,
    width = pdf_width,
    height = pdf_height,
    family = font_family,
    onefile = TRUE
  )

  for (i in seq_len(nrow(cluster_cov))) {
    draw_radar(i)
  }

  dev.off()

  message(
    "全部完成：共生成 ",
    nrow(cluster_cov),
    " 个 cluster 的雷达图。"
  )

  invisible(
    list(
      cluster_mean = cluster_cov,
      global_min_mark = global_min_mark,
      global_max_mark = global_max_mark,
      global_range_mark = global_range_mark,
      marks = marks,
      times = times
    )
  )
}


# ============================================================
# 运行
# ============================================================

radar_result <- plot_cluster_radars(
  orig_df = orig_df,
  cluster_labels = cluster_labels,
  outdir = "outplots/radar_publication",

  # 原来是1.3，会明显压缩差异
  # 改为1.05，仅保留少量外部空间
  axis_expand = 1.05,

  show_legend = FALSE
)
