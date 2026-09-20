library(dplyr)
library(fmsb)
library(scales)

# ============================================================
# Radar plots that emphasize temporal differences
#
# Normalization scheme:
# For each mark, global min-max scaling across all clusters and all time points
#
# This ensures:
# 1. Differences between time points are more obvious for the same mark
# 2. Radar plots of different clusters still use exactly the same scale
# 3. Colors, line types and fill styles of original code block 1 are preserved
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
  # 1. Input checks
  # ----------------------------------------------------------
  if (length(cluster_labels) != nrow(orig_df)) {
    stop(
      "cluster_labels has length ", length(cluster_labels),
      " but orig_df has ", nrow(orig_df),
      " rows, so the two are inconsistent."
    )
  }

  # Prefer aligning cluster_labels by name
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
  # 2. Identify coverage columns
  # ----------------------------------------------------------
  cov_cols <- grep(
    "_d[0-9]+_coverage$",
    colnames(orig_df),
    value = TRUE
  )

  if (length(cov_cols) == 0) {
    stop(
      "No columns of the form H3K27ac_d0_coverage were found; ",
      "please check the column name format."
    )
  }

  # ----------------------------------------------------------
  # 3. Extract time points
  # ----------------------------------------------------------
  times <- unique(
    sub(".*_(d[0-9]+)_coverage$", "\\1", cov_cols)
  )

  times <- times[
    order(as.numeric(sub("^d", "", times)))
  ]

  if (length(times) > 3) {
    stop(
      "Colors and line types are defined for at most 3 time points, ",
      "but detected: ",
      paste(times, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 4. Determine the marks and their ordering
  # ----------------------------------------------------------
  available_marks <- unique(
    sub("_d[0-9]+_coverage$", "", cov_cols)
  )

  marks <- mark_order[
    mark_order %in% available_marks
  ]

  if (length(marks) < 3) {
    stop("Fewer than 3 marks were detected, cannot draw a radar plot.")
  }

  missing_marks <- setdiff(
    mark_order,
    available_marks
  )

  if (length(missing_marks) > 0) {
    message(
      "The following marks were not found in the data and will be ignored: ",
      paste(missing_marks, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 5. Check that all mark × time combinations are complete
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
      "The following coverage columns are missing:\n",
      paste(missing_cols, collapse = "\n")
    )
  }

  # ----------------------------------------------------------
  # 6. Aggregate by cluster
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
  # 7. Check for missing values
  # ----------------------------------------------------------
  coverage_matrix <- as.matrix(
    cluster_cov[, expected_cols, drop = FALSE]
  )

  if (any(!is.finite(coverage_matrix))) {

    bad_cols <- expected_cols[
      colSums(!is.finite(coverage_matrix)) > 0
    ]

    stop(
      "NA, NaN or Inf found in the cluster means.\n",
      "Affected columns:\n",
      paste(bad_cols, collapse = "\n"),
      "\nTreating missing values directly as 0 is not recommended."
    )
  }

  # ----------------------------------------------------------
  # 8. Compute the global minimum and maximum for each mark
  #
  # Note:
  # The range spans all clusters and all time points
  #
  # so plots of different clusters remain mutually comparable
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

  # Global dynamic range of each mark
  global_range_mark <- global_max_mark - global_min_mark

  # Guard against division by zero when all values of a mark are identical
  zero_range_marks <- names(global_range_mark)[
    !is.finite(global_range_mark) |
      global_range_mark <= .Machine$double.eps
  ]

  if (length(zero_range_marks) > 0) {

    warning(
      "The following marks have zero global dynamic range; ",
      "the scaling denominator will be set to 1: ",
      paste(zero_range_marks, collapse = ", ")
    )

    global_range_mark[zero_range_marks] <- 1
  }

  # ----------------------------------------------------------
  # 9. Create the output directory and save raw values plus scaling parameters
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
  # 10. Plot styling
  #
  # Colors and line types of code block 1 are preserved exactly
  # ----------------------------------------------------------

    # Colorblind-friendly palette:
    # Day 0 = blue; Day 2 = orange; Day 7 = reddish purple
    time_colors <- c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#CC79A7"   # reddish purple
    )[seq_along(times)]

    # No fill colors
    fill_colors <- rep(
    NA_character_,
    length(times)
    )

    # Keep line types distinct; still distinguishable in grayscale printing
    line_types <- c(
    1,  # solid
    2,  # dashed
    3   # dotted
    )[seq_along(times)]

    # Different point shapes for different time points
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
  # 11. Plotting function for a single cluster
  # ----------------------------------------------------------
  draw_radar <- function(cluster_index) {

    cluster_name <- as.character(
      cluster_cov$cluster[cluster_index]
    )

    # Raw time × mark data matrix
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
    # Global min-max scaling
    #
    # All clusters use the same global_min/global_max
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

    # Floating-point error guard
    radar_scaled[radar_scaled < 0] <- 0
    radar_scaled[radar_scaled > 1] <- 1

    # fmsb requires the first row to be the maximum and the second the minimum
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
    # Plot margins
    #
    # Slightly wider left and right margins reduce clipping of labels such as H3K4me1
    # --------------------------------------------------------
    # --------------------------------------------------------
    # Plot margins
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
    # Radar chart
    #
    # Preserve the styling of original code block 1
    # --------------------------------------------------------
    radarchart(
    plot_data,

    # Hide radial tick labels for a cleaner plot
    axistype = 0,

    # ------------------------------------------------------
    # Data curves: no fill, lines and points stand out more
    # ------------------------------------------------------
    pcol = time_colors,
    pfcol = fill_colors,
    plwd = 3.8,
    plty = line_types,

    # Data points: slightly larger and more prominent
    pty = point_types,
    pcex = 1.35,

    # ------------------------------------------------------
    # Grid: fewer rings, markedly subdued
    # seg = 3 keeps only a few inner concentric polygons
    # ------------------------------------------------------
    cglcol = "#D9D9D9",
    cglwd = 0.7,
    seg = 2,
    cglty = 1,


    # Axis lines use a lighter gray
    axislabcol = "#666666",

    # ------------------------------------------------------
    # Marker labels: all black so they do not compete with the curves
    # ------------------------------------------------------
    vlabels = marks,
    vlcex = 1.05
    )

    # Optional legend
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
  # 12. Write a PDF and PNG for each cluster
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

    # Vector PDF
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = pdf_width,
      height = pdf_height,
      family = font_family
    )

    draw_radar(i)

    dev.off()

    # High-resolution PNG
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

    message("Saved: ", pdf_file)
    message("Saved: ", png_file)
  }

  # ----------------------------------------------------------
  # 13. Multi-page PDF
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
    "All done: generated ",
    nrow(cluster_cov),
    " cluster radar plots."
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
# Run
# ============================================================

radar_result <- plot_cluster_radars(
  orig_df = orig_df,
  cluster_labels = cluster_labels,
  outdir = "outplots/radar_publication",

  # Was 1.3, which visibly compressed the differences
  # Changed to 1.05 to leave only a little outer space
  axis_expand = 1.05,

  show_legend = FALSE
)
