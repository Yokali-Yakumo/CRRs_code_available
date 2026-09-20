
#
# Helper: plotting functions for chromHMM-state enrichment heatmaps from an enrichment table (log2 fold-enrichment and FDR).
#
library(ggplot2)
library(dplyr)
library(scales)
library(grid)

plot_chromHMM_enrichment <- function(
    enrichment.data,
    effect.col = "log2_fold_enrichment",
    fdr.col = "FDR_global",
    color.limit = 4,
    min.effect.for.star = 1,

    annotation.order = c(
      "Active promoter",
      "Active enhancer",
      "Primed enhancer",
      "Transcribed",
      "Mixed/Bivalent",
      "Polycomb repressed",
      "Heterochromatin/repeats",
      "Quiescent"
    ),

    cluster.order = NULL,

    time.order = c(
      "MSC",
      "ADI-7d",
      "ADI-14d"
    ),

    low.color = "#2E75B6",
    mid.color = "white",
    high.color = "#980000"
) {

  ## 1. Check required columns
  required.columns <- c(
    "Time",
    "Cluster",
    "Annotation",
    effect.col,
    fdr.col
  )

  missing.columns <- setdiff(
    required.columns,
    colnames(enrichment.data)
  )

  if (length(missing.columns) > 0) {
    stop(
      "The input data is missing the following columns: ",
      paste(missing.columns, collapse = ", ")
    )
  }

  if (!is.numeric(color.limit) ||
      length(color.limit) != 1 ||
      color.limit <= 0) {
    stop("color.limit must be a single numeric value greater than 0.")
  }

  ## 2. Clean character columns
  plot.data <- enrichment.data %>%
    mutate(
      Time = trimws(as.character(Time)),
      Cluster = trimws(as.character(Cluster)),
      Annotation = trimws(as.character(Annotation)),
      effect_raw = as.numeric(.data[[effect.col]]),
      FDR_plot = as.numeric(.data[[fdr.col]])
    )

  ## 3. Derive the cluster order automatically
  ## Natural sort on the trailing number of the cluster name
  if (is.null(cluster.order)) {

    cluster.values <- unique(plot.data$Cluster)

    cluster.number <- suppressWarnings(
      as.numeric(
        sub(
          pattern = ".*?([0-9]+)$",
          replacement = "\\1",
          x = cluster.values
        )
      )
    )

    if (all(!is.na(cluster.number))) {
      cluster.order <- cluster.values[
        order(cluster.number)
      ]
    } else {
      cluster.order <- sort(cluster.values)
    }
  }

  ## 4. Time-point order
  observed.times <- unique(plot.data$Time)

  time.order <- c(
    intersect(time.order, observed.times),
    setdiff(observed.times, time.order)
  )

  ## 5. Annotation order
  observed.annotations <- unique(
    plot.data$Annotation
  )

  annotation.order <- c(
    intersect(annotation.order, observed.annotations),
    setdiff(observed.annotations, annotation.order)
  )

  ## 6. Convert to factors and truncate the plotted effect sizes
  plot.data <- plot.data %>%
    mutate(
      Time = factor(
        Time,
        levels = time.order
      ),

      Cluster = factor(
        Cluster,
        levels = cluster.order
      ),

      ## rev puts the first entry of annotation.order at the top
      Annotation = factor(
        Annotation,
        levels = rev(annotation.order)
      ),

      ## Affects only the plotted colours, not the original effect sizes
      effect_plot = pmax(
        pmin(effect_raw, color.limit),
        -color.limit
      )
    )

  ## 7. Check for NAs introduced by mismatched factor levels
  if (any(is.na(plot.data$Cluster))) {
    stop(
      "NA produced after converting Cluster to a factor; check cluster.order. Original values: ",
      paste(
        unique(enrichment.data$Cluster),
        collapse = ", "
      )
    )
  }

  if (any(is.na(plot.data$Time))) {
    stop(
      "NA produced after converting Time to a factor; check time.order. Original values: ",
      paste(
        unique(enrichment.data$Time),
        collapse = ", "
      )
    )
  }

  if (any(is.na(plot.data$Annotation))) {
    stop(
      "NA produced after converting Annotation to a factor; check annotation.order. Original values: ",
      paste(
        unique(enrichment.data$Annotation),
        collapse = ", "
      )
    )
  }

  ## 8. Significance markers
  plot.data <- plot.data %>%
    mutate(
      significance = case_when(
        is.na(FDR_plot) |
          is.na(effect_raw) ~ "",

        abs(effect_raw) < min.effect.for.star ~ "",

        FDR_plot < 0.001 ~ "***",
        FDR_plot < 0.01  ~ "**",
        FDR_plot < 0.05  ~ "*",

        TRUE ~ ""
      ),

      significance_color = ifelse(
        abs(effect_plot) >= color.limit * 0.55,
        "white",
        "black"
      )
    )

  ## Use ASCII characters to avoid PDF font warnings
  legend.breaks <- c(
    -color.limit,
    -color.limit / 2,
    0,
    color.limit / 2,
    color.limit
  )

  legend.labels <- c(
    paste0("<= ", -color.limit),
    as.character(-color.limit / 2),
    "0",
    as.character(color.limit / 2),
    paste0(">= ", color.limit)
  )

  p <- ggplot(
    plot.data,
    aes(
      x = Cluster,
      y = Annotation,
      fill = effect_plot
    )
  ) +

    geom_tile(
      color = "grey85",
      linewidth = 0.35,
      width = 0.96,
      height = 0.96
    ) +

    geom_text(
      aes(
        label = significance,
        color = significance_color
      ),
      size = 3.4,
      fontface = "bold",
      show.legend = FALSE
    ) +

    scale_color_identity() +

    scale_fill_gradient2(
      low = low.color,
      mid = mid.color,
      high = high.color,
      midpoint = 0,
      limits = c(
        -color.limit,
        color.limit
      ),
      breaks = legend.breaks,
      labels = legend.labels,
      oob = scales::squish,
      name = expression(log[2]~"fold enrichment")
    ) +

    facet_grid(
      cols = vars(Time),
      scales = "fixed",
      space = "fixed"
    ) +

    scale_x_discrete(
      drop = FALSE
    ) +

    scale_y_discrete(
      drop = FALSE
    ) +

    labs(
      x = NULL,
      y = NULL,
      caption = paste0(
        "Color scale truncated at -",
        color.limit,
        " and +",
        color.limit,
        ". ",
        "* FDR < 0.05; ** FDR < 0.01; ",
        "*** FDR < 0.001. ",
        "Stars shown only when |log2 fold enrichment| >= ",
        min.effect.for.star,
        "."
      )
    ) +

    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(8.5, "cm"),
        barheight = unit(0.42, "cm"),
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      panel.grid = element_blank(),

      panel.border = element_rect(
        color = "black",
        linewidth = 0.5
      ),

      strip.background = element_rect(
        fill = "grey95",
        color = "black",
        linewidth = 0.5
      ),

      strip.text = element_text(
        face = "bold",
        size = 11,
        color = "black"
      ),

      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = "black",
        size = 9.5
      ),

      axis.text.y = element_text(
        color = "black",
        size = 9.5
      ),

      axis.ticks = element_blank(),

      legend.position = "bottom",

      legend.title = element_text(
        size = 10
      ),

      legend.text = element_text(
        size = 9
      ),

      plot.caption = element_text(
        size = 8.5,
        hjust = 0,
        color = "grey30"
      ),

      panel.spacing.x = unit(
        0.8,
        "lines"
      )
    )

  return(p)
}
