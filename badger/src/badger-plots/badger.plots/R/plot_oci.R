#' Plot Confusion matrix
#'
#' @description
#' Plot a single confusion matrix using plotly, and a provided table.
#' Return a plotly figure object.
#'
#' @details
#' This confusion matrix is plotted using a heatmap trace, along with a set of
#' annotations, to highlight counts.
#' See: https://plotly.com/r/reference/heatmap/#heatmap
#'
#' @keywords internal
#' @importFrom scales rescale col_numeric
#' @param confusion_matrix An ordinal classification index confusion matrix.
#' @param yaxis_title Title for the yaxis of the confusion matrix
#' @param xaxis_title Title for the xaxis of the confusion matrix
#' @param xaxis_title_visibility Whether or not the xaxis should be rendered
#' @param show_xaxis Unilaterally choose to display or not xaxis
#' @param ratio if TRUE, rescale values of the CM to ratios
#' @param hide_zeroes if TRUE, cells of the CM containing no observations are
#'        not annotated.
#' @param colorscale specify the colorscale for relative frequency coloring.
#' @param ticklabels vector of CM xaxis and yaxis annotation labels
#' @param tickangle angle of the CM xaxis an yaxis annotation labels..
#' @param cm_font_size font size of the cell count annotations.
#' @param tickfont_size font size of the xaxis and yaxis annotations.
#' @param axis_fontsize font size of the main confusion matrix x- and y- titles
#' @return a plotly heatmap figure object.
plot_confusion_matrix <- function(
  confusion_matrix,
  yaxis_title            = NULL,
  xaxis_title            = NULL,
  xaxis_title_visibility = FALSE,
  show_xaxis             = FALSE,
  ratio                  = TRUE,
  hide_zeroes            = TRUE,
  colorscale             = "Blues",
  ticklabels             = NULL,
  tickangle              = 0L,
  cm_font_size           = 10L,
  tickfont_size          = 10L,
  axis_fontsize          = 14L
) {

  axis_font <- list(size = axis_fontsize)

  # ---- Add default ticklabels if none were specified...
  if (length(ticklabels) == 0L || is.null(ticklabels)) {
    r_to_rel <- as.list(paste0("<b>", rownames(confusion_matrix), "</b>"))
    names(r_to_rel) <- rownames(confusion_matrix)
  } else {
    r_to_rel <- ticklabels
    names(r_to_rel) <- rownames(confusion_matrix)
  }

  # ---- Convert counts to ratio if requested
  if (ratio) {
    confusion_matrix <- prop.table(confusion_matrix, margin = 1L)
  }

  # ---- Generate colorscale for error rates.
  scaled_vals    <- unique(
    scales::rescale(c(confusion_matrix) / rowSums(confusion_matrix))
  )
  o              <- order(scaled_vals, decreasing = FALSE)
  color_palette  <- scales::col_numeric(colorscale, domain = NULL)(scaled_vals)
  heatmap_colors <- stats::setNames(
    data.frame(scaled_vals[o], color_palette[o]), NULL
  )

  # Generate heatmap with plotly
  fig <- plotly::plot_ly(
    x          = unlist(r_to_rel[colnames(confusion_matrix)]),
    y          = unlist(r_to_rel[row.names(confusion_matrix)]),
    z          = confusion_matrix / rowSums(confusion_matrix),
    type       = "heatmap",
    colorscale = heatmap_colors,
    showscale  = FALSE
  )

  # Working with ratios, ensure the sum of all rows is preserved as 1.0
  annotated_cm <- confusion_matrix
  for (i in row.names(annotated_cm)) {
    annotated_cm[i, ] <- round_preserve_sum(annotated_cm[i, ], digits = 2L)
  }
  annotated_cm <- as.data.frame(annotated_cm)

  annotated_cm$expected_r <- unlist(
    r_to_rel[as.character(annotated_cm$expected_r)]
  )
  annotated_cm$observed_r <- unlist(
    r_to_rel[as.character(annotated_cm$observed_r)]
  )

  # color frequency annotations to white if the frequency is greated than
  # 0.4 times the max.
  annotated_cm$color <- ifelse(
    as.vector(confusion_matrix / rowSums(confusion_matrix) > 0.4),
    "white",
    "black"
  )

  if (hide_zeroes) {
    annotated_cm$Freq[annotated_cm$Freq == 0.0] <- ""
  }

  annotations <- list()
  for (i in seq_len(nrow(annotated_cm))) {
    annotations[[i]] <- list(
      text      = annotated_cm$Freq[i],
      x         = annotated_cm$observed_r[i],
      y         = annotated_cm$expected_r[i],
      font      = list(color = annotated_cm$color[i], size = cm_font_size),
      showarrow = FALSE
    )
  }

  fig <- fig %>% plotly::layout(
    xaxis       = list(
      title     = list(
        text = ifelse(xaxis_title_visibility, xaxis_title, ""),
        font = axis_font
      ),
      visible   = (xaxis_title_visibility | show_xaxis),
      tickangle = tickangle,
      tickfont  = list(size = tickfont_size),
      side      = "top"
    ),
    yaxis       = list(
      title     = list(text = yaxis_title, font = axis_font),
      tickfont  = list(size = tickfont_size),
      autorange = "reversed"
    ),
    annotations = annotations
  )

  fig
}

#' Plot Classification performance of a BADGER RUN.
#'
#' @description
#' Plot a list of OCI matrices as a set of confusion matrices and a scatterplot
#' of OCI or UOC values.
#'
#' @import plotly
#' @import viridis
#' @import RColorBrewer
#' @importFrom scales rescale
#' @importFrom purrr transpose
#' @param oci_results a named list of confusion matrices and performance index,
#'        typically obtained using `badger.plots::deserialize_results`.
#' @param toolnames name of each tool found within the oci_results graph.
#' @param filename name of the exported graph, when interactively exported
#'        through plotly's embedded html exporter.
#' @param transpose Whether structural transposition should be applied on the
#'        list of confusion matrices provided with `oci_results` before plotting
#'        This has the effect of transposing the ordering when plotting cms and
#'        the scatter plot.
#' @param scatterplot_ratio ratio between the subplot containing the confusion
#'        matrices the one containing the scatterplot. higher values will have
#'        the effect of increasing the height of the scatterplot.
#' @param horizontal_margin size of the horizontal margin separating each CM.
#'        Lower values will decrease the blank space separating every CM
#' @param vertical_margin size of the vertical margin separating each CM, and
#'        the margin separating the scatterplot from the rest of the CMs.
#'        Lower values will decrease the blank space separating these subplots
#' @param axis_fontsize Font size of every axis label and annotation.
#' @param legend A named list of parameters and arguments for the legend.
#'        These arguments are directly passed on to `plotly::layout`.
#'        See: https://plotly.com/r/reference/layout/#layout-legend
#' @param cm A named list of parameters and arguments specific to the
#'        plotting of the individual confusion matrices. These arguments are
#'        directly passed on to `bader.plots::plot_confusion_matrix`. See this
#'        function for additional information regarding these arguments.
#' @param scatter A named list of parameters and arguments specific to the
#'        plottlyng of the summary scatterplot. These arguments are directly
#'        passed on to `badger.plots::make_oci_scatter_plot`. See this function
#'        for additional details regarding these arguments.
#' @param ... Additional arguments. Not used at this point, but included for
#'        for future compatibility.
#' @return a plotly figure object
#' @export
plot_oci_performance <- function(
  oci_results,
  filename            = "OCI-performance-plot",
  transpose           = FALSE,
  scatterplot_ratio   = 0.20,
  horizontal_margin   = 0.02,
  vertical_margin     = 0.02,
  axis_fontsize       = 14L,
  legend              = list(size = 10L, xpos = -0.1, ypos = -0.1),
  cm                  = list(
    condense   = FALSE,
    ratio      = FALSE,
    colorscale = "Blues",
    ticklabels = NULL,
    tickfont   = 10L,
    tickangle  = 0L,
    fontsize   = 10L,
    show_xaxis = FALSE
  ),
  scatter = list(
    dash  = c("solid"), # nolint: unnecessary_concatenation_linter.
    mode  = NULL,
    yaxis = list(
      range = c(0.38, 1.02),
      tickfont = 10L,
      dtick = 0.1
    ),
    colors = NULL
  ),
  ...
) {
  oci_cms   <- list()
  tool_figures  <- list()

  toolnames <- names(oci_results[[1L]])

  if (transpose) {
    oci_results <- purrr::transpose(oci_results)
  }

  if (is.null(scatter$mode)) {
    scatter$mode <- ifelse(transpose, "markers", "lines+markers")
  }

  if (is.null(scatter$colors)) {
    scatter$colors <- if (transpose) {
      RColorBrewer::brewer.pal(n = length(names(oci_results)), name = "Set2")
    } else {
      viridis::viridis(n = length(names(oci_results[[1L]])))
    }
  }

  # ---- If the user requests non-observable classes to be condensed, check for
  # columns which only contain 0 values in every condition / tool.
  if (cm$condense) {
    keep <- sapply(oci_results, function(x) {
      lapply(x, function(y) {
        which(colSums(y$confusion_matrix) > 0.0)
      })
    }) %>%
      unlist() %>%
      unique()
  }

  # Generate an OCI confusion matrix for each tool and coverage.
  for (coverage in names(oci_results)){
    for (tool in names(oci_results[[coverage]])){

      plot_col         <- match(tool, names(oci_results[[coverage]]))
      confusion_matrix <- oci_results[[coverage]][[tool]]$confusion_matrix
      if (cm$condense) confusion_matrix <- confusion_matrix[, keep]

      oci_cms[[tool]][[coverage]] <- plot_confusion_matrix(
        confusion_matrix       = confusion_matrix,
        yaxis_title            = paste0("<b>", toolnames[plot_col], "</b>"),
        xaxis_title            = paste0("<b>", coverage, "</b>"),
        xaxis_title_visibility = (plot_col == 1L),
        show_xaxis             = cm$show_xaxis,
        axis_fontsize          = axis_fontsize,
        ratio                  = cm$ratio,
        colorscale             = cm$colorscale,
        cm_font_size           = cm$fontsize,
        ticklabels             = cm$ticklabels,
        tickfont_size          = cm$tickfont,
        tickangle              = cm$tickangle
      )
    }
  }

  # ---- Aggregate confusion matrices into subplots. One subplot per tool / row
  for (tool in names(oci_cms)){
    tool_figures[[tool]] <- plotly::subplot(
      oci_cms[[tool]],
      shareX = FALSE,
      shareY = TRUE,
      titleX = TRUE,
      margin = horizontal_margin
    )
  }

  # ---- Create a scatter plot summarizing OCI performance values for each tool
  #      and coverage.
  scatter_figure <- make_oci_scatter_plot(
    oci_results, toolnames = toolnames,
    colors         = scatter$colors,
    dash           = scatter$dash,
    mode           = scatter$mode,
    yaxis          = scatter$yaxis,
    axis_fontsize  = axis_fontsize,
    legend         = legend
  )

  # ---- Merge all plots as a single plot.
  rows             <- length(names(tool_figures))
  cm_height_ratios <- rep(1.0 / rows, rows)

  merged_confusion_matrices <- plotly::subplot(
    tool_figures,
    nrows   = rows,
    shareX  = FALSE,
    shareY  = FALSE,
    titleX  = TRUE,
    titleY  = TRUE,
    heights = cm_height_ratios,
    margin  = vertical_margin
  )

  merged_subplot <- plotly::subplot(
    merged_confusion_matrices, scatter_figure,
    nrows  = 2L,
    shareX = FALSE,
    shareY = FALSE,
    titleX = TRUE,
    titleY = TRUE,
    heights = c(1.0 - scatterplot_ratio, scatterplot_ratio),
    margin = 0L
  )

  merged_subplot
}

#' BADGER results scatter plot.
#'
#' @description
#' Generate a scatterplot of OCI performance values, according to coverage or
#' contamination.
#'
#' @keywords internal
#' @import plotly
#' @import RColorBrewer
#' @importFrom magrittr '%>%'
#' @param oci_results a list of ordinal classification index objects.
#' @param toolnames Corresponding names of each method. used to legend the plot
#' @param dash drawing style of each line. Passed on to plotly.
#'        See: https://plotly.com/r/reference/scatter/#scatter-line-dash
#' @param mode drawing mode of the scatter plot
#'        See: https://plotly.com/r/reference/scatter/#scatter-mode
#' @param colors list of HEX colors. One for every scatter line.
#' @param yaxis named list of parameters for the yaxis. passed on to
#'        plotly::layout.
#'        See: https://plotly.com/r/reference/layout/yaxis/#layout-yaxis
#' @param axis_fontsize fontsize of axis labels
#' @param tickfont_size fontsize of axis tickfons
#' @param legend named list of parameters for the legend. passed on to
#'       plotly::layout
#'       See: https://plotly.com/r/reference/layout/#layout-legend
#' @return a plotly scatterplot figure object
make_oci_scatter_plot <- function(
  oci_results = NULL,
  toolnames = NULL,
  dash = c("solid"), # nolint: unnecessary_concatenation_linter.
  mode = "lines+markers",
  colors = RColorBrewer::brewer.pal(
    n = length(names(oci_results)),
    name = "Set2"
  ),
  yaxis = list(range = c(0.38, 1.02), dtick = 0.1, tickfont = 10L),
  axis_fontsize  = 16L,
  tickfont_size  = 14L,
  legend = list(size = 10L, xpos = -0.1, ypos = -0.1)
) {

  oci_values <- list()
  oci_errors <- list()
  axis_font <- list(size = axis_fontsize)

  for (coverage in names(oci_results)){
    for (tool in names(oci_results[[coverage]])){
      oci_values[[tool]] <- c(
        oci_values[[tool]],
        oci_results[[coverage]][[tool]]$value
      )

      oci_errors[[tool]] <- c(
        oci_errors[[tool]],
        oci_results[[coverage]][[tool]]$error
      )
    }
  }


  oci_values <- as.data.frame(oci_values, row.names = names(oci_results))

  # recycle toolnames
  dash <- rep(dash, length(toolnames))

  # check mode validity
  valid_modes <- c("lines", "markers", "lines+markers")
  if (!(mode %in% valid_modes))
    stop(
      cat("Invalid scatter mode. Valid values:",
        paste(valid_modes, sep = ", "), "\n"
      )
    )


  for (i in seq_along(toolnames)){
    line <- if (grepl("line", mode, fixed = TRUE)) {
      list(color = colors[i], dash = dash[i])
    }
    markers <- if (grepl("markers", mode, fixed = TRUE)) list(color = colors[i])

    if (i == 1L) {
      fig <- plotly::plot_ly(
        data   = oci_values,
        y      = (1.0 - oci_values[[i]]),
        x      = ~row.names(oci_values),
        type   = "scatter",
        line   = line,
        marker = markers,
        mode   = mode,
        name   = toolnames[i],
        error_y = list(
          array     = oci_errors[[i]],
          type      = "data",
          color     = "#111111",
          thickness = 1L
        )
      )
    } else {
      fig <- fig %>% plotly::add_trace(
        data   = oci_values,
        y      = (1.0 - oci_values[[i]]),
        x      = ~row.names(oci_values),
        type   = "scatter",
        line   = line,
        marker = markers,
        mode   = mode,
        name   = toolnames[i],
        error_y = list(
          array = oci_errors[[i]],
          type  = "data",
          color = "#111111"
        )
      )
    }
  }

  fig <- fig %>% plotly::layout(
    xaxis = list(
      title           = list(
        text = "<b>    </b>",
        font = axis_font
      ),
      tickfont        = axis_font,
      tickvals        = levels(factor(row.names(oci_values))),
      ticktext        = sprintf(
        "<b>%s</b>", levels(factor(row.names(oci_values)))
      ),
      range           = c(-0.5, length(row.names(oci_values)) - 0.5),
      xref            = "paper",
      autotypenumbers = "strict" # Force categorical axis
    ),
    yaxis = list(
      title = list(
        text = paste0(
          "<b>1 - ", toupper(oci_results[[1L]][[1L]]$method), "</b>"
        ),
        font = axis_font
      ),
      range   = yaxis$range,
      tickfont = list(size = yaxis$tickfont),
      dtick = yaxis$dtick
    ),
    legend = list(
      x           = legend$xpos,
      y           = legend$ypos,
      xref        = "paper",
      yref        = "paper",
      orientation = "h",
      font        = list(size = legend$size)
    )
  ) %>%
    plotly::config(
      editable    = TRUE,
      displaylogo = FALSE,
      scrollZoom  = TRUE,
      toImageButtonOptions = list(
        format   = "svg",
        filename = paste0("OCI-performance-plot")
      )
    )
}

#' Round preserve sum
#'
#' @description
#' Round a vector of numbers, while preserving their sum to 1.0
#' @keywords internal
#' @param x a vector of numeric values
#' @param digits number of kept decimal places for rounding.
#' @return a vector of rounded numeric values, whose sum equals 1.
round_preserve_sum <- function(x, digits = 0L) {
  up         <- 10L^digits
  x          <- x * up
  y          <- floor(x)
  indices    <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1L
  y / up
}