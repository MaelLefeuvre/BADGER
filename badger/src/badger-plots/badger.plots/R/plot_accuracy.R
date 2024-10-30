#' Normalized RMSD and MBE bartrace plot.
#'
#' @description
#' Loop along a named nested list of summary statistics dataframes and generate
#' a plot displaying RMSD and MBE estimates for each tool, at each coverage.
#' @importFrom plotly plot_ly layout add_annotations subplot config
#' @importFrom purrr transpose
#' @importFrom viridis viridis
#' @param rmsd_dfs a named list of dataframes, containing summary RMSD and MBE
#'        values for every studied BADGER run. Typically obtained using
#'        `badger.plots::deserialize_results`
#' @param filename name of the exported graph, when interactively exported
#'        through plotly(s embedded html exporter.)
#' @param transpose Logical specifying whether structuram transposition should
#'        be applied on the provided list of dataframes provided by `rmsd_dfs`
#'        prior to plotting. This has the effect of transposingthe ordering when
#'        plotting individual bartraces and subplots.
#' @param fixed_axis Logical specifying whether the y-axis of exery subplot
#'        should be fixed at a set range. when TRUE, the set range is selected
#'        using the maximum and minimum value found across all provided values
#' @param title Main title
#' @param rmsd A named list of parameters and arguments specific to the nRMSD
#'        subplot. These arguments are mostly passed on to plotly::layout.
#'        See: https://plotly.com/r/reference/layout
#' @param mbe A named list of parameters and arguments specific to the nMBE
#'        subplot. These arguments are mostly passed on to plotly::layout.
#'        See: https://plotly.com/r/reference/layout
#' @param mbe_plot_ratio ratio between the RMSD and MBE subplot. higher values
#'        will have the effect of increasing the width of the MBE subplot.
#' @param ticksize font sizes of the x- and y-axis ticks.
#' @param legend A named list of parameters and arguments for the legend. These
#'        arguments are directly passed on to `plotly::layout`
#'        See: https://plotly.com/r/reference/layout/#layout-legend
#' @param marker A named list of parameters and arguments for the markers.
#'        These arguments are directly passed of to `plotly::plot_ly`.
#'        See: https://plotly.com/r/reference/bar/#bar-marker
#' @param ... Additional arguments. Not used at this point, and merely included
#'        for future compatibility.
#' @return a plotly figure object
#' @export
plot_accuracy <- function(
  rmsd_dfs,
  filename         = "nRMSD-accuracy-plot",
  transpose        = FALSE,
  fixed_axis       = TRUE,
  title            = NULL,
  rmsd             = list(title = "nRMSD", dtick = 0.2, tickangle = 0L),
  mbe              = list(
    title  = "nMBE", dtick = 0.2, tickangle = 45L, scale_factor = 1L
  ),
  mbe_plot_ratio   = 0.3,
  ticksize         = list(xaxis = 16L, yaxis = 12L),
  legend           = list(
    size  = 12L,
    title = "Method"
  ),
  marker           = list(
    colors   = NULL,
    patterns = list(shape = c(""), solidity = 0.5, size = 5L)
  ),
  ...
) {
  font <- list(
    family = "Courier New, monospace",
    size   = 20L,
    color  = "#8F8F8F"
  )

  # ---- Transpose data frame if requested
  if (transpose) rmsd_dfs <- purrr::transpose(rmsd_dfs)

  # ---- Add scale factor to yaxis title if different from one.
  mbe$title <- paste0(
    "nMBE",
    ifelse(mbe$scale_factor != 1L, paste0(" (x", mbe$scale_factor, ")"), "")
  )

  rmsd_barplots <- list()

  if (is.null(marker$colors)) {
    marker$colors <- if (transpose) {
      RColorBrewer::brewer.pal(n = length(names(rmsd_dfs)), name = "Set2")
    } else {
      viridis::viridis(n = length(names(rmsd_dfs[[1L]])))
    }
  }

  for (column in names(rmsd_dfs)) {
    acc_df     <- list()
    rmsd_lines <- list()
    mbe_lines  <- list()

    # ---- recycle hatchings
    # See https://plotly.com/r/reference/bar/#bar-marker-pattern
    marker$patterns$shape <- rep(
      marker$patterns$shape, length(names(rmsd_dfs[[column]]))
    )

    for (row in names(rmsd_dfs[[column]])) {
      i <- match(column, names(rmsd_dfs))
      j <- match(row, names(rmsd_dfs[[column]]))
      show_legend <- (i == 1L)

      # ---- extract the relevant data frame of RMSE and MBE values.
      acc_df[[j]] <- rmsd_dfs[[column]][[row]]

      # ---- rescale mbe by its scale factor
      mbe_cols <- c("mbe", "mbe.upper.ci", "mbe.lower.ci")
      acc_df[[j]][, mbe_cols] <- mbe$scale_factor * acc_df[[j]][, mbe_cols]

      # ---- compute the overall average RMSD and MBE values.
      rmsd_plot <- .append_rmsd_bartrace(
        data_frame  = acc_df[[j]],
        fig         = if (j == 1L) plotly::plot_ly() else rmsd_plot,
        offsetgroup = (j - 1L),
        name        = row,
        marker      = list(
          color = marker$colors[j],
          pattern = list(
            shape    = marker$patterns$shape[j],
            solidity = marker$patterns$solidity,
            size     = marker$patterns$size
          )
        ),
        showlegend  = show_legend
      )

      mbe_plot <- .append_mbe_bartrace(
        data_frame  = acc_df[[j]],
        fig         = if (j == 1L) plotly::plot_ly() else mbe_plot,
        offsetgroup = (j - 1L),
        marker      = list(
          color   = marker$colors[j],
          pattern = list(
            shape    = marker$patterns$shape[j],
            solidity = marker$patterns$solidity,
            size     = marker$patterns$size
          )
        )
      )
    }

    # --- Compute fixed rmsd and mbe y-axis if requested
    show_ticks <- (i == length(names(rmsd_dfs)))
    if (fixed_axis) {
      rmsd_max <- 0.01 + max(
        sapply(X = rmsd_dfs, FUN = function(x) {
          sapply(X = x, FUN = function(y) {
            max(y$rmsd) + y$rmsd.upper.ci[which.max(y$rmsd)]
          })
        })
      )

      mbe_max <- 0.01 + max(
        sapply(X = rmsd_dfs, FUN = function(x) {
          sapply(X = x, function(y) {
            max(abs(y$mbe)) + y$mbe.upper.ci[which.max(abs(y$mbe))]
          })
        })
      )
      rmsd_range <- c(0.0, rmsd_max)
      mbe_range  <- c(-mbe_max, mbe_max)
      tickmode <- "linear"
      nticks <- NULL
    } else {
      rmsd_range <- mbe_range <- NULL
      rmsd$dtick <- mbe$dtick <- NULL
      tickmode <- "auto"
      nticks   <- 5L
    }


    rmsd_plot <- rmsd_plot %>% plotly::layout(
      title = title,
      xaxis = list(
        title          = "Tool",
        tickangle      = rmsd$tickangle,
        titlefont      = font,
        tickfont       = list(size = ticksize$xaxis),
        showticklabels = show_ticks
      ),
      yaxis = list(
        tickmode = tickmode,
        title     = rmsd$title,
        titlefont = font,
        dtick     = rmsd$dtick,
        tick0     = 0L,
        nticks    = nticks,
        range     = rmsd_range,
        tickfont  = list(size = ticksize$yaxis)
      ),
      legend = list(
        title = list(
          text        = paste0("<b>", legend$title, "</b>"),
          orientation = "h",
          font = list(size = legend$size)
        )
      ),
      shapes = rmsd_lines
    ) %>%
      plotly::add_annotations(
        text      = paste0("<b>", column, "</b>"),
        x         = 0.9,
        y         = 1.0,
        xref      = "paper",
        yref      = "paper",
        xanchor   = "middle",
        yanchor   = "top",
        showarrow = FALSE,
        font      = list(size = 22L)
      )


    mbe_plot <- mbe_plot %>% plotly::layout(
      legend = list(orientation = "h"),
      xaxis = list(
        tickfont       = list(size = ticksize$xaxis),
        tickangle      = mbe$tickangle,
        showticklabels = show_ticks
      ),
      yaxis = list(
        tickmode = tickmode,
        title     = "MBE (x10)", #@TODO: automate according to norm_factor,
        dtick     = mbe$dtick,
        tick0     = 0L,
        nticks    = nticks,
        range     = mbe_range,
        tickfont  = list(size = ticksize$yaxis),
        titlefont = font
      ),
      shapes = mbe_lines
    )

    rmsd_barplots[[i]] <- plotly::subplot(
      rmsd_plot, mbe_plot,
      shareY = FALSE,
      titleY = FALSE,
      widths = c(1.0 - mbe_plot_ratio, mbe_plot_ratio),
      margin = 0.03
    )
  }

  merged_rmsd_barplots <- plotly::subplot(
    rmsd_barplots,
    nrows  = length(rmsd_barplots),
    shareX = FALSE,
    titleY = TRUE
  ) %>% plotly::layout(
    annotations = list(
      plotly_yaxis_subplot(
        text    = rmsd$title,
        x       = -0.02,
        y       = 0.5,
        font    = font,
        xanchor = "right"
      ),
      plotly_yaxis_subplot(
        text    = mbe$title,
        x       = 0.68,
        y       = 0.5,
        font    = font,
        xanchor = "left"
      )
    )
  ) %>%
    plotly::config(
      editable             = FALSE,
      displaylogo          = FALSE,
      scrollZoom           = TRUE,
      toImageButtonOptions = list(
        format   = "svg",
        filename = filename
      )
    )
  merged_rmsd_barplots
}

#' Append plotly bartrace
#'
#' @description
#' Add a bartrace to an existing plotly figure (or create a new plot if needed)
#'
#' @keywords internal
#' @importFrom plotly add_trace plot_ly
#' @param fig an optional plotly figure object
#' @param data_frame a dataframe containing the traced data. Must contain:
#'        a categorical column for the xaxis (specified trough xcol)
#'        a numerical column for the yaxis (specified through ycol)
#'        two columns, containing an upper and lower length for the error bars.
#'        These can be specified with 'error_lower_col' and 'error_upper_col'
#' @param xcol column name of the categorical label to plot on the xaxis
#' @param ycol column name of the numerical values to plot on the yaxis
#' @param error_lower_col name of the column containing lower CI values for y
#' @param error_upper_col name of the column containing upper CI values for y
#' @param error_color color of the displayed error bars
#' @param error_thickness thickness of the displayed error bars
#' @param error_width whiskers width of the displayed error bars
#' @param ... additionnal optional arguments for plotly::add_trace
#' @return a plotly figure, with an added bartrace.
.append_plotly_bartrace <- function(
  fig             = plotly::plot_ly(),
  data_frame      = NULL,
  xcol            = "true_r",
  ycol            = NULL,
  error_lower_col = paste0(ycol, ".lower.ci"),
  error_upper_col = paste0(ycol, ".upper.ci"),
  error_color     = "#111111",
  error_thickness = 1L,
  error_width    = 2L,
  ...
) {
  plotly::add_trace(
    p           = fig,
    x           = data_frame[[xcol]],
    y           = data_frame[[ycol]],
    type        = "bar",
    error_y     = list(
      symmetric  = FALSE,
      array      = data_frame[[error_upper_col]],
      arrayminus = data_frame[[error_lower_col]],
      color      = error_color,
      thickness  = error_thickness,
      width      = error_width
    ),
    ...
  )
}

#' Append Root Mean Square Error bartrace
#'
#' @description
#' Wrapper function for `badger.plots::append_plotly_bartrace` to easily plot
#' Root Mean Square Deviation estimates
#'
#' @keywords internal
#' @importFrom plotly plot_ly
#' @param fig an optional plotly graph object on which the barchart is appended
#' @param data_frame a dataframe containing the traced data. Must contain
#'        columnms true_r, 'rmsd', 'rmsd.lower.ci' 'rmsd.upper.ci'. See
#'        compute_summary_statistics() function to generate such a data.frame
#'        from a kinship results file.
#' @param error_color color of the displayed error bars
#' @param error_thickness thickness of the displayed error bars
#' @param error_width whiskers width of the displayed error bars
#' @param ... additionnal optional arguments for plotly::add_trace
#' @return a plotly figure object , with an added bartrace.
.append_rmsd_bartrace <- function(
  data_frame,
  fig = plotly::plot_ly(),
  error_color     = "#111111",
  error_thickness = 1L,
  error_width    = 1L,
  ...
) {
  .append_plotly_bartrace(
    fig         = fig,
    data_frame  = data_frame,
    error_width = error_width,
    xcol        = "true_r",
    ycol        = "rmsd",
    ...
  )
}

#' Append Mean Bias Error bartrace.
#'
#' @description
#' Wrapper function for `badger.plots::append_plotly_bartrace` to easily plot
#' Mean Bias Error estimate bartraces.
#'
#' @keywords internal
#' @importFrom plotly plot_ly
#' @param fig a plotly graph object.
#' @param data_frame a dataframe. must contain the following columns:
#'                   'true_r' 'mbe', 'mbe.upper.ci', 'mbe.lower.ci'
#' @param xcol name of the column containing the expected relatedness orders
#'             within the provided dataframe
#' @param ycol name of the column containing the mean bias error estimate of
#'             each relatedness order
#' @param showlegend should the legend be displayed ?
#' @param error_color default color for error bars
#' @param error_thickness default thickness for error bars.
#' @param error_width default tail width for error bars.
#' @param ... additional arguments for plotly::add_trace
#' @return a plotly figure object , with an added bartrace.
.append_mbe_bartrace <- function(
  data_frame,
  fig             = plotly::plot_ly(),
  xcol            = "true_r",
  ycol            = "mbe",
  showlegend      = FALSE,
  error_color     = "#111111",
  error_thickness = 1L,
  error_width     = 2L,
  ...
) {
  .append_plotly_bartrace(
    fig              = fig,
    data_frame       = data_frame,
    xcol             = xcol,
    ycol             = ycol,
    showlegend       = showlegend,
    error_color      = error_color,
    error_thickness  = error_thickness,
    error_width      = error_width,
    ...
  )
}

#' Plotly Common Y-axis for subplots.
#'
#' @description
#' Generate a named list to create a common yaxis title for a given plotly
#' figure containing subplots.
#'
#' @keywords internal
#' @param text main label for the title
#' @param x x-coordinate position (in "paper" reference)
#' @param y y-coordinate position (in "paper" reference)
#' @param font named list containing various font bindings
#' @param textangle angle of the main title label (in degrees)
#' @param xanchor title x-position anchor ("left", "middle", "right")
#' @param yanchor title y-position anchor ("top", "center", "bottom")
#' @return (list) a plotly annotation in the form of a named list.
plotly_yaxis_subplot <- function(
  text      = "",
  x         = 0L,
  y         = 0L,
  font      = list(color = "#111111", size = 14L),
  textangle = 270L,
  xanchor   = "middle",
  yanchor   = "center"
) {
  list(
    text      = text,
    x         = x,
    y         = y,
    font      = font,
    textangle = textangle,
    showarrow = FALSE,
    xanchor   = xanchor,
    yanchor   = yanchor,
    xref      = "paper",
    yref      = "paper"
  )
}
