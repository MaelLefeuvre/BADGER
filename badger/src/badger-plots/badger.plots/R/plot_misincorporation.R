#' Plot Post-mortem damage misincorporation probabilities.
#'
#' @description
#' Plot misincorporation probabilities from a MapDamage `misincorporation.txt`
#' file.
#'
#' @import plotly
#' @param path Path to a misincorporation.txt file
#' @param n restrict estimation to n bases (default: 70)
#' @param title default plot title
#' @param colors default color scheme for each misincorporation type.
#' @return a plotly figure object.
plot_misincorporation <- function(
  path,
  n = 70L,
  title = paste(basename(path), "misincorporation frequency"),
  colors = list(G.A = "#008083", C.T = "#f78104")
) {
  misincorporation <- read.table(path, header = TRUE, sep = "\t")

  summed_data <- stats::aggregate(
    cbind(C, G, C.T, G.A, Total) ~ Pos + End, misincorporation, sum
  )
  summed_data <- cbind(
    summed_data, ratio = summed_data[, c("C.T", "G.A")] / summed_data$Total
  )

  summed_data$ratio.C.T <- summed_data$C.T / summed_data$C
  summed_data$ratio.G.A <- summed_data$G.A / summed_data$G

  if (!is.null(n)) {
    summed_data <- summed_data[summed_data$Pos <= n, ]
  }

  plot_5p <- plotly::plot_ly(
    data = summed_data[summed_data$End == "5p", ],
    type = "scatter",
    mode = "lines",
    line = list(shape = "splines")
  ) %>%
    plotly::add_trace(
      x = ~Pos, y = ~ratio.G.A, color = I(colors$G.A), name = "G>A"
    ) %>%
    plotly::add_trace(
      x = ~Pos, y = ~ratio.C.T, color = I(colors$C.T), name = "C>T"
    ) %>%
    plotly::layout(xaxis = list(title = ""), yaxis = list(title = ""))

  plot_3p <- plotly::plot_ly(
    data = summed_data[summed_data$End == "3p", ],
    type = "scatter",
    mode = "lines",
    line = list(shape = "splines")
  ) %>%
    plotly::add_trace(
      x = ~Pos, y = ~ratio.G.A, color = I(colors$G.A), showlegend = FALSE
    ) %>%
    plotly::add_trace(
      x = ~Pos, y = ~ratio.C.T, color = I(colors$C.T), showlegend = FALSE
    ) %>%
    plotly::layout(
      xaxis  = list(title = "", autorange = "reversed"),
      yaxis  = list(title = "", side = "right")
    )

  plotly::subplot(plot_5p, plot_3p) %>%
    plotly::layout(
      title = title,
      legend = list(title = list(text = "<b>Type</b>"), orientation = "h"),
      xaxis = list(title = ""),
      yaxis = list(title = ""),
      annotations = list(
        list(
          x         = -0.05,
          y         = 0.5,
          xref      = "paper",
          yref      = "paper",
          text      = "Frequency",
          yanchor   = "center",
          xanchor   = "right",
          showarrow = FALSE,
          font      = list(size = 16L),
          textangle = 270L
        ),
        list(
          x         = 0.5,
          xref      = "paper",
          y         = -0.05,
          yref      = "paper",
          text      = "Position (bp)",
          yanchor   = "top",
          xanchor   = "center",
          showarrow = FALSE,
          font      = list(size = 16L)
        )
      )
    ) %>%
    plotly::config(
      editable = FALSE,
      toImageButtonOptions = list(
        format   = "svg",
        filename = paste0(basename(path), "-misincorporation-plot")
      )
    )
}



#' Plot Fragment size Frequency plot.
#'
#' @description
#' Plot a read size frequency distribution from a gargammel size-frequency file.
#' @import plotly
#' @param path Path to a mapDamage-v2 or Gargammel size frequency distribution
#'        file.
#' @param color Main color of the plotted line
#' @param title Main title of the plot
#' @param showmax a logical indicating whether the maximum frequency should be
#'        highlighted on the x-axis, through a dashed line annotation
#' @param smoothing span of the kernel used during smoothing of values (loess).
#'        higher values will increase smoothing of the line, but decrease the
#'        definition.
#' @param sep field separator used within the file provided by `path`
#' @return a plotly figure object.
plot_size_frequency <- function(
  path,
  color     = "#008080",
  title     = basename(path),
  showmax   = TRUE,
  smoothing = 0.1,
  sep       = "\t"
) {
  sizefreqs <- read.table(
    file      = path,
    header    = FALSE,
    sep       = sep,
    col.names = c("length", "frequency")
  )

  smoothed <- stats::loess(frequency ~ length,
    data = sizefreqs,
    span = smoothing
  )

  sizefreqs$length    <- as.numeric(smoothed$x)
  sizefreqs$frequency <- smoothed$fitted
  sizefreqs           <- rbind(
    data.frame(length = min(sizefreqs$length) - 1.0, frequency = 0.0),
    sizefreqs
  )

  if (showmax) {
    max           <- sizefreqs[which.max(sizefreqs$frequency), ]
    minlength     <- min(sizefreqs$length)
    max$frequency <- max$frequency * 100.0
    maxlines      <- list(
      list(
        type = "line",
        y0   = 0L,
        y1   = max$frequency,
        x0   = max$length,
        x1   = max$length,
        line = list(color = "#111111", dash = "dot")
      ),
      list(
        type = "line",
        y0   = max$frequency,
        y1   = max$frequency,
        x0   = minlength,
        x1   = max$length,
        line = list(color = "#111111", dash = "dot")
      )
    )

    maxticks <- list(
      list(
        x         = max$length,
        y         = 0L,
        text      = paste("<b>", as.character(max$length), "</b>"),
        showarrow = FALSE,
        xanchor   = "right",
        yanchor   = "top",
        textangle = -45L
      )
    )
  }
  plotly::plot_ly(
    data = sizefreqs,
    x    = ~length,
    y    = ~frequency * 100.0,
    type = "scatter",
    mode = "lines",
    line = list(shape = "spline"),
    color = I(color),
    connectgaps = TRUE
  ) %>%
    plotly::layout(
      title       = title,
      xaxis       = list(title = "Length (bp)"),
      yaxis       = list(title = "Frequency (%)"),
      shapes      = maxlines,
      annotations = maxticks
    ) %>%
    plotly::config(
      editable = TRUE,
      toImageButtonOptions = list(
        format   = "svg",
        filename = paste0(basename(path), "-scatterplot")
      )
    )
}

#setwd("/data/mlefeuvre/dev/aDNA-kinship-simulations/resources/gargammel/")
#chan_size_freq <- plot_size_frequency(
#  "sizefreqs/Chan_Meso-sizefreq.txt",
#  title = "Chan Meso size frequency distribution"
#)
#ust_size_freq <- plot_size_frequency(
#  "sizefreqs/Ust_Ishim_sizefreq.size.gz",
#  title = "Chan Meso size frequency distribution", smoothing = 0.1
#)
#misincorporation_fig <- plot_misincorporation(
#  "misincorporations/Chan_meso/misincorporation.txt",
#  n = 25,
#  title = "Chan-Meso Misincorporation Frequency"
#)

#htmlwidgets::saveWidget(
#  chan_size_freq, file = "Chan-meso-sizefreq.html", selfcontained = TRUE
#)
#htmlwidgets::saveWidget(
#  ust_size_freq, file = "Ust-ishim-sizefreq.html", selfcontained = TRUE
#)
#htmlwidgets::saveWidget(
#  misincorporation_fig, file = "Chan-meso-misincorporation.html",
#  selfcontained = TRUE
#)