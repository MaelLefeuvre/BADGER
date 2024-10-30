vline <- function(x = 0L, dash = "dash", color = "black") {
  list(
    type = "line",
    y0   = 0L,
    y1   = 1L,
    yref = "paper",
    x0   = x,
    x1   = x,
    line = list(dash = dash, color = color)
  )
}



#' @import plotly
#' @import viridis
plot_density <- function(x, y) {
  df         <- data.frame(x = x, y = y)
  levels     <- sort(unique(df$y))
  histcolors <- viridis::viridis(n = length(levels), alpha = 0.5)
  linecolors <- viridis::viridis(n = length(levels), alpha = 1.0)
  fig        <- plotly::plot_ly()
  maxdens    <- 0L
  i          <- 0L
  for (lvl in levels) {
    i       <- i + 1L
    dens    <- stats::density(df$x[df$y == lvl])
    dens    <- data.frame(x = dens$x, y = dens$y)
    dens$y  <- dens$y
    maxdens <- ifelse(max(dens$y) > maxdens, max(dens$y), maxdens)

    fig <- fig %>% add_histogram(
      x = df$x[df$y == lvl],
      name     = lvl,
      alpha    = 0.5,
      marker   = list(color = histcolors[i]),
      histnorm = "probability density"
    )
    fig <- fig %>% add_lines(
      data      = dens,
      name      = lvl,
      x         = ~x,
      y         = ~y,
      line      = list(width = 3L, color = linecolors[i]),
      fill      = "tozeroy",
      fillcolor = histcolors[i]
    )
  }

  lines <- mapply(vline, x = levels, color = linecolors, SIMPLIFY = FALSE)

  fig %>% layout(barmode = "overlay", shapes = lines)
}
