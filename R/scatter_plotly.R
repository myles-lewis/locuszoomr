
#' Locus scatter plotly
#'
#' Produces a scatter plot from a 'locus' class object using plotly.
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param index_snp Specifies index SNP to be shown in a different colour and
#'   symbol. Defaults to the SNP with the lowest p-value. Set to `NULL` to not
#'   show this.
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param chromCol Colour for normal points if `LD` is `FALSE` when the locus
#'   object is made.
#' @param sigCol Colour for significant points if `LD` is `FALSE`.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param yzero Logical whether to force y axis limit to include y=0.
#' @param showLD Logical whether to show LD with colours
#' @param LD_scheme Vector of colours for plotting LD. The first colour is for
#'   SNPs which lack LD information. The next 5 colours are for r^2 or D' LD
#'   results ranging from 0 to 1 in intervals of 0.2. The final colour is for
#'   the index SNP.
#' @param marker_outline Specifies colour for outlining points.
#' @param marker_size Value for size of markers in plotly units.
#' @param recomb_col Colour for recombination rate line if recombination rate
#'   data is present. Set to `NA` to hide the line. See [link_recomb()] to add
#'   recombination rate data.
#' @return A `plotly` scatter plot.
#' @seealso [locus()] [locus_plotly()]
#' @importFrom plotly add_trace
#' @export
#' 
scatter_plotly <- function(loc,
                           index_snp = loc$index_snp,
                           pcutoff = 5e-08,
                           chromCol = 'royalblue',
                           sigCol = 'red',
                           xlab = NULL,
                           ylab = NULL,
                           yzero = (loc$yvar == "logP"),
                           showLD = TRUE,
                           LD_scheme = c('grey', 'royalblue', 'cyan2', 'green3', 
                                         'orange', 'red', 'purple'),
                           marker_outline = "black",
                           marker_size = 7,
                           recomb_col = "blue") {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  if (is.null(loc$data)) stop("No data points, only gene tracks")
  data <- loc$data
  if (is.null(xlab)) xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  if (is.null(ylab)) {
    ylab <- if (loc$yvar == "logP") "-log<sub>10</sub> P" else loc$yvar
  }
  hasLD <- "ld" %in% colnames(data)
  leg <- list()
  if (!"bg" %in% colnames(data)) {
    if (showLD & hasLD) {
      data$bg <- cut(data$ld, -1:6/5, labels = FALSE)
      data$bg[is.na(data$bg)] <- 1L
      data$bg[data[, loc$labs] == index_snp] <- 7L
      data <- data[order(data$bg), ]
      LD_scheme <- rep_len(LD_scheme, 7 - is.null(index_snp))
      data$bg <- factor(data$bg, levels = 1:7,
                        labels = c("unknown", "0.0 - 0.2", "0.2 - 0.4",
                                   "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1.0",
                                   "index"))
      leg <- list(title = list(text = "Linkage r<sup>2</sup>"))
    } else {
      showLD <- FALSE
      data$bg <- chromCol
      if (loc$yvar == "logP") data$bg[data[, loc$p] < pcutoff] <- sigCol
      data$bg[data[, loc$labs] == index_snp] <- "purple"
      LD_scheme <- c(chromCol, sigCol, "purple")
      data$bg <- factor(data$bg, levels = LD_scheme,
                        labels = c("ns", paste("P <", pcutoff), "index"))
    }
  }
  
  # scatter plotly
  recomb <- !is.null(loc$recomb) & !is.na(recomb_col)
  
  pch <- rep(21L, nrow(data))
  pch[data[, loc$labs] == index_snp] <- 23L
  if ("pch" %in% colnames(data)) pch <- data$pch
  col <- "black"
  if ("col" %in% colnames(data)) col <- data$col
  
  xlim <- loc$xrange / 1e6
  xext <- diff(xlim) * 0.05
  xlim <- xlim + c(-xext, xext)
  
  ylim <- range(data[, loc$yvar], na.rm = TRUE)
  if (yzero) ylim[1] <- min(c(0, ylim[1]))
  ydiff <- diff(ylim)
  ylim[2] <- ylim[2] + ydiff * 0.05
  ylim[1] <- if (ylim[1] != 0) ylim[1] - ydiff *0.05 else ylim[1] - ydiff *0.02
  
  hovertext <- paste0(data[, loc$labs], "<br>Chr ",
                      data[, loc$chrom], ": ", data[, loc$pos],
                      "<br>P = ", signif(data[, loc$p], 3))
  ylim2 <- c(-2, 102)
  
  if (!recomb) {
    # standard plotly
    plot_ly(x = data[, loc$pos] / 1e6, y = data[, loc$yvar],
            color = data$bg, colors = LD_scheme,
            marker = list(size = marker_size, opacity = 0.8,
                          line = list(width = 1, color = marker_outline)),
            text = hovertext,
            hoverinfo = 'text',
            type = "scattergl", mode = "markers") %>%
      plotly::layout(xaxis = list(title = xlab,
                                  ticks = "outside",
                                  range = as.list(xlim)),
                     yaxis = list(title = ylab,
                                  ticks = "outside",
                                  showline = TRUE, range = ylim),
                     legend = leg,
                     showlegend = showLD | !is.null(pcutoff)) %>%
      plotly::config(displaylogo = FALSE)
    # modeBarButtonsToRemove = c("zoom2d", "pan2d", "select2d",
    #                            "lasso2d", "zoomIn2d", "zoomOut2d", 
    #                            "autoScale2d", "toggleHover"))
  } else {
    # double y axis with recombination
    plot_ly() %>%
      add_trace(x = loc$recomb$start / 1e6, y = loc$recomb$value,
                hoverinfo = "none", colors = LD_scheme,  # colors must go here
                name = "recombination", yaxis = "y2",
                line = list(color = recomb_col),
                mode = "lines", type = "scattergl", showlegend = FALSE) %>%
      add_trace(x = data[, loc$pos] / 1e6, y = data[, loc$yvar],
                color = data$bg,
                marker = list(size = marker_size, opacity = 0.8,
                              line = list(width = 1, color = marker_outline)),
                text = hovertext,
                hoverinfo = 'text',
                type = "scattergl", mode = "markers") %>%
      plotly::layout(xaxis = list(title = xlab,
                                  ticks = "outside",
                                  range = as.list(xlim)),
                     yaxis = list(title = ylab,
                                  ticks = "outside",
                                  showline = TRUE, range = ylim),
                     yaxis2 = list(overlaying = "y", side = "right",
                                   title = "Recombination rate (%)",
                                   ticks = "outside", showgrid = FALSE,
                                   showline = TRUE,
                                   zeroline = FALSE, range = ylim2),
                     legend = c(leg, x = 1.1, y = 1),
                     showlegend = showLD | !is.null(pcutoff)) %>%
      plotly::config(displaylogo = FALSE)
  }
}
