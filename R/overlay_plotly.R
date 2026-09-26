
#' Plot overlaying eQTL and GWAS data using plotly
#'
#' Produces a scatter plot using plotly of embedded eQTL data acquired through
#' the LDlink API via [link_eqtl()] overlaid on GWAS data. As each SNP may
#' have eQTLs with multiple genes in multiple tissues, the method used is to
#' select the gene/tissue eQTL with the lowest p-value. SNPs are matched by
#' rsID.
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param gene_filter Character vector of genes to filter eQTL results.
#' @param tissue_filter Character vector of tissues to filter eQTL results.
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param scheme Vector of 3 colours if LD is not shown: 1st = normal points,
#'   2nd = colour for significant points, 3rd = index SNP(s).
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param marker_outline Specifies colour for outlining points.
#' @param marker_size Value for size of markers in plotly units.
#' @param recomb_col Colour for recombination rate line if recombination rate
#'   data is present. Set to `NA` to hide the line. See [link_recomb()] to add
#'   recombination rate data.
#' @param showlegend Logical whether to show a legend for the scatter points.
#' @param height Height in pixels (optional, defaults to automatic sizing).
#' @param webGL Logical whether to use webGL or SVG for scatter plot.
#' @returns A `plotly` scatter plot.
#' @seealso [link_eqtl()] [locus_plotly()]

overlay_plotly <- function(loc,
                           gene_filter = NULL,
                           tissue_filter = NULL,
                           pcutoff = 5e-08,
                           scheme = c('grey', 'dodgerblue', 'red'),
                           xlab = NULL,
                           ylab = NULL,
                           marker_outline = "grey",
                           marker_size = 7,
                           recomb_col = "blue",
                           showlegend = TRUE,
                           height = NULL,
                           webGL = TRUE) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  
  data <- loc$data
  if (is.null(data)) {
    return(scatter_plotly(loc, height = height))  # blank plot
  }
  
  LDX <- loc$LDexp
  if (!is.null(gene_filter)) {
    LDX <- LDX[LDX$Gene_Symbol %in% gene_filter, ]
  }
  if (!is.null(tissue_filter)) {
    LDX <- LDX[LDX$Tissue %in% tissue_filter, ]
  }
  LDX_snps <- intersect(LDX$RS_ID, data[, loc$labs])
  noEqtl <- is.null(LDX) || nrow(LDX) == 0 || length(LDX_snps) == 0
  
  xlim <- loc$xrange / 1e6
  xlim <- xlim + diff(xlim) * c(-0.01, 0.01)
  if (is.null(xlab)) xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  if (is.null(ylab)) ylab <- "-log<sub>10</sub> P"
  type <- if (webGL) "scattergl" else "scatter"
  
  recomb <- !is.null(loc$recomb) & !is.na(recomb_col)
  
  ylim <- range(data$logP, na.rm = TRUE)
  ylim[1] <- min(c(0, ylim[1]))  # yzero = TRUE
  ydiff <- diff(ylim)
  ylim[2] <- ylim[2] + ydiff * 0.05
  ylim[1] <- if (ylim[1] != 0) ylim[1] - ydiff *0.05 else ylim[1] - ydiff *0.02
  
  data$bg <- "ns"
  scheme <- scheme[1]
  symbols <- c(21L, 24L, 25L)
  sizes <- NULL
  
  leg <- list(traceorder = "reversed")
  hovertext <- paste0(data[, loc$labs], "<br>Chr ",
                      data[, loc$chrom], ": ", data[, loc$pos],
                      "<br>P = ", signif(data[, loc$p], 3))
  if (!noEqtl) {
    LDX <- min_p_by_col(LDX, "RS_ID")
    LDX <- LDX[match(LDX_snps, LDX$RS_ID), ]
    inData <- match(LDX_snps, data[, loc$labs])
    message(length(inData), " (",
            format(length(inData) / nrow(data) * 100, digits = 3), "%) eQTL SNPs")
    # colours, shapes
    LDX$sign <- sign(LDX$Effect_Size)
    LDX_by_gene <- min_p_by_col(LDX, "Gene_Symbol")
    geneset <- LDX_by_gene$Gene_Symbol
    ngene <- nrow(LDX_by_gene)
    ldx_scheme <- rainbow(ngene)
    data$bg[inData] <- LDX$Gene_Symbol
    data$bg <- factor(data$bg, levels = c("ns", geneset),
                      labels =  c("ns", geneset))
    scheme <- c(scheme[1], ldx_scheme)
    
    # beta symbols
    symbol <- rep_len("ns", nrow(data))
    symbol[inData] <- LDX$sign
    data$symbol <- factor(symbol, levels = c("ns", "1", "-1"),
                          labels = c(" ", "up", "down"))
    data$size <- 1L
    data$size[inData] <- 2L
    sizes <- c(40, 100)
    if (!webGL) sizes <- sizes/2
    
    LDX_hovertext <- paste0("<br>eQTL P = ", signif(LDX$P_value, 3),
                            "<br>eQTL beta = ", signif(LDX$Effect_Size, 3),
                            "<br>Gene: ", LDX$Gene_Symbol,
                            "<br>Tissue: ", LDX$Tissue)
    hovertext[inData] <- paste0(hovertext[inData], LDX_hovertext)
  }
  
  hline <- list(type = "line",
                line = list(width = 1, color = '#999999', dash = 'dash'),
                x0 = 0, x1 = 1, y0 = -log10(pcutoff), y1 = -log10(pcutoff),
                xref = "paper", layer = "below")
  
  if (!recomb) {
    # beta shapes
    p <- plot_ly(x = data[, loc$pos] / 1e6, y = data[, loc$yvar],
                 color = data$bg, colors = scheme,
                 symbol = data$symbol, symbols = symbols,
                 size = data$size, sizes = sizes,
                 marker = list(opacity = 0.8,
                               line = list(width = 1, color = marker_outline)),
                 text = hovertext, hoverinfo = 'text',
                 key = data[, loc$labs],
                 showlegend = showlegend,
                 source = "plotly_locus", height = height,
                 type = type, mode = "markers") %>%
      plotly::layout(xaxis = list(title = xlab,
                                  ticks = "outside",
                                  zeroline = FALSE, showgrid = FALSE,
                                  range = as.list(xlim)),
                     yaxis = list(title = ylab,
                                  ticks = "outside",
                                  fixedrange = TRUE,
                                  showline = TRUE,
                                  range = ylim),
                     shapes = hline, legend = leg, dragmode = "pan")
  } else {
    # double y axis with recombination
    ylim2 <- c(-2, 102)
    
    # beta shapes
    p <- plot_ly(source = "plotly_locus", height = height) %>%
      # scatter plot
      add_trace(x = data[, loc$pos] / 1e6, y = data[, loc$yvar],
                color = data$bg,
                symbol = data$symbol,
                size = data$size,
                colors = scheme,  # colors, symbols, sizes must go here
                symbols = symbols, sizes = sizes,
                marker = list(opacity = 0.8,
                              line = list(width = 1, color = marker_outline)),
                text = hovertext, hoverinfo = 'text', key = data[, loc$labs],
                showlegend = showlegend,
                type = type, mode = "markers") %>%
      # recombination line
      add_trace(x = loc$recomb$start / 1e6, y = loc$recomb$value,
                hoverinfo = "none",
                name = "recombination", yaxis = "y2",
                line = list(color = recomb_col, width = 1.5),
                mode = "lines", type = type, showlegend = FALSE) %>%
      plotly::layout(xaxis = list(title = xlab,
                                  ticks = "outside",
                                  zeroline = FALSE,
                                  range = as.list(xlim)),
                     yaxis = list(title = ylab,
                                  ticks = "outside", showgrid = FALSE,
                                  showline = TRUE, fixedrange = TRUE,
                                  range = ylim),
                     yaxis2 = list(overlaying = "y", side = "right",
                                   title = "Recombination rate (%)",
                                   ticks = "outside", showgrid = FALSE,
                                   showline = TRUE, fixedrange = TRUE,
                                   zeroline = FALSE, range = ylim2),
                     shapes = hline,
                     legend = c(leg, x = 1.1, y = 1), showlegend = TRUE)
  }
  
  p %>%
    plotly::config(displaylogo = FALSE,
                   modeBarButtonsToRemove = c("select2d", "lasso2d",
                                              "autoScale2d", "resetScale2d",
                                              "hoverClosest", "hoverCompare"),
                   toImageButtonOptions = list(format = "svg"))
}
