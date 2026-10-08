
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
#' @param eqtl_scheme Vector of colours for eQTL genes, which can be named.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param marker_outline Specifies colour for outlining points.
#' @param marker_size Value for size of markers in plotly units.
#' @param recomb_col Colour for recombination rate line if recombination rate
#'   data is present. Set to `NA` to hide the line. See [link_recomb()] to add
#'   recombination rate data.
#' @param showlegend Logical whether to show a legend for the scatter points.
#' @param show_annot Logical whether to show an annotation of how many eQTL SNPs
#'   were retrieved from LDlink and how many are shown.
#' @param height Height in pixels (optional, defaults to automatic sizing).
#' @param webGL Logical whether to use webGL or SVG for scatter plot.
#' @returns A `plotly` scatter plot.
#' @seealso [link_eqtl()] [locus_plotly()]

overlay_plotly <- function(loc,
                           gene_filter = NULL,
                           tissue_filter = NULL,
                           pcutoff = 5e-08,
                           eqtl_scheme = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           marker_outline = "grey",
                           marker_size = 7,
                           recomb_col = "blue",
                           showlegend = TRUE,
                           show_annot = TRUE,
                           height = NULL,
                           webGL = TRUE) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  
  data <- loc$data
  if (is.null(data)) {
    return(scatter_plotly(loc, height = height))  # blank plot
  }
  
  LDX <- loc$LDexp
  nsnp <- length(unique(LDX$RS_ID))
  if (!is.null(gene_filter)) {
    LDX <- LDX[LDX$Gene_Symbol %in% gene_filter, ]
  }
  if (!is.null(tissue_filter)) {
    LDX <- LDX[LDX$Tissue %in% tissue_filter, ]
  }
  nsnp[2] <- length(unique(LDX$RS_ID))  # filtered
  LDX_snps <- intersect(LDX$RS_ID, data[, loc$labs])
  nsnp[3] <- length(LDX_snps)
  nsnp[4] <- length(unique(LDX$RS_ID[LDX$pos < loc$xrange[1] |
                                       LDX$pos > loc$xrange[2]]))
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
  scheme <- "grey"
  symbols <- c(21L, 24L, 25L)
  sizes <- NULL
  inData <- FALSE
  ns <- TRUE
  
  leg <- list(traceorder = "reversed")
  hovertext <- paste0(data[, loc$labs], "<br>Chr ",
                      data[, loc$chrom], ": ", data[, loc$pos],
                      "<br>P = ", signif(data[, loc$p], 3))
  annot <- geneset <- NULL
  if (!noEqtl) {
    gtab <- tapply(LDX$Gene_Symbol, LDX$RS_ID, function(x) length(unique(x)))
    tisstab <- tapply(LDX$Tissue, LDX$RS_ID, function(x) length(unique(x)))
    LDX <- min_p_by_col(LDX, "RS_ID")
    LDX <- LDX[match(LDX_snps, LDX$RS_ID), ]
    LDX$ngene <- gtab[LDX$RS_ID] -1
    LDX$ntissue <- tisstab[LDX$RS_ID] -1
    inData <- match(LDX_snps, data[, loc$labs])
    ns <- -inData
    if (!show_annot) {
      message(nsnp[3], " / ", nsnp[1], " eQTL SNPs shown")
    }
    # colours, shapes
    LDX$sign <- sign(LDX$Effect_Size)
    LDX_by_gene <- min_p_by_col(LDX, "Gene_Symbol")
    geneset <- LDX_by_gene$Gene_Symbol
    ngene <- nrow(LDX_by_gene)
    ldx_scheme <- if (is.null(eqtl_scheme)) {
      rainbow(ngene)
    } else {
      if (is.null(names(eqtl_scheme))) eqtl_scheme else eqtl_scheme[geneset]
    }
    data$bg[inData] <- LDX$Gene_Symbol
    data$bg <- factor(data$bg, levels = c("ns", geneset),
                      labels =  c("ns", geneset))
    scheme <- c("grey", ldx_scheme)
    names(scheme) <- NULL
    
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
    w <- LDX$ngene > 0
    LDX_hovertext[w] <- paste0(LDX_hovertext[w], "<br>+ ",
                               plural(LDX$ngene[w], "gene(s)"))
    w <- LDX$ntissue > 0
    LDX_hovertext[w] <- paste0(LDX_hovertext[w], "<br>+ ",
                               plural(LDX$ntissue[w], "tissue(s)"))
    hovertext[inData] <- paste0(hovertext[inData], LDX_hovertext)
    # annotation
    if (show_annot) {
      msg <- paste0(nsnp[1], " eQTL SNPs<br>",
                    (if (nsnp[1] != nsnp[2]) paste0(nsnp[2], " filtered<br>")),
                    nsnp[3], " shown",
                    (if (nsnp[4] > 0) paste0("<br>", nsnp[4], " outside window")))
      annot <- list(x = 0.01, y = 1,
                    text = msg, font = list(size = 10),
                    bgcolor = "rgba(255, 255, 255, 0.9)",
                    xref = "paper", yref = "paper", align = "left",
                    yanchor = "top", showarrow = FALSE)
    }
  }
  
  hline <- list(type = "line",
                line = list(width = 1, color = '#999999', dash = 'dash'),
                x0 = 0, x1 = 1, y0 = -log10(pcutoff), y1 = -log10(pcutoff),
                xref = "paper", layer = "below")
  
  if (!recomb) {
    # beta shapes
    p <- plot_ly(x = data[ns, loc$pos] / 1e6, y = data[ns, loc$yvar],
                 color = data$bg[ns], colors = scheme,
                 symbol = data$symbol[ns], symbols = symbols,
                 marker = list(opacity = 0.5, size = 6.5,
                               line = list(width = 1, color = marker_outline)),
                 text = hovertext[ns], hoverinfo = 'text',
                 key = data[ns, loc$labs],
                 showlegend = showlegend,
                 source = "plotly_locus", height = height,
                 type = type, mode = "markers") %>%
      add_trace(x = data[inData, loc$pos] / 1e6, y = data[inData, loc$yvar],
                color = data$bg[inData], colors = scheme,
                symbol = data$symbol[inData], symbols = symbols,
                marker = list(opacity = 0.8, size = 11,
                              line = list(width = 1, color = marker_outline)),
                text = hovertext[inData], hoverinfo = 'text',
                key = data[inData, loc$labs],
                showlegend = showlegend,
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
                     annotations = annot,
                     shapes = hline, legend = leg, dragmode = "pan")
  } else {
    # double y axis with recombination
    ylim2 <- c(-2, 102)
    
    # beta shapes
    p <- plot_ly(source = "plotly_locus", height = height) %>%
      add_trace(x = data[ns, loc$pos] / 1e6, y = data[ns, loc$yvar],
                color = data$bg[ns], colors = scheme,
                symbol = data$symbol[ns], symbols = symbols,
                marker = list(opacity = 0.5, size = 6.5,
                              line = list(width = 1, color = marker_outline)),
                text = hovertext[ns], hoverinfo = 'text',
                key = data[ns, loc$labs],
                showlegend = showlegend,
                type = type, mode = "markers") %>%
      add_trace(x = data[inData, loc$pos] / 1e6, y = data[inData, loc$yvar],
                color = data$bg[inData], colors = scheme,
                symbol = data$symbol[inData], symbols = symbols,
                marker = list(opacity = 0.8, size = 11,
                              line = list(width = 1, color = marker_outline)),
                text = hovertext[inData], hoverinfo = 'text',
                key = data[inData, loc$labs],
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
                     annotations = annot,
                     shapes = hline,
                     legend = c(leg, x = 1.1, y = 1), showlegend = TRUE)
  }
  
  attr(p, "eqtl_genes") <- geneset
  
  p %>%
    plotly::config(displaylogo = FALSE,
                   modeBarButtonsToRemove = c("select2d", "lasso2d",
                                              "autoScale2d", "resetScale2d",
                                              "hoverClosest", "hoverCompare"),
                   toImageButtonOptions = list(format = "svg"))
}


plural <- function(n, text) {
  is1 <- n == 1
  textout <- rep_len(text, length(n))
  textout[is1] <- gsub("\\(s\\)", "", text)
  textout[!is1] <- gsub("\\(|\\)", "", text)
  paste(n, textout)
}
