
#' Locus eQTL scatter plotly
#'
#' Produces a scatter plot using plotly of embedded eQTL data acquired through
#' the LDlink API via [link_eqtl()].
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param gene_filter Character vector of genes to filter eQTL results.
#' @param tissue_filter Character vector of tissues to filter eQTL results.
#' @param ... Optional arguments passed to `plot_ly()`.
#' @returns A `plotly` scatter plot.
#' @seealso [link_eqtl()] [locus_plotly()]

eqtl_plotly <- function(loc,
                        gene_filter = NULL,
                        tissue_filter = NULL, ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  
  xlim <- loc$xrange / 1e6
  xlim <- xlim + diff(xlim) * c(-0.01, 0.01)
  xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  
  LDX <- loc$LDexp
  if (!is.null(gene_filter)) {
    LDX <- LDX[LDX$Gene_Symbol %in% gene_filter, ]
  }
  if (!is.null(tissue_filter)) {
    LDX <- LDX[LDX$Tissue %in% tissue_filter, ]
  }
  if (is.null(LDX) || nrow(LDX) == 0) {
    # blank plotly
    LDX <- data.frame(matrix(nrow = 0, ncol = 2))
    colnames(LDX) <- c("pos", "logP")
    p <- plot_ly(LDX,
                 x = ~pos, y = ~logP,
                 showlegend = FALSE,
                 source = "plotly_locus",
                 type = "scattergl", mode = "markers", ...) %>%
      plotly::layout(xaxis = list(title = xlab,
                                  ticks = "outside",
                                  zeroline = FALSE, showgrid = FALSE,
                                  range = as.list(xlim)),
                     yaxis = list(title = "",
                                  showticklabels = FALSE,
                                  zeroline = FALSE, showgrid = FALSE),
                     annotations = list(x = 0.5, y = 0.5, text = "No eQTL data",
                                        xref = "paper", yref = "paper",
                                        showarrow = FALSE)) %>%
      plotly::config(displaylogo = FALSE,
                     modeBarButtonsToRemove = c("select2d", "lasso2d",
                                                "autoScale2d", "resetScale2d",
                                                "hoverClosest", "hoverCompare"),
                     toImageButtonOptions = list(format = "svg"))
    return(p)
  }
  
  ylab <- "eQTL -log<sub>10</sub> P"
  ylim <- range(LDX$logP, na.rm = TRUE)
  ylim <- ylim + diff(ylim) * c(-0.01, 0.01)
  LDX$sign <- factor(sign(LDX$Effect_Size), levels = c(1, -1),
                     labels = c("up", "down"))
  symbols <- c(24L, 25L)
  ngene <- length(unique(LDX$Gene_Symbol))
  scheme <- rainbow(ngene)
  
  hovertext <- paste0("eQTL: ", LDX$RS_ID, "<br>Chr ",
                      loc$seqname, ": ", LDX$pos,
                      "<br>P = ", signif(LDX$P_value, 3),
                      "<br>Gene: ", LDX$Gene_Symbol,
                      "<br>Tissue: ", LDX$Tissue)
  
  p <- plot_ly(x = LDX$pos / 1e6, y = LDX$logP,
               color = LDX$Gene_Symbol, colors = scheme,
               symbol = LDX$sign, symbols = symbols,
               marker = list(size = 9, opacity = 0.8,
                             line = list(width = 0.5, color = "black")),
               text = hovertext, hoverinfo = 'text',
               key = LDX$RS_ID,
               showlegend = FALSE,
               source = "plotly_locus",
               type = "scattergl", mode = "markers", ...) %>%
    plotly::layout(xaxis = list(title = xlab,
                                ticks = "outside",
                                zeroline = FALSE, showgrid = FALSE,
                                range = as.list(xlim)),
                   yaxis = list(title = ylab,
                                ticks = "outside",
                                fixedrange = TRUE,
                                showline = TRUE,
                                range = ylim),
                   dragmode = "pan") %>%
    plotly::config(displaylogo = FALSE,
                   modeBarButtonsToRemove = c("select2d", "lasso2d",
                                              "autoScale2d", "resetScale2d",
                                              "hoverClosest", "hoverCompare"),
                   toImageButtonOptions = list(format = "svg"))
  p
}
