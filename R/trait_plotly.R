
#' Locus scatter plotly for GWAS Catologue traits
#'
#' Produces a scatter plot using plotly of embedded GWAS Catalogue trait data
#' acquired through the LDlink API via [link_trait()].
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param ... Optional arguments passed to `plot_ly()`.
#' @returns A `plotly` scatter plot.
#' @seealso [link_trait()] [locus_plotly()]

trait_plotly <- function(loc,
                         ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  
  xlim <- loc$xrange / 1e6
  xlim <- xlim + diff(xlim) * c(-0.01, 0.01)
  xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  LDT <- loc$LDtrait
  
  if (is.null(LDT) || nrow(LDT) == 0) {
    # blank plotly
    LDT <- data.frame(matrix(nrow = 0, ncol = 2))
    colnames(LDT) <- c("pos", "logP")
    p <- plot_ly(LDT,
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
                     annotations = list(x = 0.5, y = 0.5, text = "No trait data",
                                        xref = "paper", yref = "paper",
                                        showarrow = FALSE)) %>%
      plotly::config(displaylogo = FALSE,
                     modeBarButtonsToRemove = c("select2d", "lasso2d",
                                                "autoScale2d", "resetScale2d",
                                                "hoverClosest", "hoverCompare"),
                     toImageButtonOptions = list(format = "svg"))
    return(p)
  }
  
  ylab <- "GWAS trait -log<sub>10</sub> P"
  ylim <- range(LDT$logP, na.rm = TRUE)
  ylim <- ylim + diff(ylim) * c(-0.01, 0.01)
  
  # 27/10/26
  # LDtrait GWAS catalog API columns `Effect_size_95_CI` and `Beta_or_OR`
  # are the wrong way round!
  
  hovertext <- paste0("GWAS Catologue: ", LDT$RS_Number,
                      "<br>Chr ", loc$seqname, ": ", LDT$pos,
                      "<br>P = ", signif(LDT$P_value, 3),
                      "<br>RAF = ", LDT$Risk_Allele,
                      "<br>OR/beta = ", LDT$Effect_size_95_CI,
                      "<br>CI = ", LDT$Beta_or_OR,
                      "<br>Trait: ", LDT$GWAS_Trait)
  
  p <- plot_ly(x = LDT$pos / 1e6, y = LDT$logP,
               marker = list(size = 9, opacity = 0.8,
                             line = list(width = 0.5, color = "black")),
               text = hovertext, hoverinfo = 'text',
               key = LDT$RS_Number,
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
