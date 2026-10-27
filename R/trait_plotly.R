
#' Locus scatter plotly for GWAS Catologue traits
#'
#' Produces a scatter plot using plotly of embedded GWAS Catalogue trait data
#' acquired through the LDlink API via [link_trait()].
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param trait_filter Character vector of traits to filter results.
#' @param ... Optional arguments passed to `plot_ly()`.
#' @returns A `plotly` scatter plot.
#' @seealso [link_trait()] [locus_plotly()]

trait_plotly <- function(loc,
                         trait_filter = NULL,
                         ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  
  xlim <- loc$xrange / 1e6
  xlim <- xlim + diff(xlim) * c(-0.01, 0.01)
  xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  LDT <- loc$LDtrait
  if (!is.null(trait_filter)) {
    LDT <- LDT[LDT$GWAS_Trait %in% trait_filter, ]
  }
  traits <- sort(unique(LDT$GWAS_Trait))
  
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
  ylim <- ylim + diff(ylim) * c(-0.05, 0.05)
  
  traits <- sort(unique(LDT$GWAS_Trait))
  if (length(traits) > 1) {
    trait_group <- similarGroups(traits)
    scheme <- rainbow(length(unique(trait_group)))
    LDT$colour <- scheme[trait_group[LDT$GWAS_Trait]]
  } else {
    LDT$colour <- "red"
  }
  
  # 27/10/26
  # LDtrait GWAS catalog API columns `Effect_size_95_CI` and `Beta_or_OR`
  # are the wrong way round!
  
  hovertext <- paste0("GWAS Catalogue: ", LDT$RS_Number,
                      "<br>Chr ", loc$seqname, ": ", LDT$pos,
                      "<br>P = ", signif(LDT$P_value, 3),
                      "<br>RAF = ", LDT$Risk_Allele,
                      "<br>OR/beta = ", LDT$Effect_size_95_CI,
                      "<br>CI = ", LDT$Beta_or_OR,
                      "<br>Trait: ", LDT$GWAS_Trait)
  
  p <- plot_ly(x = LDT$pos / 1e6, y = LDT$logP,
               marker = list(opacity = 0.8,
                             color = I(LDT$colour),
                             line = list(width = 1, color = "black")),
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


#' @importFrom stats as.dist cutree hclust
similarGroups <- function(x, thresh = 0.7, method = "jw"){
  ox <- x
  # remove stop words
  x <- gsub("-", " ", tolower(x))
  x <- gsub("\\(|\\)|,", "", x)
  stop_words <- c("to", "in", "and", "or", "at", "of", "ratio", "concentration",
                  "ukb", "data", "field", "fraction", "absolute", "percentage", "total",
                  "fraction", "maximum", "minimum", "mean", "inv norm transformed")
  x <- removeWords(stop_words, x)
  x <- gsub('[0-9]+', '', x)
  
  if (!requireNamespace("stringdist")) stop("stringdist package required")
  simmat <- stringdist::stringsimmatrix(x, x, method)
  hc <- hclust(as.dist(1 - simmat), method = "complete")
  grp <- cutree(hc, h = 1 - thresh)
  setNames(grp, ox)
}


removeWords <- function(words, x) {
  p <- paste0("(?:\\s+|^)\\S*(?<!\\w)(?:", paste(words, collapse = "|"), ")(?!\\w)\\S*")
  trimws(gsub(p, "", x, perl = TRUE))
}
