
#' Locus scatter plotly for GWAS Catologue traits
#'
#' Produces a scatter plot using plotly of embedded GWAS Catalogue trait data
#' acquired through the LDlink API via [link_trait()].
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param trait_filter Character vector of traits to filter results.
#' @param category_filter Character vector of GWAS Catalog trait categories to
#'   filter results. See [link_trait()].
#' @param collapse Logical whether to collapse duplicate associations of the same
#'   SNP and trait, keeping the smallest P value.
#' @param r2_min Optional minimum r^2 with the index SNP: associations below it
#'   or without r^2 are shown in grey.
#' @param ld_scheme Optional vector of 7 colours to colour associations by r^2
#'   with the index SNP instead of by trait, as `LD_scheme` in
#'   [scatter_plotly()]; markers then have a fixed size rather than one
#'   scaled by r^2.
#' @param show_ld Whether to use LD for colour associations with the index SNP.
#' @param size_by_r2 Logical whether to scale marker size by r^2 with the index
#'   SNP.
#' @param ... Optional arguments passed to `plot_ly()`.
#' @returns A `plotly` scatter plot.
#' @seealso [link_trait()] [locus_plotly()]
#' @export

trait_plotly <- function(loc,
                         trait_filter = NULL,
                         category_filter = NULL,
                         collapse = TRUE,
                         r2_min = NULL,
                         ld_scheme = NULL,
                         show_ld = "ld" %in% colnames(loc$data),
                         size_by_r2 = TRUE,
                         ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  
  xlim <- loc$xrange / 1e6
  xlim <- xlim + diff(xlim) * c(-0.01, 0.01)
  xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  LDT <- loc$LDtrait
  if (!is.null(trait_filter)) {
    LDT <- LDT[LDT$GWAS_Trait %in% trait_filter, ]
  }
  tc <- loc$trait_category
  if (!is.null(tc)) {
    LDT$category <- tc$category[match(LDT$GWAS_Trait, tc$trait)]
    if (!is.null(category_filter)) {
      LDT <- LDT[LDT$category %in% category_filter, ]
    }
  }
  if (collapse && !is.null(LDT)) LDT <- collapse_traits(LDT)
  traits <- sort(unique(LDT$GWAS_Trait))
  hasLD <- "R2" %in% colnames(LDT)
  r2_cut <- hasLD && !is.null(r2_min) && r2_min > 0
  
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
  
  if (hasLD && show_ld) {
    if (is.null(ld_scheme)) ld_scheme <- eval(formals(scatter_plotly)$LD_scheme)
    # same bins as scatter_plotly(); strongest LD drawn on top
    bin <- cut(LDT$R2, -1:6/5, labels = FALSE)
    bin[LDT$R2 %in% 0] <- 2L
    bin[is.na(bin)] <- 1L
    bin[LDT$RS_Number %in% loc$index_snp] <- 7L
    hex <- rgb(t(col2rgb(rep_len(ld_scheme, 7))), maxColorValue = 255)
    LDT$colour <- hex[bin]
    LDT <- LDT[order(bin), ]
  } else if (!is.null(tc) && length(unique(LDT$category)) > 1) {
    # colour by trait categories
    scheme <- setNames(rainbow(16, s = 0.9), gcat_categories)
    LDT$colour <- scheme[LDT$category]
    LDT$colour[is.na(LDT$colour)] <- "#CCCCCC"
  } else if (length(traits) > 1) {
    # cluster based on text similarity
    trait_group <- similarGroups(traits)
    scheme <- rainbow(length(unique(trait_group)))
    LDT$colour <- scheme[trait_group[LDT$GWAS_Trait]]
  } else {
    LDT$colour <- "red"
  }
  
  outline <- "black"
  if (r2_cut) {
    grey <- which(is.na(LDT$R2) | LDT$R2 < r2_min)
    if (length(grey) > 0) {
      LDT$colour[grey] <- "#CCCCCC"
      outline <- rep_len("black", nrow(LDT))
      outline[grey] <- "#999999"
    }
  }
  size <- 7
  if (hasLD && size_by_r2) {
    size <- rep_len(6, nrow(LDT))
    ok <- !is.na(LDT$R2)
    size[ok] <- 6 + 8 * LDT$R2[ok]
  }
  
  # 27/10/26
  # LDtrait GWAS catalog API columns `Effect_size_95_CI` and `Beta_or_OR`
  # are the wrong way round!
  
  hovertext <- paste0("GWAS Catalog: ", LDT$RS_Number,
                      "<br>Chr ", loc$seqname, ": ", LDT$pos,
                      "<br>P = ", signif(LDT$P_value, 3),
                      (if (hasLD) {
                        paste0("<br>r<sup>2</sup> = ", signif(LDT$R2, 3))
                      }),
                      "<br>RAF = ", LDT$Risk_Allele,
                      "<br>OR/beta = ", LDT$Effect_Size_95_CI,
                      "<br>CI = ", LDT$Beta_or_OR,
                      "<br>Trait: ", strwrap2(LDT$GWAS_Trait))
  if (!is.null(tc)) {
    hovertext <- paste0(hovertext, "<br>Category: ", LDT$category)
  }
  n_assoc <- if (collapse) {
    w <- LDT$n_assoc > 1
    hovertext[w] <- paste0(hovertext[w], "<br>Associations: ", LDT$n_assoc[w])
  }
  
  p <- plot_ly(x = LDT$pos / 1e6, y = LDT$logP,
               marker = list(opacity = 0.8,
                             size = size,
                             color = I(LDT$colour),
                             line = list(width = 1, color = outline)),
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
  
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    return(setNames(seq_along(ox), ox))
  }
  simmat <- stringdist::stringsimmatrix(x, x, method)
  hc <- hclust(as.dist(1 - simmat), method = "complete")
  grp <- cutree(hc, h = 1 - thresh)
  setNames(grp, ox)
}


removeWords <- function(words, x) {
  p <- paste0("(?:\\s+|^)\\S*(?<!\\w)(?:", paste(words, collapse = "|"), ")(?!\\w)\\S*")
  trimws(gsub(p, "", x, perl = TRUE))
}


strwrap2 <- function(x, w = 40) {
  sw <- strwrap(x, w, simplify = FALSE)
  vapply(sw, paste, collapse = "<br>   ", character(1))
}
