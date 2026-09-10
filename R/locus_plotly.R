
#' Locus plotly
#' 
#' Genomic locus plot similar to locuszoom, using plotly.
#' 
#' @details 
#' This is an R/plotly version of locuszoom for exploring regional Manhattan
#' plots of gene loci. Use [locus()] first to generate an object of class
#' 'locus' for plotting. This references a selected Ensembl database for
#' annotating genes and exons. Hover over the points or gene tracks to reveal
#' more information.
#' 
#' @param loc Object of class 'locus' to use for plot. See [locus()].
#' @param loc2 Optional 2nd 'locus' object to be layered underneath the 1st
#'   scatter plot.
#' @param heights Vector controlling relative height of each panel on 0-1 scale.
#'   Alternatively a vector of length 2 of height in pixels passed to
#'   `scatter_plotly()` and `genetrack_ly()`.
#' @param filter_gene_name Vector of gene names to display.
#' @param filter_gene_biotype Vector of gene biotypes to be filtered. Use
#' [ensembldb::listGenebiotypes()] to display possible biotypes. For example, 
#' `ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)`
#' @param cex.text Font size for gene text.
#' @param italics Logical whether gene text is in italics.
#' @param gene_col Colour for gene lines.
#' @param exon_col Fill colour for exons.
#' @param exon_border Border line colour outlining exons (or genes if
#'   `showExons` is `FALSE`). Set to `NA` for no border.
#' @param showExons Logical whether to show exons or simply show whole gene as a
#'   rectangle. If `showExons = FALSE` colours are specified by `exon_border`
#'   for rectangle border and `gene_col` for the fill colour.
#' @param maxrows Specifies maximum number of rows to display in gene 
#' annotation panel.
#' @param width Width of plotly plot in pixels which is purely used to prevent
#'   overlapping text for gene names.
#' @param xlab Title for x axis. Defaults to chromosome `seqname` specified 
#' in `locus`.
#' @param ylab Title for y axis, or a vector of 2 titles for each y axis if
#'   `loc2` is provided.
#' @param prioritise Vector of genes to be placed first in the gene tracks.
#' @param blanks Controls handling of genes with blank names: `"fill"` replaces
#'   blank gene symbols with ensembl gene ids. `"hide"` completely hides genes
#'   which are missing gene symbols. `"show"` shows gene lines but no label
#'   (hovertext is still available).
#' @param beta Optional column name for beta coefficient to display upward
#'   triangles for positive beta and downward triangles for negative beta
#'   (significant SNPs only). If `loc2` is supplied, then a vector can be used
#'   to specify different beta columns in `loc` and `loc2`; use `NA` to indicate
#'   no `beta`.
#' @param ... Optional arguments passed to [scatter_plotly()] to control the
#'   scatter plot.
#' @returns A 'plotly' plotting object showing a scatter plot above gene tracks.
#' @seealso [locus()] [genetrack_ly()] [scatter_plotly()]
#' @examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = "IRF5", flank = c(7e4, 2e5), LD = "r2",
#'              ens_db = "EnsDb.Hsapiens.v75")
#' locus_plotly(loc)
#' }
#' @export

locus_plotly <- function(loc,
                         loc2 = NULL,
                         heights = c(0.6, 0.4),
                         filter_gene_name = NULL,
                         filter_gene_biotype = NULL,
                         cex.text = 0.7,
                         italics = FALSE,
                         gene_col = ifelse(showExons, 'blue4', 'skyblue'),
                         exon_col = 'blue4',
                         exon_border = 'blue4',
                         showExons = TRUE,
                         maxrows = 8,
                         width = 600,
                         xlab = NULL,
                         ylab = NULL,
                         prioritise = NULL,
                         blanks = "show",
                         beta = NULL,
                         ...) {
  if (!is.null(loc2) && length(heights) == 2) {
    heights <- c(0.375, 0.375, 0.25)
  }
  pheights <- NULL
  if (any(heights > 1)) {
    pheights <- heights
    pheights[length(pheights)] <- sum(heights)
    heights <- heights / sum(heights)
  }
  
  g <- genetrack_ly(loc, filter_gene_name, filter_gene_biotype, cex.text, 
                    italics, gene_col, exon_col, exon_border, showExons, 
                    maxrows, width, xlab, prioritise, blanks,
                    height = pheights[length(pheights)])
  p <- scatter_plotly(loc, xlab = xlab, ylab = ylab[1], height = pheights[1],
                      beta = beta[1], ...)
  if (!is.null(loc2)) {
    if (!is.null(beta)) beta <- rep_len(beta, 2)
    p2 <- scatter_plotly(loc2, xlab = xlab, ylab = ylab[2],
                         height = pheights[2], showlegend = FALSE, beta = beta[2])
    pp <- plotly::subplot(p, p2, g, shareX = TRUE, nrows = 3, heights = heights,
                           titleY = TRUE, margin = c(0, 0, 0, 0.02))
    return(remap_overlaying_yaxes(pp))
  }
  
  plotly::subplot(p, g, shareX = TRUE, nrows = 2, heights = heights,
                  titleY = TRUE, margin = 0)
}


# fix double y axis with >1 scatter_plotly
# from Tom Willis
remap_overlaying_yaxes <- function(p) {
  lay <- p$x$layout
  nms <- grep("^yaxis[0-9]*$", names(lay), value = TRUE)
  base <- nms[vapply(nms, function(n) is.null(lay[[n]]$overlaying), logical(1))]
  for (n in setdiff(nms, base)) {
    dom <- lay[[n]]$domain
    if (is.null(dom)) next
    hit <- base[vapply(base,
                       function(b) isTRUE(all.equal(lay[[b]]$domain, dom)),
                       logical(1))]
    if (length(hit) == 1L) {
      p$x$layout[[n]]$overlaying <- sub("^yaxis", "y", hit)
    }
  }
  p
}
