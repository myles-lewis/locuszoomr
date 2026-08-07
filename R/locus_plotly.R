
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
#' @param prioritise Vector of genes to be placed first in the gene tracks.
#' @param blanks Controls handling of genes with blank names: `"fill"` replaces
#'   blank gene symbols with ensembl gene ids. `"hide"` completely hides genes
#'   which are missing gene symbols. `"show"` shows gene lines but no label
#'   (hovertext is still available).
#' @param dynamic Logical whether to re-pack the gene tracks in the browser as
#'   the user zooms or pans. Requires JavaScript, so set to `FALSE` when
#'   exporting to a static image.
#' @param scrollZoom Logical whether the mouse wheel zooms the plot. Defaults to
#'   `FALSE`, matching plotly's own default for cartesian plots, because an
#'   enabled plot captures the wheel and prevents the reader scrolling past it
#'   in a vignette or R Markdown document. Only the shared x axis is affected,
#'   since both y axes are fixed, so the -log10 P scale cannot be distorted.
#'   Combines with `dynamic`: scrolling re-packs the gene tracks just as
#'   dragging does.
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

locus_plotly <- function(loc, heights = c(0.6, 0.4),
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
                         prioritise = NULL,
                         blanks = "show",
                         dynamic = TRUE,
                         scrollZoom = FALSE,
                         ...) {
  pheights <- NULL
  if (any(heights > 1)) {
    pheights <- heights
    pheights[2] <- sum(heights)
    heights <- heights / sum(heights)
  }
  
  g <- genetrack_ly(loc, filter_gene_name, filter_gene_biotype, cex.text,
                    italics, gene_col, exon_col, exon_border, showExons,
                    maxrows, width, xlab, prioritise, blanks,
                    height = pheights[2], dynamic = FALSE)
  p <- scatter_plotly(loc, xlab = xlab, height = pheights[1], ...)

  gt <- attr(g, "genetrack_data")

  sp <- plotly::subplot(p, g, shareX = TRUE, nrows = 2, heights = heights,
                        titleY = TRUE, margin = 0)

  out <- if (!dynamic || is.null(gt)) {
    sp
  } else {
    # See the matching comment in genetrack_ly() (R/genetrack_ly.R): the R
    # side of add_genetrack_relayout() can stop() (e.g. if the gene track
    # traces can't be resolved after subplot()), which, left unguarded, would
    # turn a working locus_plotly(loc) call into a hard error since
    # dynamic = TRUE is the default. Warn and fall back to the static
    # subplot `sp` instead.
    tryCatch(
      add_genetrack_relayout(sp, gt$TX, gt$EX, gt$cfg),
      error = function(e) {
        warning("dynamic gene track disabled: ", conditionMessage(e),
                call. = FALSE)
        sp
      })
  }

  # Applied to the finished object rather than to the constituent panels, so
  # it reaches the combined widget on both the dynamic and static paths.
  # plotly::config() merges, so the modebar settings set upstream survive.
  if (scrollZoom) out <- plotly::config(out, scrollZoom = TRUE)
  out
}
