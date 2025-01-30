
#' Locus plot using ggplot2
#'
#' Genomic locus plot similar to locuszoom.
#'
#' Arguments to control plotting of the gene tracks are passed onto
#' [gg_genetracks()] and for the scatter plot are passed via `...` to
#' [gg_scatter()]. See the documentation for each of these functions for
#' details.
#'
#' @param loc Object of class 'locus' to use for plot. See [locus()].
#' @param heights Vector supplying the ratio of top to bottom plot.
#' @param filter_gene_name Vector of gene names to display.
#' @param filter_gene_biotype Vector of gene biotypes to be filtered. Use
#' [ensembldb::listGenebiotypes()] to display possible biotypes. For example, 
#' `ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)`
#' @param border Logical whether a bounding box is plotted.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.lab Specifies font size for axis titles.
#' @param cex.text Font size for gene text.
#' @param gene_col Colour for gene lines.
#' @param exon_col Fill colour for exons.
#' @param exon_border Border line colour outlining exons (or genes if
#'   `showExons` is `FALSE`). Set to `NA` for no border.
#' @param showExons Logical whether to show exons or simply show whole gene as a
#'   rectangle. If `showExons = FALSE` colours are specified by `exon_border`
#'   for rectangle border and `gene_col` for the fill colour.
#' @param maxrows Specifies maximum number of rows to display in gene annotation
#'   panel.
#' @param text_pos Character value of either 'top' or 'left' specifying
#'   placement of gene name labels.
#' @param italics Logical whether gene text is in italics.
#' @param xticks Logical whether x axis ticks and numbers are plotted.
#' @param xlab Title for x axis. Defaults to chromosome `seqname` specified in
#'   `locus`.
#' @param highlight Vector of genes to highlight.
#' @param highlight_col Single colour or vector of colours for highlighted
#'   genes.
#' @param blanks Controls handling of genes with blank names: `"fill"` replaces
#'   blank gene symbols with ensembl gene ids. `"hide"` hides genes which are
#'   missing gene symbols.
#' @param ... Additional arguments passed to [gg_scatter()] to control
#'   the scatter plot, e.g. `pcutoff`, `scheme`, `recomb_offset` etc.
#' @return Returns a ggplot2 plot containing a scatter plot with genetracks
#'   underneath.
#' @seealso [gg_scatter()] [gg_genetracks()]
#' @examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
#'              ens_db = "EnsDb.Hsapiens.v75")
#' locus_ggplot(loc)
#' }
#' @importFrom cowplot plot_grid
#' @export

locus_ggplot <- function(loc, heights = c(3, 2),
                         filter_gene_name = NULL,
                         filter_gene_biotype = NULL,
                         border = FALSE,
                         cex.axis = 1,
                         cex.lab = 1,
                         cex.text = 0.7,
                         gene_col = ifelse(showExons, 'blue4', 'skyblue'),
                         exon_col = 'blue4',
                         exon_border = 'blue4',
                         showExons = TRUE,
                         maxrows = 12,
                         text_pos = 'top',
                         italics = FALSE,
                         xticks = "top",
                         xlab = NULL,
                         highlight = NULL,
                         highlight_col = "red",
                         blanks = "fill",
                         ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  if (is.null(loc$data)) stop("No SNPs/data points")
  p <- gg_scatter(loc,
                  cex.axis = cex.axis,
                  cex.lab = cex.lab,
                  xlab = xlab,
                  xticks = (xticks == "top"),
                  border = border, ...)
  g <- gg_genetracks(loc,
                     filter_gene_name, filter_gene_biotype,
                     border,
                     cex.axis, cex.lab, cex.text,
                     gene_col, exon_col, exon_border,
                     showExons,
                     maxrows, text_pos, italics,
                     xticks = (xticks != "top"), xlab,
                     highlight, highlight_col,
                     blanks)

  plot_grid(p, g, nrow = 2, rel_heights = heights, align = "v")
}
