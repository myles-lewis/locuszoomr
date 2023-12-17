
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
#' @param filter_gene_name Vector of gene names to display.
#' @param filter_gene_biotype Vector of gene biotypes to be filtered. Use
#' [ensembldb::listGenebiotypes()] to display possible biotypes. For example, 
#' `ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)`
#' @param cex.text Font size for gene text.
#' @param maxrows Specifies maximum number of rows to display in gene 
#' annotation panel.
#' @param xlab Title for x axis. Defaults to chromosome `seqname` specified 
#' in `locus`.
#' @param gene_col Colour for gene lines.
#' @param exon_col Fill colour for exons.
#' @param exon_border Border line colour outlining exons. Set to `NA` for no 
#' border.
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
                         gene_col = 'blue4',
                         exon_col = 'blue4',
                         exon_border = 'blue4',
                         maxrows = 8,
                         xlab = NULL,
                         ...) {
  g <- genetrack_ly(loc, filter_gene_name, filter_gene_biotype, cex.text, 
                    gene_col, exon_col, exon_border, maxrows, xlab)
  p <- scatter_plotly(loc, xlab = xlab, ...)
  
  plotly::subplot(p, g, shareX = TRUE, nrows = 2, heights = heights,
                  titleY = TRUE)
}
