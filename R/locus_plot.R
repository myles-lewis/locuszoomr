
#' Locus plot
#' 
#' Genomic locus plot similar to locuszoom.
#' 
#' @details 
#' This is an R version of locuszoom for generating publication ready Manhattan 
#' plots of gene loci. It references Ensembl databases for annotating genes 
#' and exons. Use [locus()] first to generate an object of class 'locus' for
#' plotting. LDlink web server can be queried using function [link_LD()] to
#' retrieve linkage disequilibrium (LD) information on the index SNP.
#' 
#' @param x Object of class 'locus' to use for plot. See [locus()].
#' @param filter_gene_name Vector of gene names to display.
#' @param filter_gene_biotype Vector of gene biotypes to be filtered. Use
#' [ensembldb::listGenebiotypes()] to display possible biotypes. For example, 
#' `ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)`
#' @param xlab x axis title.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.text Font size for gene text.
#' @param use_layout Logical whether `graphics::layout` is called. Default
#'   `TRUE` is for a standard single plot. Set to `FALSE` if a more complex
#'   layout with multiple plots is required e.g. using [multi_layout()].
#' @param heights Ratio of top to bottom plot. See [layout].
#' @param maxrows Specifies maximum number of rows to display in gene 
#' annotation panel.
#' @param xticks Character value of either 'top' or 'bottom' specifying 
#' whether x axis ticks and numbers are plotted on top or bottom plot window.
#' @param border Logical whether a bounding box is plotted around upper and 
#' lower plots.
#' @param gene_col Colour for genes and exons.
#' @param exon_col Fill colour for exons.
#' @param exon_border Border line colour outlining exons. Set to `NA` for no 
#' border.
#' @param text_pos Character value of either 'top' or 'left' specifying 
#' placement of gene name labels.
#' @param ... Other arguments passed to [scatter_plot()] to control the scatter
#'   plot.
#' @return No return value.
#' @seealso [locus()] [scatter_plot()] [genetracks()]
#' @examples
#' library(EnsDb.Hsapiens.v75)
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5)
#' locus_plot(loc)
#' loc2 <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5)
#' locus_plot(loc2)
#' @export

locus_plot <- function(x,
                        filter_gene_name = NULL,
                        filter_gene_biotype = NULL,
                        xlab = NULL,
                        cex.axis = 1,
                        cex.text = 0.7,
                        use_layout = TRUE,
                        heights = c(3, 2),
                        maxrows = 7,
                        xticks = 'bottom',
                        border = FALSE,
                        gene_col = 'blue4',
                        exon_col = 'blue4',
                        exon_border = 'blue4',
                        text_pos = 'top', ...) {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  
  if (use_layout) {
    op0 <- set_layers(1, heights, rev = TRUE)
    on.exit(par(op0), add = TRUE)
  }
  # lower panel gene tracks at locus
  genetracks(x, filter_gene_name, filter_gene_biotype,
             border, cex.axis, cex.text, gene_col, exon_col, exon_border,
             maxrows, text_pos, xticks = (xticks == 'bottom'),
             xlab = if (xticks == 'bottom') xlab else "")
  
  # upper panel plot points
  scatter_plot(x, xticks = (xticks == 'top'),
               border = border, cex.axis = cex.axis, xlab = xlab, ...)
}