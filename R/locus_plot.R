
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
#' Arguments to control plotting of the gene tracks are passed onto
#' [genetracks()] and for the scatter plot are passed via `...` to
#' [scatter_plot()]. See the documentation for each of these functions for
#' details.
#' 
#' @param loc Object of class 'locus' to use for plot. See [locus()].
#' @param filter_gene_name Vector of gene names to display.
#' @param filter_gene_biotype Vector of gene biotypes to be filtered. Use
#' [ensembldb::listGenebiotypes()] to display possible biotypes. For example, 
#' `ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)`
#' @param xlab x axis title.
#' @param cex Specifies size for points.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.lab Specifies font size for axis titles.
#' @param cex.text Font size for gene text.
#' @param use_layout Logical whether `graphics::layout` is called. Default
#'   `TRUE` is for a standard single plot. Set to `FALSE` if a more complex
#'   layout with multiple plots is required e.g. using [multi_layout()].
#' @param heights Ratio of top to bottom plot. See [layout].
#' @param showExons Logical whether to show exons or simply show whole gene as a
#'   rectangle
#' @param maxrows Specifies maximum number of rows to display in gene 
#' annotation panel.
#' @param xticks Character value of either 'top' or 'bottom' specifying 
#' whether x axis ticks and numbers are plotted on top or bottom plot window.
#' @param border Logical whether a bounding box is plotted around upper and 
#' lower plots.
#' @param gene_col Colour for gene lines.
#' @param exon_col Fill colour for exons.
#' @param exon_border Border line colour outlining exons (or genes if
#'   `showExons` is `FALSE`). Set to `NA` for no border.
#' @param text_pos Character value of either 'top' or 'left' specifying
#'   placement of gene name labels.
#' @param italics Logical whether gene text is in italics.
#' @param highlight Vector of genes to highlight.
#' @param highlight_col Single colour or vector of colours for highlighted
#'   genes.
#' @param blanks Controls handling of genes with blank names: `"fill"` replaces
#'   blank gene symbols with ensembl gene ids. `"hide"` hides genes which are
#'   missing gene symbols.
#' @param recomb_col Colour for recombination rate line if recombination rate
#'   data is present. Set to `NA` to hide the line. See [link_recomb()] to add
#'   recombination rate data.
#' @param ... Other arguments passed to [scatter_plot()] and [plot()] to control
#'   the scatter plot, e.g. `ylab`, `main`, etc.
#' @return No return value.
#' @seealso [locus()] [scatter_plot()] [genetracks()]
#' @examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5,
#'              ens_db = "EnsDb.Hsapiens.v75")
#' locus_plot(loc)
#' 
#' ## Use embedded LD information in column `r2`
#' loc2 <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
#'               ens_db = "EnsDb.Hsapiens.v75")
#' ## Add label for index SNP
#' locus_plot(loc2, labels = "index")
#' }
#' @export

locus_plot <- function(loc,
                       filter_gene_name = NULL,
                       filter_gene_biotype = NULL,
                       xlab = NULL,
                       cex = 1,
                       cex.axis = 0.9,
                       cex.lab = 1,
                       cex.text = 0.7,
                       use_layout = TRUE,
                       heights = c(3, 2),
                       showExons = TRUE,
                       maxrows = 7,
                       xticks = 'bottom',
                       border = FALSE,
                       gene_col = ifelse(showExons, 'blue4', 'skyblue'),
                       exon_col = 'blue4',
                       exon_border = 'blue4',
                       text_pos = 'top',
                       italics = FALSE,
                       highlight = NULL,
                       highlight_col = "red",
                       blanks = 'fill',
                       recomb_col = "blue", ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  if (is.null(loc$data)) stop("No SNPs/data points")
  
  if (use_layout) {
    op0 <- set_layers(1, heights, rev = TRUE)
    on.exit(par(op0), add = TRUE)
  }
  
  # lower panel gene tracks at locus
  genetracks(loc, filter_gene_name, filter_gene_biotype,
             border, cex.axis, cex.lab, cex.text, gene_col, exon_col, exon_border,
             showExons, maxrows, text_pos, italics,
             xticks = (xticks == 'bottom'),
             xlab = if (xticks == 'bottom') xlab else "",
             highlight = highlight, highlight_col = highlight_col,
             blanks = blanks, showRecomb = !is.na(recomb_col))
  
  # upper panel plot points
  scatter_plot(loc, xticks = (xticks == 'top'),
               border = border, xlab = xlab,
               cex = cex, cex.axis = cex.axis, cex.lab = cex.lab,
               recomb_col = recomb_col, ...)
}
