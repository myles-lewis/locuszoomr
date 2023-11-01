
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
#' @param ... Optional arguments passed to [genetrack_ly()] to control the gene
#'   tracks.
#' @returns A 'plotly' plotting object showing a scatter plot above gene tracks.
#' @seealso [genetrack_ly()] [locus()]
#' @examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' library(EnsDb.Hsapiens.v75)
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = "IRF5", flank = c(7e4, 2e5), LD = "r2",
#'              ens_db = "EnsDb.Hsapiens.v75")
#' locus_plotly(loc)
#' }
#' @export

locus_plotly <- function(loc, heights = c(0.6, 0.4), ...) {
  g <- genetrack_ly(loc, ...)
  p <- scatter_plotly(loc)
  
  plotly::subplot(p, g, shareX = T, nrows = 2, heights = c(0.6, 0.4),
                  titleY = TRUE)
}
