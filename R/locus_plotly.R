
#' Locus plotly
#' 
#' Genomic locus plot similar to locuszoom, using plotly.
#' 
#' @details 
#' This is an R/plotly version of locuszoom for generating publication ready
#' Manhattan plots of gene loci. It references Ensembl databases for annotating
#' genes and exons. Use [locus()] first to generate an object of class 'locus'
#' for plotting. LDlink web server can be queried using function [link_LD()] to
#' retrieve linkage disequilibrium (LD) information on the index SNP.
#' 
#' @param loc Object of class 'locus' to use for plot. See [locus()].
#' @param heights Vector controlling relative height of each panel on 0-1 scale. 
#' @param ... Optional arguments passed to [genetrack_ly()] to control the gene
#'   tracks.
#' @export

locus_plotly <- function(loc, heights = c(0.6, 0.4), ...) {
  g <- genetrack_ly(loc, ...)
  p <- scatter_plotly(loc)
  
  plotly::subplot(p, g, shareX = T, nrows = 2, heights = c(0.6, 0.4),
                  titleY = TRUE)
}
