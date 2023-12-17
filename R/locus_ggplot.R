
#' Locus plot using ggplot2
#'
#' Genomic locus plot similar to locuszoom.
#'
#' @param loc Object of class 'locus' to use for plot. See [locus()].
#' @param heights Vector supplying the ratio of top to bottom plot.
#' @param index_snp Specifies index SNP to be shown in a different colour and
#'   symbol. Defaults to the SNP with the lowest p-value. Set to `NULL` to not
#'   show this.
#' @param xticks Character value of either 'top' or 'bottom' specifying whether
#'   x axis ticks and numbers are plotted on top or bottom plot window.
#' @param border Logical whether a bounding box is plotted around upper and
#'   lower plots.
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param scheme Vector of 3 colors if LD is not shown: 1st = normal points, 2nd
#'   = colour for significant points, 3rd = index SNP.
#' @param size Specifies size for points.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.lab Specifies font size for axis titles.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param xticks Logical whether x axis numbers and axis title are plotted.
#' @param border Logical whether a bounding box is plotted around top and bottom
#'   plots.
#' @param showLD Logical whether to show LD with colours
#' @param LD_scheme Vector of colours for plotting LD. The first colour is for
#'   SNPs which lack LD information. The next 5 colours are for r2 or D' LD
#'   results ranging from 0 to 1 in intervals of 0.2. The final colour is for
#'   the index SNP.
#' @param recomb_col Colour for recombination rate line if recombination rate
#'   data is present. Set to NA to hide the line. See [link_recomb()] to add
#'   recombination rate data.
#' @param legend_pos Position of legend e.g. "topleft", "topright" or ggplot2
#'   settings. Set to `NULL` to hide legend.
#' @param ... Additional arguments passed to [gg_genetracks()] to control
#'   colours of gene tracks etc.
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
                         index_snp = loc$index_snp,
                         pcutoff = 5e-08,
                         scheme = c('royalblue', 'red', 'purple'),
                         size = 2,
                         cex.axis = 1,
                         cex.lab = 1,
                         xlab = NULL,
                         ylab = NULL,
                         xticks = "top",
                         border = FALSE,
                         showLD = TRUE,
                         LD_scheme = c('grey', 'royalblue', 'cyan2', 'green3', 
                                       'orange', 'red', 'purple'),
                         recomb_col = "blue",
                         legend_pos = 'topleft',
                         ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  if (is.null(loc$data)) stop("No data points, only gene tracks")
  p <- gg_scatter(loc,
                  index_snp = index_snp,
                  pcutoff = pcutoff,
                  scheme = scheme,
                  size = size,
                  cex.axis = cex.axis,
                  cex.lab = cex.lab,
                  xlab = xlab,
                  ylab = ylab,
                  xticks = (xticks == "top"),
                  border = border,
                  showLD = showLD,
                  LD_scheme = LD_scheme,
                  recomb_col = recomb_col,
                  legend_pos = legend_pos)
  g <- gg_genetracks(loc, xticks = (xticks != "top"),
                     border = border, xlab = xlab,
                     cex.axis = cex.axis,
                     cex.lab = cex.lab, ...)
  
  plot_grid(p, g, nrow = 2, rel_heights = heights, align = "v")
}
