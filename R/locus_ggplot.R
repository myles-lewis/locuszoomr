
#' Locus plot using ggplot2
#'
#' Genomic locus plot similar to locuszoom.
#'
#' @param loc Object of class 'locus' to use for plot. See [locus()].
#' @param height A unit object specifying height of the lower gene track.
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
#' @param legend_pos Position of legend e.g. "topleft", "topright" or ggplot2
#'   settings. Set to `NULL` to hide legend.
#' @param draw Logical whether to call [grid.draw()] to draw the plot.
#' @param ... Additional arguments passed to [gg_genetracks()] to control
#'   colours of gene tracks etc.
#' @return Returns a ggplot2 plot containing a scatter plot with genetracks
#'   underneath.
#' @seealso [gg_scatter()] [gg_genetracks()]
#' @examples
#' library(EnsDb.Hsapiens.v75)
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
#'              ens_db = "EnsDb.Hsapiens.v75")
#' locus_ggplot(loc)
#' @importFrom ggplot2 ggplotGrob find_panel
#' @importFrom gtable gtable_add_rows gtable_add_grob
#' @importFrom grid grid.newpage grid.draw
#' @export

locus_ggplot <- function(loc, height = unit(5, "cm"),
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
                         legend_pos = 'topleft',
                         draw = TRUE,
                         ...) {
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
                  legend_pos = legend_pos)
  g <- ggplotGrob(p)
  panels_extent <- g %>% find_panel()
  pg <- g %>%
    gtable_add_rows(heights = height) %>%
    gtable_add_grob(gg_genetracks(loc, xticks = (xticks != "top"),
                                  border = border, xlab = xlab,
                                  cex.axis = cex.axis,
                                  cex.lab = cex.lab, draw = FALSE, ...),
                    t = -1, b = -1,
                    l = panels_extent$l, r = panels_extent$l +1)
  if (draw) {
    grid.newpage()
    grid.draw(pg)
  }
  invisible(pg)
}
