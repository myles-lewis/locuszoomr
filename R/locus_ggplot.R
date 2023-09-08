
#' Locus plot using ggplot2
#'
#' Genomic locus plot similar to locuszoom.
#'
#' @param loc Object of class 'locus' to use for plot. See [locus()].
#' @param height A unit object specifying height of the lower gene track.
#' @param xticks Character value of either 'top' or 'bottom' specifying whether
#'   x axis ticks and numbers are plotted on top or bottom plot window.
#' @seealso [gg_scatter()] [gg_genetracks()]
#' @importFrom ggplot2 ggplotGrob find_panel
#' @importFrom gtable gtable_add_rows gtable_add_grob
#' @importFrom grid grid.newpage grid.draw
#' @export

locus_ggplot <- function(loc, height = unit(5, "cm"), xticks = "top") {
  p <- gg_scatter(loc, xticks = (xticks == "top"))
  g <- ggplotGrob(p)
  panels_extent <- g %>% find_panel()
  g %>%
    gtable_add_rows(heights = height) %>%
    gtable_add_grob(gg_genetracks(loc, xticks = (xticks != "top")),
                    t = -1, b = -1,
                    l = panels_extent$l, r = panels_extent$l +1) %>%
    {grid.newpage(); grid.draw(.)}
}
