
#' Add gene tracks to a ggplot2 plot
#'
#' Adds gene tracks to an existing ggplot2 plot.
#'
#' @param p ggplot2 plot object
#' @param loc Object of class 'locus' to use for plot. See [locus()].
#' @param height A unit object specifying height of the lower gene track.
#' @param ... Additional arguments passed to [gg_genetracks()] to control
#'   colours of gene tracks etc.
#' @return A gtable plotting object.
#' @details
#' The combined plot is outputted to the current device. A gtable plotting object is
#' returned invisibly. This can be plotted using [grid.draw()].
#' @seealso [gg_scatter()] [gg_genetracks()]
#' @importFrom ggplot2 ggplotGrob find_panel
#' @importFrom gtable gtable_add_rows gtable_add_grob
#' @importFrom grid grid.newpage grid.draw
#' @export

gg_addgenes <- function(p, loc,
                        height = unit(5, "cm"),
                        ...) {
  g <- ggplotGrob(p)
  panels_extent <- g %>% find_panel()
  pg <- g %>%
    gtable_add_rows(heights = height) %>%
    gtable_add_grob(gg_genetracks(loc, ...),
                    t = -1, b = -1,
                    l = panels_extent$l, r = panels_extent$l +1)
  grid.newpage()
  grid.draw(pg)
  invisible(pg)
}
