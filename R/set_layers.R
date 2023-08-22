
#' Set up a column of multiple plots
#'
#' Uses [layout()] to set up multiple locus plots aligned in a column.
#'
#' @param n Number of plots (not including gene tracks on bottom)
#' @param heights Vector of length `nrow + 1` specifying height for plots with
#'   a gene track on the bottom
#' @param rev Logical whether to reverse plotting order and plot from bottom to 
#'   top
#' @return Sets [layout()] to enable multiple plots aligned in a column. The
#'   gene track is assumed to be positioned on the bottom. Returns `par()`
#'   invisibly so that layout can be reset to default at the end of plotting.
#' @seealso [layout()]
#' @export

set_layers <- function(n = 1,
                       heights = c(rep(3, n), 2),
                       rev = FALSE) {
  op <- par(no.readonly = TRUE)
  s <- if (rev) rev(seq_len(n +1)) else seq_len(n +1)
  mat <- matrix(s)
  graphics::layout(mat, heights = heights)
  invisible(op)
}
