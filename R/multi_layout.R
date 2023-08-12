
#' Layout multiple locus plots
#'
#' Uses [layout()] to set up a pages of multiple locus plots.
#'
#' @param nrow Number of rows of plots
#' @param ncol Number of columns of plots
#' @param heights Vector of length 2 specifying height for plot and gene tracks
#' @return Sets [layout()] to enable multiple plots per page. Returns `par()` so
#'   that layout can be reset to default at the end of plotting.
#' @export

multi_layout <- function(nrow = 1, ncol = 1, heights = c(3, 2)) {
  op <- par(no.readonly = TRUE)
  xr <- seq_len(ncol) * 2
  x2 <- c(xr, xr-1)
  x3 <- rep(x2, nrow) + (rep(seq_len(nrow), each = length(x2))-1)*length(x2)
  mat <- matrix(x3, nrow = nrow * 2, ncol = ncol, byrow = TRUE)
  graphics::layout(mat, heights = rep(heights, nrow))
  op
}
