
#' Layout multiple locus plots
#'
#' Produces pages with multiple locus plots on.
#'
#' @param plots Either an 'expression' to be evaluated which is a series of
#'   calls to [locus_plot()] or similar plotting functions, or a list of 'locus'
#'   class objects which are plotted in sequence.
#' @param nrow Number of rows of plots
#' @param ncol Number of columns of plots
#' @param heights Vector of length 2 specifying height for plot and gene tracks
#' @param legend_pos A keyword either "topleft" or "topright" or `NULL` to hide
#'   the legend. Not invoked if `plots` is an expression. The legend is only
#'   shown on one plot on each page.
#' @param ... Optional arguments passed to [locus_plot()] if `plots` contains a
#'   list
#' @return No return value.
#' @seealso [locus_plot()]
#' @examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' library(EnsDb.Hsapiens.v75)
#' data(SLE_gwas_sub)
#' genes <- c("STAT4", "UBE2L3", "IRF5")
#' loclist <- lapply(genes, locus,
#'                   data = SLE_gwas_sub,
#'                   ens_db = "EnsDb.Hsapiens.v75",
#'                   LD = "r2")
#' ## produce 3 locus plots, one on each page
#' multi_layout(loclist)
#' 
#' ## place 3 locus plots in a row on a single page
#' multi_layout(loclist, ncol = 3)
#' 
#' ## full control
#' loc <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5, LD = "r2",
#'              ens_db = "EnsDb.Hsapiens.v75")
#' loc2 <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
#'               ens_db = "EnsDb.Hsapiens.v75")
#' loc3 <- locus(SLE_gwas_sub, gene = 'UBE2L3', LD = "r2",
#'               ens_db = "EnsDb.Hsapiens.v75")
#' multi_layout(ncol = 3,
#'              plots = {
#'                locus_plot(loc, use_layout = FALSE, legend_pos = 'topleft')
#'                locus_plot(loc2, use_layout = FALSE, legend_pos = NULL)
#'                locus_plot(loc3, use_layout = FALSE, legend_pos = NULL)
#'              })
#' 
#' }
#' @export

multi_layout <- function(plots,
                         nrow = 1, ncol = 1, heights = c(3, 2),
                         legend_pos = "topleft",
                         ...) {
  op <- par(no.readonly = TRUE)
  xr <- seq_len(ncol) * 2
  x2 <- c(xr, xr-1)
  x3 <- rep(x2, nrow) + (rep(seq_len(nrow), each = length(x2)) -1) * length(x2)
  mat <- matrix(x3, nrow = nrow * 2, ncol = ncol, byrow = TRUE)
  graphics::layout(mat, heights = rep(heights, nrow))
  on.exit(par(op))
  if (is.list(plots)) {
    npage <- nrow * ncol
    lpos <- if (!is.null(legend_pos) && legend_pos == "topleft") {1
      } else if (nrow == 1) 0 else ncol
    lapply(seq_along(plots), function(i) {
      leg <- NULL
      if (i %% npage == lpos) leg <- legend_pos
      locus_plot(plots[[i]], use_layout = FALSE, legend_pos = leg, ...)
    })
  } else plots
  invisible()
}
