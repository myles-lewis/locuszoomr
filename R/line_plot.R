
#' Locus line plot
#'
#' Produces a line plot from a 'locus' class object. Intended for use with
#' [set_layers()].
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param cex.axis Specifies font size for axis numbering.
#' @param xticks Logical whether x axis numbers and axis title are plotted.
#' @param border Logical whether a bounding box is plotted around upper and
#'   lower plots.
#' @param align Logical whether set [par()] to align the plot.
#' @param ... Other arguments passed to [plot()] for the scatter plot.
#' @return No return value. Produces a scatter plot using base graphics.
#' @seealso [locus()] [set_layers()] [scatter_plot()]
#' @export
#' 
line_plot <- function(loc,
                      pcutoff = 5e-08,
                      xlab = NULL,
                      ylab = expression("-log"[10] ~ "P"),
                      cex.axis = 1,
                      xticks = FALSE,
                      border = FALSE,
                      align = TRUE, ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  if (is.null(loc$data)) stop("No data points, only gene tracks")
  data <- loc$data
  if (is.null(xlab)) xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  
  # line plot
  if (align) {
    op <- par(mar = c(ifelse(xticks, 3, 0.1), 3.5, 2, 1.5))
    on.exit(par(op))
  }
  
  if (!is.null(pcutoff)) {
    abl <- quote(abline(h = -log10(pcutoff), col = 'darkgrey', lty = 2))
  } else abl <- NULL
  
  new.args <- list(...)
  plot.args <- list(x = data[, loc$pos], y = data$logP,
                    type = "l",
                    las = 1, font.main = 1,
                    xlim = loc$xrange,
                    xlab = if (xticks) xlab else "",
                    ylab = ylab,
                    bty = if (border) 'o' else 'l',
                    cex.axis = cex.axis,
                    xaxt = 'n',
                    tcl = -0.3, 
                    mgp = c(1.7, 0.5, 0),
                    panel.first = abl)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  
  if (xticks) {
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis,
         mgp = c(1.7, 0.4, 0), tcl = -0.3)
  } else if (!border) {
    axis(1, at = axTicks(1), labels = FALSE, tcl = -0.3)
  }
}
