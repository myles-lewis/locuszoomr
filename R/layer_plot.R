
#' Set up a column of multiple plots
#'
#' Uses [layout()] to set up multiple locus plots aligned in a column.
#'
#' @param nrow Number of rows of locus plots
#' @param heights Vector of length `nrow + 1` specifying height for plots with
#'   a gene track on the bottom
#' @return Sets [layout()] to enable multiple plots aligned in a column. The
#'   gene track is assumed to be positioned on the bottom. Returns `par()`
#'   invisibly so that layout can be reset to default at the end of plotting.
#' @seealso [layout()]
#' @export

set_layers <- function(nrow = 1,
                       heights = c(rep(3, nrow-1), 2)) {
  op <- par(no.readonly = TRUE)
  mat <- matrix(seq_len(nrow +1))
  graphics::layout(mat, heights = heights)
  invisible(op)
}


#' Locus scatter plot
#'
#' Produces a scatter plot from a 'locus' class object. Intended for use with
#' [set_layers()].
#'
#' @param x Object of class 'locus' to use for plot. See [locus].
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param chromCols Colour for normal points if `LD` is `FALSE` when the locus
#'   object is made.
#' @param sigCol Colour for significant points if `LD` is `FALSE`.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param cex.axis Specifies font size for axis numbering.
#' @param xticks Logical whether x axis numbers and axis title are plotted.
#' @param border Logical whether a bounding box is plotted around upper and
#'   lower plots.
#' @param LDcols Vector of colours for plotting LD. The first colour is for SNPs
#'   which lack LD information. The next 5 colours are for r2 or D' LD results
#'   ranging from 0 to 1 in intervals of 0.2. The final colour is for the index
#'   SNP.
#' @param legend_pos Position of legend. See [legend()]. Set to `NULL` to hide
#'   legend.
#' @param align Logical whether to set [par()] to align the plot.
#' @param ... Other arguments passed to [plot()] for the scatter plot.
#' @return No return value. Produces a scatter plot using base graphics.
#' @seealso [locus()] [set_layers()]
#' @export
#' 
scatter_plot <- function(x,
                         pcutoff = 5e-08,
                         chromCols = 'royalblue',
                         sigCol = 'red',
                         xlab = NULL,
                         ylab = expression("-log"[10] ~ "P"),
                         cex.axis = 1,
                         xticks = FALSE,
                         border = FALSE,
                         showLD = TRUE,
                         LDcols = c('grey', 'royalblue', 'cyan2', 'green3', 
                                    'orange', 'red', 'purple'),
                         legend_pos = 'topleft',
                         align = TRUE, ...) {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  data <- x$data
  TX <- x$TX
  EX <- x$EX
  if (is.null(xlab)) xlab <- paste("Chromosome", x$seqname, "(Mb)")
  hasLD <- "ld" %in% colnames(data)
  if (showLD & hasLD) {
    data$col <- cut(data$ld, -1:6/5, labels = FALSE)
    data$col[is.na(data$col)] <- 1L
    data$col[which.max(data$logP)] <- 7L
    data <- data[order(data$col), ]
    LDcols <- rep_len(LDcols, 7)
    data$col <- LDcols[data$col]
  } else {
    data$col <- chromCols
    data$col[data[, x$p] < pcutoff] <- sigCol
  }
  
  # scatter plot
  if (align) {
    op <- par(tcl = -0.25,
              mar = c(ifelse(xticks, 3, 0.1), 4, 2, 1.5),
              mgp = c(1.7, 0.5, 0))
    on.exit(par(op))
  }
  
  new.args <- list(...)
  plot.args <- list(x = data[, x$pos], y = data$logP,
               pch = 21, bg = data$col,
               las = 1, font.main = 1,
               xlim = x$xrange,
               xlab = if (xticks) xlab else "",
               ylab = ylab,
               bty = if (border) 'o' else 'l',
               cex.axis = cex.axis,
               xaxt = 'n',
               panel.first = {
                 if (!is.null(pcutoff)) {
                   abline(h = -log10(pcutoff), col = 'darkgrey', lty = 2)
                 }
               })
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  
  # plot(data[, x$pos], data$logP,
  #      pch = 21, bg = data$col,
  #      xlim = x$xrange,
  #      xlab = if (xticks) xlab else "",
  #      ylab = ylab,
  #      bty = if (border) 'o' else 'l',
  #      cex.axis = cex.axis,
  #      xaxt = 'n',
  #      panel.first = {
  #        if (!is.null(pcutoff)) {
  #          abline(h = -log10(pcutoff), col = 'darkgrey', lty = 2)
  #        }
  #      }, ...)
  if (xticks) {
    par(mgp = c(1.6, 0.3, 0))
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis)
  } else {
    axis(1, at = axTicks(1), labels = FALSE)
  }
  if (!is.null(legend_pos)) {
    if (showLD & hasLD) {
      legend(legend_pos,
             legend = c('Index SNP',
                        expression({0.8 < r^2} <= "1.0"),
                        expression({0.6 < r^2} <= 0.8),
                        expression({0.4 < r^2} <= 0.6),
                        expression({0.2 < r^2} <= 0.4),
                        expression({"0.0" < r^2} <= 0.2),
                        expression("No" ~ r^2 ~ "data")),
             pch = 21, col = 'black', pt.bg = rev(LDcols), bty = 'n', cex = 0.7)
    }
  }
}

