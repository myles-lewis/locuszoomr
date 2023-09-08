
#' Locus scatter plot
#'
#' Produces a scatter plot from a 'locus' class object. Intended for use with
#' [set_layers()].
#'
#' @param x Object of class 'locus' to use for plot. See [locus].
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param chromCol Colour for normal points if `LD` is `FALSE` when the locus
#'   object is made.
#' @param sigCol Colour for significant points if `LD` is `FALSE`.
#' @param cex Specifies size for points.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.lab Specifies font size for axis titles.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param xticks Logical whether x axis numbers and axis title are plotted.
#' @param border Logical whether a bounding box is plotted around upper and
#'   lower plots.
#' @param showLD Logical whether to show LD with colours
#' @param LD_scheme Vector of colours for plotting LD. The first colour is for SNPs
#'   which lack LD information. The next 5 colours are for r2 or D' LD results
#'   ranging from 0 to 1 in intervals of 0.2. The final colour is for the index
#'   SNP.
#' @param legend_pos Position of legend. See [legend()]. Set to `NULL` to hide
#'   legend.
#' @param add Logical whether to add points to an existing plot or generate a
#'   new plot.
#' @param align Logical whether to set [par()] to align the plot.
#' @param ... Other arguments passed to [plot()] for the scatter plot.
#' @return No return value. Produces a scatter plot using base graphics.
#' @seealso [locus()] [set_layers()]
#' @importFrom graphics par legend
#' @export
#' 
scatter_plot <- function(x,
                         pcutoff = 5e-08,
                         chromCol = 'royalblue',
                         sigCol = 'red',
                         cex = 1,
                         cex.axis = 1,
                         cex.lab = 1,
                         xlab = NULL,
                         ylab = expression("-log"[10] ~ "P"),
                         xticks = FALSE,
                         border = FALSE,
                         showLD = TRUE,
                         LD_scheme = c('grey', 'royalblue', 'cyan2', 'green3', 
                                       'orange', 'red', 'purple'),
                         legend_pos = 'topleft',
                         add = FALSE,
                         align = TRUE, ...) {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  data <- x$data
  if (is.null(xlab)) xlab <- paste("Chromosome", x$seqname, "(Mb)")
  hasLD <- "ld" %in% colnames(data)
  if (!"bg" %in% colnames(data)) {
    if (showLD & hasLD) {
      data$bg <- cut(data$ld, -1:6/5, labels = FALSE)
      data$bg[is.na(data$bg)] <- 1L
      data$bg[which.max(data$logP)] <- 7L
      data <- data[order(data$bg), ]
      LD_scheme <- rep_len(LD_scheme, 7)
      data$bg <- LD_scheme[data$bg]
    } else {
      data$bg <- chromCol
      data$bg[data[, x$p] < pcutoff] <- sigCol
    }
  }
  
  # scatter plot
  if (align) {
    op <- par(mar = c(ifelse(xticks, 3, 0.1), 3.5, 2, 1.5))
    on.exit(par(op))
  }
  
  if (!is.null(pcutoff)) {
    abl <- quote(abline(h = -log10(pcutoff), col = 'darkgrey', lty = 2))
  } else abl <- NULL
  
  pch <- 21
  if ("pch" %in% colnames(data)) pch <- data$pch
  col <- "black"
  if ("col" %in% colnames(data)) col <- data$col
  
  new.args <- list(...)
  if (add) {
    plot.args <- list(x = data[, x$pos], y = data$logP,
                      pch = pch, bg = data$bg, cex = cex)
    if (length(new.args)) plot.args[names(new.args)] <- new.args
    do.call("points", plot.args)
    return()
  }
  plot.args <- list(x = data[, x$pos], y = data$logP,
               pch = pch, bg = data$bg, col = col,
               las = 1, font.main = 1,
               cex = cex, cex.axis = cex.axis, cex.lab = cex.lab,
               xlim = x$xrange,
               ylim = c(0, max(data$logP, na.rm = TRUE)),
               xlab = if (xticks) xlab else "",
               ylab = ylab,
               bty = if (border) 'o' else 'l',
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
             pch = 21, col = 'black', pt.bg = rev(LD_scheme), bty = 'n', cex = 0.8)
    }
  }
}

