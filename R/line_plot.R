
line_plot <- function(x,
                      pcutoff = 5e-08,
                      xlab = NULL,
                      ylab = expression("-log"[10] ~ "P"),
                      cex.axis = 0.8,
                      cex.text = 0.7,
                      xticks = FALSE,
                      border = FALSE,
                      align = TRUE, ...) {
  data <- x$data
  TX <- x$TX
  EX <- x$EX
  if (is.null(xlab)) xlab <- paste("Chromosome", x$seqname, "(Mb)")
  
  # line plot
  if (align) {
    op <- par(tcl = -0.25, las = 1, font.main = 1,
              mar = c(ifelse(xticks, 3, 0.1), 4, 2, 1.5),
              mgp = c(1.7, 0.5, 0))
    on.exit(par(op))
  }
  
  plot(data[, x$pos], data$logP,
       type = "l",
       xlim = x$xrange,
       xlab = if (xticks == 'top') xlab else "",
       ylab = ylab,
       bty = if (border) 'o' else 'l',
       cex.axis = cex.axis,
       xaxt = 'n',
       panel.first = {
         if (!is.null(pcutoff)) {
           abline(h = -log10(pcutoff), col = 'darkgrey', lty = 2)
         }
       }, ...)
  if (xticks) {
    par(mgp = c(1.6, 0.3, 0))
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis)
  } else {
    axis(1, at = axTicks(1), labels = FALSE)
  }
}
