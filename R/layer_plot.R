





set_layers <- function(nrow = 1,
                       heights = c(rep(3, nrow), 2)) {
  op <- par(no.readonly = TRUE)
  mat <- matrix(seq_len(nrow +1))
  graphics::layout(mat, heights = heights)
  invisible(op)
}

scatter_plot <- function(x,
                         pcutoff = 5e-08,
                         chromCols = 'royalblue',
                         sigCol = 'red',
                         xlab = NULL,
                         ylab = expression("-log"[10] ~ "P"),
                         cex.axis = 0.8,
                         cex.text = 0.7,
                         xticks = FALSE,
                         border = FALSE,
                         LDcols = c('grey', 'royalblue', 'cyan2', 'green3', 
                                    'orange', 'red', 'purple'),
                         legend_pos = 'topleft',
                         align = TRUE, ...) {
  data <- x$data
  TX <- x$TX
  EX <- x$EX
  if (is.null(xlab)) xlab <- paste("Chromosome", x$seqname, "(Mb)")
  LD <- "ld" %in% colnames(data)
  if (LD) {
    data$col <- cut(data$ld, -1:6/5, labels = FALSE)
    data$col[is.na(data$col)] <- 1L
    data$col[which.max(data$logP)] <- 7L
    data <- data[order(data$col), ]
    data$col <- LDcols[data$col]
  } else {
    data$col <- chromCols
    data$col[data[, x$p] < pcutoff] <- sigCol
  }
  
  # scatter plot
  if (align) {
    op <- par(tcl = -0.25, las = 1, font.main = 1,
              mar = c(ifelse(xticks, 3, 0.1), 4, 2, 1.5),
              mgp = c(1.7, 0.5, 0))
    on.exit(par(op))
  }
  
  plot(data[, x$pos], data$logP,
       pch = 21, bg = data$col,
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
  if (!is.null(legend_pos)) {
    if (LD) {
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

