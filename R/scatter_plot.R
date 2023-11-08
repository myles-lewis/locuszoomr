
#' Locus scatter plot
#'
#' Produces a scatter plot from a 'locus' class object. Intended for use with
#' [set_layers()].
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param index_snp Specifies index SNP to be shown in a different colour and
#'   symbol. Defaults to the SNP with the lowest p-value. Set to `NULL` to not
#'   show this.
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
#' @param yzero Logical whether to force y axis limit to include y=0.
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
#' @param labels Character vector of SNP or genomic feature IDs to label. The
#'   value "index" selects the highest point or index SNP as defined when
#'   [locus()] is called. Set to `NULL` to remove all labels.
#' @param label_x Value or vector for position of label as percentage of x axis
#'   scale.
#' @param label_y Value or vector for position of label as percentage of y axis
#'   scale.
#' @param add Logical whether to add points to an existing plot or generate a
#'   new plot.
#' @param align Logical whether to set [par()] to align the plot.
#' @param ... Other arguments passed to [plot()] to control the scatter plot
#'   e.g. `main`, `ylim` etc.
#' @return No return value. Produces a scatter plot using base graphics.
#' @seealso [locus()] [set_layers()]
#' @importFrom graphics par legend
#' @export
#' 
scatter_plot <- function(loc,
                         index_snp = loc$index_snp,
                         pcutoff = 5e-08,
                         chromCol = 'royalblue',
                         sigCol = 'red',
                         cex = 1,
                         cex.axis = 0.9,
                         cex.lab = 1,
                         xlab = NULL,
                         ylab = NULL,
                         yzero = (loc$yvar == "logP"),
                         xticks = TRUE,
                         border = FALSE,
                         showLD = TRUE,
                         LD_scheme = c('grey', 'royalblue', 'cyan2', 'green3', 
                                       'orange', 'red', 'purple'),
                         legend_pos = 'topleft',
                         labels = NULL,
                         label_x = 4, label_y = 4,
                         add = FALSE,
                         align = TRUE, ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  data <- loc$data
  if (is.null(xlab)) xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  if (is.null(ylab)) {
    ylab <- if (loc$yvar == "logP") expression("-log"[10] ~ "P") else loc$yvar
  }
  hasLD <- "ld" %in% colnames(data)
  if (!"bg" %in% colnames(data)) {
    if (showLD & hasLD) {
      data$bg <- cut(data$ld, -1:6/5, labels = FALSE)
      data$bg[is.na(data$bg)] <- 1L
      data$bg[data[, loc$labs] == index_snp] <- 7L
      data <- data[order(data$bg), ]
      LD_scheme <- rep_len(LD_scheme, 7)
      data$bg <- LD_scheme[data$bg]
    } else {
      data$bg <- chromCol
      if (loc$yvar == "logP") data$bg[data[, loc$p] < pcutoff] <- sigCol
      data$bg[data[, loc$labs] == index_snp] <- "purple"
    }
  }
  
  # scatter plot
  if (align) {
    op <- par(mar = c(ifelse(xticks, 3, 0.1), 3.5, 2, 1.5))
    on.exit(par(op))
  }
  
  if (loc$yvar == "logP" & !is.null(pcutoff)) {
    abl <- quote(abline(h = -log10(pcutoff), col = 'darkgrey', lty = 2))
  } else abl <- NULL
  
  pch <- rep(21L, nrow(data))
  pch[data[, loc$labs] == index_snp] <- 23L
  if ("pch" %in% colnames(data)) pch <- data$pch
  col <- "black"
  if ("col" %in% colnames(data)) col <- data$col
  
  new.args <- list(...)
  if (add) {
    plot.args <- list(x = data[, loc$pos], y = data[, loc$yvar],
                      pch = pch, bg = data$bg, cex = cex)
    if (length(new.args)) plot.args[names(new.args)] <- new.args
    return(do.call("points", plot.args))
  }
  ylim <- range(data[, loc$yvar], na.rm = TRUE)
  if (yzero) ylim[1] <- min(c(0, ylim[1]))
  plot.args <- list(x = data[, loc$pos], y = data[, loc$yvar],
               pch = pch, bg = data$bg, col = col,
               las = 1, font.main = 1,
               cex = cex, cex.axis = cex.axis, cex.lab = cex.lab,
               xlim = loc$xrange,
               ylim = ylim,
               xlab = if (xticks) xlab else "",
               ylab = ylab,
               bty = if (border) 'o' else 'l',
               xaxt = 'n',
               tcl = -0.3, 
               mgp = c(1.7, 0.5, 0),
               panel.first = abl)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  
  # add labels
  if (!is.null(labels)) {
    i <- grep("index", labels, ignore.case = TRUE)
    if (i) labels[i] <- index_snp
    ind <- match(labels, data[, loc$labs])
    if (any(is.na(ind))) {
      message("label ", paste(labels[is.na(ind)], collapse = ", "),
              " not found")
    }
    lx <- data[ind, loc$pos]
    ly <- data[ind, loc$yvar]
    labs <- data[ind, loc$labs]
    add_labels(lx, ly, labs, label_x, label_y, cex = cex.axis *0.95)
  }
  
  if (xticks) {
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis,
         mgp = c(1.7, 0.4, 0), tcl = -0.3)
  } else if (!border) {
    axis(1, at = axTicks(1), labels = FALSE, tcl = -0.3)
  }
  if (!is.null(legend_pos)) {
    if (showLD & hasLD) {
      legend(legend_pos,
             legend = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4",
                        "0.0 - 0.2"),
             title = expression({r^2}), y.intersp = 0.96,
             pch = 21, col = 'black', pt.bg = rev(LD_scheme[-c(1, 7)]), 
             bty = 'n', cex = 0.8)
    }
  }
}


add_labels <- function(lx, ly, labs, label_x, label_y, cex = 1) {
  label_x <- rep_len(label_x, length(lx))
  label_y <- rep_len(label_y, length(ly))
  dx <- diff(par("usr")[1:2]) * label_x /100
  dy <- diff(par("usr")[3:4]) * label_y /100
  dlines(lx, ly, dx, dy, xpd = NA)
  adj1 <- -sign(dx) *0.56+0.5
  adj2 <- -sign(dy) +0.5
  adj2[abs(label_x) > abs(label_y)] <- 0.5
  adj1[abs(label_x) < abs(label_y)] <- 0.5
  if (length(unique(adj1)) == 1 & length(unique(adj2)) == 1) {
    # unique adj
    adj <- c(adj1[1], adj2[1])
    text(lx + dx, ly + dy, labs,
         adj = adj, cex = cex, xpd = NA)
  } else {
    # varying adj
    adj <- cbind(adj1, adj2)
    for (i in seq_along(labs)) {
      text(lx[i] + dx[i], ly[i] + dy[i], labs[i],
           adj = adj[i,], cex = cex, xpd = NA)
    }
  }
}


dlines <- function(x, y, dx, dy, ...) {
  mx <- cbind(x, x + dx, NA)
  my <- cbind(y, y + dy, NA)
  xs <- as.vector(t(mx))
  ys <- as.vector(t(my))
  lines(xs, ys, ...)
}
