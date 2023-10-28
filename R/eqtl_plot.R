
#' Locus eQTL plot
#'
#' Produces a plot of eQTL data embedded in a 'locus' class object. Intended for
#' use with [set_layers()].
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param tissue GTex tissue in which eQTL has been measured
#' @param eqtl_gene Gene showing eQTL effect
#' @param scheme Character string specifying palette for effect size showing
#'   up/downregulation eQTL using [grDevices::hcl.colors]. Alternatively a
#'   vector of 6 colours.
#' @param col Outline point colour. `NA` for no outlines.
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param cex.axis Specifies font size for axis numbering.
#' @param xticks Logical whether x axis numbers and axis title are plotted.
#' @param border Logical whether a bounding box is plotted around upper and
#'   lower plots.
#' @param add Logical whether to add points to an existing plot or generate a
#'   new plot.
#' @param align Logical whether set [par()] to align the plot.
#' @param legend_pos Character value specifying legend position. See [legend()].
#' @param ... Other arguments passed to [plot()] for the scatter plot.
#' @return No return value. Produces a scatter plot using base graphics.
#' @seealso [locus()] [set_layers()] [scatter_plot()]
#' @importFrom graphics points
#' @export
#' 
eqtl_plot <- function(loc,
                      tissue = "Whole Blood",
                      eqtl_gene = loc$gene,
                      scheme = "RdYlBu",
                      col = NA,
                      pcutoff = NULL,
                      xlab = NULL,
                      ylab = expression("-log"[10] ~ "P"),
                      cex.axis = 0.9,
                      xticks = TRUE,
                      border = FALSE,
                      add = FALSE,
                      align = TRUE, 
                      legend_pos = "topright", ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  
  if (!"LDexp" %in% names(loc)) stop("Contains no eQTL data")
  data <- loc$LDexp
  data <- data[data$Tissue == tissue & data$Gene_Symbol == eqtl_gene, ]
  if (nrow(data) == 0) stop("No data")
  
  data$pos <- gsub(".*:", "", data$Position_grch37)  # remove up to ':'
  data$pos <- as.numeric(data$pos)
  data$logP <- -log10(data$P_value)
  
  # fix effect allele not being minor allele
  data$Effect_Allele_Freq <- gsub(".*=", "", data$Effect_Allele_Freq)
  data$Effect_Allele_Freq <- as.numeric(data$Effect_Allele_Freq)
  swap <- !is.na(data$Effect_Allele_Freq) & data$Effect_Allele_Freq > 0.5
  data$Effect_Size[swap] <- -data$Effect_Size[swap]
  data$pch <- -sign(data$Effect_Size) / 2 + 24.5
  equp <- sign(data$Effect_Size) == 1
  if (is.character(scheme)) {
    scheme <- hcl.colors(9, scheme)[-c(4:6)]
  }
  up_cols <- rev(scheme[1:3])
  down_cols <- scheme[4:6]
  ecol <- cut(abs(data$Effect_Size), breaks = 3)
  data$bg[equp] <- up_cols[ecol[equp]]
  data$bg[!equp] <- down_cols[ecol[!equp]]
  labs <- levels(ecol)
  cutlev <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
                  upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  cutlev <- signif(cutlev, 2)
  
  if (is.null(xlab)) xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  
  # scatter plot
  if (align) {
    op <- par(mar = c(ifelse(xticks, 3, 0.1), 3.5, 2, 1.5))
    on.exit(par(op))
  }
  
  if (!is.null(pcutoff)) {
    abl <- quote(abline(h = -log10(pcutoff), col = 'darkgrey', lty = 2))
  } else abl <- NULL
  
  new.args <- list(...)
  if (add) {
    plot.args <- list(x = data$pos, y = data$logP,
           pch = data$pch, bg = data$bg, col = col)
    if (length(new.args)) plot.args[names(new.args)] <- new.args
    do.call("points", plot.args)
    return()
  }
  plot.args <- list(x = data$pos, y = data$logP,
                    pch = data$pch, bg = data$bg, col = col,
                    las = 1, font.main = 1,
                    xlim = loc$xrange,
                    ylim = c(0, max(data$logP, na.rm = TRUE)),
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
  if (!is.null(legend_pos)) {
    legtext <- c(rev(paste(cutlev[,1], cutlev[,2], sep=" : ")),
                 paste(-cutlev[,2], -cutlev[,1], sep=" : "))
    legend(legend_pos, legend = legtext, pch = rep(c(24, 25), each=3),
           col = col, pt.bg = c(rev(up_cols), down_cols),
           title = "eQTL effect",
           bty = "n", cex = 0.85, pt.cex = 1, y.intersp = 0.96)
  }
}
