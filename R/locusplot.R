
#' Locus plot
#' 
#' Genomic locus plot similar to locuszoom.
#' 
#' @details 
#' This is an R version of locuszoom for generating publication ready Manhattan 
#' plots of gene loci. It references Ensembl databases for annotating genes 
#' and exons. It queries LDlink to retrieve linkage disequilibrium (LD) 
#' information on the index SNP. Use [locus()] first to generate an object of 
#' class 'locus' for plotting.
#' 
#' @param x Object of class 'locus' to use for plot. See [locus].
#' @param filter_gene_name Vector of gene names to display.
#' @param filter_gene_biotype Vector of gene biotypes to be filtered. Use
#' [ensembldb::listGenebiotypes()] to display possible biotypes. For example, 
#' `ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)`
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set 
#' to `NULL` to disable.
#' @param chromCols Colour for normal points if `LD` is `FALSE` when the locus 
#' object is made.
#' @param sigCol Colour for significant points if `LD` is `FALSE`.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.text Font size for gene text.
#' @param use_layout Logical whether `graphics::layout` is called. Default
#'   `TRUE` is for a standard single plot. Set to `FALSE` if a more complex
#'   layout with multiple plots is required.
#' @param heights Ratio of top to bottom plot. See [layout].
#' @param maxrows Specifies maximum number of rows to display in gene 
#' annotation panel.
#' @param xticks Character value of either 'top' or 'bottom' specifying 
#' whether x axis ticks and numbers are plotted on top or bottom plot window.
#' @param border Logical whether a bounding box is plotted around upper and 
#' lower plots.
#' @param LDcols Vector of colours for plotting LD. The first colour is for SNPs
#'   which lack LD information. The next 5 colours are for r2 or D' LD results
#'   ranging from 0 to 1 in intervals of 0.2. The final colour is for the index
#'   SNP.
#' @param gene_col Colour for genes and exons.
#' @param exon_col Fill colour for exons.
#' @param exon_border Border line colour outlining exons. Set to `NA` for no 
#' border.
#' @param text_pos Character value of either 'top' or 'left' specifying 
#' placement of gene name labels.
#' @param legend_pos Position of legend. See [legend()]. Set to `NULL` to hide 
#' legend.
#' @param ... Other arguments passed to [plot()] for the scatter plot.
#' @return No return value.
#' @examples
#' library(EnsDb.Hsapiens.v75)
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5, LD = FALSE)
#' plot(loc)
#' loc2 <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5, LD = FALSE)
#' plot(loc2)
#' @import EnsDb.Hsapiens.v75
#' @importFrom graphics axTicks axis layout par strwidth abline legend
#' @export

plot.locus <- function(x, ...,
                       filter_gene_name = NULL,
                       filter_gene_biotype = NULL,
                       pcutoff = 5e-08,
                       chromCols = 'royalblue',
                       sigCol = 'red',
                       xlab = NULL, ylab = expression("-log"[10] ~ "P"),
                       cex.axis = 0.8,
                       cex.text = 0.7,
                       use_layout = TRUE,
                       heights = c(3, 2),
                       maxrows = 7,
                       xticks = 'bottom',
                       border = FALSE,
                       LDcols = c('grey', 'royalblue', 'cyan2', 'green3', 
                                  'orange', 'red', 'purple'),
                       gene_col = 'blue4',
                       exon_col = 'blue4',
                       exon_border = 'blue4',
                       text_pos = 'top',
                       legend_pos = 'topleft') {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
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
  
  if (use_layout) {
    op0 <- par(no.readonly = TRUE)
    on.exit(par(op0), add = TRUE)
    graphics::layout(matrix(2:1, nrow = 2), heights = heights)
  }
  # lower locus plot
  genetracks(x, filter_gene_name, filter_gene_biotype,
             border, cex.axis, cex.text, gene_col, exon_col, exon_border,
             maxrows, text_pos, xticks = (xticks == 'bottom'),
             xlab = if (xticks == 'bottom') xlab else "")
  
  # scatter plot
  op <- par(mar = c(ifelse(xticks == 'top', 3, 0.1), 4, 2, 1.5),
      mgp = c(1.7, 0.5, 0))
  on.exit(par(op), add = TRUE)
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
  if (xticks == 'top') {
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
