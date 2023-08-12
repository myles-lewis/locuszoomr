
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
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    graphics::layout(matrix(2:1, nrow = 2), heights = heights)
  }
  # lower locus plot
  par(tcl = -0.3, las = 1, font.main = 1,
      mgp = c(1.8, 0.5, 0), 
      mar = c(ifelse(xticks == 'bottom', 4, 2), 4, 0.25, 1.5))
  genetracks(x, filter_gene_name, filter_gene_biotype,
             border, cex.axis, cex.text, gene_col, exon_col, exon_border,
             maxrows, text_pos, xticks = (xticks == 'bottom'),
             xlab = if (xticks == 'bottom') xlab else "")
  
  # scatter plot
  par(mar = c(ifelse(xticks == 'top', 3, 0), 4, 2, 1.5))
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
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis)
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

#' Plot gene tracks
#' 
#' Plot gene annotation tracks from `ensembldb` data.
#' 
#' @details This function is called by plot.locus(). It can used to plot the 
#' gene annotation tracks on their own. It uses base graphics, so layout() can 
#' be used to position adjacent plots above or below.
#' @param locus Object of class 'locus' generated by [locus()].
#' @param filter_gene_name Vector of gene names to display.
#' @param filter_gene_biotype Vector of gene biotypes to be filtered. Use
#' [ensembldb::listGenebiotypes()] to display possible biotypes. For example, 
#' `ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)`
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.text Font size for gene text.
#' @param maxrows Specifies maximum number of rows to display in gene 
#' annotation panel.
#' @param xticks Logical whether x axis ticks and numbers are plotted.
#' @param xlab Title for x axis. Defaults to chromosome `seqname` specified 
#' in `locus`.
#' @param border Logical whether a bounding box is plotted.
#' @param gene_col Colour for gene lines.
#' @param exon_col Fill colour for exons.
#' @param exon_border Border line colour outlining exons. Set to `NA` for no 
#' border.
#' @param text_pos Character value of either 'top' or 'left' specifying 
#' placement of gene name labels.
#' @return No return value.
#' @examples
#' library(EnsDb.Hsapiens.v75)
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5, LD = FALSE)
#' genetracks(loc)
#' 
#' ## Limit the number of tracks
#' genetracks(loc, maxrows = 4)
#' 
#' ## Filter by gene biotype
#' genetracks(loc, filter_gene_biotype = 'protein_coding')
#' 
#' ## Customise colours
#' genetracks(loc, gene_col = 'grey', exon_col = 'orange',
#'            exon_border = 'darkgrey')
#' 
#' @import EnsDb.Hsapiens.v75
#' @importFrom BiocGenerics start end
#' @importFrom graphics axTicks axis lines rect text plot.new
#' @export

genetracks <- function(locus,
                       filter_gene_name = NULL,
                       filter_gene_biotype = NULL,
                       border = FALSE,
                       cex.axis = 0.8,
                       cex.text = 0.7,
                       gene_col = 'blue4',
                       exon_col = 'blue4',
                       exon_border = 'blue4',
                       maxrows = NULL,
                       text_pos = 'top',
                       xticks = TRUE,
                       xlab = NULL) {
  if (!inherits(locus, "locus")) stop("Object of class 'locus' required")
  TX <- locus$TX
  EX <- locus$EX
  xrange <- locus$xrange
  if (!is.null(filter_gene_name)) {
    TX <- TX[TX$gene_name %in% filter_gene_name, ]
  }
  if (!is.null(filter_gene_biotype)) {
    TX <- TX[TX$gene_biotype %in% filter_gene_biotype, ]
  }
  if (nrow(TX) == 0) {
    message('No genes to plot')
    return(plot.new())
  }
  TX <- mapRow(TX, xlim = xrange, cex.text = cex.text, text_pos = text_pos)
  maxrows <- if (is.null(maxrows)) max(TX$row) else min(c(max(TX$row), maxrows))
  TX <- TX[TX$row <= maxrows, ]
  if (is.null(xlab)) xlab <- paste("Chromosome", locus$seqname, "(Mb)")
  
  plot(NA, xlim = xrange,
       ylim = c(-maxrows - 0.3, -0.3), 
       bty = if (border) 'o' else 'n',
       yaxt = 'n', xaxt = 'n',
       xlab = if (xticks) xlab else "",
       ylab = "",
       xaxt = 'n')
  if (xticks) {
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis)
  }
  exheight <- switch(text_pos, "top" = 0.15, "left" = 0.3)
  for (i in seq_len(nrow(TX))) {
    lines(TX[i, c('start', 'end')], rep(-TX[i, 'row'], 2),
          col = gene_col, lwd = 1, lend = 1)
    e <- EX[EX$gene_id == TX$gene_id[i], ]
    exstart <- start(e)
    exend <- end(e)
    rect(exstart, -TX[i, 'row'] - exheight, exend, -TX[i, 'row'] + exheight,
         col = exon_col, border = exon_border, lwd = 0.5, lend = 2, ljoin = 1)
  }
  if (text_pos == "top") {
    tfilter <- which(TX$tmin > (xrange[1] - diff(xrange) * 0.04) & 
              (TX$tmax < xrange[2] + diff(xrange) * 0.04))
    for (i in tfilter) {
      text(TX$mean[i], -TX[i, 'row'] + 0.45,
           labels = if (TX$strand[i] == "+") {
             bquote(.(TX$gene_name[i]) * symbol("\256"))
           } else {     
             bquote(symbol("\254") * .(TX$gene_name[i]))
           }, cex = cex.text, xpd = NA)
    }
  } else if (text_pos == "left") {
    tfilter <- if (border) {
      which(TX$tmin > xrange[1])
    } else seq_len(nrow(TX))
    for (i in tfilter) {
      text(max(c(TX$start[i], xrange[1] - diff(xrange) * 0.04)), -TX[i, 'row'],
           labels = if (TX$strand[i] == "+") {
             bquote(.(TX$gene_name[i]) * symbol("\256"))
           } else {     
             bquote(symbol("\254") * .(TX$gene_name[i]))
           }, cex = cex.text, pos = 2, xpd = NA)
    }
  }
  
}


# map genes into rows without overlap
mapRow <- function(TX, gap = diff(xlim) * 0.02, cex.text = 0.7, 
                   xlim = range(TX[, c('start', 'end')]),
                   text_pos = 'top') {
  gw <- strwidth(paste0("--", TX$gene_name), units = "inch", 
                 cex = cex.text) * diff(xlim) / par("pin")[1]
  TX$mean <- rowMeans(TX[, c('start', 'end')])
  if (text_pos == 'top') {
    TX$tmin <- TX$mean - gw / 2
    TX$tmax <- TX$mean + gw / 2
  } else if (text_pos == 'left') {
    TX$tmin <- TX$start - gw - gap
    TX$tmax <- TX$end
  }
  TX$min <- apply(TX[, c('start', 'end', 'tmin')], 1, min) - gap / 2
  TX$max <- apply(TX[, c('start', 'end', 'tmax')], 1, max) + gap / 2
  TX$row <- 0
  j <- 1
  while (any(TX$row == 0)) {
    xset <- which(TX$row == 0)
    for (i in xset) {
      # overlap detection
      if (!any(TX$min[i] < TX$max[TX$row == j] &
               TX$max[i] > TX$min[TX$row == j])) {
        TX$row[i] <- j
      }
    }
    j <- j + 1
  }
  TX
}

# use memoise to reduce calls to LDlink API
mem_LDmatrix <- memoise(LDmatrix)
mem_LDexpress <- memoise(LDexpress)

#' @importFrom grDevices col2rgb rgb

col2hex <- function(cname) {
  colMat <- col2rgb(cname)
  rgb(red = colMat[1, ]/255, green = colMat[2, ]/255, blue = colMat[3, ]/255)
}

