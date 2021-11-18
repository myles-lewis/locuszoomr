#' Create locus object for plotting
#' 
#' Creates object of class 'locus' for genomic locus plot similar to locuszoom.
#' 
#' @details 
#' This is an R version of locuszoom for generating publication ready Manhattan 
#' plots of gene loci. It references ensembl databases for annotating genes 
#' and exons. It queries LDlink to retrieve LD information on the index SNP.
#' 
#' @param data Dataset (data.frame or data.table) to use for plot.
#' @param xrange Vector of genomic position range for the x axis.
#' @param seqname Specifies which chromosome to plot.
#' @param gene Optional specifies which gene to view. Either `xrange` with 
#' `seqname`, or `gene` must be specified.
#' @param flank Vector of how much flanking region left and right of the gene 
#' to show. 
#' @param ens_version Specifies which ensembl database to query for gene and 
#' exon positions. See `ensembldb` Bioconductor package.
#' @param chrom Determines which column in `data` contains chromosome 
#' information.
#' @param pos Determines which column in `data` contains position information.
#' @param p Determines which column in `data` contains SNP p values.
#' @param labs Determines which column in `data` contains SNP rs IDs.
#' @param index_snp Specifies the index SNP for displaying linkage 
#' disequilibrium (LD). If not specifiied, the SNP with the lowest P value is 
#' selected.
#' @param LD Logical whether LD is plotted. Queries 1000 genomes via LDlinkR 
#' package. See [LDlinkR]. Results are cached using the `memoise` package, so
#' that if exactly the same locus is requested the system does not repeatedly 
#' call the API.
#' @param pop A 1000 Genomes Project population, (e.g. YRI or CEU), multiple 
#' allowed, default = "CEU". Passed to [LDlinkR::LDmatrix].
#' @param r2d Either "r2" for LD r^2 or "d" for LD D', default = "r2". Passed 
#' to [LDlinkR::LDmatrix].
#' @param LDtoken Peronsal access token for accessing 1000 genomes LD data via 
#' LDlink API. See [LDlinkR].
#' @return Returns an object of class 'locus' ready for plotting, containing 
#' the subset of GWAS data to be plotted, 
#' chromosome and genomic position range, 
#' ensembl database version number, 
#' column names for chromosome, position, SNP ID, p value
#' locus gene information from ensembl and
#' locus exon information from ensembl.
#' @importFrom ensembldb genes exons
#' @importFrom BiocGenerics start end
#' @importFrom LDlinkR LDmatrix
#' @importFrom AnnotationFilter GeneNameFilter AnnotationFilterList 
#' SeqNameFilter GeneIdFilter
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom memoise memoise
#' @export

locus <- function(data, xrange = NULL, seqname = NULL,
                      gene = NULL, flank = 5e4,
                      ens_version = "EnsDb.Hsapiens.v75",
                      chrom = 'chrom', pos = 'pos', p = 'p',
                      labs = 'rsid',
                      index_snp = NULL,
                      LD = TRUE,
                      pop = "CEU",
                      r2d = "r2",
                      LDtoken = "") {
  require(ens_version, character.only = TRUE)
  edb <- get(ens_version)
  if (!is.null(gene)) {
    locus <- genes(edb, filter = GeneNameFilter(gene))
    xrange <- c(start(locus) - flank, end(locus) + flank)
    seqname <- names(seqlengths(locus))
  }
  if (is.null(xrange) | is.null(seqname)) stop('No locus specified')
  data <- data[data[, chrom] == seqname &
                 data[, pos] > xrange[1] & data[, pos] < xrange[2], ]
  data$logP <- -log10(data[, p])
  data <- as.data.frame(data)
  if (LD) {
    if (LDtoken == "") stop("LDtoken is missing")
    rslist <- data[, labs]
    if (length(rslist) > 1000) {
      rslist <- rslist[order(data$logP, decreasing = TRUE)[1:1000]]
    }
    cat(paste("Obtaining LD on", length(rslist), "SNPs"))
    ldm <- mem_LDmatrix(rslist, pop = pop, r2d = r2d, token = LDtoken)
    if (is.null(index_snp)) index_snp <- data[which.max(data$logP), labs]
    ld <- ldm[, index_snp]
    data$ld <- ld[match(data[, labs], ldm$RS_number)]
  }
  
  TX <- ensembldb::genes(edb, filter = AnnotationFilterList(
    SeqNameFilter(c(1:22, "X", "Y")),
    GeneIdFilter("ENSG", "startsWith")))
  TX <- data.frame(TX)
  TX <- TX[! is.na(TX$start), ]
  TX <- TX[TX$seqnames == seqname, ]
  TX <- TX[TX$end > xrange[1], ]
  TX <- TX[TX$start < xrange[2], ]
  EX <- ensembldb::exons(edb, filter = GeneIdFilter(TX$gene_id))
  
  loc <- list(seqname = seqname, xrange = xrange,
              ens_version = ens_version,
              chrom = chrom, pos = pos, p = p, labs = labs,
              data = data, TX = TX, EX = EX)
  class(loc) <- "locus"
  loc
}

#' Locus plot
#' 
#' Genomic locus plot similar to locuszoom.
#' 
#' @details 
#' This is an R version of locuszoom for generating publication ready Manhattan 
#' plots of gene loci. It references ensembl databases for annotating genes 
#' and exons. It queries LDlink to retrieve LD information on the index SNP.
#' 
#' @param x Object of class 'locus' to use for plot. See [locus].
#' @param pcutoff Cut-off for p value significance. Defaults to 5e-08.
#' @param chromCols Colour for normal points if `LD` is `FALSE` when the locus 
#' object is made.
#' @param sigCol Colour for significant points if `LD` is `FALSE`.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.text Font size for gene text.
#' @param heights Ratio of top to bottom plot. See [layout].
#' @param maxrows Specifies maximum nunber of rows to display in gene 
#' annotation panel.
#' @param xticks Character value of either 'top' or 'bottom' specifying 
#' whether x axis ticks and numbers are plotted on top or bottom plot window.
#' @param border Logical whether a bounding box is plotted around upper and 
#' lower plots.
#' @param LDcols Vector of colours for plotting LD. The first colour is for SNPs 
#' which lack LD information. The next 5 colours are for r2 or D' LD results 
#' ranging from 0 to 1 in intervals of 0.2. The final colour is for the index 
#' SNP.
#' @param genecol Colour for genes and exons.
#' @param legend_pos Position of legend. See [legend()]. Set to `NULL` to hide 
#' legend.
#' @param ... Other arguments passed to [plot()] for the scatter plot.
#' @return No return value.
#' @importFrom graphics axTicks axis layout par strwidth abline legend
#' @export

plot.locus <- function(x, ...,
                      pcutoff = 5e-08,
                      chromCols = 'royalblue',
                      sigCol = 'red',
                      xlab = NULL, ylab = expression("-log"[10] ~ "P"),
                      cex.axis = 0.8,
                      cex.text = 0.7,
                      heights = c(3, 2),
                      maxrows = 7,
                      xticks = 'bottom',
                      border = FALSE,
                      LDcols = c('grey', 'royalblue', 'cyan2', 'green3', 'orange', 'red', 
                                 'purple'),
                      genecol = 'blue4',
                      legend_pos = 'topleft') {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  data <- x$data
  TX <- x$TX
  EX <- x$EX
  if (is.null(xlab)) xlab <- paste("Chromosome", x$seqname, "(Mb)")
  LD <- "ld" %in% colnames(data)
  if (LD) {
    data$col <- LDcols[cut(data$ld, -1:6/5, labels = FALSE)]
    data$col[is.na(data$col)] <- LDcols[1]
    data$col[which.max(data$logP)] <- LDcols[7]
  } else {
    data$col <- chromCols
    data$col[data[, x$p] < pcutoff] <- sigCol
  }
  # TX <- mapRow(TX, xlim = x$xrange, cex.text = cex.text)
  # maxrows <- if (is.null(maxrows)) max(TX$row) else min(c(max(TX$row), maxrows))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  layout(matrix(2:1, nrow = 2), heights = heights)
  
  # lower locus plot
  par(tcl = -0.3, las = 1, font.main = 1,
      mgp = c(1.8, 0.5, 0), mar = c(ifelse(xticks == 'bottom', 4, 2), 4, 1, 2))
  genetracks(x, border, cex.axis, cex.text, genecol, 
             maxrows, xticks = (xticks == 'bottom'))
  
  # scatter plot
  par(mar = c(ifelse(xticks == 'top', 3, 0), 4, 2, 2))
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
#' Plot gene annotation tracks from ensembldb data.
#' 
#' @details This function is called by plot.locus(). It can used to plot the 
#' gene annotation tracks on their own. It uses base graphics, so layout() can 
#' be used to position adjacent plots above or below.
#' @param locus Object of class 'locus' generated by [locus()].
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.text Font size for gene text.
#' @param maxrows Specifies maximum nunber of rows to display in gene 
#' annotation panel.
#' @param xticks Logical whether x axis ticks and numbers are plotted.
#' @param xlab Title for x axis. Defaults to chromosome `seqname` specified 
#' in `locus`.
#' @param border Logical whether a bounding box is plotted.
#' @param genecol Colour for genes and exons.
#' @return No return value.
#' @importFrom BiocGenerics start end
#' @importFrom graphics axTicks axis lines rect text
#' @export

genetracks <- function(locus,
                       border = FALSE,
                       cex.axis = 0.8,
                       cex.text = 0.7,
                       genecol = 'blue4', 
                       maxrows = 7,
                       xticks = TRUE,
                       xlab = NULL) {
  if (!inherits(locus, "locus")) stop("Object of class 'locus' required")
  TX <- locus$TX
  EX <- locus$EX
  xrange <- locus$xrange
  TX <- mapRow(TX, xlim = xrange, cex.text = cex.text)
  maxrows <- if (is.null(maxrows)) max(TX$row) else min(c(max(TX$row), maxrows))
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
  for (i in 1:nrow(TX)) {
    lines(TX[i, c('start', 'end')], rep(-TX[i, 'row'], 2), lwd = 1, lend = 1)
    e <- EX[EX$gene_id == TX$gene_id[i], ]
    exstart <- start(e)
    exend <- end(e)
    rect(exstart, -TX[i, 'row'] - 0.15, exend, -TX[i, 'row'] + 0.15,
         col = genecol, border = genecol, lwd = 0.5, lend = 2, ljoin = 1)
  }
  tfilter <- TX$tmin > xrange[1] & TX$tmax < xrange[2]
  for (i in which(tfilter)) {
    text(TX$mean[i], -TX[i, 'row'] + 0.45,
         labels = if (TX$strand[i] == "+") {
           bquote(.(TX$gene_name[i]) ~ symbol("\256"))
         } else {     
           bquote(symbol("\254") ~ .(TX$gene_name[i]))
         }, cex = cex.text)
  }
}

# map genes into rows without overlap
mapRow <- function(TX, gap = 2e3, cex.text = 0.7, 
                   xlim = range(TX[, c('start', 'end')])) {
  gw <- strwidth(paste("--", TX$gene_name), units = "inch", cex = cex.text) * diff(xlim) / par("pin")[1]
  TX$mean <- rowMeans(TX[, c('start', 'end')])
  TX$tmin <- TX$mean - gw / 2
  TX$tmax <- TX$mean + gw / 2
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
