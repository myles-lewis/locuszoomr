#' Locus plot
#' 
#' Genomic locus plot similar to locuszoom.
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
#' @param pcutoff Cut-off for p value significance. Defaults to 5e-08.
#' @param chromCols Colour for normal points if `LD` is `FALSE`.
#' @param sigCol Colour for significant points if `LD` is `FALSE`.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param cex.axis Specifies font size for axis numbering.
#' @param heights Ratio of top to bottom plot. See [layout].
#' @param maxrows Specifies maximum nunber of rows to display in gene 
#' annotation panel.
#' @param xticks Character value specifying whether x axis ticks and numbers.
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
#' @param border Logical whether a bounding box is plotted around upper and 
#' lower plots.
#' @param LDcols Vector of colours for plotting LD. The first colour is for SNPs 
#' which lack LD information. The next 5 colours are for r2 or D' LD results 
#' ranging from 0 to 1 in intervals of 0.2. The final colour is for the index 
#' SNP.
#' @param ... Other arguments passed to [plot()] for the scatter plot.
#' @return Returns a list containing the subset of data plotted, chromosome and
#' genomic position range.
#' @importFrom ensembldb genes exons
#' @importFrom BiocGenerics start end
#' @importFrom LDlinkR LDmatrix
#' @importFrom AnnotationFilter GeneNameFilter AnnotationFilterList 
#' SeqNameFilter GeneIdFilter
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom graphics axTicks axis layout lines par rect strwidth text
#' @importFrom memoise memoise
#' @export

locusplot <- function(data, xrange = NULL, seqname = NULL,
                      gene = NULL, flank = 5e4,
                      ens_version = "EnsDb.Hsapiens.v75",
                      chrom = 'chrom', pos = 'pos', p = 'p',
                      labs = 'rsid',
                      pcutoff = 5e-08,
                      chromCols = 'royalblue',
                      sigCol = 'red',
                      xlab = NULL, ylab = '-log10 P value',
                      cex.axis = 0.8,
                      heights = c(3, 2),
                      maxrows = 7,
                      xticks = 'bottom',
                      index_snp = NULL,
                      LD = TRUE,
                      pop = "CEU",
                      r2d = "r2",
                      LDtoken = "",
                      border = FALSE,
                      LDcols = c('grey', 'blue', 'cyan', 'green3', 'orange', 'red', 
                                 'purple'),
                      ...) {
  args <- list(...)
  require(ens_version, character.only = TRUE)
  edb <- get(ens_version)
  if (!is.null(gene)) {
    locus <- genes(edb, filter = GeneNameFilter(gene))
    xrange <- c(start(locus) - flank, end(locus) + flank)
    seqname <- names(seqlengths(locus))
  }
  if (is.null(xrange) | is.null(seqname)) stop('No locus specified')
  if (is.null(xlab)) xlab <- paste("Chromosome", seqname, "(Mb)")
  data <- data[data[, chrom] == seqname &
                 data[, pos] > xrange[1] & data[, pos] < xrange[2], ]
  data$logP <- -log10(data[, p])
  data <- as.data.frame(data)
  if (LD) {
    rslist <- data[, labs]
    if (length(rslist) > 1000) {
      rslist <- rslist[order(data$logP, decreasing = TRUE)[1:1000]]
    }
    cat(paste("Obtaining LD on", length(rslist), "SNPs"))
    ldm <- mem_LDmatrix(rslist, pop = pop, r2d = r2d, token = LDtoken)
    if (is.null(index_snp)) index_snp <- data[which.max(data$logP), labs]
    ld <- ldm[, index_snp]
    data$ld <- ld[match(data[, labs], ldm$RS_number)]
    data$col <- LDcols[cut(data$ld, -1:6/5, labels = FALSE)]
    data$col[is.na(data$col)] <- LDcols[1]
    data$col[which.max(data$logP)] <- LDcols[7]
  } else {
    data$col <- chromCols
    data$col[data[, p] < pcutoff] <- sigCol
  }
  
  TX <- genes(edb, filter = AnnotationFilterList(
    SeqNameFilter(c(1:22, "X", "Y")),
    GeneIdFilter("ENSG", "startsWith")))
  TX <- data.frame(TX)
  TX <- TX[! is.na(TX$start), ]
  TX <- TX[TX$seqnames == seqname, ]
  TX <- TX[TX$end > xrange[1], ]
  TX <- TX[TX$start < xrange[2], ]
  TX <- mapRow(TX)
  maxrows <- if (is.null(maxrows)) max(TX$row) else min(c(max(TX$row), maxrows))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  layout(matrix(2:1, nrow = 2), heights = heights)
  
  # lower locus plot
  par(tcl = -0.3, las = 1, font.main = 1,
          mgp = c(1.8, 0.5, 0), mar = c(4, 4, 1, 2))
  plot(NA, xlim = xrange,
       ylim = c(-maxrows - 0.3, -0.3), 
       bty = if (border) 'o' else 'n',
       yaxt = 'n', xaxt = 'n',
       xlab = if (xticks == 'bottom') xlab else "",
       ylab = "",
       xaxt = 'n')
  if (xticks == 'bottom') {
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis)
  }
  for (i in 1:nrow(TX)) {
    lines(TX[i, c('start', 'end')], rep(-TX[i, 'row'], 2), lwd = 2)
    e <- exons(edb, filter = GeneIdFilter(TX$gene_id[i]))
    exstart <- start(e)
    exend <- end(e)
    for (j in 1:length(e)) {
      rect(exstart[j], -TX[i, 'row'] - 0.1, exend[j], -TX[i, 'row'] + 0.1,
           col = 'black', border = NA)
    }
  }
  
  tfilter <- TX$tmin > xrange[1] & TX$tmax < xrange[2]
  for (i in which(tfilter)) {
    text(TX$mean[i], -TX[i, 'row'] + 0.4,
         labels = if (TX$strand[i] == "+") {
           bquote(.(TX$gene_name[i]) ~ symbol("\256"))
         } else {     
           bquote(symbol("\254") ~ .(TX$gene_name[i]))
         }, cex = 0.7)
  }
  
  # scatter plot
  par(mar = c(ifelse(xticks == 'top', 3, 0), 4, 2, 2))
  plot(data[, pos], data$logP,
       pch = 21, bg = data$col,
       xlim = xrange,
       xlab = if (xticks == 'top') xlab else "",
       ylab = ylab,
       bty = if (border) 'o' else 'l',
       cex.axis = cex.axis,
       xaxt = 'n', ...)
  if (xticks == 'top') {
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis)
  }
  invisible(list(data = data, seqname = seqname, xrange = xrange))
}

# map genes into rows without overlap
mapRow <- function(TX, gap = 2e3, cex.text = 0.7, 
                   xlim = range(TX[, c('start', 'end')])) {
  gw <- strwidth(TX$gene_name, units = "inch", cex = cex.text) * diff(xlim) / 12
  TX$mean <- rowMeans(TX[, c('start', 'end')])
  TX$tmin <- TX$mean - gw
  TX$tmax <- TX$mean + gw
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
