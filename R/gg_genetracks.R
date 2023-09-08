

#' @importFrom grid 
#' @export

gg_genetracks <- function(locus,
                          filter_gene_name = NULL,
                          filter_gene_biotype = NULL,
                          border = FALSE,
                          cex.axis = 1,
                          cex.lab = 1,
                          cex.text = 0.7,
                          gene_col = ifelse(showExons, 'blue4', 'skyblue'),
                          exon_col = 'blue4',
                          exon_border = 'blue4',
                          showExons = TRUE,
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
    return()
  }
  if (is.null(xlab)) xlab <- paste("Chromosome", locus$seqname, "(Mb)")
  
  pos <- TX$strand == "+"
  TX$gene_name[pos] <- paste0(TX$gene_name[pos], sprintf("\u2192"))
  TX$gene_name[!pos] <- paste0(sprintf("\u2190"), TX$gene_name[!pos])
  
  TX <- mapRow(TX, xlim = xrange, cex.text = cex.text, text_pos = text_pos)
  maxrows <- if (is.null(maxrows)) max(TX$row) else min(c(max(TX$row), maxrows))
  if (max(TX$row) > maxrows) message(max(TX$row), " tracks needed to show all genes")
  TX <- TX[TX$row <= maxrows, ]
  
  ylim <- c(-maxrows - 0.3, -0.3)
  xrange <- xrange / 1e6
  TX$start <- TX$start / 1e6
  TX$end <- TX$end / 1e6
  TX[, c("mean", "tmin", "tmax", "min", "max")] <- TX[, c("mean", "tmin", "tmax", "min", "max")] / 1e6
  
  exheight <- switch(text_pos, "top" = 0.15, "left" = 0.3)
  
  vp <- viewport(x = 0,
                 y = if (xticks) unit(4, "lines") else 0,
                 width = unit(1, "npc"),
                 height = unit(1, "npc") - unit(5, "lines"),
                 just = c("left", "bottom"),
                 xscale = xrange + c(-0.04, 0.04) * diff(xrange),
                 yscale = ylim
                 )
  
  pushViewport(vp)
  if (showExons) {
    for (i in seq_len(nrow(TX))) {
      grid.lines(unit(TX[i, c('start', 'end')], "native"),
                 unit(rep(-TX[i, 'row'], 2), "native"),
                 gp = gpar(col = gene_col, lwd = 1.5, lineend = "butt"))
      e <- EX[EX$gene_id == TX$gene_id[i], ]
      exstart <- start(e) / 1e6
      exwidth <- end(e) / 1e6 - exstart
      grid.rect(x = unit(exstart, "native"),
                y = unit(-TX[i, 'row'] - exheight, "native"),
                width = unit(exwidth, "native"),
                height = unit(2 * exheight, "native"),
                just = c("left", "bottom"),
                gp = gpar(fill = exon_col, col = exon_border,
                          lwd = 0.5, lineend = "square", linejoin = "mitre"))
    }
  } else {
    grid.rect(x = unit(TX$start, "native"),
              y = unit(-TX[, 'row'] - exheight, "native"),
              width = unit(TX$end - TX$start, "native"),
              height = unit(exheight*2, "native"),
              just = c("left", "bottom"),
              gp = gpar(fill = gene_col, col = exon_border,
                        lineend = "square", linejoin = "mitre"))
  }
  
  if (text_pos == "top") {
    tfilter <- which(TX$tmin > (xrange[1] - diff(xrange) * 0.04) & 
                       (TX$tmax < xrange[2] + diff(xrange) * 0.04))
    grid.text(TX$gene_name[tfilter],
              x = unit(TX$mean[tfilter], "native"),
              y = unit(-TX$row[tfilter] + 0.45, "native"),
              gp = gpar(cex = cex.text))
  } else if (text_pos == "left") {
    tfilter <- if (border) {
      which(TX$tmin > xrange[1])
    } else seq_len(nrow(TX))
    grid.text(TX$gene_name[tfilter],
              x = unit(pmax(TX$start[tfilter],
                            xrange[1] - diff(xrange) * 0.04) - diff(xrange) * 0.01,
                       "native"),
              y = unit(-TX$row[tfilter], "native"),
              just = "right",
              gp = gpar(cex = cex.text))
  }
  
  if (border) grid.rect()
  if (xticks) {
    grid.text(xlab, y = unit(-3, "lines"), gp = gpar(fontsize = 12))
    grid.xaxis()
  }
  
  popViewport()
}
