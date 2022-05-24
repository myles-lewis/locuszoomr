

#' @import EnsDb.Hsapiens.v75
#' @importFrom graphics axTicks axis layout par strwidth abline legend
#' @importFrom grDevices adjustcolor hcl.colors
#' @export

eQTLplot <- function(x, ...,
                     filter_gene_name = NULL,
                     filter_gene_biotype = NULL,
                     pcutoff = 5e-08,
                     col = 'black',
                     up_palette = "OrRd",
                     down_palette = "Blues 3",
                     alpha = 0.5,
                     tissue = "Whole Blood",
                     eqtl_gene = x$gene,
                     xlab = NULL, ylab = expression("-log"[10] ~ "P"),
                     cex.axis = 0.8,
                     cex.text = 0.7,
                     use_layout = TRUE,
                     heights = c(3, 2),
                     maxrows = 7,
                     xticks = 'bottom',
                     border = FALSE,
                     gene_col = 'blue4',
                     exon_col = 'blue4',
                     exon_border = 'blue4',
                     text_pos = 'top',
                     legend_pos = 'topleft') {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  if (!"LDexp" %in% names(x)) stop("Missing eQTL data")
  if (is.null(eqtl_gene)) stop("eqtl_gene not specified")
  data <- x$data
  TX <- x$TX
  EX <- x$EX
  if (is.null(xlab)) xlab <- paste("Chromosome", x$seqname, "(Mb)")
  LDX <- x$LDexp[x$LDexp$Tissue == tissue & x$LDexp$Gene_Symbol == eqtl_gene, ]
  ind <- match(data[, x$labs], LDX$RS_ID)
  data$eqtl_effect <- NA
  data$eqtl_effect <- LDX$Effect_Size[ind]
  data$eqtl_p <- LDX$P_value[ind]
  data$eqtl_effect_allele <- LDX$Effect_Allele[ind]
  # gwas allele and eqtl effect allele are the same
  mismatch <- which(data$eqtl_effect_allele != data$effect_allele)
  which_rev <- data$other_allele[mismatch] == data$eqtl_effect_allele[mismatch]
  rev_effect <- mismatch[which_rev]
  mismatch <- mismatch[!which_rev]
  data$eqtl_effect[rev_effect] <- -data$eqtl_effect[rev_effect]
  data$eqtl_effect[mismatch] <- NA
  up_cols <- hcl.colors(8, up_palette, rev = TRUE)[-c(1,2)]
  down_cols <- hcl.colors(8, down_palette, rev = TRUE)[-c(1,2)]
  ecol <- cut(abs(data$eqtl_effect), breaks=6)
  data$col <- 'black'
  data$col <- adjustcolor(data$col, alpha.f = alpha)
  eqind <- !is.na(data$eqtl_effect)
  equp <- eqind & sign(data$eqtl_effect) == -1
  eqdown <- eqind & sign(data$eqtl_effect) == 1
  data$col[equp] <- down_cols[ecol[equp]]
  data$col[eqdown] <- up_cols[ecol[eqdown]]
  data$pch <- 21
  data$pch[eqind] <- 24.5 - sign(data$eqtl_effect[eqind]) / 2
  
  if (use_layout) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    graphics::layout(matrix(2:1, nrow = 2), heights = heights)
  }
  # lower locus plot
  par(tcl = -0.3, las = 1, font.main = 1,
      mgp = c(1.8, 0.5, 0), 
      mar = c(ifelse(xticks == 'bottom', 4, 2), 4, 0.25, 2))
  genetracks(x, filter_gene_name, filter_gene_biotype,
             border, cex.axis, cex.text, gene_col, exon_col, exon_border,
             maxrows, text_pos, xticks = (xticks == 'bottom'),
             xlab = if (xticks == 'bottom') xlab else "")
  
  # scatter plot
  par(mar = c(ifelse(xticks == 'top', 3, 0), 4, 2, 2))
  plot(data[, x$pos], data$logP,
       pch = data$pch, col = data$col, bg = data$col,
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
  # if (!is.null(legend_pos)) {
  #     legend(legend_pos,
  #            legend = c('Index SNP',
  #                       expression({0.8 < r^2} <= "1.0"),
  #                       expression({0.6 < r^2} <= 0.8),
  #                       expression({0.4 < r^2} <= 0.6),
  #                       expression({0.2 < r^2} <= 0.4),
  #                       expression({"0.0" < r^2} <= 0.2),
  #                       expression("No" ~ r^2 ~ "data")),
  #            pch = 21, col = 'black', pt.bg = rev(LDcols), bty = 'n', cex = 0.7)
  # }
}
