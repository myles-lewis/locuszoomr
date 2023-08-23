
#' Plot overlaying eQTL and GWAS data
#'
#' Experimental plotting function for overlaying eQTL data from GTEx on top of
#' GWAS results. y axis shows the -log10 p-value for the GWAS result.
#' Significant eQTL for the specified gene are overlaid using colours and
#' symbols.
#' 
#' @param x Object of class 'locus' to use for plot. See [locus].
#' @param filter_gene_name Vector of gene names to display.
#' @param filter_gene_biotype Vector of gene biotypes to be filtered. Use
#' [ensembldb::listGenebiotypes()] to display possible biotypes. For example, 
#' `ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)`
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set 
#' to `NULL` to disable.
#' @param col Colour of points for SNPs which do not have eQTLs.
#' @param outline_col Colour of symbol outlines. Default is `NA` for no outlines.
#' @param up_palette Character string specifying palette for upregulation eQTL
#'   using [grDevices::hcl.colors]
#' @param down_palette Character string specifying palette for downregulation
#'   eQTL using [grDevices::hcl.colors]
#' @param alpha Alpha opacity for points
#' @param tissue GTex tissue in which eQTL has been measured
#' @param eqtl_gene Gene showing eQTL effect
#' @param xlab x axis title
#' @param ylab y axis title
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
#' @param gene_col Colour for genes and exons.
#' @param exon_col Fill colour for exons.
#' @param exon_border Border line colour outlining exons. Set to `NA` for no 
#' border.
#' @param text_pos Character value of either 'top' or 'left' specifying 
#' placement of gene name labels.
#' @param legend_pos Position of legend. See [legend()]. Set to `NULL` to hide 
#' legend.
#' @param ... Other arguments passed to [plot()] for the scatter plot.
#' @return No return value. Produces a plot using base graphics.
#' 
#' @importFrom graphics axTicks axis layout par strwidth abline legend
#' @importFrom grDevices adjustcolor hcl.colors
#' @export

overlay_plot <- function(x, ...,
                     filter_gene_name = NULL,
                     filter_gene_biotype = NULL,
                     pcutoff = 5e-08,
                     col = 'black',
                     outline_col = NA,
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
  if (is.null(xlab)) xlab <- paste("Chromosome", x$seqname, "(Mb)")
  LDX <- x$LDexp[x$LDexp$Tissue == tissue & x$LDexp$Gene_Symbol == eqtl_gene, ]
  ind <- match(data[, x$labs], LDX$RS_ID)
  data$col <- 'black'
  data$col <- adjustcolor(data$col, alpha.f = alpha)
  data$pch <- 21
  if (all(is.na(ind))) {
    message("No significant eQTL")
  } else {
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
    eqind <- !is.na(data$eqtl_effect)
    equp <- eqind & sign(data$eqtl_effect) == -1
    eqdown <- eqind & sign(data$eqtl_effect) == 1
    data$col[equp] <- down_cols[ecol[equp]]
    data$col[eqdown] <- up_cols[ecol[eqdown]]
    data$pch[eqind] <- 24.5 - sign(data$eqtl_effect[eqind]) / 2
  }
  
  if (use_layout) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    graphics::layout(matrix(2:1, nrow = 2), heights = heights)
  }
  # lower locus plot
  genetracks(x, filter_gene_name, filter_gene_biotype,
             border, cex.axis, cex.text, gene_col, exon_col, exon_border,
             maxrows, text_pos, xticks = (xticks == 'bottom'),
             xlab = if (xticks == 'bottom') xlab else "")
  
  # scatter plot
  op <- par(mar = c(ifelse(xticks == 'top', 3, 0), 4, 2, 2))
  on.exit(par(op), add = TRUE)
  plot(data[, x$pos], data$logP,
       pch = data$pch, col = outline_col, bg = data$col,
       tcl = -0.3, las = 1, font.main = 1,
       mgp = c(1.7, 0.5, 0),
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
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis,
         tcl = -0.3, mgp = c(1.7, 0.4, 0))
  } else if (!border) {
    axis(1, at = axTicks(1), labels = FALSE, tcl = -0.3)
  }
}
