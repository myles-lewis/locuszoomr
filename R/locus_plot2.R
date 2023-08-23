

locus_plot2 <- function(x,
                        filter_gene_name = NULL,
                        filter_gene_biotype = NULL,
                        xlab = NULL,
                        cex.axis = 1,
                        cex.text = 0.7,
                        use_layout = TRUE,
                        maxrows = 7,
                        xticks = 'bottom',
                        border = FALSE,
                        gene_col = 'blue4',
                        exon_col = 'blue4',
                        exon_border = 'blue4',
                        text_pos = 'top', ...) {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  
  if (use_layout) {
    op0 <- set_layers(rev = TRUE)
    on.exit(par(op0), add = TRUE)
  }
  # lower locus plot
  genetracks(x, filter_gene_name, filter_gene_biotype,
             border, cex.axis, cex.text, gene_col, exon_col, exon_border,
             maxrows, text_pos, xticks = (xticks == 'bottom'),
             xlab = if (xticks == 'bottom') xlab else "")
  
  scatter_plot(x, xticks = (xticks == 'top'),
               border = border, cex.axis = cex.axis, xlab = xlab, ...)
}
