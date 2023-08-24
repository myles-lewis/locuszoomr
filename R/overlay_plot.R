
#' Plot overlaying eQTL and GWAS data
#'
#' Experimental plotting function for overlaying eQTL data from GTEx on top of
#' GWAS results. y axis shows the -log10 p-value for the GWAS result.
#' Significant eQTL for the specified gene are overlaid using colours and
#' symbols.
#' 
#' @param x Object of class 'locus' to use for plot. See [locus].
#' @param base_col Colour of points for SNPs which do not have eQTLs.
#' @param up_palette Character string specifying palette for upregulation eQTL
#'   using [grDevices::hcl.colors]
#' @param down_palette Character string specifying palette for downregulation
#'   eQTL using [grDevices::hcl.colors]
#' @param alpha Alpha opacity for points
#' @param tissue GTex tissue in which eQTL has been measured
#' @param eqtl_gene Gene showing eQTL effect
#' @param ... Other arguments passed to [locus_plot()] for the locus plot.
#' @return No return value. Produces a plot using base graphics.
#' 
#' @importFrom grDevices adjustcolor hcl.colors
#' @export

overlay_plot <- function(x,
                         base_col = 'black',
                         up_palette = "OrRd",
                         down_palette = "Blues 3",
                         alpha = 0.5,
                         tissue = "Whole Blood",
                         eqtl_gene = x$gene,
                         ...) {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  if (!"LDexp" %in% names(x)) stop("Missing eQTL data")
  if (is.null(eqtl_gene)) stop("eqtl_gene not specified")
  
  x$data$bg <- adjustcolor(base_col, alpha.f = alpha)
  x$data$pch <- 21
  
  LDX <- x$LDexp[x$LDexp$Tissue == tissue & x$LDexp$Gene_Symbol == eqtl_gene, ]
  # match by rsid
  ind <- match(x$data[, x$labs], LDX$RS_ID)
  message(sum(!is.na(ind)), "/", nrow(LDX), " matched eQTL SNPs (total ", 
          nrow(x$data), ")")
  
  if (all(is.na(ind))) {
    message("No significant eQTL")
  } else {
    x$data$eqtl_effect <- NA
    x$data$eqtl_effect <- LDX$Effect_Size[ind]
    x$data$eqtl_p <- LDX$P_value[ind]
    x$data$eqtl_effect_allele <- LDX$Effect_Allele[ind]
    # gwas allele and eqtl effect allele are the same
    mismatch <- which(x$data$eqtl_effect_allele != x$data$effect_allele)
    which_rev <- x$data$other_allele[mismatch] == x$data$eqtl_effect_allele[mismatch]
    rev_effect <- mismatch[which_rev]
    mismatch <- mismatch[!which_rev]
    x$data$eqtl_effect[rev_effect] <- -x$data$eqtl_effect[rev_effect]
    x$data$eqtl_effect[mismatch] <- NA
    up_cols <- hcl.colors(8, up_palette, rev = TRUE)[-c(1,2)]
    down_cols <- hcl.colors(8, down_palette, rev = TRUE)[-c(1,2)]
    ecol <- cut(abs(x$data$eqtl_effect), breaks=6)
    eqind <- !is.na(x$data$eqtl_effect)
    equp <- eqind & sign(x$data$eqtl_effect) == -1
    eqdown <- eqind & sign(x$data$eqtl_effect) == 1
    x$data$bg[equp] <- down_cols[ecol[equp]]
    x$data$bg[eqdown] <- up_cols[ecol[eqdown]]
    x$data$pch[eqind] <- 24.5 - sign(x$data$eqtl_effect[eqind]) / 2
  }
  
  locus_plot(x, col = NA, showLD = FALSE, ...)
}
