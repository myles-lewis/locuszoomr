
#' Plot overlaying eQTL and GWAS data
#'
#' Experimental plotting function for overlaying eQTL data from GTEx on top of
#' GWAS results. y axis shows the -log10 p-value for the GWAS result.
#' Significant eQTL for the specified gene are overlaid using colours and
#' symbols.
#' 
#' @param x Object of class 'locus' to use for plot. See [locus()].
#' @param base_col Colour of points for SNPs which do not have eQTLs.
#' @param alpha Alpha opacity for non-eQTL points
#' @param scheme Character string specifying palette for effect size showing
#'   up/downregulation eQTL using [grDevices::hcl.colors]. Alternatively a
#'   vector of 6 colours.
#' @param tissue GTex tissue in which eQTL has been measured
#' @param eqtl_gene Gene showing eQTL effect
#' @param legend_pos Character value specifying legend position. See [legend()].
#' @param ... Other arguments passed to [locus_plot()] for the locus plot.
#' @return No return value. Produces a plot using base graphics.
#' 
#' @importFrom grDevices adjustcolor hcl.colors
#' @export

overlay_plot <- function(x,
                         base_col = 'black',
                         alpha = 0.5,
                         scheme = "RdYlBu",
                         tissue = "Whole Blood",
                         eqtl_gene = x$gene,
                         legend_pos = "topright",
                         ...) {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  if (!"LDexp" %in% names(x)) stop("Missing eQTL data")
  if (is.null(eqtl_gene)) stop("eqtl_gene not specified")
  
  x$data$bg <- adjustcolor(base_col, alpha.f = alpha)
  x$data$pch <- 21
  # x$data$col <- NA
  
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
    if (length(scheme) == 1) {
      scheme <- hcl.colors(9, scheme)[-c(4:6)]
    }
    up_cols <- rev(scheme[1:3])
    down_cols <- scheme[4:6]
    ecol <- cut(abs(x$data$eqtl_effect), breaks = 3)
    eqind <- !is.na(x$data$eqtl_effect)
    eqdown <- eqind & sign(x$data$eqtl_effect) == -1
    equp <- eqind & sign(x$data$eqtl_effect) == 1
    x$data$bg[equp] <- up_cols[ecol[equp]]
    x$data$bg[eqdown] <- down_cols[ecol[eqdown]]
    x$data$pch[eqind] <- 24.5 - sign(x$data$eqtl_effect[eqind]) / 2
    # x$data$col[eqind] <- "black"
    x$data <- x$data[order(x$data$pch), ]
    pcex <- rep_len(0.9, nrow(x$data))
    pcex[x$data$pch != 21] <- 1.1
    labs <- levels(ecol)
    cutlev <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
                    upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
    cutlev <- signif(cutlev, 2)
  }
  
  if (!is.null(legend_pos)) {
    legtext <- c(rev(paste(cutlev[,1], cutlev[,2], sep=" : ")),
                 paste(-cutlev[,2], -cutlev[,1], sep=" : "))
    legendFUN <- substitute(legend(lpos,
                                   legend = legtext,
                                   pch = rep(c(24, 25), each=3),
                                   pt.bg = cols, col = NA,
                                   title = "eQTL effect",
                                   bty = "n", cex = 0.85, pt.cex = 1,
                                   y.intersp = 0.95),
                            list(lpos = legend_pos, legtext = legtext,
                                 cols = c(rev(up_cols), down_cols)))
  } else legendFUN <- NULL
  
  locus_plot(x, cex = pcex, col = NA, showLD = FALSE, panel.last = legendFUN, ...)
}
