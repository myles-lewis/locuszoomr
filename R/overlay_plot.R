
#' Plot overlaying eQTL and GWAS data
#'
#' Experimental plotting function for overlaying eQTL data from GTEx on top of
#' GWAS results. y axis shows the -log10 p-value for the GWAS result.
#' Significant eQTL for the specified gene are overlaid using colours and
#' symbols.
#' 
#' @param loc Object of class 'locus' to use for plot. See [locus()].
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

overlay_plot <- function(loc,
                         base_col = 'black',
                         alpha = 0.5,
                         scheme = "RdYlBu",
                         tissue = "Whole Blood",
                         eqtl_gene = loc$gene,
                         legend_pos = "topright",
                         ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  if (!"LDexp" %in% names(loc)) stop("Missing eQTL data")
  if (is.null(eqtl_gene)) stop("eqtl_gene not specified")
  
  loc$data$bg <- adjustcolor(base_col, alpha.f = alpha)
  loc$data$pch <- 21
  # loc$data$col <- NA
  
  LDX <- loc$LDexp[loc$LDexp$Tissue == tissue & loc$LDexp$Gene_Symbol == eqtl_gene, ]
  # match by rsid
  ind <- match(loc$data[, loc$labs], LDX$RS_ID)
  message(sum(!is.na(ind)), "/", nrow(LDX), " matched eQTL SNPs (total ", 
          nrow(loc$data), ")")
  
  if (all(is.na(ind))) {
    message("No significant eQTL")
  } else {
    loc$data$eqtl_effect <- NA
    loc$data$eqtl_effect <- LDX$Effect_Size[ind]
    loc$data$eqtl_p <- LDX$P_value[ind]
    loc$data$eqtl_effect_allele <- LDX$Effect_Allele[ind]
    # gwas allele and eqtl effect allele are the same
    mismatch <- which(loc$data$eqtl_effect_allele != loc$data$effect_allele)
    which_rev <- loc$data$other_allele[mismatch] == loc$data$eqtl_effect_allele[mismatch]
    rev_effect <- mismatch[which_rev]
    mismatch <- mismatch[!which_rev]
    loc$data$eqtl_effect[rev_effect] <- -loc$data$eqtl_effect[rev_effect]
    loc$data$eqtl_effect[mismatch] <- NA
    if (length(scheme) == 1) {
      scheme <- hcl.colors(9, scheme)[-c(4:6)]
    }
    up_cols <- rev(scheme[1:3])
    down_cols <- scheme[4:6]
    ecol <- cut(abs(loc$data$eqtl_effect), breaks = 3)
    eqind <- !is.na(loc$data$eqtl_effect)
    eqdown <- eqind & sign(loc$data$eqtl_effect) == -1
    equp <- eqind & sign(loc$data$eqtl_effect) == 1
    loc$data$bg[equp] <- up_cols[ecol[equp]]
    loc$data$bg[eqdown] <- down_cols[ecol[eqdown]]
    loc$data$pch[eqind] <- 24.5 - sign(loc$data$eqtl_effect[eqind]) / 2
    # loc$data$col[eqind] <- "black"
    loc$data <- loc$data[order(loc$data$pch), ]
    pcex <- rep_len(0.9, nrow(loc$data))
    pcex[loc$data$pch != 21] <- 1.1
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
                                   y.intersp = 0.96),
                            list(lpos = legend_pos, legtext = legtext,
                                 cols = c(rev(up_cols), down_cols)))
  } else legendFUN <- NULL
  
  locus_plot(loc, cex = pcex, col = NA, showLD = FALSE, panel.last = legendFUN, ...)
}
