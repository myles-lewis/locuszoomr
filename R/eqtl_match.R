# Count how many fetched eQTL SNPs can be shown on the locus plot
# @param loc 'locus' object with `LDexp` from [link_eqtl()]
# @param data Full GWAS dataframe
# @param chrom_val Chromosome value as stored in `data`
# @param chrom,labs Column names in `data` for chromosome and SNP id
# @param gene_filter,tissue_filter As in [overlay_plotly()]
# @return Named integer vector: total unique eQTL SNPs, shown in the window,
#   in the GWAS but outside the window, and not in the GWAS at all. `NULL` if
#   there is no eQTL data.

eqtl_match_counts <- function(loc, data, chrom_val, chrom, labs,
                              gene_filter = NULL, tissue_filter = NULL) {
  LDX <- loc$LDexp
  if (!is.null(gene_filter)) LDX <- LDX[LDX$Gene_Symbol %in% gene_filter, ]
  if (!is.null(tissue_filter)) LDX <- LDX[LDX$Tissue %in% tissue_filter, ]
  if (is.null(LDX) || nrow(LDX) == 0) return(NULL)
  snps <- unique(LDX$RS_ID)
  shown <- snps %in% loc$data[, labs]
  in_gwas <- shown | snps %in% data[which(data[, chrom] == chrom_val), labs]
  c(total = length(snps), shown = sum(shown),
    outside = sum(in_gwas & !shown), absent = sum(!in_gwas))
}


eqtl_match_msg <- function(counts, trait = NULL) {
  if (is.null(counts)) return(NULL)
  msg <- paste0(counts["shown"], " of ", counts["total"], " shown")
  if (counts["outside"] > 0) {
    msg <- paste0(msg, ", ", counts["outside"], " outside window")
  }
  if (counts["absent"] > 0) {
    msg <- paste0(msg, ", ", counts["absent"], " not in GWAS")
  }
  paste0(if (!is.null(trait)) paste0(trait, " "), "eQTL SNPs: ", msg)
}
