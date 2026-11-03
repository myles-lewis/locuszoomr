
#' Collapse duplicate GWAS trait associations
#'
#' Keeps the row with the smallest P value for each SNP-trait pair and records
#' the number of associations collapsed in column `n_assoc`.
#'
#' @param LDT Dataframe of LDtrait results.
#' @return Dataframe with one row per SNP-trait pair.
#' @noRd

collapse_traits <- function(LDT) {
  if (nrow(LDT) == 0) return(cbind(LDT, n_assoc = integer(0)))
  key <- paste(LDT$RS_Number, LDT$GWAS_Trait, sep = "\t")
  n <- table(key)
  o <- order(LDT$P_value)
  keep <- sort(o[!duplicated(key[o])])
  LDT <- LDT[keep, ]
  LDT$n_assoc <- as.integer(n[key[keep]])
  LDT
}
