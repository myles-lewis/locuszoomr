
#' Table of GWAS Catalog trait associations
#'
#' Builds the table shown under [zoom()]'s Table tab from [link_trait()]
#' results, applying the same filters and collapsing as [trait_plotly()].
#'
#' @param loc List or 'locus' object with element `LDtrait` and optionally
#'   `trait_category`.
#' @param trait_filter,category_filter,collapse See [trait_plotly()].
#' @return Dataframe sorted by P value, or `NULL` if there is no trait data.
#' @noRd

trait_table <- function(loc, trait_filter = NULL, category_filter = NULL,
                        collapse = TRUE) {
  LDT <- loc$LDtrait
  if (is.null(LDT)) return(NULL)
  tc <- loc$trait_category
  if (!is.null(trait_filter)) {
    LDT <- LDT[LDT$GWAS_Trait %in% trait_filter, ]
  }
  if (!is.null(category_filter) && !is.null(tc)) {
    LDT <- LDT[LDT$GWAS_Trait %in% tc$trait[tc$category %in% category_filter], ]
  }
  if (collapse) {
    key <- paste(LDT$RS_Number, LDT$GWAS_Trait, sep = "\t")
    pmid <- tapply(LDT$PMID, key, function(x) paste(unique(x), collapse = ", "))
    LDT <- collapse_traits(LDT)
    LDT$PMID <- as.character(
      pmid[paste(LDT$RS_Number, LDT$GWAS_Trait, sep = "\t")])
  }
  tab <- data.frame(SNP = LDT$RS_Number, Position = LDT$pos,
                    r2 = LDT$R2, Trait = LDT$GWAS_Trait,
                    check.names = FALSE)
  if (!is.null(tc)) {
    cats <- tapply(tc$category, tc$trait, paste, collapse = "; ")
    tab$Category <- as.character(cats[LDT$GWAS_Trait])
  }
  tab$P <- LDT$P_value
  # LDtrait's Effect_Size_95_CI and Beta_or_OR columns are swapped
  tab$`OR/beta` <- LDT$Effect_Size_95_CI
  tab$CI <- LDT$Beta_or_OR
  tab$RAF <- LDT$Risk_Allele
  if (collapse) tab$Associations <- LDT$n_assoc
  tab$PMID <- LDT$PMID
  tab$`In data` <- LDT$indata
  tab <- tab[order(tab$P), ]
  rownames(tab) <- NULL
  tab
}

#' Table of LDlink eQTL associations
#'
#' @param loc List or 'locus' object with element `LDexp` from [link_eqtl()].
#' @param gene_filter,tissue_filter See [overlay_plotly()].
#' @return Dataframe sorted by P value, or `NULL` if there is no eQTL data.
#' @noRd

eqtl_table <- function(loc, gene_filter = NULL, tissue_filter = NULL) {
  LDX <- loc$LDexp
  if (is.null(LDX)) return(NULL)
  if (!is.null(gene_filter)) LDX <- LDX[LDX$Gene_Symbol %in% gene_filter, ]
  if (!is.null(tissue_filter)) LDX <- LDX[LDX$Tissue %in% tissue_filter, ]
  tab <- data.frame(SNP = LDX$RS_ID, Position = LDX$pos, r2 = LDX$R2,
                    Gene = LDX$Gene_Symbol, Tissue = LDX$Tissue,
                    `Effect allele` = LDX$Effect_Allele,
                    `Effect size` = LDX$Effect_Size, P = LDX$P_value,
                    `In data` = LDX$inData, check.names = FALSE)
  tab <- tab[order(tab$P), ]
  rownames(tab) <- NULL
  tab
}

# comma-separated PMIDs -> PubMed links; non-digits dropped so output is safe HTML
pubmed_links <- function(x) {
  vapply(strsplit(ifelse(is.na(x), "", x), ",\\s*"), function(ids) {
    ids <- gsub("[^0-9]", "", ids)
    ids <- ids[nzchar(ids)]
    if (!length(ids)) return("")
    paste0('<a href="https://pubmed.ncbi.nlm.nih.gov/', ids,
           '/" target="_blank">', ids, '</a>', collapse = ", ")
  }, character(1))
}
