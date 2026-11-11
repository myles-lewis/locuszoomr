
#' Prepare a local copy of the GWAS Catalog
#'
#' Reads the GWAS Catalog associations download, keeps single-SNP associations
#' and the columns used by the GWAS Catalog track, and optionally lifts
#' positions over to GRCh37 and adds trait categories. Save the result with
#' `saveRDS()` and pass it, or the file path, to [zoom()] as `gwas_catalog`.
#'
#' @details The associations file ("All associations v1.0.2 - with added
#'   ontology annotations", unzipped) and the trait mappings file are at
#'   <https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/>; the
#'   hg38ToHg19 chain file (uncompressed) is at
#'   <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/>.
#'
#' @param file Path to the GWAS Catalog associations file (tsv), e.g.
#'   `gwas-catalog-download-associations-alt-full.tsv`.
#' @param chain Optional path to a UCSC hg38ToHg19 chain file (uncompressed),
#'   to add GRCh37 positions.
#' @param mappings Optional path or URL of the GWAS Catalog trait mappings file,
#'   to add trait categories.
#' @return List of class 'gwas_catalog' with elements `assoc` (one row per
#'   association), `category` (trait-category pairs, or `NULL`) and `date` (the
#'   latest date an association was added to the Catalog).
#' @seealso [link_gwas_catalog()] [zoom()]
#' @examples
#' \dontrun{
#' g <- prep_gwas_catalog("gwas-catalog-download-associations-alt-full.tsv",
#'                        chain = "hg38ToHg19.over.chain",
#'                        mappings = "gwas-efo-trait-mappings.tsv")
#' saveRDS(g, "gwas_catalog.rds")
#' }
#' @export

prep_gwas_catalog <- function(file, chain = NULL, mappings = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' must be installed to use this feature",
         call. = FALSE)
  }
  cols <- c("DATE ADDED TO CATALOG", "PUBMEDID", "DISEASE/TRAIT", "CHR_ID", "CHR_POS", "SNPS",
            "RISK ALLELE FREQUENCY", "P-VALUE", "PVALUE_MLOG", "OR or BETA",
            "95% CI (TEXT)")
  d <- data.table::fread(file, select = cols, sep = "\t", quote = "",
                         colClasses = "character", data.table = FALSE)
  date <- max(as.Date(d$`DATE ADDED TO CATALOG`), na.rm = TRUE)
  d$`DATE ADDED TO CATALOG` <- NULL
  d <- unique(d)
  d <- d[grepl("^rs[0-9]+$", d$SNPS) & grepl("^[0-9]+$", d$CHR_POS), ]
  # compact types (expanded per window by link_gwas_catalog()); rs numbers
  # exceed the integer range
  num <- function(x) suppressWarnings(as.numeric(x))
  assoc <- data.frame(rs = num(sub("^rs", "", d$SNPS)),
                      GWAS_Trait = factor(d$`DISEASE/TRAIT`),
                      chr = factor(d$CHR_ID),
                      pos38 = as.integer(d$CHR_POS),
                      P_value = num(d$`P-VALUE`),
                      logP = num(d$PVALUE_MLOG),
                      Risk_Allele = num(d$`RISK ALLELE FREQUENCY`),
                      # named to match LDtrait, whose two columns are swapped
                      Effect_Size_95_CI = num(d$`OR or BETA`),
                      Beta_or_OR = factor(d$`95% CI (TEXT)`),
                      PMID = as.integer(d$PUBMEDID))
  assoc$pos37 <- NA_integer_
  if (!is.null(chain)) {
    ch <- rtracklayer::import.chain(chain)
    gr <- GenomicRanges::GRanges(paste0("chr", assoc$chr),
                                 IRanges::IRanges(assoc$pos38, width = 1))
    assoc$pos37 <- as.integer(lifted_pos(rtracklayer::liftOver(gr, ch), assoc$chr))
  }
  assoc <- assoc[order(assoc$chr, assoc$pos38), ]
  rownames(assoc) <- NULL
  category <- NULL
  if (!is.null(mappings)) {
    category <- trait_category(levels(assoc$GWAS_Trait),
                               read_trait_mappings(mappings))
  }
  structure(list(assoc = assoc, category = category, date = date),
            class = "gwas_catalog")
}

# start positions from liftOver(); NA unless 1:1 and on the same chromosome
lifted_pos <- function(lo, chr) {
  pos <- rep(NA_real_, length(lo))
  ok <- lengths(lo) == 1
  u <- unlist(lo[ok])
  same <- as.character(GenomicRanges::seqnames(u)) == paste0("chr", chr[ok])
  pos[which(ok)[same]] <- GenomicRanges::start(u)[same]
  pos
}

#' Add local GWAS Catalog associations to a locus
#'
#' Fills `loc$LDtrait` (and `loc$trait_category`) from a local GWAS Catalog with
#' every association in the locus window, so [trait_plotly()] can show it. If
#' `loc$data` has an `ld` column (from [link_LD()]), it is used as the r^2 of
#' associations at SNPs in the data.
#'
#' @param loc Object of class 'locus'.
#' @param gcat Object from [prep_gwas_catalog()].
#' @return The 'locus' object with `LDtrait` and `trait_category` added.
#' @seealso [link_trait()] for the LDlink alternative.
#' @export

link_gwas_catalog <- function(loc, gcat) {
  a <- gcat$assoc
  build38 <- grepl("38", loc$genome)
  pos <- if (build38) a$pos38 else a$pos37
  if (all(is.na(pos))) {
    stop("Local GWAS Catalog has no ", if (build38) "GRCh38" else "GRCh37",
         " positions; prepare it with a liftOver chain", call. = FALSE)
  }
  w <- which(a$chr == sub("^chr", "", loc$seqname) & !is.na(pos) &
               pos >= loc$xrange[1] & pos <= loc$xrange[2])
  LDT <- a[w, ]
  LDT <- data.frame(RS_Number = sprintf("rs%.0f", LDT$rs),
                    GWAS_Trait = as.character(LDT$GWAS_Trait),
                    chr = as.character(LDT$chr),
                    pos = pos[w],
                    P_value = LDT$P_value, logP = LDT$logP,
                    Risk_Allele = LDT$Risk_Allele,
                    Effect_Size_95_CI = LDT$Effect_Size_95_CI,
                    Beta_or_OR = as.character(LDT$Beta_or_OR),
                    PMID = as.character(LDT$PMID))
  # P can underflow to 0: cap the axis, and order so collapse_traits() keeps
  # the strongest of tied P = 0 rows
  LDT$logP <- pmin(LDT$logP, 300)
  LDT <- LDT[order(LDT$pos, -a$logP[w]), ]
  ind <- match(LDT$RS_Number, loc$data[, loc$labs])
  LDT$indata <- !is.na(ind)
  if ("ld" %in% colnames(loc$data)) LDT$R2 <- loc$data$ld[ind]
  rownames(LDT) <- NULL
  loc$LDtrait <- LDT
  if (!is.null(gcat$category)) {
    loc$trait_category <- gcat$category[gcat$category$trait %in% LDT$GWAS_Trait, ]
  }
  loc
}
