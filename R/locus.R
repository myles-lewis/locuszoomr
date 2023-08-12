
#' Create locus object for plotting
#' 
#' Creates object of class 'locus' for genomic locus plot similar to 
#' `locuszoom`.
#' 
#' @details 
#' This is an R version of `locuszoom` (http://locuszoom.org) for generating
#' publication ready Manhattan plots of gene loci. It references Ensembl
#' databases using the `ensembldb` Bioconductor package framework for annotating
#' genes and exons. It queries LDlink (https://ldlink.nci.nih.gov/) via the
#' [LDlinkR] package to retrieve linkage disequilibrium (LD) information on a
#' reference SNP.
#' 
#' @param data Dataset (data.frame or data.table) to use for plot.
#' @param xrange Vector of genomic position range for the x axis.
#' @param seqname Specifies which chromosome to plot.
#' @param gene Optional specifies which gene to view. Either `xrange` with 
#' `seqname`, or `gene` must be specified.
#' @param flank Single value or vector with 2 values for how much flanking 
#' region left and right of the gene to show.
#' @param ens_version Specifies which Ensembl database to query for gene and 
#' exon positions. See `ensembldb` Bioconductor package.
#' @param chrom Determines which column in `data` contains chromosome 
#' information. Automatically looks for PLINK headings.
#' @param pos Determines which column in `data` contains position information. 
#' Automatically looks for PLINK headings.
#' @param p Determines which column in `data` contains SNP p values. 
#' Automatically looks for PLINK headings.
#' @param labs Determines which column in `data` contains SNP rs IDs.
#' Automatically looks for PLINK headings.
#' @param index_snp Specifies the index SNP for displaying linkage 
#' disequilibrium (LD). If not specified, the SNP with the lowest P value is 
#' selected.
#' @param LD Logical or character. If logical specifies whether LD is plotted by
#'   querying 1000 Genomes via `LDlinkR` package. See [LDlinkR::LDmatrix].
#'   Results are cached using the `memoise` package, so that if exactly the same
#'   locus is requested the system does not repeatedly call the API. If
#'   set to a character value, this determines which column in `data` contains
#'   LD information.
#' @param eQTL Logical whether to obtain eQTL data. Queries GTEx eQTL data via
#'   [LDlinkR::LDexpress] using the SNP specified by `index_snp`
#' @param pop A 1000 Genomes Project population, (e.g. YRI or CEU), multiple 
#' allowed, default = "CEU". Passed to [LDlinkR::LDmatrix].
#' @param r2d Either "r2" for LD r^2 or "d" for LD D', default = "r2". Passed 
#' to [LDlinkR::LDmatrix].
#' @param LDtoken Personal access token for accessing 1000 Genomes LD data via 
#' LDlink API. See [LDlinkR].
#' @return Returns an object of class 'locus' ready for plotting, containing 
#' the subset of GWAS data to be plotted, 
#' chromosome and genomic position range, 
#' Ensembl database version number, 
#' column names for chromosome, position, SNP ID, p value
#' locus gene information from Ensembl and
#' locus exon information from Ensembl.
#' @examples
#' library(EnsDb.Hsapiens.v75)
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5, LD = FALSE)
#' plot(loc)
#' loc2 <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5, LD = FALSE)
#' plot(loc2)
#' @import EnsDb.Hsapiens.v75
#' @importFrom ensembldb genes exons
#' @importFrom BiocGenerics start end
#' @importFrom LDlinkR LDmatrix LDexpress
#' @importFrom AnnotationFilter GeneNameFilter AnnotationFilterList 
#' SeqNameFilter GeneIdFilter
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom memoise memoise
#' @export

locus <- function(data, xrange = NULL, seqname = NULL,
                  gene = NULL, flank = 5e4,
                  ens_version = "EnsDb.Hsapiens.v75",
                  chrom = 'chrom', pos = 'pos', p = 'p',
                  labs = 'rsid',
                  index_snp = NULL,
                  LD = FALSE,
                  eQTL = FALSE,
                  pop = "CEU",
                  r2d = "r2",
                  LDtoken = "") {
  if (!ens_version %in% (.packages())) {
    stop("Ensembl database not loaded. Try: library(", ens_version, ")",
         call. = FALSE)
  }
  edb <- get(ens_version)
  flank <- rep_len(flank, 2)
  if (!is.null(gene)) {
    locus <- genes(edb, filter = GeneNameFilter(gene))
    xrange <- c(start(locus) - flank[1], end(locus) + flank[2])
    seqname <- names(seqlengths(locus))
  }
  if (is.null(xrange) | is.null(seqname)) stop('No locus specified')
  # PLINK headings
  if (!chrom %in% colnames(data)) {
    if ("CHR" %in% colnames(data)) {chrom <- "CHR"
    } else stop("Column specified by `chrom` not found in `data`")
  }
  if (!pos %in% colnames(data)) {
    if ("BP" %in% colnames(data)) {pos <- "BP"
    } else stop("Column specified by `pos` not found in `data`")
  }
  if (!p %in% colnames(data)) {
    if ("P" %in% colnames(data)) {p <- "P"
    } else stop("Column specified by `p` not found in `data`")
  }
  if (!labs %in% colnames(data)) {
    if ("SNP" %in% colnames(data)) {labs <- "SNP"
    } else stop("Column specified by `labs` not found in `data`")
  }
  data <- data[data[, chrom] == seqname &
                 data[, pos] > xrange[1] & data[, pos] < xrange[2], ]
  data$logP <- -log10(data[, p])
  data <- as.data.frame(data)
  if (is.null(index_snp)) index_snp <- data[which.max(data$logP), labs]
  if (is.character(LD)) {
    colnames(data)[which(colnames(data) == LD)] <- "ld"
  } else if (LD) {
    if (LDtoken == "") stop("LDtoken is missing")
    rslist <- data[, labs]
    if (length(rslist) > 1000) {
      rslist <- rslist[order(data$logP, decreasing = TRUE)[seq_len(1000)]]
    }
    message("Obtaining LD on ", length(rslist), " SNPs", appendLF = FALSE)
    ldm <- mem_LDmatrix(rslist, pop = pop, r2d = r2d, token = LDtoken)
    ld <- ldm[, index_snp]
    data$ld <- ld[match(data[, labs], ldm$RS_number)]
  }
  if (eQTL) {
    if (LDtoken == "") stop("LDtoken is missing")
    LDexp <- mem_LDexpress(snps = index_snp, pop = pop, r2d = r2d, 
                           token = LDtoken)
    for (i in c("R2", "D'", "Effect_Size", "P_value")) {
      LDexp[, i] <- as.numeric(LDexp[, i])
    }
    LDexp$Effect_Allele <- gsub("=.*", "", LDexp$Effect_Allele_Freq)
  }
  
  TX <- ensembldb::genes(edb, filter = AnnotationFilterList(
    SeqNameFilter(c(seq_len(22), "X", "Y")),
    GeneIdFilter("ENSG", "startsWith")))
  TX <- data.frame(TX)
  TX <- TX[! is.na(TX$start), ]
  TX <- TX[TX$seqnames == seqname, ]
  TX <- TX[TX$end > xrange[1], ]
  TX <- TX[TX$start < xrange[2], ]
  EX <- ensembldb::exons(edb, filter = GeneIdFilter(TX$gene_id))
  
  loc <- list(seqname = seqname, xrange = xrange, gene = gene,
              ens_version = ens_version,
              chrom = chrom, pos = pos, p = p, labs = labs,
              data = data, TX = TX, EX = EX)
  if (eQTL) loc$LDexp <- LDexp
  class(loc) <- "locus"
  loc
}
