
#' Create locus object for plotting
#' 
#' Creates object of class 'locus' for genomic locus plot similar to 
#' `locuszoom`.
#' 
#' @details 
#' This is an R version of `locuszoom` (http://locuszoom.org) for generating
#' publication ready Manhattan plots of gene loci. It references Ensembl
#' databases using the `ensembldb` Bioconductor package framework for annotating
#' genes and exons in the locus.
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
#' information. If `NULL` tries to autodetect the column.
#' @param pos Determines which column in `data` contains position information. 
#' If `NULL` tries to autodetect the column.
#' @param p Determines which column in `data` contains SNP p values. 
#' If `NULL` tries to autodetect the column.
#' @param labs Determines which column in `data` contains SNP rs IDs.
#' If `NULL` tries to autodetect the column.
#' @param index_snp Specifies the index SNP. If not specified, the SNP with the
#'   lowest P value is selected.
#' @param LD Optional character value to specify which column in `data` contains
#'   LD information.
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
#' summary(loc)
#' locus_plot(loc)
#' loc2 <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5, LD = FALSE)
#' locus_plot(loc2)
#' @importFrom ensembldb genes exons
#' @importFrom BiocGenerics start end
#' @importFrom AnnotationFilter GeneNameFilter AnnotationFilterList 
#' SeqNameFilter GeneIdFilter
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom memoise memoise
#' @export

locus <- function(data, xrange = NULL, seqname = NULL,
                  gene = NULL, flank = 5e4,
                  ens_version = "EnsDb.Hsapiens.v75",
                  chrom = NULL, pos = NULL, p = NULL,
                  labs = NULL,
                  index_snp = NULL,
                  LD = NULL) {
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
  message("Chromosome ", seqname, ", position ", xrange[1], " to ", xrange[2])
  
  # autodetect headings
  if (is.null(chrom)) {
    w <- grep("chr", colnames(data), ignore.case = TRUE)
    if (length(w) == 1) {
      chrom <- colnames(data)[w]
    } else stop("cannot autodetect chromosome column")
  }
  if (is.null(pos)) {
    w <- grep("pos", colnames(data), ignore.case = TRUE)
    if (length(w) == 1) {
      pos <- colnames(data)[w]
    } else stop("cannot autodetect SNP position column")
  }
  if (is.null(p)) {
    if ("p" %in% colnames(data)) {
      p <- "p"
    } else {
      w <- grep("^p?val", colnames(data), ignore.case = TRUE)
      if (length(w) == 1) {
        p <- colnames(data)[w]
      } else stop("cannot autodetect p-value column")
    }
  }
  if (is.null(labs)) {
    w <- grep("rs?id|SNP", colnames(data), ignore.case = TRUE)
    if (length(w) == 1) {
      labs <- colnames(data)[w]
    } else stop("cannot autodetect SNP id column")
  }
  
  # check headings
  if (!chrom %in% colnames(data)) {
    stop("Column specified by `chrom` not found in `data`")}
  if (!pos %in% colnames(data)) {
    stop("Column specified by `pos` not found in `data`")}
  if (!p %in% colnames(data)) {
    stop("Column specified by `p` not found in `data`")}
  if (!labs %in% colnames(data)) {
    stop("Column specified by `labs` not found in `data`")}
  
  data <- data[data[, chrom] == seqname &
                 data[, pos] > xrange[1] & data[, pos] < xrange[2], ]
  data$logP <- -log10(data[, p])
  data <- as.data.frame(data)
  if (is.null(index_snp)) index_snp <- data[which.max(data$logP), labs]
  if (is.character(LD)) {
    colnames(data)[which(colnames(data) == LD)] <- "ld"
  }
  message(nrow(data), " SNPs/datapoints")
  
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
              index_snp = index_snp,
              data = data, TX = TX, EX = EX)
  class(loc) <- "locus"
  loc
}


#' @export
summary.locus <- function(object, ...) {
  cat("Gene", object$gene, "\n")
  cat(paste0("Chromosome ", object$seqname, ", position ",
             format(object$xrange[1], big.mark=","), " to ",
             format(object$xrange[2], big.mark=","), "\n"))
  cat(nrow(object$data), "SNPs/datapoints\n")
  cat(nrow(object$TX), "gene transcripts\n")
  tb <- sort(c(table(object$TX$gene_biotype)), decreasing = TRUE)
  cat(paste(tb, names(tb), collapse = ", "), "\n")
}

