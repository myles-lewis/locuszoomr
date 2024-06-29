
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
#' @param gene Optional character value specifying which gene to view. Either
#'   `gene`, or `xrange` plus `seqname`, or `index_snp` must be specified.
#' @param data Dataset (data.frame or data.table) to use for plot. If
#'   unspecified or `NULL`, gene track information alone is returned.
#' @param xrange Optional vector of genomic position range for the x axis.
#' @param seqname Optional, specifies which chromosome to plot.
#' @param flank Single value or vector with 2 values for how much flanking 
#' region left and right of the gene to show. Defaults to 100kb.
#' @param fix_window Optional alternative to `flank`, which allows users to
#'   specify a fixed genomic window centred on the specified gene. Both `flank`
#'   and `fix_window` cannot be specified simultaneously.
#' @param ens_db Either a character string which specifies which Ensembl
#'   database package (version 86 and earlier for Homo sapiens) to query for
#'   gene and exon positions (see `ensembldb` Bioconductor package). Or an
#'   `ensembldb` object which can be obtained from the AnnotationHub database.
#'   See the vignette and the `AnnotationHub` Bioconductor package for how to
#'   create this object.
#' @param chrom Determines which column in `data` contains chromosome 
#' information. If `NULL` tries to autodetect the column.
#' @param pos Determines which column in `data` contains position information. 
#' If `NULL` tries to autodetect the column.
#' @param p Determines which column in `data` contains SNP p-values. 
#' If `NULL` tries to autodetect the column.
#' @param yvar Specifies column in `data` for plotting on the y axis as an
#'   alternative to specifying p-values. Both `p` and `yvar` cannot be specified
#'   simultaneously.
#' @param labs Determines which column in `data` contains SNP rs IDs.
#' If `NULL` tries to autodetect the column.
#' @param index_snp Specifies the index SNP. If not specified, the SNP with the
#'   lowest P value is selected. Can be used to specify locus region instead of
#'   specifying `gene`, or `seqname` and `xrange`.
#' @param LD Optional character value to specify which column in `data` contains
#'   LD information.
#' @return Returns an object of class 'locus' ready for plotting, containing:
#' \item{seqname}{chromosome value}
#' \item{xrange}{vector of genomic position range}
#' \item{gene}{gene name}
#' \item{ens_db}{Ensembl or AnnotationHub database version}
#' \item{chrom}{column name in `data` containing chromosome information}
#' \item{pos}{column name in `data` containing position}
#' \item{p}{column name in `data` containing p-value}
#' \item{yvar}{column name in `data` to be plotted on y axis as alternative to 
#' `p`}
#' \item{labs}{column name in `data` containing SNP IDs}
#' \item{index_snp}{id of the most significant SNP}
#' \item{data}{the subset of GWAS data to be plotted}
#' \item{TX}{dataframe of transcript annotations}
#' \item{EX}{`GRanges` object of exon annotations}
#' If `data` is `NULL` when `locus()` is called then gene track information
#' alone is returned.
#' @seealso [locus_plot()] [locus_ggplot()] [locus_plotly()]
#' @examples
#' ## Bioconductor package EnsDb.Hsapiens.v75 is needed for these examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5,
#'              ens_db = "EnsDb.Hsapiens.v75")
#' summary(loc)
#' locus_plot(loc)
#' loc2 <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5,
#'               ens_db = "EnsDb.Hsapiens.v75")
#' locus_plot(loc2)
#' }
#' @importFrom ensembldb genes exons
#' @importFrom BiocGenerics start end
#' @importFrom AnnotationFilter GeneNameFilter AnnotationFilterList 
#' SeqNameFilter GeneIdFilter TxStartFilter TxEndFilter ExonStartFilter 
#' ExonEndFilter
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom memoise memoise
#' @export

locus <- function(gene = NULL,
                  data = NULL,
                  xrange = NULL, seqname = NULL,
                  flank = NULL, fix_window = NULL,
                  ens_db,
                  chrom = NULL, pos = NULL, p = NULL, yvar = NULL,
                  labs = NULL,
                  index_snp = NULL,
                  LD = NULL) {
  if (is.character(ens_db)) {
    if (!ens_db %in% (.packages())) {
      stop("Ensembl database not loaded. Try: library(", ens_db, ")",
           call. = FALSE)
    }
    edb <- get(ens_db)
  } else edb <- ens_db
  if (!is.null(flank) & !is.null(fix_window))
    stop("both `flank` and `fix_window` cannot be specified at the same time")
  if (is.null(flank)) flank <- 1e5
  flank <- rep_len(flank, 2)
  
  if (!is.null(gene)) {
    locus <- genes(edb, filter = AnnotationFilterList(
      GeneNameFilter(gene),
      SeqNameFilter(c(1:22, 'X', 'Y'))
    ))
    
    if(length(locus) > 1) {
      message(sprintf('Identified %d genes matching name \'%s\', taking first',
                      length(locus), gene))
      locus <- locus[1]
    }
    seqname <- names(seqlengths(locus))
    if (is.null(fix_window)) {
      xrange <- c(start(locus) - flank[1], end(locus) + flank[2])
    } else {
      m <- mean(c(start(locus), end(locus)))
      xrange <- as.integer(c(m - fix_window/2, m + fix_window/2))
    }
    xrange[xrange < 0] <- 0
  }
  
  if (!is.null(data)) {
    # autodetect headings
    dc <- detect_cols(data, chrom, pos, p, labs, yvar)
    chrom <- dc$chrom
    pos <- dc$pos
    p <- dc$p
    labs <- dc$labs
    
    if (!is.null(index_snp) & is.null(gene) & is.null(seqname) & is.null(xrange)) {
      # region based on index SNP
      if (!index_snp %in% data[, labs])
        stop("SNP specified by `index_snp` not found")
      ind <- which(data[, labs] == index_snp)
      if (length(ind) > 1) message("SNP found more than once")
      seqname <- data[ind[1], chrom]
      snp_pos <- data[ind[1], pos]
      xrange <- if (is.null(fix_window)) {
        c(snp_pos - flank[1], snp_pos + flank[2])
      } else {
        as.integer(c(snp_pos - fix_window/2, snp_pos + fix_window/2))
      }
      xrange[xrange < 0] <- 0
    }
  }
  
  if (is.null(xrange) | is.null(seqname)) stop('No locus specified')
  msg <- paste0("chromosome ", seqname, ", position ", xrange[1], " to ",
                xrange[2])
  if (!is.null(gene)) msg <- paste(gene, msg, sep = ", ")
  if (!is.null(index_snp)) msg <- paste(index_snp, msg, sep = ", ")
  message(msg)
  
  if (!is.null(data)) {
    data <- data[which(data[, chrom] == seqname), ]
    data <- data[which(data[, pos] > xrange[1] & data[, pos] < xrange[2]), ]
    # smallest floating point
    data[data[, p] < 5e-324, p] <- 5e-324
    if (is.null(yvar)) {
      data$logP <- -log10(data[, p])
      yvar <- "logP"
    }
    data <- as.data.frame(data)

    if (nrow(data) == 0) {
      message("Locus contains no SNPs/datapoints")
      data <- NULL
    } else {
      message(nrow(data), " SNPs/datapoints")
      if (is.null(index_snp)) index_snp <- data[which.max(data[, yvar]), labs]
      if (is.character(LD)) {
        colnames(data)[which(colnames(data) == LD)] <- "ld"
      }
    }
  }
  
  seqname <- gsub("chr|[[:punct:]]", "", seqname, ignore.case = TRUE)
  if (!seqname %in% c(1:22, "X", "Y")) 
    warning("`seqname` refers to a non-conventional chromosome")
  TX <- ensembldb::genes(edb, filter = AnnotationFilterList(
    SeqNameFilter(seqname),
    TxStartFilter(xrange[2], condition = "<"),
    TxEndFilter(xrange[1], condition = ">"),
    GeneIdFilter("ENSG", "startsWith")))
  TX <- data.frame(TX)
  TX <- TX[! is.na(TX$start), ]
  TX <- TX[!duplicated(TX$gene_id), ]
  
  if (nrow(TX) == 0) {
    message("No gene transcripts")
    # Creating empty exons object here in suitable format
    EX <- ensembldb::exons(edb, filter = AnnotationFilterList(
      SeqNameFilter(seqname),
      ExonStartFilter(xrange[2], condition = "<"),
      ExonEndFilter(xrange[1], condition = ">"),
      GeneIdFilter("ENSG", "startsWith")))
  } else {
    EX <- ensembldb::exons(edb, filter = GeneIdFilter(TX$gene_id))
  }

  loc <- list(seqname = seqname, xrange = xrange, gene = gene,
              ens_db = ens_db,
              chrom = chrom, pos = pos, p = p, yvar = yvar, labs = labs,
              index_snp = index_snp,
              data = data, TX = TX, EX = EX)
  class(loc) <- "locus"
  loc
}


#' @export
summary.locus <- function(object, ...) {
  if (!is.null(object$gene)) {
    cat("Gene", object$gene, "\n")
  } else if (!is.null(object$index_snp)) {
    cat("Index SNP", object$index_snp, "\n")
  }
  cat(paste0("Chromosome ", object$seqname, ", position ",
             format(object$xrange[1], big.mark=","), " to ",
             format(object$xrange[2], big.mark=","), "\n"))
  nr <- nrow(object$data)
  if (is.null(nr)) nr <- 0
  cat(nr, "SNPs/datapoints\n")
  cat(nrow(object$TX), "gene transcripts\n")
  if (nrow(object$TX) > 0) {
    tb <- sort(c(table(object$TX$gene_biotype)), decreasing = TRUE)
    cat(paste(tb, names(tb), collapse = ", "), "\n")
  }
}


detect_cols <- function(data, chrom, pos, p, labs = NULL, yvar = NULL) {
  # autodetect headings
  if (is.null(chrom)) {
    w <- grep("chr", colnames(data), ignore.case = TRUE)
    if (length(w) == 1) {
      chrom <- colnames(data)[w]
    } else stop("unable to autodetect chromosome column")
  }
  if (is.null(pos)) {
    w <- grep("pos", colnames(data), ignore.case = TRUE)
    if (length(w) == 1) {
      pos <- colnames(data)[w]
    } else stop("unable to autodetect SNP position column")
  }
  if (!is.null(p) && !is.null(yvar)) stop("cannot specify both `p` and `yvar`")
  if (is.null(p) && is.null(yvar)) {
    if ("p" %in% colnames(data)) {
      p <- "p"
    } else {
      w <- grep("^p?val", colnames(data), ignore.case = TRUE)
      if (length(w) == 1) {
        p <- colnames(data)[w]
      } else stop("unable to autodetect p-value column")
    }
  }
  if (is.null(labs)) {
    w <- grep("rs?id", colnames(data), ignore.case = TRUE)
    if (length(w) > 1) stop("unable to autodetect SNP id column")
    if (length(w) == 0) {
      w <- grep("SNP", colnames(data), ignore.case = TRUE)
    }
    if (length(w) == 1) {
      labs <- colnames(data)[w]
    } else stop("unable to autodetect SNP id column")
  }
  
  # check headings
  if (!chrom %in% colnames(data)) {
    stop("Column specified by `chrom` not found in `data`")}
  if (!pos %in% colnames(data)) {
    stop("Column specified by `pos` not found in `data`")}
  if (is.null(yvar) && !p %in% colnames(data)) {
    stop("Column specified by `p` not found in `data`")}
  if (!labs %in% colnames(data)) {
    stop("Column specified by `labs` not found in `data`")}
  if (!is.null(yvar)) {
    if (!yvar %in% colnames(data)) {
      stop("Column specified by `yvar` not found in `data`")
    }
  }
  
  list(chrom = chrom, pos = pos, p = p, labs = labs)
}
