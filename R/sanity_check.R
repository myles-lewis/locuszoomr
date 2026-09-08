
#' Sanity checks on GWAS dataframes
#' 
#' Performs sanity checks on GWAS dataframes to check similar chromosome
#' content, overlap in SNP ids and alignment of SNPs by chromosome/position.
#' 
#' @param data Dataframe of GWAS results with columns for chromosome, position,
#'   p value and SNP rs IDs. Data.tables are coerced to dataframe.
#' @param data2 Dataframe of 2nd GWAS results for comparison.
#' @param chrom Vector determines which column in `data` and `data2` contains
#'   chromosome information. If `NULL` or `NA` tries to autodetect the columns.
#'   The 1st element refers to `data` and 2nd element refers to `data2`.
#' @param pos Determines which column in `data` and `data2` contains position
#'   information. See `chrom`.
#' @param p Determines which column in `data` and `data2` contains SNP
#'   p-values. See `chrom`.
#' @param labs Determines which column in `data` and `data2` contains SNP rs
#'   IDs. See `chrom`.
#' @returns No return value. Prints diagnostic messages.
#' @export

sanity_check <- function(data, data2,
                         chrom = NULL, pos = NULL, p = NULL, labs = NULL) {
  dat_name <- deparse(substitute(data))
  dat2_name <- deparse(substitute(data2))
  data <- data.frame(data)
  data2 <- data.frame(data2)
  
  # autodetect headings
  dc <- detect_cols(data, chrom[1], pos[1], p[1], labs[1])
  if (!is.null(chrom)) chrom <- rep_len(chrom, 2)
  if (!is.null(pos)) pos <- rep_len(pos, 2)
  if (!is.null(p)) p <- rep_len(p, 2)
  if (!is.null(labs)) labs <- rep_len(labs, 2)
  dc2 <- detect_cols(data2, chrom[2], pos[2], p[2], labs[2])
  
  chrom <- c(dc$chrom, dc2$chrom)
  pos <- c(dc$pos, dc2$pos)
  p <- c(dc$p, dc2$p)
  labs <- c(dc$labs, dc2$labs)
  
  # check chromosomes
  chr_set <- list()
  chr_set[[1]] <- unique(data[, chrom[1]])
  chr_set[[2]] <- unique(data2[, chrom[2]])
  if (!setequal(chr_set[[1]], chr_set[[2]])) {
    message("Difference in chromosomes")
    s <- setdiff(chr_set[[1]], chr_set[[2]])
    if (length(s)) {
      message("  Chromosome ", paste(s, collapse = ", "),
              " found in ", dat_name, ", but not ", dat2_name)
    }
    s <- setdiff(chr_set[[2]], chr_set[[1]])
    if (length(s)) {
      message("  Chromosome ", paste(s, collapse = ", "),
              " found in ", dat2_name, ", but not ", dat_name)
    }
  }
  
  # check SNP ids
  snp1 <- unique(data[, labs[1]])
  snp2 <- unique(data2[, labs[2]])
  overlap <- intersect(snp1, snp2)
  message(length(snp1), " unique SNP ids in ", dat_name)
  message(length(snp2), " unique SNP ids in ", dat2_name)
  message(length(overlap), " SNP ids overlap")
  message("  ", format(length(overlap) / length(snp1) * 100, digits = 3), "% of ",
          dat_name, " and ",
          format(length(overlap) / length(snp2) * 100, digits = 3), "% of ",
          dat2_name, " overlap")
  
  # check SNP alignment
  cpos1 <- paste0(data[, chrom[1]], ":", data[, pos[1]])
  cpos2 <- paste0(data2[, chrom[2]], ":", data2[, pos[2]])
  overlap <- intersect(cpos1, cpos2)
  message(length(overlap), " SNPs align by chrom/position")
  message("  ", format(length(overlap) / length(cpos1) * 100, digits = 3), "% of ",
          dat_name, " and ",
          format(length(overlap) / length(cpos2) * 100, digits = 3), "% of ",
          dat2_name, " align")
}
