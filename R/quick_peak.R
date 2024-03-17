
#' Fast peak finder in GWAS data
#' 
#' Simple but fast function for finding peaks in genome-wide association study
#' (GWAS) data based on setting a minimum distance between peaks.
#' 
#' @param data GWAS dataset (data.frame or data.table)
#' @param npeaks Number of peaks to find. If set to `NA`, algorithm finds all
#'   distinct peaks separated from one another by region size specified by
#'   `span`.
#' @param p_cutoff Specifies cut-off for p-value significance above which
#'   p-values are ignored.
#' @param span Minimum genomic distance between peaks (default 1 Mb)
#' @param min_points Minimum number of p-value significant points which must lie
#'   within the span of a peak. This removes peaks with single or only a few low
#'   p-value SNPs. To disable set `min_points` to 1 or less.
#' @param chrom Determines which column in `data` contains chromosome
#'   information. If `NULL` tries to autodetect the column.
#' @param pos Determines which column in `data` contains position information.
#'   If `NULL` tries to autodetect the column.
#' @param p Determines which column in `data` contains SNP p-values. If `NULL`
#'   tries to autodetect the column.
#' @details
#' This function is designed for speed. SNP p-values are filtered to only those
#' which are significant as specified by `p_cutoff`. Each peak is identified as
#' the SNP with the lowest p-value and then SNPs in proximity to each peak
#' within the distance specified by `span` are removed. Regions such as the HLA
#' whose peaks may well be broader than `span` may produce multiple entries.
#' @returns Vector of row indices
#' @export

quick_peak <- function(data, npeaks = NA, p_cutoff = 5e-08, span = 1e6,
                       min_points = 2,
                       chrom = NULL, pos = NULL, p = NULL) {
  start <- Sys.time()
  # autodetect column headings
  dc <- detect_cols(data, chrom, pos, p)
  chrom <- dc$chrom
  pos <- dc$pos
  p <- dc$p
  
  if (is.na(npeaks)) npeaks <- Inf
  
  i <- 2
  index <- order(data[, p])
  index <- index[data[index, p] < p_cutoff]
  pks <- index[1]
  del <- which(abs(data[pks, pos] - data[index, pos]) <= span &
                 data[pks, chrom] == data[index, chrom])[-1]
  if (length(del) > 0) index <- index[-del]
  while (length(pks) < npeaks & i <= length(index)) {
    if (all(abs(data[index[i], pos] - data[pks, pos]) > span |
            data[index[i], chrom] != data[pks, chrom])) {
      del <- which(abs(data[index[i], pos] - data[index, pos]) <= span & 
                     data[index[i], chrom] == data[index, chrom])
      if (length(del) >= min_points) {
        pks <- c(pks, index[i])
        del <- del[-1]
        if (length(del) > 0) index <- index[-del]
      }
    }
    i <- i + 1
  }
  end <- Sys.time()
  message(length(pks), " peaks found (", format(end - start, digits = 3), ")")
  if (!is.infinite(npeaks) & length(pks) < npeaks)
    message("lower p_cutoff to find more peaks")
  pks
}
