
#' Download and cache the GWAS Catalog
#'
#' Downloads the GWAS Catalog associations and trait mappings and the UCSC
#' hg38ToHg19 liftOver chain, prepares them with [prep_gwas_catalog()] and
#' saves the result in the locuszoomr cache directory, from where
#' `zoom(..., gwas_catalog = TRUE)` loads it. Run it again to update the cached
#' copy; the Catalog is updated weekly.
#'
#' @param dir Directory to save the prepared Catalog in. Defaults to the
#'   locuszoomr user cache directory, see [tools::R_user_dir()].
#' @param ask Logical whether to ask before downloading.
#' @return Invisibly, the path of the saved file, or `NULL` if the user
#'   declined.
#' @details About 90 MB is downloaded; the saved file is about 25 MB.
#'   Preparation takes about 20 seconds.
#' @seealso [prep_gwas_catalog()] [zoom()]
#' @examples
#' \dontrun{
#' download_gwas_catalog()
#' zoom(SLE_gwas_sub, ens_db = "EnsDb.Hsapiens.v75", gwas_catalog = TRUE)
#' }
#' @export

download_gwas_catalog <- function(dir = gwas_catalog_dir(),
                                  ask = interactive()) {
  if (!has_data_table()) {
    stop("Package 'data.table' must be installed to use this feature",
         call. = FALSE)
  }
  if (ask && !isTRUE(gcat_ask(paste0(
    "Download the GWAS Catalog (about 90 MB) and save a prepared copy ",
    "(about 25 MB) in ", dir, "?")))) {
    return(invisible(NULL))
  }
  tmp <- tempfile("gwas_catalog")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  src <- file.path(tmp, c(associations = "associations.zip",
                          mappings = "mappings.tsv", chain = "chain.gz"))
  names(src) <- names(gwas_catalog_urls)
  for (i in names(gwas_catalog_urls)) {
    message("Downloading ", basename(gwas_catalog_urls[[i]]))
    gcat_fetch(gwas_catalog_urls[[i]], src[[i]])
  }
  tsv <- utils::unzip(src[["associations"]], exdir = tmp)
  tsv <- tsv[grepl("\\.tsv$", tsv)][1]
  if (is.na(tsv)) {
    stop("no .tsv in ", basename(gwas_catalog_urls[["associations"]]),
         call. = FALSE)
  }
  chain <- gunzip_file(src[["chain"]], file.path(tmp, "chain"))
  message("Preparing GWAS Catalog")
  g <- prep_gwas_catalog(tsv, chain = chain, mappings = src[["mappings"]])
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, "gwas_catalog.rds")
  # write then rename, so an interrupted download keeps the previous copy
  part <- tempfile("gwas_catalog", tmpdir = dir, fileext = ".part")
  on.exit(unlink(part), add = TRUE)
  gcat_save(g, part)
  if (!file.rename(part, f)) stop("could not write ", f, call. = FALSE)
  message("Saved GWAS Catalog (associations added up to ", g$date, ") to ", f)
  invisible(f)
}

#' Load the cached GWAS Catalog
#'
#' Loads the GWAS Catalog saved by [download_gwas_catalog()].
#'
#' @param dir Directory the Catalog was saved in.
#' @return Object of class 'gwas_catalog', or `NULL` with a message if there is
#'   no cached copy.
#' @seealso [download_gwas_catalog()]
#' @export

load_gwas_catalog <- function(dir = gwas_catalog_dir()) {
  f <- file.path(dir, "gwas_catalog.rds")
  if (!file.exists(f)) {
    message("No downloaded GWAS Catalog found; run download_gwas_catalog() ",
            "to get one")
    return(NULL)
  }
  g <- readRDS(f)
  message("GWAS Catalog: associations added up to ", g$date,
          " (update with download_gwas_catalog())")
  g
}

gwas_catalog_dir <- function() {
  if (r_version() < "4.0.0") {
    stop("The GWAS Catalog cache needs R >= 4.0; pass `dir` instead",
         call. = FALSE)
  }
  tools::R_user_dir("locuszoomr", "cache")
}

r_version <- function() getRversion()

gwas_catalog_urls <- c(
  associations = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip",
  mappings = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-efo-trait-mappings.tsv",
  chain = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz")

gcat_fetch <- function(url, dest) {
  old <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old))
  if (utils::download.file(url, dest, mode = "wb", quiet = TRUE) != 0) {
    stop("download of ", url, " failed", call. = FALSE)
  }
}

gcat_save <- function(x, f) saveRDS(x, f)

has_data_table <- function() requireNamespace("data.table", quietly = TRUE)

# zoom()'s gwas_catalog argument: NULL/FALSE, TRUE (cached copy), an RDS path
# or a prepared object
as_gwas_catalog <- function(x) {
  if (is.null(x) || isFALSE(x)) return(NULL)
  if (isTRUE(x)) return(load_gwas_catalog())
  if (is.character(x)) x <- readRDS(x)
  if (!inherits(x, "gwas_catalog")) {
    stop("`gwas_catalog` must be TRUE, or an object (or RDS file) from ",
         "prep_gwas_catalog()", call. = FALSE)
  }
  x
}

gcat_ask <- function(msg) utils::askYesNo(msg)

gunzip_file <- function(src, dest) {
  con <- gzfile(src, "rb")
  out <- file(dest, "wb")
  on.exit({close(con); close(out)})
  repeat {
    b <- readBin(con, raw(), 1e6)
    if (length(b) == 0) break
    writeBin(b, out)
  }
  dest
}

gcat_categories <- c(
  "Biological process",
  "Body measurement",
  "Cancer",
  "Cardiovascular disease",
  "Cardiovascular measurement",
  "Digestive system disorder",
  "Hematological measurement",
  "Immune system disorder",
  "Inflammatory measurement",
  "Lipid or lipoprotein measurement",
  "Metabolic disorder",
  "Neurological disorder",
  "Other disease",
  "Other measurement",
  "Other trait",
  "Unmapped"
)
