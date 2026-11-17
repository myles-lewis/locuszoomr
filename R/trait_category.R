
#' Read GWAS Catalog trait mappings
#'
#' Reads the GWAS Catalog file mapping reported traits to EFO parent
#' categories (e.g. "Immune system disorder", "Cancer").
#'
#' @param file Path or URL of the GWAS Catalog trait mappings file.
#' @return Dataframe with columns `trait` and `category`, one row per distinct
#'   trait-category pair.
#' @noRd

read_trait_mappings <- function(
    file = "https://www.ebi.ac.uk/gwas/api/search/downloads/trait_mappings") {
  message("Downloading GWAS Catalog mappings", appendLF = FALSE)
  start <- Sys.time()
  txt <- readLines(file, encoding = "UTF-8", warn = FALSE)
  fields <- strsplit(txt, "\t", fixed = TRUE)
  cols <- match(c("Disease trait", "Parent term"), fields[[1]])
  fields <- fields[-1]
  map <- data.frame(trait = vapply(fields, `[`, character(1), cols[1]),
                    category = vapply(fields, `[`, character(1), cols[2]))
  map <- unique(map)
  rownames(map) <- NULL
  end <- Sys.time()
  message(" (", format(end - start, digits = 3),")")
  map
}

mem_trait_mappings <- memoise(read_trait_mappings)

#' Map GWAS traits to GWAS Catalog categories
#'
#' Traits with several parent categories get one row per category; traits not
#' in the mappings get category "Unmapped".
#'
#' @param traits Character vector of reported GWAS traits.
#' @param map Dataframe from `read_trait_mappings()`.
#' @return Dataframe with columns `trait` and `category`.
#' @noRd

trait_category <- function(traits, map) {
  traits <- unique(traits)
  tc <- map[map$trait %in% traits, ]
  miss <- setdiff(traits, tc$trait)
  if (length(miss)) {
    tc <- rbind(tc, data.frame(trait = miss, category = "Unmapped"))
  }
  tc$category[tc$category == "NR"] <- "Unmapped"
  rownames(tc) <- NULL
  tc
}

# fetch categories for link_trait(); NA with a message if download fails
link_category <- function(traits) {
  map <- try(mem_trait_mappings(), silent = TRUE)
  if (inherits(map, "try-error")) {
    message("GWAS Catalog trait categories unavailable: ",
            trimws(attr(map, "condition")$message))
    return(NA)
  }
  trait_category(traits, map)
}
