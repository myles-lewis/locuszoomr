
#' R code for a static figure of a zoom() view
#'
#' Mirrors what `zoom()`'s locus plot draws (see `output$locus` in zoom.R):
#' the same `locus()`, `link_*()` calls and filters, then base-graphics panels.
#'
#' @param state List from `figure_state()`.
#' @param backend Plotting system; only `"base"` so far.
#' @return The code as a single string.
#' @noRd

figure_code <- function(state, backend = "base") {
  backend <- match.arg(backend, "base")
  s <- state
  ocall <- c("# Original call and objects needed",
             code_call("zoom", lapply(s$call.[-1], deparse)))
  pre <- character(0)
  # name of a user object that was originally loaded
  obj <- function(name, what) {
    pre <<- c(pre, sprintf('# `%s` = %s', name, what))
    name
  }
  libs <- "library(locuszoomr)"
  ens <- ens_db_source(s$call.$ens_db, s$ens_db)
  ens_code <- if (isTRUE(ens$package)) {
    libs <- c(libs, sprintf("library(%s)", ens$name))
    code_value(ens$name)
  } else {
    obj(ens$name, "ensembl database")
  }
  linked <- !is.null(s$link_type) && !is.null(s$link_snp) &&
    isTRUE(s$link_snp_in_view)
  lt <- if (linked) s$link_type else ""
  # write the token as the user passed it to zoom(), never its value
  token <- if (is.null(s$call.$ld_token)) 'Sys.getenv("LDLINK_TOKEN")' else
    paste(deparse(s$call.$ld_token, backtick = TRUE), collapse = " ")
  
  # data
  out <- character(0)
  add <- function(...) out <<- c(out, ...)
  locus_call <- function(data, i, assign) {
    args <- list(data, seqname = code_value(s$chr),
                 xrange = sprintf("c(%s, %s)", format(s$xrange[1], scientific = FALSE),
                                  format(s$xrange[2], scientific = FALSE)),
                 ens_db = ens_code,
                 chrom = code_value(s$chrom[i]), pos = code_value(s$pos[i]),
                 p = code_value(s$p[i]), labs = code_value(s$labs[i]))
    # as the app, loc2 only takes the link SNP for LD and eQTL
    if (linked && (i == 1 || lt %in% c("LD", "eQTL"))) {
      args$index_snp <- code_value(s$link_snp)
    }
    code_call("locus", args, assign)
  }
  locs <- "loc"
  add(locus_call(obj(arg_name(s$call.$data), "dataset 1"), 1, "loc"))
  if (isTRUE(s$man2)) {
    add(locus_call(obj(arg_name(s$call.$data2), "dataset 2"), 2, "loc2"))
    locs <- c(locs, "loc2")
  }
  each <- function(fun, ...) {
    for (l in locs) add(code_call(fun, c(list(l), list(...)), l))
  }
  if (isTRUE(s$recomb_on)) {
    each("link_recomb",
         recomb = obj(arg_name(s$call.$recomb), "recombination rate object"))
  }
  if (lt == "LD") each("link_LD", token = token, pop = code_value(s$ld_pop))
  if (lt == "eQTL") each("link_eqtl", token = token)
  if (lt == "GWAS trait") {
    add(code_call("link_trait", list("loc", token = token, category = "TRUE"), "loc"))
  }
  qtl <- lt == "Open Targets" && !isFALSE(s$qtl_found)
  if (qtl) add(code_call("link_ot", list("loc"), "loc"))
  if (lt == "Open Targets" && !qtl) {
    add("# no Open Targets credible sets were found for this SNP")
  }
  gcat_on <- isTRUE(s$catalog_on) && !is.null(s$gwas_catalog) && lt != "GWAS trait"
  if (gcat_on) {
    g <- s$gwas_catalog
    gname <- switch(g$how,
                    cache = {add("gcat <- load_gwas_catalog()"); "gcat"},
                    path = {add(sprintf("gcat <- readRDS(%s)", code_value(g$value))); "gcat"},
                    object = obj(g$value, "gwas_catalog object"))
    add(sprintf("loc <- link_gwas_catalog(loc, %s)", gname))
  }
  if (!is.null(s$eqtl_gene) && !is.null(s$select_gene) && s$select_gene != "all") {
    add(sprintf("loc$data <- loc$data[loc$data[, %s] >= %s | loc$data[, %s] == %s, ]",
                code_value(s$p[1]), code_value(s$pcutoff),
                code_value(s$eqtl_gene), code_value(s$select_gene)))
  }

  # one extra track: Open Targets takes the GWAS Catalog's place
  track <- NULL
  if (qtl) {
    track <- code_call("qtl_plot", c(list("loc"), code_args(
      type_filter = s$qtl_type_filter, gene_filter = s$qtl_gene_filter)))
  } else if (lt == "GWAS trait" || gcat_on) {
    track <- code_call("trait_plot", c(list("loc"), code_args(
      trait_filter = s$trait_filter, category_filter = s$category_filter,
      collapse = if (!isTRUE(s$collapse)) FALSE,
      r2_min = if (lt == "LD" && isTRUE(s$r2_min > 0)) s$r2_min),
      if (gcat_on && isTRUE(s$ld_colour) && lt == "LD") {
        list(ld_scheme = "eval(formals(scatter_plot)$LD_scheme)")
      },
      if (gcat_on) list(size_by_r2 = "FALSE")))
  }

  # plot
  plot_out <- character(0)
  padd <- function(...) plot_out <<- c(plot_out, ...)
  pcut <- if (!identical(s$pcutoff, 5e-8)) s$pcutoff
  gene_args <- c(code_args(filter_gene_biotype = s$biotype),
                 list(maxrows = code_value(s$maxrows)))
  overlay <- lt == "eQTL"
  if (is.null(track) && !isTRUE(s$man2) && !overlay) {
    padd(code_call("locus_plot", c(list("loc"), gene_args, code_args(
      pcutoff = pcut, eqtl_gene = s$eqtl_gene, beta = s$beta))))
  } else {
    k <- length(locs) + !is.null(track)
    padd(sprintf("oldpar <- set_layers(%d)", k))
    for (i in seq_along(locs)) {
      ylab <- if (isTRUE(s$man2) && length(s$traits) >= i) {
        list(ylab = sprintf("expression(%s ~ -log[10] ~ P)", code_value(s$traits[i])))
      }
      if (overlay) {
        padd(code_call("overlay_scatter", c(list(locs[i]), code_args(
          gene_filter = s$eqtl_gene_filter, tissue_filter = s$eqtl_tissue_filter,
          xticks = FALSE, pcutoff = pcut), ylab)))
      } else {
        padd(code_call("scatter_plot", c(list(locs[i]), code_args(
          xticks = FALSE, pcutoff = pcut,
          eqtl_gene = if (i == 1) s$eqtl_gene,
          beta = if (length(s$beta) >= i) s$beta[i]),
          r2_min = if (s$r2_min > 0) s$r2_min, ylab)))
      }
    }
    if (!is.null(track)) padd(track)
    padd(code_call("genetracks", c(list("loc"), gene_args)), "par(oldpar)")
  }

  notes <- character(0)
  if (!is.null(s$eqtl_gene)) {
    notes <- "# eQTL gene colours are chosen by scatter_plot() and may differ from the app"
  }
  if (overlay && isTRUE(s$man2)) {
    notes <- c(notes, "# eQTL gene colours are not matched between the two panels as in the app")
  }
  paste(c(ocall, pre, "", "# R code for the current zoom() view", libs, "",
          out, "", notes, plot_out), collapse = "\n")
}

# a value as R code; vector elements are joined with a soft break ("\001")
# that code_call() turns into ", " or a wrapped line
code_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (length(x) > 1 && is.atomic(x) && is.null(names(x))) {
    el <- vapply(x, function(e) deparse(e, width.cutoff = 500L), character(1))
    return(paste0("c(", paste(el, collapse = ",\001"), ")"))
  }
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

# named, non-NULL arguments as code
code_args <- function(...) {
  a <- list(...)
  a <- a[!vapply(a, is.null, logical(1))]
  lapply(a, code_value)
}

# fun(args) as code, optionally assigned, wrapped at 80 characters; args is a
# list of code strings, named for named arguments
code_call <- function(fun, args, assign = NULL, width = 80L) {
  prefix <- paste0(if (!is.null(assign)) paste(assign, "<- "), fun, "(")
  indent <- strrep(" ", nchar(prefix))
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  pieces <- ifelse(nzchar(nms), paste(nms, "=", unlist(args)), unlist(args))
  lines <- prefix
  wrapped <- FALSE
  for (i in seq_along(pieces)) {
    el <- strsplit(pieces[i], "\001", fixed = TRUE)[[1]]
    sep <- if (i > 1) ", " else ""
    one <- paste(el, collapse = " ")
    cur <- lines[length(lines)]
    if (i > 1 && (wrapped || nchar(cur) + nchar(sep) + nchar(one) > width - 2)) {
      lines[length(lines)] <- paste0(cur, ",")
      lines <- c(lines, indent)
      sep <- ""
    }
    lines[length(lines)] <- paste0(lines[length(lines)], sep)
    # vector elements: fill lines, continuing under the first element
    cont <- strrep(" ", nchar(lines[length(lines)]) + regexpr("(", el[1], fixed = TRUE)[1])
    wrapped <- FALSE
    for (j in seq_along(el)) {
      e <- el[j]
      cur <- lines[length(lines)]
      if (j > 1 && nchar(cur) + nchar(e) + 3 > width) {
        lines <- c(lines, paste0(cont, e))
        wrapped <- TRUE
      } else {
        lines[length(lines)] <- paste0(cur, if (j > 1) " ", e)
      }
    }
  }
  lines[length(lines)] <- paste0(lines[length(lines)], ")")
  paste(lines, collapse = "\n")
}

# how zoom() was given its arguments, for the code it writes
arg_name <- function(expr) {
  # if (!is.symbol(expr) || identical(expr, as.name("."))) return(NULL)
  deparse(expr, backtick = TRUE)
}

catalog_source <- function(expr, value) {
  if (is.null(value) || isFALSE(value)) return(NULL)
  if (isTRUE(value)) return(list(how = "cache"))
  if (is.character(value)) return(list(how = "path", value = value))
  if (!is.null(arg_name(expr))) return(list(how = "object", value = arg_name(expr)))
  list(how = "unknown")
}

ens_db_source <- function(expr, value) {
  if (is.character(value)) return(list(name = value, package = TRUE))
  if (!is.null(arg_name(expr))) return(list(name = arg_name(expr), package = FALSE))
  NULL
}
