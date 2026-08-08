# Client-side gene track re-layout.
#
# The gene annotation panel is packed once by mapRow() at render time. These
# helpers ship the gene and exon tables into the widget so that JavaScript can
# re-pack them against the visible window whenever the user zooms or pans.


#' Build the gene track payload passed to JavaScript
#'
#' @param TX Data frame of transcripts, already ordered by [mapRow()] and with
#'   `start`/`end` in Mb.
#' @param EX Data frame of exons with `gene_id`, `start`, `end` in Mb.
#' @param cfg Named list of render settings. `italics` is required; the caller
#'   adds the remaining display settings.
#' @return A list with elements `tx`, `ex` and `cfg`. `tx$hover` carries the
#'   per-gene hover string in the same format as [genetrack_ly()]'s line
#'   trace, so the JS re-layout handler can keep the line trace's `text`
#'   array in sync with its `x`/`y` arrays after a re-pack.
#' @noRd
genetrack_payload <- function(TX, EX, cfg) {
  nm <- TX$gene_name
  if (isTRUE(cfg$italics)) nm <- paste0("<i>", nm, "</i>")
  pos <- as.character(TX$strand) == "+"
  label <- ifelse(pos, paste0(nm, "&#8594;"), paste0("&#8592;", nm))
  label[TX$gene_name == ""] <- ""

  # Matches the hovertext built in genetrack_ly() (R/genetrack_ly.R:188-193)
  # exactly, including that TX$fullname contributes nothing when absent.
  hover <- paste0(TX$gene_name,
                   TX$fullname,
                   "<br>Gene ID: ", TX$gene_id,
                   "<br>Biotype: ", TX$gene_biotype,
                   "<br>Start: ", TX$start * 1e6,
                   "<br>End: ", TX$end * 1e6)

  tx <- data.frame(
    id       = as.character(TX$gene_id),
    name     = as.character(TX$gene_name),
    label    = label,
    hover    = hover,
    start    = as.numeric(TX$start),
    end      = as.numeric(TX$end),
    strand   = as.character(TX$strand),
    priority = seq_len(nrow(TX)),
    stringsAsFactors = FALSE
  )

  gene_idx <- match(as.character(EX$gene_id), tx$id) - 1L
  keep <- !is.na(gene_idx)
  ex <- data.frame(
    gene_idx = as.integer(gene_idx[keep]),
    start    = as.numeric(EX$start[keep]),
    end      = as.numeric(EX$end[keep]),
    stringsAsFactors = FALSE
  )

  list(tx = tx, ex = ex, cfg = cfg)
}


#' Locate the gene track traces and shapes within a built plotly object
#'
#' `subplot()` renames axes and merges every panel's shapes into one array, so
#' positions cannot be assumed. The two gene track traces carry `meta` tags; the
#' exon shapes are then identified by sharing those traces' y axis reference.
#'
#' Note: `meta` is set as a plain (non-`I()`) value inside a data-bound
#' `add_segments()`/`add_text()` pipeline, so `plotly_build()` recycles it to
#' match each trace's row count rather than keeping it scalar; `add_segments()`
#' additionally interleaves `NA` at its line-break separators. Matching is
#' therefore done against any non-`NA` element rather than requiring a scalar.
#'
#' @param built A plotly object that has been through [plotly::plotly_build()].
#' @return A list of 0-based indices and axis identifiers.
#' @noRd
resolve_genetrack_idx <- function(built) {
  has_tag <- function(d, tag) {
    m <- d$meta
    if (is.null(m) || !is.character(m)) return(FALSE)
    any(m == tag, na.rm = TRUE)
  }

  li <- which(vapply(built$x$data, has_tag, logical(1),
                      tag = "locuszoomr_genetrack_lines"))
  ti <- which(vapply(built$x$data, has_tag, logical(1),
                      tag = "locuszoomr_genetrack_labels"))
  if (length(li) != 1L || length(ti) != 1L) {
    stop("could not locate gene track traces in the built plotly object",
         call. = FALSE)
  }

  xax <- built$x$data[[li]]$xaxis
  yax <- built$x$data[[li]]$yaxis
  if (is.null(xax)) xax <- "x"
  if (is.null(yax)) yax <- "y"

  shapes <- built$x$layout$shapes
  shp <- if (length(shapes) == 0L) integer(0) else {
    which(vapply(shapes, function(s) identical(s$yref, yax), logical(1)))
  }

  # I() forces jsonlite/htmlwidgets to serialise shapeIdx as a JSON array even
  # when it has exactly one element or zero elements; without it auto_unbox
  # collapses a length-1 integer vector to a bare number, which breaks the JS
  # side's `.indexOf()` call (TypeError: idxs.indexOf is not a function).
  list(lineTrace  = as.integer(li - 1L),
       labelTrace = as.integer(ti - 1L),
       shapeIdx   = I(as.integer(shp - 1L)),
       xaxis      = xax,
       yaxis      = yax)
}


#' Attach the client-side re-layout handler to a plotly object
#'
#' @param p A plotly object. It is built internally so trace and shape indices
#'   can be resolved.
#' @param TX,EX Transcript and exon data frames in Mb, `TX` already ordered by
#'   [mapRow()].
#' @param cfg Named list of render settings. Every current call site routes
#'   `cfg` through [genetrack_cfg()], which always populates `showExons`,
#'   `geneCol`, `exonCol` and `exonBorder`, so in practice this function never
#'   receives them missing. The `TRUE`/`col2hex("skyblue")`/
#'   `col2hex("blue4")`/`col2hex("blue4")` fallbacks below are defensive only
#'   — not a documented part of the expected input — in case a future or
#'   external caller supplies an incomplete `cfg`. `geneCol` is only read
#'   client-side when `showExons` is `FALSE`; `exonCol` only when it is
#'   `TRUE`; `exonBorder` is read in both cases.
#' @return The built plotly object with an `onRender` handler attached.
#' @importFrom htmlwidgets onRender
#' @noRd
add_genetrack_relayout <- function(p, TX, EX, cfg) {
  built <- plotly::plotly_build(p)
  if (is.null(cfg$showExons)) cfg$showExons <- TRUE
  if (is.null(cfg$geneCol)) cfg$geneCol <- col2hex("skyblue")
  if (is.null(cfg$exonCol)) cfg$exonCol <- col2hex("blue4")
  if (is.null(cfg$exonBorder)) cfg$exonBorder <- col2hex("blue4")
  payload <- genetrack_payload(TX, EX, cfg)
  payload$idx <- resolve_genetrack_idx(built)

  js <- paste(readLines(
    system.file("js", "genetrack-relayout.js", package = "locuszoomr"),
    warn = FALSE), collapse = "\n")

  htmlwidgets::onRender(
    built,
    sprintf("function(el, x, data) {\n%s\nLZR.attach(el, data);\n}", js),
    data = payload)
}


#' Assemble the render settings passed to the client
#'
#' @param italics Logical, whether gene names are italicised.
#' @param cex.text Font size multiplier.
#' @param maxrows Resolved row capacity of the panel.
#' @param showExons Logical, whether exons are drawn individually.
#' @param gene_col,exon_col,exon_border Colours.
#' @return A named list.
#' @noRd
genetrack_cfg <- function(italics, cex.text, maxrows, showExons, gene_col,
                          exon_col, exon_border) {
  list(italics    = italics,
       fontSizePx = 14 * cex.text,
       maxrows    = maxrows,
       showExons  = showExons,
       gapFrac    = 0.02,
       textPos    = "top",
       geneCol    = gene_col,
       exonCol    = exon_col,
       exonBorder = exon_border)
}
