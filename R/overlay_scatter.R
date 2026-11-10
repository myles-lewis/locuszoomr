
#' Locus scatter plot with an eQTL overlay
#'
#' Static (base graphics) version of the eQTL overlay in [zoom()] and
#' [overlay_plotly()]: a [scatter_plot()] panel of the GWAS data in which SNPs
#' that are eQTLs in `loc$LDexp` (from [link_eqtl()]) are filled with their
#' gene's colour and drawn as up or down triangles by the sign of the eQTL
#' effect. Unlike [overlay_plot()], it draws a single panel, so it can be
#' stacked with other panels using [set_layers()].
#'
#' @param loc Object of class 'locus' with element `LDexp` from [link_eqtl()].
#' @param gene_filter Character vector of eQTL genes to show.
#' @param tissue_filter Character vector of eQTL tissues to show.
#' @param eqtl_scheme Optional vector of gene colours, named by gene or in
#'   order of each gene's smallest eQTL P value. Defaults to `rainbow()`.
#' @param legend_pos Position of the gene legend, see [legend()], or `NULL` for
#'   none.
#' @param ... Optional arguments passed to [scatter_plot()], e.g. `xticks`,
#'   `pcutoff`, `ylab`.
#' @return Invisibly, the eQTL associations plotted (smallest P per SNP) with
#'   columns `colour` and `pch`.
#' @seealso [link_eqtl()] [overlay_plotly()] [set_layers()]
#' @examples
#' \dontrun{
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = "IRF5", flank = 1e5,
#'              ens_db = "EnsDb.Hsapiens.v75")
#' loc <- link_eqtl(loc, token = Sys.getenv("LDLINK_TOKEN"))
#' oldpar <- set_layers(1)
#' overlay_scatter(loc, xticks = FALSE)
#' genetracks(loc)
#' par(oldpar)
#' }
#' @export

overlay_scatter <- function(loc,
                            gene_filter = NULL,
                            tissue_filter = NULL,
                            eqtl_scheme = NULL,
                            legend_pos = "topright",
                            ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  LDX <- loc$LDexp
  if (is.null(LDX)) LDX <- data.frame(RS_ID = character(0),
                                      Gene_Symbol = character(0),
                                      Tissue = character(0),
                                      P_value = numeric(0),
                                      Effect_Size = numeric(0))
  if (!is.null(gene_filter)) LDX <- LDX[LDX$Gene_Symbol %in% gene_filter, ]
  if (!is.null(tissue_filter)) LDX <- LDX[LDX$Tissue %in% tissue_filter, ]
  # as overlay_plotly(): smallest P per SNP, SNPs in the data only
  LDX <- min_p_by_col(LDX, "RS_ID")
  LDX <- LDX[LDX$RS_ID %in% loc$data[, loc$labs], ]
  genes <- min_p_by_col(LDX, "Gene_Symbol")$Gene_Symbol
  scheme <- if (is.null(eqtl_scheme)) {
    rainbow(length(genes))
  } else if (is.null(names(eqtl_scheme))) {
    rep_len(eqtl_scheme, length(genes))
  } else unname(eqtl_scheme[genes])
  LDX$colour <- scheme[match(LDX$Gene_Symbol, genes)]
  LDX$pch <- ifelse(LDX$Effect_Size < 0, 25, 24)
  rownames(LDX) <- NULL

  if (!is.null(loc$data)) {
    d <- loc$data
    d$bg <- "grey"
    d$pch <- 21
    d$col <- "grey40"
    d$cex <- 0.8
    i <- match(LDX$RS_ID, d[, loc$labs])
    d$bg[i] <- LDX$colour
    d$pch[i] <- LDX$pch
    d$col[i] <- "black"
    d$cex[i] <- 1.2
    # eQTL SNPs drawn last, on top
    loc$data <- d[order(seq_len(nrow(d)) %in% i), ]
  }
  leg <- NULL
  if (!is.null(legend_pos) && length(genes) > 0) {
    leg <- function() {
      legend(legend_pos, legend = c(genes, "up", "down"),
             pch = c(rep(21, length(genes)), 2, 6),
             pt.bg = c(scheme, NA, NA), col = "black", title = "eQTL gene",
             bty = "n", cex = 0.8, y.intersp = 0.96)
    }
  }
  # the legend goes in plot()'s panel.last, before scatter_plot() restores
  # the margins
  scatter_plot(loc, showLD = FALSE, legend_pos = NULL,
               panel.last = if (!is.null(leg)) as.call(list(leg)), ...)
  invisible(LDX)
}
