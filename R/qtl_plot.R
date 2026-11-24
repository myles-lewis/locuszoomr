
#' Locus scatter plot for Open Targets QTL credible sets
#'
#' Static (base graphics) version of the Open Targets track in [zoom()]: one
#' point per credible set from [link_ot()] at its lead variant, coloured by QTL
#' type and sized by the posterior inclusion probability (PIP) of the index SNP.
#' Stack it between [scatter_plot()] and [genetracks()] using [set_layers()].
#'
#' @param loc Object of class 'locus' with element `OTqtl` from [link_ot()].
#' @param type_filter Character vector of QTL types to show.
#' @param gene_filter Character vector of QTL genes to show.
#' @param legend_pos Position of the QTL type legend, see [legend()], or `NULL`
#'   for none.
#' @inheritParams trait_plot
#' @return Invisibly, the credible sets plotted with columns `colour` and
#'   `cex`.
#' @seealso [link_ot()] [set_layers()]
#' @examples
#' \dontrun{
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = "IRF5", flank = 1e5,
#'              ens_db = "EnsDb.Hsapiens.v75")
#' loc <- link_ot(loc)
#' oldpar <- set_layers(2)
#' scatter_plot(loc, xticks = FALSE)
#' qtl_plot(loc)
#' genetracks(loc)
#' par(oldpar)
#' }
#' @export

qtl_plot <- function(loc,
                     type_filter = NULL,
                     gene_filter = NULL,
                     legend_pos = "topright",
                     xticks = FALSE,
                     border = FALSE,
                     showRecomb = TRUE,
                     cex = 1,
                     cex.axis = 0.9,
                     cex.lab = 1,
                     xlab = NULL,
                     ylab = NULL,
                     align = TRUE, ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  q <- loc$OTqtl
  if (is.null(q)) q <- data.frame(type = character(0), gene = character(0),
                                  pos = numeric(0), logP = numeric(0),
                                  pip = numeric(0))
  if (!is.null(type_filter)) q <- q[q$type %in% type_filter, ]
  if (!is.null(gene_filter)) q <- q[q$gene %in% gene_filter, ]
  q$colour <- unname(qtl_scheme[q$type])
  q$colour[is.na(q$colour)] <- "grey"
  q$cex <- (6 + 14 * sqrt(ifelse(is.na(q$pip), 0, q$pip))) / 7 * cex
  if (is.null(ylab)) ylab <- expression("QTL -log"[10] ~ "P")
  leg <- NULL
  if (!is.null(legend_pos) && nrow(q) > 0) {
    types <- sort(unique(q$type))
    cols <- qtl_scheme[types]
    cols[is.na(cols)] <- "grey"
    leg <- function() {
      legend(legend_pos, legend = types, pch = 21, pt.bg = cols, col = "black",
             bty = "n", cex = 0.8, y.intersp = 0.96)
    }
  }
  track_panel(loc, q$pos, q$logP, bg = q$colour, col = "black", cex = q$cex,
              blank = "No QTL data", xticks = xticks, border = border,
              showRecomb = showRecomb, cex.axis = cex.axis, cex.lab = cex.lab,
              xlab = xlab, ylab = ylab, align = align, panel_last = leg, ...)
  invisible(q)
}
