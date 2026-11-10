
#' Locus scatter plot for GWAS Catalog traits
#'
#' Static (base graphics) version of [trait_plotly()]: a panel of GWAS Catalog
#' associations from [link_trait()] or [link_gwas_catalog()], to stack between
#' [scatter_plot()] and [genetracks()] using [set_layers()]. Points are coloured
#' and sized as in [trait_plotly()].
#'
#' @param loc Object of class 'locus' with element `LDtrait`.
#' @param trait_filter,category_filter,collapse,r2_min,ld_scheme,size_by_r2 See
#'   [trait_plotly()].
#' @param xticks Logical whether x axis numbers and title are plotted.
#' @param border Logical whether a bounding box is plotted.
#' @param showRecomb Logical whether to leave room for a recombination axis on
#'   the right, as [scatter_plot()] does when `loc` has recombination data.
#' @param cex,cex.axis,cex.lab Point, axis number and axis title sizes.
#' @param xlab,ylab Axis titles.
#' @param align Logical whether to set the margins to line up with
#'   [scatter_plot()] and [genetracks()].
#' @param ... Optional arguments passed to [plot()].
#' @return Invisibly, the associations plotted with columns `colour`,
#'   `outline` and `size` (marker size in plotly pixels).
#' @seealso [trait_plotly()] [set_layers()]
#' @examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = "IRF5", flank = 1e5,
#'              ens_db = "EnsDb.Hsapiens.v75")
#' \dontrun{
#' loc <- link_trait(loc, token = Sys.getenv("LDLINK_TOKEN"))
#' oldpar <- set_layers(2)
#' scatter_plot(loc, xticks = FALSE)
#' trait_plot(loc)
#' genetracks(loc)
#' par(oldpar)
#' }
#' }
#' @export

trait_plot <- function(loc,
                       trait_filter = NULL,
                       category_filter = NULL,
                       collapse = TRUE,
                       r2_min = NULL,
                       ld_scheme = NULL,
                       size_by_r2 = TRUE,
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
  LDT <- trait_colours(loc, trait_filter, category_filter, collapse, r2_min,
                       ld_scheme, size_by_r2)
  if (is.null(ylab)) ylab <- expression("GWAS trait -log"[10] ~ "P")
  track_panel(loc, LDT$pos, LDT$logP, bg = LDT$colour, col = LDT$outline,
              cex = LDT$size / 7 * cex, blank = "No trait data",
              xticks = xticks, border = border, showRecomb = showRecomb,
              cex.axis = cex.axis, cex.lab = cex.lab, xlab = xlab, ylab = ylab,
              align = align, ...)
  invisible(LDT)
}

# GWAS Catalog associations as trait_plotly() draws them: filtered, collapsed,
# in drawing order, with marker colour, outline and size (plotly px)
trait_colours <- function(loc, trait_filter = NULL, category_filter = NULL,
                          collapse = TRUE, r2_min = NULL, ld_scheme = NULL,
                          size_by_r2 = TRUE) {
  LDT <- loc$LDtrait
  if (is.null(LDT)) {
    return(data.frame(pos = numeric(0), logP = numeric(0), colour = character(0),
                      outline = character(0), size = numeric(0)))
  }
  if (!is.null(trait_filter)) LDT <- LDT[LDT$GWAS_Trait %in% trait_filter, ]
  if (!is.null(category_filter) && !is.null(loc$trait_category)) {
    tc <- loc$trait_category
    LDT <- LDT[LDT$GWAS_Trait %in% tc$trait[tc$category %in% category_filter], ]
  }
  if (collapse) LDT <- collapse_traits(LDT)
  if (nrow(LDT) == 0) {
    return(cbind(LDT, colour = character(0), outline = character(0),
                 size = numeric(0)))
  }
  traits <- sort(unique(LDT$GWAS_Trait))
  if (length(traits) > 1) {
    trait_group <- similarGroups(traits)
    scheme <- rainbow(length(unique(trait_group)))
    LDT$colour <- scheme[trait_group[LDT$GWAS_Trait]]
  } else {
    LDT$colour <- "red"
  }
  if (!is.null(ld_scheme)) {
    bin <- cut(LDT$R2, -1:6/5, labels = FALSE)
    bin[LDT$R2 %in% 0] <- 2L
    bin[is.na(bin)] <- 1L
    bin[LDT$RS_Number %in% loc$index_snp] <- 7L
    hex <- rgb(t(col2rgb(rep_len(ld_scheme, 7))), maxColorValue = 255)
    LDT$colour <- hex[bin]
    LDT <- LDT[order(bin), ]
  }
  LDT$outline <- "black"
  if (!is.null(r2_min) && r2_min > 0) {
    grey <- is.na(LDT$R2) | LDT$R2 < r2_min
    LDT$colour[grey] <- "#CCCCCC"
    LDT$outline[grey] <- "#999999"
  }
  LDT$size <- if (size_by_r2 && is.null(ld_scheme)) {
    6 + 8 * ifelse(is.na(LDT$R2), 0, LDT$R2)
  } else 7
  LDT
}

# panel margins lining up with scatter_plot() and genetracks()
track_mar <- function(xticks, recomb) {
  c(ifelse(xticks, 3, 0.1), 3.5, 0.25, ifelse(recomb, 3.5, 1.5))
}

# scatter panel for the static tracks, styled as scatter_plot()
track_panel <- function(loc, x, y, bg, col, cex, blank, xticks, border,
                        showRecomb, cex.axis, cex.lab, xlab, ylab, align,
                        pch = 21, panel_last = NULL, ...) {
  recomb <- !is.null(loc$recomb) && showRecomb
  if (align) {
    op <- par(mar = track_mar(xticks, recomb))
    on.exit(par(op))
  }
  if (is.null(xlab)) xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  empty <- !any(is.finite(y))
  ylim <- if (empty) c(0, 1) else range(y[is.finite(y)])
  ylim <- ylim + diff(ylim) * c(-0.05, 0.05)
  if (diff(ylim) == 0) ylim <- ylim + c(-1, 1)
  plot.args <- list(x = if (empty) NA else x, y = if (empty) NA else y,
                    pch = pch, bg = bg, col = col, cex = cex,
                    las = 1, cex.axis = cex.axis, cex.lab = cex.lab,
                    xlim = loc$xrange, ylim = ylim,
                    xlab = if (xticks) xlab else "",
                    ylab = if (empty) "" else ylab,
                    yaxt = if (empty) "n" else "s",
                    bty = if (border | recomb) "o" else "l",
                    xaxt = "n", tcl = -0.3, mgp = c(1.7, 0.5, 0))
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  if (empty) text(mean(loc$xrange), 0.5, blank, adj = 0.5, cex = cex.axis)
  if (xticks) {
    axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, cex.axis = cex.axis,
         mgp = c(1.7, 0.4, 0), tcl = -0.3)
  } else if (!border) {
    axis(1, at = axTicks(1), labels = FALSE, tcl = -0.3)
  }
  # e.g. a legend, drawn before the margins are restored
  if (!is.null(panel_last)) panel_last()
}
