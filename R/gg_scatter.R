
#' Locus scatter plot using ggplot2
#'
#' Produces a scatter plot from a 'locus' class object (without gene tracks).
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param index_snp Specifies index SNP to be shown in a different colour and
#'   symbol. Defaults to the SNP with the lowest p-value. Set to `NULL` to not
#'   show this.
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param scheme Vector of 3 colors if LD is not shown: 1st = normal points, 2nd
#'   = colour for significant points, 3rd = index SNP.
#' @param size Specifies size for points.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.lab Specifies font size for axis titles.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param yzero Logical whether to force y axis limit to include y=0.
#' @param xticks Logical whether x axis numbers and axis title are plotted.
#' @param border Logical whether a bounding box is plotted around the plot.
#' @param showLD Logical whether to show LD with colours
#' @param LD_scheme Vector of colours for plotting LD. The first colour is for SNPs
#'   which lack LD information. The next 5 colours are for r2 or D' LD results
#'   ranging from 0 to 1 in intervals of 0.2. The final colour is for the index
#'   SNP.
#' @param legend_pos Position of legend. Set to `NULL` to hide legend.
#' @return Returns a ggplot2 plot.
#' @seealso [locus()] [gg_addgenes()]
#' @examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' library(EnsDb.Hsapiens.v75)
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
#'              ens_db = "EnsDb.Hsapiens.v75")
#' gg_scatter(loc)
#' }
#' @importFrom ggplot2 ggplot geom_point xlim ylim labs theme_classic theme
#'  scale_fill_manual scale_color_manual aes guide_legend element_text
#'  element_blank element_rect unit geom_hline
#' @importFrom rlang .data
#' @export
#' 
gg_scatter <- function(loc,
                       index_snp = loc$index_snp,
                       pcutoff = 5e-08,
                       scheme = c('royalblue', 'red', 'purple'),
                       size = 2,
                       cex.axis = 1,
                       cex.lab = 1,
                       xlab = NULL,
                       ylab = NULL,
                       yzero = TRUE,
                       xticks = TRUE,
                       border = FALSE,
                       showLD = TRUE,
                       LD_scheme = c('grey', 'royalblue', 'cyan2', 'green3', 
                                     'orange', 'red', 'purple'),
                       legend_pos = 'topleft') {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  data <- loc$data
  if (is.null(xlab) & xticks) xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  if (is.null(ylab)) {
    ylab <- if (loc$yvar == "logP") expression("-log"[10] ~ "P") else loc$yvar
  }
  hasLD <- "ld" %in% colnames(data)
  if (!"bg" %in% colnames(data)) {
    if (showLD & hasLD) {
      data$bg <- cut(data$ld, -1:6/5, labels = FALSE)
      data$bg[is.na(data$bg)] <- 1L
      data$bg[data[, loc$labs] == index_snp] <- 7L
      data$bg <- factor(data$bg)
      data <- data[order(data$bg), ]
      scheme <- rep_len(LD_scheme, 7)
      if (is.null(index_snp)) scheme <- scheme[1:6]
    } else {
      data$bg <- scheme[1]
      if (loc$yvar == "logP") data$bg[data[, loc$p] < pcutoff] <- scheme[2]
      data$bg[data[, loc$labs] == index_snp] <- scheme[3]
      data$bg <- factor(data$bg, levels = scheme)
    }
  }
  
  # scatter plot
  data[, loc$pos] <- data[, loc$pos] / 1e6
  if (!"col" %in% colnames(data)) data$col <- "black"
  data$col <- as.factor(data$col)
  # if (!"pch" %in% colnames(data)) data$pch <- 21
  # data$pch <- as.factor(data$pch)
  
  legend.justification <- NULL
  legend_labels <- NULL
  if (!is.null(legend_pos)) {
    if (legend_pos == "topleft") {
      legend.justification <- c(0, 1)
      legend.position <- c(0.01, 0.99)
    } else if (legend_pos == "topright") {
      legend.justification <- c(1, 1)
      legend.position <- c(0.99, 0.99)
    } else {
      legend.position <- legend_pos
    }
    if (showLD & hasLD) {
      legend_labels <- rev(c("Index SNP", "0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4",
                             "0.0 - 0.2", "NA"))
      if (is.null(index_snp)) legend_labels <- legend_labels[1:6]
    } else legend.position = "none"
  } else legend.position = "none"
  yrange <- range(data[, loc$yvar], na.rm = TRUE)
  yrange[1] <- if (yzero) min(c(0, yrange[1]))
  ycut <- -log10(pcutoff)
  
  p <- ggplot(data, aes(x = .data[[loc$pos]], y = .data[[loc$yvar]], color = .data$col,
                   fill = .data$bg)) +
    (if (loc$yvar == "logP" & !is.null(pcutoff) & ycut >= yrange[1] & ycut <= yrange[2]) {
      geom_hline(yintercept = ycut,
                 colour = "grey", linetype = "dashed")
    }
    ) +
    geom_point(shape = 21, size = size) +
    scale_fill_manual(breaks = levels(data$bg), values = scheme,
                      guide = guide_legend(reverse = TRUE),
                      labels = legend_labels, name = expression({r^2})) +
    scale_color_manual(breaks = levels(data$col), values = levels(data$col),
                       guide = "none") +
    # scale_shape_manual(breaks = levels(data$pch), values = levels(data$pch)) +
    xlim(loc$xrange[1] / 1e6, loc$xrange[2] / 1e6) + ylim(yrange[1], NA) +
    labs(x = xlab, y = ylab) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 10 * cex.axis),
          axis.title = element_text(size = 10 * cex.lab),
          legend.justification = legend.justification,
          legend.position = legend.position,
          legend.title.align = 0.5,
          legend.text.align = 0,
          legend.key.size = unit(1, 'lines'),
          legend.spacing.y = unit(0, 'lines')) +
    if (!xticks) theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
  if (border) p <- p + theme(panel.border = element_rect(colour = "black", fill = NA))
  p
}

