
#' Locus scatter plot using ggplot2
#'
#' Produces a scatter plot from a 'locus' class object. Intended for use with
#' [set_layers()].
#'
#' @param x Object of class 'locus' to use for plot. See [locus].
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param chromCol Colour for normal points if `LD` is `FALSE` when the locus
#'   object is made.
#' @param sigCol Colour for significant points if `LD` is `FALSE`.
#' @param size Specifies size for points.
#' @param cex.axis Specifies font size for axis numbering.
#' @param cex.lab Specifies font size for axis titles.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param xticks Logical whether x axis numbers and axis title are plotted.
#' @param border Logical whether a bounding box is plotted around the plot.
#' @param showLD Logical whether to show LD with colours
#' @param LD_scheme Vector of colours for plotting LD. The first colour is for SNPs
#'   which lack LD information. The next 5 colours are for r2 or D' LD results
#'   ranging from 0 to 1 in intervals of 0.2. The final colour is for the index
#'   SNP.
#' @param legend_pos Position of legend. See [legend()]. Set to `NULL` to hide
#'   legend.
#' @return Returns a ggplot2 plot.
#' @seealso [locus()] [set_layers()]
#' @importFrom ggplot2 ggplot geom_point labs theme_classic theme
#'  scale_fill_manual scale_color_manual aes guide_legend element_text
#'  element_blank element_rect unit
#' @importFrom rlang .data
#' @export
#' 
gg_scatter <- function(x,
                       pcutoff = 5e-08,
                       chromCol = 'royalblue',
                       sigCol = 'red',
                       size = 2,
                       cex.axis = 1,
                       cex.lab = 1,
                       xlab = NULL,
                       ylab = expression("-log"[10] ~ "P"),
                       xticks = FALSE,
                       border = FALSE,
                       showLD = TRUE,
                       LD_scheme = c('grey', 'royalblue', 'cyan2', 'green3', 
                                     'orange', 'red', 'purple'),
                       legend_pos = 'topleft') {
  if (!inherits(x, "locus")) stop("Object of class 'locus' required")
  data <- x$data
  if (is.null(xlab) & xticks) xlab <- paste("Chromosome", x$seqname, "(Mb)")
  hasLD <- "ld" %in% colnames(data)
  if (!"bg" %in% colnames(data)) {
    if (showLD & hasLD) {
      data$bg <- cut(data$ld, -1:6/5, labels = FALSE)
      data$bg[is.na(data$bg)] <- 1L
      data$bg[which.max(data$logP)] <- 7L
      data$bg <- factor(data$bg)
      data <- data[order(data$bg), ]
      scheme <- rep_len(LD_scheme, 7)
    } else {
      data$bg <- chromCol
      data$bg[data[, x$p] < pcutoff] <- sigCol
      scheme <- c(chromCol, sigCol)
    }
  }
  
  # scatter plot
  data[, x$pos] <- data[, x$pos] / 1e6
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
      legend.position <- c(1, 0.99)
    } else {
      legend.position <- legend_pos
    }
    if (showLD & hasLD) {
      legend_labels <- rev(c('Index SNP',
                             expression({0.8 < r^2} <= "1.0"),
                             expression({0.6 < r^2} <= 0.8),
                             expression({0.4 < r^2} <= 0.6),
                             expression({0.2 < r^2} <= 0.4),
                             expression({"0.0" < r^2} <= 0.2),
                             expression("No" ~ r^2 ~ "data")))
    } else legend.position = "none"
  } else legend.position = "none"
  
  p <- ggplot(data, aes(x = .data[[x$pos]], y = .data$logP, color = .data$col,
                   fill = .data$bg)) +
    geom_point(shape = 21, size = size) +
    scale_fill_manual(breaks = levels(data$bg), values = scheme,
                      guide = guide_legend(reverse = TRUE),
                      labels = legend_labels) +
    scale_color_manual(breaks = levels(data$col), values = levels(data$col),
                       guide = "none") +
    # scale_shape_manual(breaks = levels(data$pch), values = levels(data$pch)) +
    labs(x = xlab, y = ylab) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 10 * cex.axis),
          axis.title = element_text(size = 10 * cex.lab),
          legend.title = element_blank(),
          legend.justification = legend.justification,
          legend.position = legend.position,
          legend.text.align = 0,
          legend.key.size = unit(1, 'lines'),
          legend.spacing.y = unit(0, 'lines')) +
    if (!xticks) theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
  if (border) p <- p + theme(panel.border = element_rect(colour = "black", fill = NA))
  p
}

