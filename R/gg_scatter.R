
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
#' @param scheme Vector of 3 colours if LD is not shown: 1st = normal points,
#'   2nd = colour for significant points, 3rd = index SNP.
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
#' @param recomb_col Colour for recombination rate line if recombination rate
#'   data is present. Set to NA to hide the line. See [link_recomb()] to add
#'   recombination rate data.
#' @param legend_pos Position of legend. Set to `NULL` to hide legend.
#' @param labels Character vector of SNP or genomic feature IDs to label. The
#'   value "index" selects the highest point or index SNP as defined when
#'   [locus()] is called. Set to `NULL` to remove all labels.
#' @param eqtl_gene Optional column name in `loc$data` for colouring eQTL genes.
#' @param beta Optional column name for beta coefficient to display upward
#'   triangles for positive beta and downward triangles for negative beta
#'   (significant SNPs only).
#' @param shape Optional column name in `loc$data` for controlling shapes.
#'   `beta` and `shape` cannot both be set. This column is expected to be a factor.
#' @param shape_values Vector of shape values which match levels of the column
#'   specified by `shape`. This vector is passed to
#'   `ggplot2::scale_shape_manual()` as the argument `values`. See [points()]
#'   for a list of shapes and the numbers they map to.
#' @param ... Optional arguments passed to `geom_text_repel()` to configure
#'   label drawing.
#' @return Returns a ggplot2 plot.
#' @details
#' If recombination rate data is included in the locus object following a call
#' to [link_recomb()], this is plotted as an additional line with a secondary y
#' axis. In the base graphics version the line is placed under the scatter
#' points, but this is not possible with ggplot2 as the secondary y axis data
#' must be plotted on top of the primary scatter point data.
#' 
#' @seealso [locus()] [gg_addgenes()]
#' @examples
#' if(require(EnsDb.Hsapiens.v75)) {
#' data(SLE_gwas_sub)
#' loc <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
#'              ens_db = "EnsDb.Hsapiens.v75")
#' gg_scatter(loc)
#' }
#' @importFrom ggplot2 ggplot geom_point xlim ylim labs theme_classic theme
#'  scale_fill_manual scale_color_manual aes guide_legend element_text
#'  element_blank element_rect unit geom_hline scale_y_continuous sec_axis
#'  geom_line scale_shape_manual guides
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr bind_rows
#' @importFrom rlang .data
#' @importFrom zoo na.approx
#' @export
#' 
gg_scatter <- function(loc,
                       index_snp = loc$index_snp,
                       pcutoff = 5e-08,
                       scheme = c('grey', 'dodgerblue', 'red'),
                       size = 2,
                       cex.axis = 1,
                       cex.lab = 1,
                       xlab = NULL,
                       ylab = NULL,
                       yzero = (loc$yvar == "logP"),
                       xticks = TRUE,
                       border = FALSE,
                       showLD = TRUE,
                       LD_scheme = c('grey', 'royalblue', 'cyan2', 'green3', 
                                     'orange', 'red', 'purple'),
                       recomb_col = "blue",
                       legend_pos = 'topleft',
                       labels = NULL,
                       eqtl_gene = NULL,
                       beta = NULL,
                       shape = NULL,
                       shape_values = c(21, 24, 25), ...) {
  if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
  if (is.null(loc$data)) stop("No data points, only gene tracks")
  
  .call <- match.call()
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
      data$bg[data[, loc$labs] %in% index_snp] <- 7L
      data$bg <- factor(data$bg, levels = 1:7)
      data <- data[order(data$bg), ]
      scheme <- rep_len(LD_scheme, 7)
      if (is.null(index_snp)) {
        scheme <- scheme[1:6]
        data$bg <- factor(data$bg, levels = 1:6)
      }
    } else if (!is.null(eqtl_gene)) {
      # eqtl gene colours
      bg <- data[, eqtl_gene]
      bg[data[, loc$p] > pcutoff] <- "ns"
      bg <- relevel(factor(bg, levels = unique(bg)), "ns")
      if (is.null(.call$scheme)) scheme <- eqtl_scheme(nlevels(bg))
      data$bg <- bg
    } else {
      data$bg <- scheme[1]
      if (loc$yvar == "logP") data$bg[data[, loc$p] < pcutoff] <- scheme[2]
      data$bg[data[, loc$labs] %in% index_snp] <- scheme[3]
      data$bg <- factor(data$bg, levels = scheme)
    }
  }
  
  # scatter plot
  if (!"col" %in% colnames(data)) data$col <- "black"
  data$col <- as.factor(data$col)
  
  # shapes
  if (!is.null(shape)) {
    if (!is.null(beta)) stop("cannot set both `shape` and `beta`")
    if (!shape %in% colnames(data)) stop("incorrect column name for `shape`")
    shape_breaks <- shape_labels <- levels(data[, shape])
  }
  if (!is.null(beta)) {
    # beta symbols
    data[, beta] <- signif(data[, beta], 3)
    symbol <- as.character(sign(data[, beta]))
    ind <- data[, loc$p] > pcutoff
    symbol[ind] <- "ns"
    data$.beta <- factor(symbol, levels = c("ns", "1", "-1"),
                         labels = c("ns", "up", "down"))
    shape <- ".beta"
    shape_breaks <- c("ns", "up", "down")
    shape_labels <- c("ns", expression({beta > 0}), expression({beta < 0}))
  }
  
  # legend
  legend.justification <- NULL
  legend_labels <- legend_title <- NULL
  legend.position <- "none"
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
      legend_title <- expression({r^2})
      legend_labels <- rev(c("Index SNP", "0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6",
                             "0.2 - 0.4", "0.0 - 0.2", "NA"))
      if (is.null(index_snp)) legend_labels <- legend_labels[1:6]
    } else if (!is.null(eqtl_gene)) {
      legend_labels <- levels(bg)
    } else if (is.null(beta) & is.null(shape)) legend.position <- "none"
  }
  
  yrange <- range(data[, loc$yvar], na.rm = TRUE)
  if (yzero) yrange[1] <- min(c(0, yrange[1]))
  ycut <- -log10(pcutoff)
  
  # recombination line
  recomb <- !is.null(loc$recomb) & !is.na(recomb_col)
  if (recomb) {
    df <- loc$recomb[, c("start", "value")]
    colnames(df) <- c(loc$pos, "recomb")
    data <- dplyr::bind_rows(data, df)
    data <- data[order(data[, loc$pos]), ]
    data$recomb <- zoo::na.approx(data$recomb, data[, loc$pos], na.rm = FALSE)
  }
  data[, loc$pos] <- data[, loc$pos] / 1e6

  # add labels
  if (!is.null(labels)) {
    i <- grep("index", labels, ignore.case = TRUE)
    if (length(i) > 0) {
      if (length(index_snp) == 1) {
        labels[i] <- index_snp
      } else {
        labels <- labels[-i]
        labels <- c(index_snp, labels)
      }
    }
    text_label_ind <- match(labels, data[, loc$labs])
    if (any(is.na(text_label_ind))) {
      message("label ", paste(labels[is.na(text_label_ind)], collapse = ", "),
              " not found")
    }
  }
  ind <- data[, loc$labs] %in% index_snp

  if (!recomb) {
    if (is.null(shape)) {
      # standard plot
      p <- ggplot(data[!ind, ], aes(x = .data[[loc$pos]], y = .data[[loc$yvar]],
                                    color = .data$col, fill = .data$bg)) +
        (if (loc$yvar == "logP" & !is.null(pcutoff) &
             ycut >= yrange[1] & ycut <= yrange[2]) {
          geom_hline(yintercept = ycut,
                     colour = "grey", linetype = "dashed")
        }) +
        geom_point(shape = 21, size = size) +
        # index SNP
        (if (any(ind)) {
          geom_point(data = data[ind, ],
                     aes(y = .data[[loc$yvar]], color = .data$col,
                         fill = .data$bg),
                     shape = 23, size = size)
        })
    } else {
      # shapes or beta triangles
      p <- ggplot(data, aes(x = .data[[loc$pos]], y = .data[[loc$yvar]],
                            color = .data$col, fill = .data$bg,
                            shape = .data[[shape]])) +
        (if (loc$yvar == "logP" & !is.null(pcutoff) &
             ycut >= yrange[1] & ycut <= yrange[2]) {
          geom_hline(yintercept = ycut,
                     colour = "grey", linetype = "dashed")
        }) +
        geom_point(size = size) +
        scale_shape_manual(values = shape_values, name = NULL,
                           breaks = shape_breaks,
                           labels = shape_labels) +
        (if (showLD & hasLD) {
          guides(fill = guide_legend(override.aes = list(shape = 21),
                                     reverse = TRUE, order = 1))
        } else {
          guides(fill = "none")
        })
    }
    p <- p +
      scale_fill_manual(breaks = levels(data$bg), values = scheme,
                        guide = guide_legend(reverse = TRUE),
                        labels = legend_labels, name = legend_title) +
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
            legend.key.size = unit(0.9, 'lines'),
            legend.spacing.y = unit(0, 'lines')) +
      if (!xticks) theme(axis.text.x=element_blank(),
                         axis.ticks.x=element_blank())
  } else {
    # recombination plot with dual y axis
    ymult <- 100 / diff(yrange)
    if (is.null(shape)) {
      # standard plot
      p <- ggplot(data[!ind, ], aes(x = .data[[loc$pos]])) +
        (if (loc$yvar == "logP" & !is.null(pcutoff) &
             ycut >= yrange[1] & ycut <= yrange[2]) {
          geom_hline(yintercept = ycut,
                     colour = "grey", linetype = "dashed")
        }) +
        geom_point(aes(y = .data[[loc$yvar]], color = .data$col,
                       fill = .data$bg), shape = 21, size = size, na.rm = TRUE) +
        # index SNP
        (if (any(ind)) {
          geom_point(data = data[ind, ],
                     aes(y = .data[[loc$yvar]], color = .data$col,
                         fill = .data$bg), shape = 23, size = size, na.rm = TRUE)
        })
    } else {
      # shapes or beta triangles
      p <- ggplot(data, aes(x = .data[[loc$pos]])) +
        (if (loc$yvar == "logP" & !is.null(pcutoff) &
             ycut >= yrange[1] & ycut <= yrange[2]) {
          geom_hline(yintercept = ycut,
                     colour = "grey", linetype = "dashed")
        }) +
        geom_point(aes(y = .data[[loc$yvar]], color = .data$col, fill = .data$bg,
                       shape = .data[[shape]]), size = size, na.rm = TRUE) +
        scale_shape_manual(values = shape_values, name = NULL,
                           breaks = shape_breaks,
                           labels = shape_labels) +
        (if (showLD & hasLD) {
          guides(fill = guide_legend(override.aes = list(shape = 21),
                                     reverse = TRUE, order = 1))
        } else {
          guides(fill = "none")
        })
    }
    p <- p +
      scale_fill_manual(breaks = levels(data$bg), values = scheme,
                        guide = guide_legend(reverse = TRUE),
                        labels = legend_labels, name = legend_title) +
      scale_color_manual(breaks = levels(data$col), values = levels(data$col),
                         guide = "none") +
      geom_line(aes(y = .data$recomb / ymult + yrange[1]), color = recomb_col,
                na.rm = TRUE) +
      scale_y_continuous(name = ylab,
                         sec.axis = sec_axis(~(. - yrange[1]) * ymult,
                                             name = "Recombination rate (%)")) +
      xlim(loc$xrange[1] / 1e6, loc$xrange[2] / 1e6) +
      xlab(xlab) +
      theme_classic() +
      theme(axis.text = element_text(colour = "black", size = 10 * cex.axis),
            axis.title = element_text(size = 10 * cex.lab),
            legend.justification = legend.justification,
            legend.position = legend.position,
            legend.title.align = 0.5,
            legend.text.align = 0,
            legend.key.size = unit(0.9, 'lines'),
            legend.spacing.y = unit(0, 'lines')) +
      if (!xticks) theme(axis.text.x=element_blank(),
                         axis.ticks.x=element_blank())
  }
  
  if (!is.null(labels)) {
    p <- p +
      geom_text_repel(data = data[text_label_ind, ],
                      mapping = aes(x = .data[[loc$pos]], y = .data[[loc$yvar]],
                                    label = .data[[loc$labs]]),
                      point.size = size, ...)
  }

  if (border | recomb) {
    p <- p + theme(panel.border = element_rect(colour = "black", fill = NA))
  }
  p
}

