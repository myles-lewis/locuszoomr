
manhattan <- function(data,
                      chrom = NULL, pos = NULL, p = NULL, labs = NULL,
                      pcutoff = 5e-08,
                      chromGap = NULL,
                      chromCols = c('royalblue', 'skyblue'),
                      sigCol = 'red',
                      npoints = 1e6) {
  # autodetect headings
  dc <- detect_cols(data, chrom, pos, p, labs)
  chrom <- dc$chrom
  pos <- dc$pos
  p <- dc$p
  labs <- dc$labs
  
  if (!is.na(npoints) & nrow(data) > npoints) {
    index <- order(data[, p])
    if (npoints <= 1e5) {
      data <- data[index[seq_len(npoints)], ]
    } else {
      # thin points near x axis
      nplotly <- 1e5
      s1 <- seq_len(nplotly)
      s2len <- nrow(data) - nplotly
      s2 <- round(seq_len(npoints - nplotly) * s2len / (npoints - nplotly)) + nplotly
      data <- data[index[unique(c(s1, s2))], ]
    }
  }
  
  data$logP <- -log10(data[, p])
  chrom_list <- mixedsort(unique(data[, chrom]), na.last = NA)
  chrom_list <- as.character(chrom_list)
  
  data[, chrom] <- factor(data[, chrom], levels = chrom_list)
  if (length(chrom_list) == 1) {
    data$genome_pos <- data[, pos]  # single chrom
  } else {
    maxpos <- tapply(data[, pos], data[, chrom], max, na.rm = TRUE)
    maxpos <- maxpos[chrom_list]  # reorder
    minpos <- tapply(data[, pos], data[, chrom], min, na.rm = TRUE)
    minpos <- minpos[chrom_list]  # reorder
    # calculate gap
    if (is.null(chromGap)) {
      chromGap <- sum(maxpos - minpos) / length(chrom_list) / 4.15
    }
    chrom_cumsum <- c(0, cumsum(maxpos - minpos + chromGap))
    chrom_cumsum2 <- chrom_cumsum - c(minpos, 0)
    chrom_cumsum <- chrom_cumsum[1:length(maxpos)]
    chrom_cumsum2 <- chrom_cumsum2[1:length(maxpos)]
    data$genome_pos <- data[, pos] + chrom_cumsum2[as.numeric(data[, chrom])]
  }
  data <- data[order(data$genome_pos), ]
  data$col <- ((as.numeric(data[, chrom]) - 1) %% length(chromCols)) + 1
  colScheme <- chromCols
  if (!is.na(sigCol)) {
    data$col[data[, p] < pcutoff] <- length(chromCols) + 1
    colScheme <- c(chromCols, sigCol)
  }
  if (length(chrom_list) > 1) {
    xticks <- list(at = chrom_cumsum + 0.5 * (maxpos - minpos), 
                   labels = levels(data[, chrom]))
  } else xticks <- NULL
  
  yrange <- range(data$logP, na.rm = TRUE)
  
  ret <- list(data = data, xticks = xticks, pcutoff = pcutoff,
              chrom_list = chrom_list, labs = labs, yrange = yrange)
  class(ret) <- "manhattan"
  ret
}


plotly_manhattan <- function(obj,
                             scheme = c('royalblue', 'skyblue', 'red'),
                             xlab = "Chromosome",
                             ylab = "-log<sub>10</sub> P",
                             pcutline = NULL,
                             source = "plotly_manh") {
  
  df <- obj$data
  df$col <- as.factor(df$col)
  labs <- obj$labs
  scheme <- scheme[as.numeric(levels(df$col))]
  if (is.null(obj$xticks)) {
    # single chrom
    df$genome_pos <- df$genome_pos / 1e6
    if (xlab == "Chromosome") xlab <- paste(xlab, obj$chrom_list, "(Mb)")
  }
  xr <- range(df$genome_pos, na.rm = TRUE)
  xr <- xr + diff(xr) * c(-0.01, 0.01)
  yr <- range(df$logP, na.rm = TRUE)
  yr <- yr + diff(yr) * c(-0.05, 0.05)
  
  hline <- if (!is.null(pcutline)) {
    list(type = "line",
         line = list(width = 1, color = '#AAAAAA', dash = 'dash'),
         x0 = 0, x1 = 1, y0 = -log10(pcutline), y1 = -log10(pcutline),
         xref = "paper", layer = "below")
  } else NULL
  xlayout <- list(range = xr, title = xlab, ticks = "outside",
                  zeroline = FALSE, showline = TRUE, showgrid = FALSE)
  if (!is.null(obj$xticks)) {
    xlayout <- c(xlayout, list(tickvals = obj$xticks$at,
                               ticktext = obj$xticks$labels))
  }
  
  plot_ly(data = df, x = ~genome_pos, y = ~logP,
          color = ~col, colors = scheme,
          marker = list(size = 4, opacity = 0.8),
          text = as.formula(paste0('~', labs)),
          hoverinfo = 'text', key = as.formula(paste0('~', labs)),
          showlegend = FALSE,
          type = "scattergl", mode = "markers",
          source = source) %>%
    plotly::layout(xaxis = xlayout,
                   yaxis = list(range = yr,
                                title = ylab,
                                ticks = "outside",
                                zeroline = FALSE, showline = TRUE),
                   shapes = hline) %>%
    config(displayModeBar = FALSE)
}


seg2line <- function(x, xend) {
  m <- rbind(x, xend, NA)
  as.vector(m)
}


unique_snps <- function(data, labs, append) {
  snps <- data[, labs]
  dups <- which(duplicated(snps))
  if (length(dups) > 0) {
    message("Duplicated SNPs found")
    append_col <- data[dups, append]
    # only append X, Y, or genes (not pure numbers)
    ok <- is.na(suppressWarnings(as.numeric(append_col)))
    snps2 <- snps[dups]
    snps2[ok] <- paste(snps2[ok], append_col[ok], sep = ".")
    snps[dups] <- make.unique(snps2)
  }
  snps
}


#' @importFrom ensembldb listColumns
fullGeneNames <- function(edb, AnnotationDb) {
  # check ens_db cols for 'description' first
  if ("description" %in% listColumns(edb) | is.null(AnnotationDb)) return(NULL)
  
  if (!requireNamespace(AnnotationDb)) {
    stop("Gene annotation database '", AnnotationDb, "' is not installed")
  }
  if (is.character(AnnotationDb)) {
    AnnotationDb <- eval(str2lang(paste0(AnnotationDb, "::", AnnotationDb)))
  }
  alias <- AnnotationDbi::keys(AnnotationDb, "ALIAS")
  suppressMessages(
    AnnotationDbi::mapIds(AnnotationDb, alias,
                          "GENENAME", "ALIAS", multiVals = 'first')
  )
}


# Full genename lookup, returns hovertext
expandGenes <- function(TX, fullnames) {
  genelist <- TX$gene_name
  if ("description" %in% colnames(TX)) {
    # check ensembldb first
    out <- gsub(" \\[[^][]*]", "", TX$description)
  } else {
    if (is.null(fullnames)) return(genelist)
    out <- fullnames[genelist]
  }
  bad <- is.na(out) | out == "NULL" | out == ""
  out <- paste0("<br>", out)
  out[bad] <- ""
  out
}


eqtl_colours <- function(sigdat, chrom, pos, eqtl_gene, eqtl_scheme) {
  message("Setting eQTL colours")
  sigdat <- sigdat[order(sigdat[, chrom], sigdat[, pos]), ]
  eqtl_set <- unique(sigdat[, eqtl_gene])
  message(length(eqtl_set), " eQTL genes")
  setNames(rep_len(eqtl_scheme, length(eqtl_set)), eqtl_set)
}


suppress_warnings <- function(expr, pattern) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(pattern, conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
