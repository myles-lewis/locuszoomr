
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
  
  chrom_list <- as.character(mixedsort(unique(data[, chrom]), na.last = NA))
  data[, chrom] <- factor(data[, chrom], levels = chrom_list)
  
  if (length(chrom_list) == 1) {
    data$genome_pos <- data[, pos]  # single chrom
    lim <- NULL
  } else {
    maxpos <- tapply(data[, pos], data[, chrom], max, na.rm = TRUE)
    maxpos <- maxpos[chrom_list]  # reorder
    minpos <- tapply(data[, pos], data[, chrom], min, na.rm = TRUE)
    minpos <- minpos[chrom_list]  # reorder
    lim <- matrix(c(minpos, maxpos), ncol = 2,
                  dimnames = list(chrom_list, c("min", "max")))
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
  
  # thin points
  if (!is.na(npoints) && nrow(data) > npoints) {
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
  
  data <- data[order(data$genome_pos), ]
  data$col <- ((as.numeric(data[, chrom]) - 1) %% length(chromCols)) + 1
  colScheme <- chromCols
  if (!is.na(sigCol)) {
    data$col[data[, p] < pcutoff] <- length(chromCols) + 1
    colScheme <- c(chromCols, sigCol)
  }
  xticks <- chrom_range <- NULL
  if (length(chrom_list) > 1) {
    chrom_widths <- maxpos - minpos
    xticks <- list(at = chrom_cumsum + 0.5 * chrom_widths, 
                   labels = levels(data[, chrom]))
    chrom_range <- matrix(c(chrom_cumsum, chrom_cumsum + chrom_widths),
                          ncol = 2, dimnames = list(chrom_list, NULL))
  }
  
  data$logP <- -log10(data[, p])
  yrange <- range(data$logP, na.rm = TRUE)
  
  ret <- list(data = data, xticks = xticks, chrom_range = chrom_range,
              chrom_lim = lim, pcutoff = pcutoff, chrom_list = chrom_list,
              labs = labs, yrange = yrange, chrom = chrom, pos = pos)
  class(ret) <- "manhattan"
  ret
}


plotly_manhattan <- function(obj,
                             scheme = c('royalblue', 'skyblue', 'red'),
                             xlab = "Chromosome",
                             ylab = "-log<sub>10</sub> P",
                             pcutline = NULL,
                             xlim = NULL,
                             source = "plotly_manh") {
  
  df <- obj$data
  df$col <- as.factor(df$col)
  labs <- obj$labs
  scheme <- scheme[as.numeric(levels(df$col))]
  if (is.null(obj$xticks)) {
    # single chrom
    df$genome_pos <- df$genome_pos / 1e6
    if (!is.null(xlim)) xlim <- xlim / 1e6
    if (xlab == "Chromosome") xlab <- paste(xlab, obj$chrom_list, "(Mb)")
  }
  xr <- if (is.null(xlim)) range(df$genome_pos, na.rm = TRUE) else xlim
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


#' @importFrom stats complete.cases

complete_data <- function(data, chrom, pos, p) {
  ok <- complete.cases(data[, c(chrom, pos, p)])
  if (!all(ok)) {
    message(sum(!ok), " rows with incomplete data")
    data <- data[which(ok), ]
  }
  data
}


# x is a list of manhattans
align_chrom_lim <- function(x) {
  chr_mins <- lapply(x, function(i) i$chrom_lim[, "min"])
  chr_maxs <- lapply(x, function(i) i$chrom_lim[, "max"])
  full_set <- unique(unlist(lapply(x, function(i) i$chrom_list)))
  t(vapply(full_set, function(i) {
    c(min(unlist(lapply(chr_mins, function(x) x[i])), na.rm = TRUE),
      max(unlist(lapply(chr_maxs, function(x) x[i])), na.rm = TRUE))
  }, numeric(2)))
}


man_lim <- function(x) {
  mrange <- vapply(x, function(i) range(i$data$genome_pos, na.rm = TRUE),
                   numeric(2))
  c(min(mrange[1, ]), max(mrange[2, ]))
}


realign_manhat <- function(m, chrom_lim) {
  minpos <- chrom_lim[m$chrom_list, 1]
  maxpos <- chrom_lim[m$chrom_list, 2]
  chromGap <- sum(maxpos - minpos) / length(m$chrom_list) / 4.15
  
  chrom_cumsum <- c(0, cumsum(maxpos - minpos + chromGap))
  chrom_cumsum2 <- chrom_cumsum - c(minpos, 0)
  chrom_cumsum <- chrom_cumsum[1:length(maxpos)]
  chrom_cumsum2 <- chrom_cumsum2[1:length(maxpos)]
  m$data$genome_pos <- m$data[, m$pos] + chrom_cumsum2[as.numeric(m$data[, m$chrom])]
  
  chrom_widths <- maxpos - minpos
  m$xticks <- list(at = chrom_cumsum + 0.5 * chrom_widths, 
                   labels = levels(m$data[, m$chrom]))
  m$chrom_range <- matrix(c(chrom_cumsum, chrom_cumsum + chrom_widths),
                          ncol = 2, dimnames = list(m$chrom_list, NULL))
  m
}


align_manhats <- function(x) {
  chrom_lim <- align_chrom_lim(x)
  manhats <- lapply(x, realign_manhat, chrom_lim)
  list(manhats = manhats, chrom_lim = chrom_lim, man_lim = man_lim(manhats))
}


# Find minimum p-value per SNP/gene in eQTL results
# @param res Dataframe of LDlink eQTL output
# @param col Column e.g. SNP or gene
# @return Dataframe containing the minimum p-value for each SNP/gene 

min_p_by_col <- function(res, col) {
  ord <- order(res$P_value)
  res_col <- res[ord, col]
  dup <- duplicated(res_col)
  res[ord[!dup], ]
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
