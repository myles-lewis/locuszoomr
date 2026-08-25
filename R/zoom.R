
#' Zoom browser to explore GWAS/eQTL results
#' 
#' Interactive genome browser to explore GWAS/eQTL results using a shiny
#' interface.
#' 
#' @details 
#' This launches a shiny app to explore the GWAS/eQTL results through visualising 
#' the Manhattan plot and exploring regional Manhattan plots of gene loci 
#' through selecting points or searching SNPs/genes. 
#' @param data Dataframe of GWAS results with columns for chromosome, position,
#'   p value and SNP rs IDs. Data.tables are coerced to dataframe.
#' @param ens_db Either a character string which specifies which Ensembl
#'   database package (version 86 and earlier for Homo sapiens) to query for
#'   gene and exon positions (see `ensembldb` Bioconductor package). Or an
#'   `ensembldb` object which can be obtained from the AnnotationHub database.
#'   See the vignette and the `AnnotationHub` Bioconductor package for how to
#'   create this object.
#' @param chrom Determines which column in `data` contains chromosome
#'   information If `NULL` tries to autodetect the column.
#' @param pos Determines which column in `data` contains position information.
#'   If `NULL` tries to autodetect the column.
#' @param p Determines which column in `data` contains SNP p-values. If `NULL`
#'   tries to autodetect the column.
#' @param labs Determines which column in `data` contains SNP rs IDs. If `NULL`
#'   tries to autodetect the column.
#' @param scheme Vector of 3 colours: 1st = normal points, 2nd = colour for
#'   significant points, 3rd = index SNP(s).
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param eqtl_gene Determines which column in `data` contains eQTL genes.
#' @param eqtl_beta Optional column name for beta coefficient to display upward
#'   triangles for positive beta and downward triangles for negative beta
#'   (significant SNPs only).
#' @param eqtl_scheme Colour scheme for eQTL genes.
#' @param add_hover Optional vector of column names in 'data' to add to the
#'   plotly hover text for scatter points.
#' @param mh_points Number of points to display in manhattan plot. Default is
#'   `1e5`.
#' @param recomb Optional `GRanges` class object of recombination data.
#' @param AnnotationDb An `AnnotationDb` gene annotation database, specified
#'   either as a character string or as an `AnnotationDb` class object, used to
#'   obtain expanded gene names. The ensembl database specified in `ens_db` is
#'   queried first. Set to `NULL` to disable this feature.
#' @returns No return value. Opens an interactive shiny window.
#' @importFrom plotly plotlyOutput renderPlotly event_data config plotlyProxy
#' @importFrom plotly plotlyProxyInvoke layout
#' @importFrom shiny fluidPage tabsetPanel tabPanel fluidRow column actionButton 
#' @importFrom shiny icon uiOutput checkboxInput textOutput splitLayout req
#' @importFrom shiny textInput conditionalPanel h5 runApp debounce isolate
#' @importFrom shiny renderUI reactiveValues reactive observe observeEvent radioButtons
#' @importFrom shiny reactiveVal validate need renderText updateTextInput outputOptions
#' @importFrom shinyFeedback useShinyFeedback hideFeedback showFeedback
#' @importFrom shinyWidgets pickerInput pickerOptions dropdown
#' @importFrom shinycssloaders withSpinner
#' @importFrom htmltools tags br
#' @importFrom DT datatable formatSignif
#' @importFrom gtools mixedsort
#' @importFrom stats as.formula setNames
#' @export

zoom <- function(data, ens_db,
                 chrom = NULL, pos = NULL, p = NULL, labs = NULL,
                 scheme = c('royalblue', 'skyblue', 'red'),
                 pcutoff = 5e-8,
                 eqtl_gene = NULL,
                 eqtl_beta = NULL,
                 eqtl_scheme = c("#FF0000", "#00FFFF", "#FF9000", "#0080FF", "#FFFF00",
                                 "#0000FF", "#80DD00", "#8000FF", "#009900", "#FF00FF"),
                 add_hover = NULL,
                 mh_points = 1e5,
                 recomb = NULL,
                 AnnotationDb = "org.Hs.eg.db") {
  data <- data.frame(data)
  # autodetect headings
  dc <- detect_cols(data, chrom, pos, p, labs)
  chrom <- dc$chrom
  pos <- dc$pos
  p <- dc$p
  labs <- dc$labs
  if (is.null(eqtl_gene)) {
    data[, labs] <- unique_snps(data, labs, chrom)
  } else {
    data[, labs] <- unique_snps(data, labs, eqtl_gene)
  }
  
  message("Generating Manhattan plot")
  chr_set <- unique(data[, chrom])
  if (is.character(ens_db)) {
    if (!ens_db %in% (.packages())) {
      stop("Ensembl database not loaded. Try: library(", ens_db, ")",
           call. = FALSE)
    }
    edb <- get(ens_db)
  } else edb <- ens_db
  
  gene_db <- genes(edb, filter = AnnotationFilterList(
    SeqNameFilter(c(1:22, 'X', 'Y'))))
  gene_set <- unique(gene_db$gene_name)
  biotypes <- sort(unique(gene_db$gene_biotype))
  
  # lookup table for full length gene names using org.Hs.eg.db
  fullnames <- fullGeneNames(edb, AnnotationDb)
   
  if (!is.null(eqtl_gene)) {
    eqtl_colour <- eqtl_colours(data[data[, p] < pcutoff, ], chrom, pos,
                                eqtl_gene, eqtl_scheme)
  }
  
  # apply min_p_snp to data for manhat?
  # smallest floating point
  data[which(data[, p] < 5e-324), p] <- 5e-324
  manhat <- manhattan(data, chrom, pos, p, labs, pcutoff = pcutoff,
                      npoints = mh_points)
  yrange <- range(manhat$data$logP, na.rm = TRUE)
  ymax <- yrange[2] + diff(yrange) * 0.05
  
  js <- '$(document).on("keyup", function(e) {
          if(e.key === "Enter") {
            Shiny.onInputChange("enter", Math.random());
          }
        });'
  
  # https://shiny.posit.co/r/articles/build/packaging-javascript/
  
  ui <- fluidPage(
    tags$script(js),
    # 3 plotly scattergl figures gives error "too many active WebGL contexts"
    # see https://plotly.com/python/webgl-vs-svg/
    tags$script(src = "https://unpkg.com/virtual-webgl@1.0.6/src/virtual-webgl.js"),
    useShinyFeedback(),
    tabsetPanel(
      tabPanel("Plot",
               fluidRow(
                 column(11,
                        withSpinner(
                          plotlyOutput("manhattan", width = "85vw", height = "300px"),
                          type = 8, size = 0.7)
                 ),
                 column(1,
                        br(),
                        actionButton("m_zoomin", NULL, icon = icon("magnifying-glass-plus")),
                        actionButton("m_zoomout", NULL, icon = icon("magnifying-glass-minus"))
                 )),
               fluidRow(
                 column(12,
                        uiOutput("ui_chrom")
                 )
               ),
               fluidRow(
                 column(3,
                        checkboxInput("show_chrom", "show chromosome")
                 )
               ),
               fluidRow(
                 column(4,
                        actionButton("left2", NULL, icon = icon("angles-left")),
                        actionButton("left", NULL, icon = icon("angle-left")),
                        actionButton("right", NULL, icon = icon("angle-right")),
                        actionButton("right2", NULL, icon = icon("angles-right")),
                        actionButton("zoomin", NULL, icon = icon("magnifying-glass-plus")),
                        actionButton("zoomout", NULL, icon = icon("magnifying-glass-minus"))
                        ),
                 column(3,
                        textOutput("pos"),
                        align = "centre", style='margin-top:7px;'),
                 column(4,
                        splitLayout(
                          textInput("tex", NULL, placeholder = "chr:start-end, rs or gene",
                                    width = "100%"),
                          actionButton("text_go", NULL, icon = icon("magnifying-glass"),
                                       class = "btn-success"),
                          cellWidths = c("75%", "25%")
                        )),
                 column(1,
                        dropdown(
                          (if (!is.null(recomb)) {
                            checkboxInput("recomb", "show recombination rate", value = TRUE)
                          } else NULL),
                          checkboxInput("alltracks", "show all gene tracks"),
                          pickerInput("biotype", h5("Select gene biotypes"),
                                      choices = biotypes, selected = biotypes,
                                      multiple = TRUE,
                                      options = pickerOptions(actionsBox = TRUE,
                                                              selectedTextFormat = 'count > 1')),
                          (if (!is.null(eqtl_gene)) {
                            uiOutput("ui_genes")
                          } else NULL),
                          right = TRUE, icon = icon("gear")
                        ))
                 ),
                 fluidRow(
                   column(12,
                          plotlyOutput("locus", width = "95vw", height = 600),
                          br(), br()
                   )
                 )
      ),
      tabPanel("Table",
               fluidRow(
                 column(12, br(), DT::dataTableOutput("table"))))
    )
  )
  
  server <- function(input, output, session) {
    
    output$manhattan <- renderPlotly({
      plotly_manhattan(manhat, labs, pcutline = NULL) %>%
        config(displayModeBar = FALSE)
    })
    
    output$ui_chrom <- renderUI({
      req(input$show_chrom, coords$chr)
      fluidRow(
        column(11,
               withSpinner(
                 plotlyOutput("chrom", width = "85vw", height = "220px"),
                 type = 8, size = 0.7)
        ),
        column(1,
               br(),
               actionButton("chr_zoomin", NULL, icon = icon("magnifying-glass-plus")),
               actionButton("chr_zoomout", NULL, icon = icon("magnifying-glass-minus"))
        )
      )
    })
    
    output$chrom <- renderPlotly({
      req(coords$chr)
      chr_manhat <- manhattan(data[which(data[, chrom] == coords$chr), ],
                              chrom, pos, p, labs, pcutoff = pcutoff,
                              npoints = 1e5)
      chr <- suppressWarnings(as.numeric(coords$chr))
      if ((!is.na(chr) && chr %% 2 == 0 || coords$chr == "Y")) {
        scheme[1] <- scheme[2]
      }
      yr <- range(chr_manhat$data$logP, na.rm = TRUE)
      isolate(chr_y$range <- yr)
      isolate(chr_y$max <- yr[2])
      isolate(xr <- coords$xrange)
      
      plotly_manhattan(chr_manhat, labs, scheme = scheme,
                       source = "plotly_chrom") %>%
        layout(margin = list(t = 5),
               shapes = list(
                 list(type = "rect",
                      line = list(width = 1, color = "#00CD00"),
                      x0 = xr[1] / 1e6,
                      x1 = xr[2] / 1e6, y0 = 0, y1 = 1,
                      xref = "x", yref = "paper", layer = "below"))
               ) %>%
        config(displayModeBar = FALSE)
    })
    
    coords <- reactiveValues(chr = NULL, xrange = NULL)
    
    # hide picker at start
    output$coords_ok <- reactive({!is.null(coords$chr)})
    outputOptions(output, "coords_ok", suspendWhenHidden = FALSE)
    
    observe({
      s <- event_data("plotly_click", source = "plotly_manh")
      req(s)
      w <- which(data[, labs] == s$key)
      if (length(w) > 0) {
        coords$chr <- data[w[1], chrom]
        coords$xrange <- data[w[1], pos] + c(-5e5, 5e5)
      }
    })
    
    observe({
      s <- event_data("plotly_click", source = "plotly_chrom")
      req(s)
      w <- which(data[, labs] == s$key)
      if (length(w) > 0) {
        coords$chr <- data[w[1], chrom]
        coords$xrange <- data[w[1], pos] + c(-5e5, 5e5)
      }
    })
    
    # zoom manhattan y axis
    m_ylim <- reactiveValues(max = yrange[2])
    
    observeEvent(input$m_zoomin, {
      m_ylim$max <- pmax(m_ylim$max * 0.88, 5)
      yr <- c(yrange[1], m_ylim$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = "-log<sub>10</sub> P",
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    observeEvent(input$m_zoomout, {
      m_ylim$max <- pmin(m_ylim$max / 0.88, yrange[2])
      yr <- c(yrange[1], m_ylim$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = "-log<sub>10</sub> P",
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    # zoom chrom y axis
    chr_y <- reactiveValues(max = 0, range = c(0, 0))
    
    observeEvent(input$chr_zoomin, {
      chr_y$max <- pmax(chr_y$max * 0.88, 5)
      yr <- c(chr_y$range[1], chr_y$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("chrom", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = "-log<sub>10</sub> P",
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    observeEvent(input$chr_zoomout, {
      chr_y$max <- pmin(chr_y$max / 0.88, chr_y$range[2])
      yr <- c(chr_y$range[1], chr_y$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("chrom", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = "-log<sub>10</sub> P",
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    input_biotype <- reactive({input$biotype}) %>% debounce(2000)
    
    loc <- reactiveValues(i = NULL)
    ntrace <- reactiveVal()
    genes <- reactiveValues(x = NULL)
    
    output$locus <- renderPlotly({
      req(coords$chr %in% chr_set, coords$xrange)
      loc1 <- locus(data = data, xrange = coords$xrange,
                     seqname = coords$chr, ens_db = ens_db,
                     chrom = chrom, pos = pos, p = p, labs = labs)
      validate(need(loc1$data, "Locus contains no SNPs/datapoints"))
      validate(need(nrow(loc1$data) < 1.5e5, "Too many datapoints. Zoom in."))
      if (!is.null(recomb) && input$recomb) {
        loc1 <- link_recomb(loc1, recomb = recomb)
      }
      loc1$TX$fullname <- expandGenes(loc1$TX, fullnames)
      loc$i <- loc1
      
      # req(nrow(loc1$data) > 0)
      if (is.null(eqtl_gene) & is.null(eqtl_beta)) {
        nt <- (sum(loc1$data[, p] < pcutoff) > 1) + 2L
      } else {
        # eqtl
        ind <- loc1$data[, p] < pcutoff
        eqtls <- loc1$data[ind, eqtl_gene]
        genes$x <- unique(eqtls)
        if (!is.null(eqtl_gene) & is.null(eqtl_beta)) {
          # eqtl genes only
          ngene <- length(unique(eqtls))
          nt <- ngene + 1L
        } else if (!is.null(eqtl_gene) & !is.null(eqtl_beta)) {
          # eqtl + beta
          sgn <- sign(loc1$data[ind, eqtl_beta])
          eb <- paste0(eqtls, sgn)
          nt <- length(unique(eb)) + 1L
        }
      }
      nt <- nt + (!is.null(recomb) && input$recomb)
      ntrace(nt)
      
      if (!is.null(eqtl_gene)) {
        genes1 <- unique(eqtls)
        locscheme <- unname(c('grey', eqtl_colour[genes1]))
        if (!is.null(input$select_gene) && input$select_gene != "all") {
          # filter gene
          req(input$select_gene %in% unique(eqtls))  # stops double plot
          ok <- !ind | loc1$data[, eqtl_gene] == input$select_gene
          loc1$data <- loc1$data[ok, ]
          locscheme <- unname(c('grey', eqtl_colour[input$select_gene]))
        }
      } else locscheme <- c('grey', 'dodgerblue', 'red')
      
      isolate(width <- loc_width())
      isolate(biotype <- input_biotype())
      h <- c(0.6, 0.4)
      maxrows <- 8
      if (input$alltracks) {
        cex.width <- 0.7 * par("pin")[1] * 80 / (width - 250)
        tryTX <- mapRow(loc1$TX, xlim = loc1$xrange, cex.text = cex.width,
                        blanks = "show")
        needrow <- pmax(max(tryTX$row, na.rm = TRUE), 8)
        h <- c(360, 20 * needrow + 80)
        maxrows <- NULL
      }
      hideFeedback("tex")
      locus_plotly(loc1, h, filter_gene_biotype = biotype, pcutoff = pcutoff,
                   width = width, eqtl_gene = eqtl_gene, beta = eqtl_beta,
                   add_hover = add_hover, scheme = locscheme, maxrows = maxrows)
    })
    
    output$ui_genes <- renderUI({
      # req(length(genes$x) > 1)
      g <- c("all", genes$x)
      isolate(ig <- input$select_gene)
      if (length(ig) == 0 || !ig %in% genes$x) ig <- "all"
      conditionalPanel("output.coords_ok",
                       radioButtons("select_gene", h5("eQTL genes"), 
                                    choices = g, selected = ig)
      )
    })
    
    outputOptions(output, "ui_genes", suspendWhenHidden = FALSE)
    
    observeEvent(input$left2, {
      dif <- diff(coords$xrange)
      coords$xrange <- pmax(coords$xrange - dif, 0)
    })
    
    observeEvent(input$right2, {
      dif <- diff(coords$xrange)
      coords$xrange <- coords$xrange + dif
    })
    
    observeEvent(input$left, {
      dif <- round(diff(coords$xrange) / 2)
      coords$xrange <- pmax(coords$xrange - dif, 0)
    })
    
    observeEvent(input$right, {
      dif <- round(diff(coords$xrange) / 2)
      coords$xrange <- coords$xrange + dif
    })
    
    observeEvent(input$zoomin, {
      dif <- round(diff(coords$xrange) / 4)
      coords$xrange <- coords$xrange + c(dif, -dif)
    })
    
    observeEvent(input$zoomout, {
      dif <- round(diff(coords$xrange) / 2)
      coords$xrange <- pmax(coords$xrange + c(-dif, dif), 0)
    })
    
    output$pos <- renderText({
      req(coords$chr %in% chr_set, coords$xrange)
      paste0("chr ", coords$chr, ": ", coords$xrange[1], " - ",
             coords$xrange[2])
    })
    
    # parse text box
    observeEvent(c(input$text_go, input$enter), {
      hideFeedback("tex")
      req(input$tex)
      chr <- NULL
      tex <- input$tex
      tex <- gsub(" ", "", tex)
      if (grepl(":", tex) && grepl("-", tex)) {
        # chr & range
        tex <- gsub("chr", "", tex, ignore.case = TRUE)
        ss <- strsplit(tex, ":")[[1]]
        chr <- ss[1]
        xr <- as.integer(strsplit(ss[2], "-")[[1]])
      } else if (grepl(":", tex)) {
        # single position
        tex <- gsub("chr", "", tex, ignore.case = TRUE)
        ss <- strsplit(tex, ":")[[1]]
        chr <- ss[1]
        xr <- as.integer(ss[2]) + c(-5e5, 5e5)
      } else if (any(w <- which(toupper(gene_set) == toupper(input$tex)))) {
        gene <- gene_set[w]
        if (input$tex != gene) updateTextInput(session, "tex", value = gene)
        loc <- genes(edb, filter = AnnotationFilterList(
          GeneNameFilter(gene),
          SeqNameFilter(c(1:22, 'X', 'Y'))))
        if (length(loc) > 1) loc <- loc[1]
        chr <- names(seqlengths(loc))
        m <- mean(c(start(loc), end(loc)))
        xr <- as.integer(c(m - 5e5, m + 5e5))
      } else if (grepl("^rs", input$tex)) {
        w <- which(data[, labs] == input$tex)
        if (length(w) > 0) {
          chr <- data[w[1], chrom]
          xr <- data[w[1], pos] + c(-5e5, 5e5)
        } else {
          showFeedback("tex", "SNP not found")
          return()
        }
      } else {
        if (nchar(tex) > 1) showFeedback("tex", "not found")
        return()
      }
      xr <- as.integer(pmax(xr, 0))
      
      if (chr %in% chr_set) {
        coords$chr <- chr
        if (any(is.na(xr))) {
          showFeedback("tex", "invalid entry")
          return()
        }
        coords$xrange <- xr
        hideFeedback("tex")
      } else {
        showFeedback("tex", "not present")
      }
    })
    
    # Table tab
    output$table <- DT::renderDataTable({
      cols <- colnames(data)[sapply(data, class) == "numeric"]
      datatable(data) %>% formatSignif(cols, digits = 3)
    })
    
    # detect change to x axis range
    observeEvent(event_data("plotly_relayout", source = "plotly_locus"), {
      req(coords$chr %in% chr_set, coords$xrange)
      s <- event_data("plotly_relayout", source = "plotly_locus")
      req(c("xaxis.range[0]", "xaxis.range[1]") %in% names(s))
      xr <- c(s$`xaxis.range[0]`, s$`xaxis.range[1]`)
      coords$xrange <- as.integer(xr * 1e6)
    })
    
    loc_width <- reactiveVal(600)
    
    observe({
      loc_width(session$clientData$output_locus_width)
    })
    
    # redo gene tracks only
    observeEvent(c(loc_width(), input_biotype()), {
      req(loc$i)
      gt <- genetrack_ly(loc$i, filter_gene_biotype = input_biotype(),
                         width = loc_width(), blanks = "show", plot = FALSE)
      req(nrow(gt$TX) != 0)
      TX <- gt$TX
      EX <- gt$EX
      lx <- seg2line(TX$start, TX$end)
      ly <- seg2line(-TX$row, -TX$row)
      hovertext <- paste0(TX$gene_name,
                          TX$fullname,
                          "<br>Gene ID: ", TX$gene_id,
                          "<br>Biotype: ", TX$gene_biotype,
                          "<br>Start: ", TX$start * 1e6,
                          "<br>End: ", TX$end * 1e6)
      ht <- seg2line(hovertext, hovertext)
      exon_col <- exon_border <- "#00008B"
      yref <- if (is.null(recomb) || !input$recomb) "y2" else "y3"
      shapes <- lapply(seq_len(nrow(EX)), function(i) {
        list(type = "rect", fillcolor = exon_col, line = list(color = exon_border,
                                                              width = 0.5),
             x0 = EX$start[i], x1 = EX$end[i], xref = "x",
             y0 = -EX$row[i] - 0.15, y1 = -EX$row[i] + 0.15, yref = yref)
      })
      ok <- !is.na(TX$gene_name2)
      
      plotlyProxy("locus", session) %>%
        plotlyProxyInvoke("restyle",
                          list(x = list(lx), y = list(ly), text = list(ht),
                               hoverinfo = "text"),
                          list(ntrace())) %>%
        plotlyProxyInvoke("update",
                          list(x = list(TX$tx[ok]), y = list(TX$ty[ok]),
                               text = list(TX$gene_name2[ok]), hoverinfo = "none"),
                          list(shapes = shapes),
                          list(ntrace() + 1L))
    })
    
    # chrom highlight
    observeEvent(coords$xrange, {
      req(input$show_chrom, coords$chr)
      plotlyProxy("chrom", session) %>%
        plotlyProxyInvoke("relayout",
                          list(shapes = list(
                            list(type = "rect",
                                 line = list(width = 1, color = "#00CD00"),
                                 x0 = coords$xrange[1] / 1e6,
                                 x1 = coords$xrange[2] / 1e6, y0 = 0, y1 = 1,
                                 xref = "x", yref = "paper", layer = "below")
                          )))
    })
    
  }
  
  runApp(list(ui = ui, server = server))
}


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
  
  ret <- list(data = data, xticks = xticks, pcutoff = pcutoff,
              chrom_list = chrom_list)
  class(ret) <- "manhattan"
  ret
}


plotly_manhattan <- function(obj,
                             labs,
                             scheme = c('royalblue', 'skyblue', 'red'),
                             xlab = "Chromosome",
                             pcutline = NULL,
                             source = "plotly_manh") {
  
  df <- obj$data
  df$col <- as.factor(df$col)
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
                                title = "-log<sub>10</sub> P",
                                ticks = "outside",
                                zeroline = FALSE, showline = TRUE),
                   shapes = hline)
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
    snps[dups] <- make.unique(paste(snps[dups], data[dups, append], sep = "."))
  }
  snps
}


#' @importFrom ensembldb listColumns
fullGeneNames <- function(edb, AnnotationDb) {
  # check ens_db cols for 'description' first
  if (!"description" %in% listColumns(edb) | is.null(AnnotationDb)) return(NULL)
  
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
