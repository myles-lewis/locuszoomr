
#' Zoom browser to explore GWAS/eQTL results
#' 
#' Interactive genome browser to explore GWAS/eQTL results using a shiny
#' interface.
#' 
#' @details 
#' This launches a shiny app to explore the GWAS/eQTL results through visualising 
#' the Manhattan plot and exploring regional Manhattan plots of gene loci 
#' through selecting points or searching SNPs/genes.
#' 
#' Linkage disequilibrium data can be pulled from LDlink API within the app
#' using the 'Get LD' button in the Settings dropdown. These buttons only appear
#' once a user provides a LDlink API token via the `ld_token` argument. Each API
#' request takes around 5-10 secs. 'Get LD' pins the current index SNP as the
#' reference variant and colours points by r^2 with it. The reference stays
#' pinned while you pan and zoom, so the colouring keeps its meaning and repeat
#' queries are served from the `memoise` cache rather than the API.
#' @param data Dataframe of GWAS results with columns for chromosome, position,
#'   p value and SNP rs IDs. Data.tables are coerced to dataframe.
#' @param data2 Optional 2nd dataframe of GWAS results for comparison.
#' @param traits Vector of trait names for identifying `data` and `data2`
#'   datasets.
#' @param ens_db Either a character string which specifies which Ensembl
#'   database package (version 86 and earlier for Homo sapiens) to query for
#'   gene and exon positions (see `ensembldb` Bioconductor package). Or an
#'   `ensembldb` object which can be obtained from the AnnotationHub database.
#'   See the vignette and the `AnnotationHub` Bioconductor package for how to
#'   create this object.
#' @param chrom Determines which column in `data` contains chromosome
#'   information. If `NULL` or `NA` tries to autodetect the column. If `data2`
#'   is provided, this is a vector where 1st element refers to `data` and 2nd
#'   element refers to `data2`.
#' @param pos Determines which column in `data` (and optionally `data2`)
#'   contains position information. See `chrom`.
#' @param p Determines which column in `data` (and optionally `data2`) contains
#'   SNP p-values. See `chrom`.
#' @param labs Determines which column in `data` (and optionally `data2`)
#'   contains SNP rs IDs. See `chrom`.
#' @param scheme Vector of 3 colours for main Manhattan plot: 1st, 2nd colours
#'   for alternating chromosomes, 3rd colour for significant points.
#' @param scheme2 Vector of colours for 2nd Manhattan plot.
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
#' @param ld_token Personal access token for the LDlink API, available from
#'   <https://ldlink.nih.gov/?tab=apiaccess>. See `LDlinkR` package
#'   documentation and [link_LD()]. When empty the LD controls are hidden. LD
#'   information can only be requested if `eqtl_gene` is left as `NULL`. LD is
#'   fetched on demand, not automatically, by pressing the "Get LD" button.
#' @param ld_pop 1000 Genomes population used for LD. Defaults to `"EUR"`. See
#'   `LDlinkR::LDproxy()` for the available codes.
#' @param seq_filter Vector of acceptable chromosomes. Used to restrict queries
#'   to standard chromosome assembly.
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
#' @importFrom shiny showNotification removeNotification HTML downloadButton downloadHandler
#' @importFrom shinyFeedback useShinyFeedback hideFeedback showFeedback
#' @importFrom shinyWidgets pickerInput pickerOptions dropdown
#' @importFrom shinycssloaders withSpinner
#' @importFrom htmltools tags br
#' @importFrom DT datatable formatSignif
#' @importFrom gtools mixedsort
#' @importFrom stats as.formula setNames
#' @importFrom grDevices dev.off pdf
#' @export

zoom <- function(data, ens_db,
                 chrom = NULL, pos = NULL, p = NULL, labs = NULL,
                 data2 = NULL,
                 traits = NULL,
                 scheme = c('royalblue', 'skyblue', 'red'),
                 scheme2 = c("#33a02c", "#b2df8a", "purple"),
                 pcutoff = 5e-8,
                 eqtl_gene = NULL,
                 eqtl_beta = NULL,
                 eqtl_scheme = c("#FF0000", "#00FFFF", "#FF9000", "#0080FF", "#FFFF00",
                                 "#0000FF", "#80DD00", "#8000FF", "#009900", "#FF00FF"),
                 add_hover = NULL,
                 mh_points = 1e5,
                 recomb = NULL,
                 ld_token = Sys.getenv("LDLINK_TOKEN"),
                 ld_pop = "EUR",
                 seq_filter = c(1:22, 'X', 'Y'),
                 AnnotationDb = "org.Hs.eg.db") {
  dat_name <- deparse(substitute(data))
  data <- data.frame(data)
  # autodetect headings
  dc <- detect_cols(data, chrom[1], pos[1], p[1], labs[1])
  
  # autodetect headings gwas2
  dc2 <- NULL
  man2 <- FALSE
  if (!is.null(data2)) {
    dat2_name <- deparse(substitute(data2))
    data2 <- data.frame(data2)
    if (!is.null(chrom)) chrom <- rep_len(chrom, 2)
    if (!is.null(pos)) pos <- rep_len(pos, 2)
    if (!is.null(p)) p <- rep_len(p, 2)
    if (!is.null(labs)) labs <- rep_len(labs, 2)
    dc2 <- detect_cols(data2, chrom[2], pos[2], p[2], labs[2])
    if (is.null(traits)) traits <- c(dat_name, dat2_name)
    man2 <- TRUE
  }
  
  chrom <- c(dc$chrom, dc2$chrom)
  pos <- c(dc$pos, dc2$pos)
  p <- c(dc$p, dc2$p)
  labs <- c(dc$labs, dc2$labs)
  
  message("Generating Manhattan plot", (if (man2) " 1"))
  if (is.null(eqtl_gene)) {
    data[, labs[1]] <- unique_snps(data, labs[1], chrom[1])
  } else {
    data[, labs[1]] <- unique_snps(data, labs[1], eqtl_gene)
  }
  # currently eqtl_gene can only apply to data1
  
  chr_set <- list()
  chr_set[[1]] <- unique(data[, chrom[1]])
  if (is.character(ens_db)) {
    if (!ens_db %in% (.packages())) {
      stop("Ensembl database not loaded. Try: library(", ens_db, ")",
           call. = FALSE)
    }
    edb <- get(ens_db)
  } else edb <- ens_db
  
  gene_db <- genes(edb, filter = AnnotationFilterList(
    SeqNameFilter(seq_filter)))
  gene_set <- unique(gene_db$gene_name)
  biotypes <- sort(unique(gene_db$gene_biotype))
  
  # lookup table for full length gene names using org.Hs.eg.db
  fullnames <- fullGeneNames(edb, AnnotationDb)
   
  if (!is.null(eqtl_gene)) {
    eqtl_colour <- eqtl_colours(data[data[, p] < pcutoff, ], chrom, pos,
                                eqtl_gene, eqtl_scheme)
  }
  
  show_ld <- nzchar(ld_token) && is.null(eqtl_gene)
  
  # apply min_p_snp to data for manhat?
  # smallest floating point
  data[which(data[, p[1]] < 5e-324), p[1]] <- 5e-324
  manhat <- manhattan(data, chrom[1], pos[1], p[1], labs[1], pcutoff = pcutoff,
                      npoints = mh_points)
  man_ylab <- "-log<sub>10</sub> P"
    
  if (man2) {
    message("Generating Manhattan plot 2")
    data2[, labs[2]] <- unique_snps(data2, labs[2], chrom[2])
    data2[which(data2[, p[2]] < 5e-324), p[2]] <- 5e-324
    chr_set[[2]] <- unique(data2[, chrom[2]])
    manhat2 <- manhattan(data2, chrom[2], pos[2], p[2], labs[2], pcutoff = pcutoff,
                         npoints = mh_points)
    man_ylab <- paste(traits, man_ylab)
  }
  
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
    tags$head(
      tags$style(HTML("
      .shiny-output-error-validation {
        display: flex;
        align-items: center;
      }"))
    ),
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
               (if (man2) {
                 fluidRow(
                   column(11,
                          withSpinner(
                            plotlyOutput("manhattan2", width = "85vw", height = "300px"),
                            type = 8, size = 0.7)
                   ),
                   column(1,
                          br(),
                          actionButton("m_zoomin2", NULL, icon = icon("magnifying-glass-plus")),
                          actionButton("m_zoomout2", NULL, icon = icon("magnifying-glass-minus"))
                   )
                 )
               }),
               fluidRow(
                 column(12,
                        conditionalPanel('input.show_chrom & output.coords_ok',
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
                        )
                 )
               ),
               (if (man2) {
                 fluidRow(
                   column(12,
                          conditionalPanel('input.show_chrom & output.coords_ok',
                                           fluidRow(
                                             column(11,
                                                    withSpinner(
                                                      plotlyOutput("chrom2", width = "85vw", height = "220px"),
                                                      type = 8, size = 0.7)
                                             ),
                                             column(1,
                                                    br(),
                                                    actionButton("chr_zoomin2", NULL, icon = icon("magnifying-glass-plus")),
                                                    actionButton("chr_zoomout2", NULL, icon = icon("magnifying-glass-minus"))
                                             )
                                           )
                          )
                   )
                 )
               }),
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
                        actionButton("zoomout", NULL, icon = icon("magnifying-glass-minus")),
                        uiOutput("save_ui", inline = T)
                        ),
                 column(3,
                        textOutput("pos"),
                        textOutput("ld_status"),
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
                          (if (show_ld) {
                            fluidRow(
                              column(12,
                                     h5("Linkage disequilibrium"),
                                     actionButton("ld_get", "Get LD", icon = icon("circle-nodes"),
                                                  class = "btn-primary btn-sm"),
                                     actionButton("ld_clear", "Clear",
                                                  class = "btn-default btn-sm")
                              ))
                          } else NULL),
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
                          plotlyOutput("locus", width = "95vw", height = 624)
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
      plotly_manhattan(manhat, ylab = man_ylab[1], pcutline = NULL,
                       scheme = scheme)
    })
    
    output$chrom <- renderPlotly({
      req(coords$chr)
      validate(need(coords$chr %in% chr_set[[1]], "No data for this chromosome"))
      chr_manhat <- manhattan(data[which(data[, chrom[1]] == coords$chr), ],
                              chrom[1], pos[1], p[1], labs[1], pcutoff = pcutoff,
                              npoints = 1e5)
      chr <- suppressWarnings(as.numeric(coords$chr))
      if ((!is.na(chr) && chr %% 2 == 0 || coords$chr == "Y")) {
        scheme[1] <- scheme[2]
      }
      yr <- range(chr_manhat$data$logP, na.rm = TRUE)
      isolate(chr_y$range <- yr)
      isolate(chr_y$max <- yr[2])
      isolate(xr <- coords$xrange)
      
      plotly_manhattan(chr_manhat, scheme = scheme, ylab = man_ylab[1],
                       source = "plotly_chrom") %>%
        layout(margin = list(t = 5),
               shapes = list(
                 list(type = "rect",
                      line = list(width = 1, color = "#00CD00"),
                      x0 = xr[1] / 1e6,
                      x1 = xr[2] / 1e6, y0 = 0, y1 = 1,
                      xref = "x", yref = "paper", layer = "below")))
    })
    
    coords <- reactiveValues(chr = NULL, xrange = NULL)
    
    # hide picker at start
    output$coords_ok <- reactive({!is.null(coords$chr)})
    outputOptions(output, "coords_ok", suspendWhenHidden = FALSE)
    
    observe({
      s <- event_data("plotly_click", source = "plotly_manh")
      req(s)
      w <- which(data[, labs[1]] == s$key)
      if (length(w) > 0) {
        coords$chr <- data[w[1], chrom[1]]
        xr <- data[w[1], pos[1]] + c(-5e5, 5e5)
        if (xr[1] < 0) xr <- c(0, 1e6)
        coords$xrange <- xr
      }
    })
    
    observe({
      s <- event_data("plotly_click", source = "plotly_chrom")
      req(s)
      w <- which(data[, labs[1]] == s$key)
      if (length(w) > 0) {
        coords$chr <- data[w[1], chrom[1]]
        xr <- data[w[1], pos[1]] + c(-5e5, 5e5)
        if (xr[1] < 0) xr <- c(0, 1e6)
        coords$xrange <- xr
      }
    })
    
    # zoom manhattan y axis
    m_ylim <- reactiveValues(max = manhat$yrange[2])
    
    observeEvent(input$m_zoomin, {
      m_ylim$max <- pmax(m_ylim$max * 0.88, 5)
      yr <- c(manhat$yrange[1], m_ylim$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = man_ylab[1],
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    observeEvent(input$m_zoomout, {
      m_ylim$max <- pmin(m_ylim$max / 0.88, manhat$yrange[2])
      yr <- c(manhat$yrange[1], m_ylim$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = man_ylab[1],
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
                                            title = man_ylab[1],
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
                                            title = man_ylab[1],
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    # 2nd gwas
    if (man2) {
      output$manhattan2 <- renderPlotly({
        req(man2)
        plotly_manhattan(manhat2, ylab = man_ylab[2], pcutline = NULL,
                         scheme = scheme2,
                         source = "plotly_manh2")
      })
      
      # chromosome plotly 2
      output$chrom2 <- renderPlotly({
        req(man2, coords$chr)
        validate(need(coords$chr %in% chr_set[[2]], "No data for this chromosome"))
        chr_manhat2 <- manhattan(data2[which(data2[, chrom[2]] == coords$chr), ],
                                 chrom[2], pos[2], p[2], labs[2], pcutoff = pcutoff,
                                 npoints = 1e5)
        chr <- suppressWarnings(as.numeric(coords$chr))
        if ((!is.na(chr) && chr %% 2 == 0 || coords$chr == "Y")) {
          scheme2[1] <- scheme2[2]
        }
        yr <- range(chr_manhat2$data$logP, na.rm = TRUE)
        isolate(chr_y2$range <- yr)
        isolate(chr_y2$max <- yr[2])
        isolate(xr <- coords$xrange)
        
        plotly_manhattan(chr_manhat2, scheme = scheme2, ylab = man_ylab[2],
                         source = "plotly_chrom2") %>%
          layout(margin = list(t = 5),
                 shapes = list(
                   list(type = "rect",
                        line = list(width = 1, color = "red"),
                        x0 = xr[1] / 1e6,
                        x1 = xr[2] / 1e6, y0 = 0, y1 = 1,
                        xref = "x", yref = "paper", layer = "below")))
      })
      
      # 2nd manhattan click
      observe({
        s <- event_data("plotly_click", source = "plotly_manh2")
        req(s)
        w <- which(data2[, labs[2]] == s$key)
        if (length(w) > 0) {
          coords$chr <- data2[w[1], chrom[2]]
          xr <- data2[w[1], pos[2]] + c(-5e5, 5e5)
          if (xr[1] < 0) xr <- c(0, 1e6)
          coords$xrange <- xr
        }
      })
      
      # 2nd chrom click
      observe({
        s <- event_data("plotly_click", source = "plotly_chrom2")
        req(s)
        w <- which(data2[, labs[2]] == s$key)
        if (length(w) > 0) {
          coords$chr <- data2[w[1], chrom[2]]
          xr <- data2[w[1], pos[2]] + c(-5e5, 5e5)
          if (xr[1] < 0) xr <- c(0, 1e6)
          coords$xrange <- xr
        }
      })
      
      # zoom manhattan2 y axis
      m_ylim2 <- reactiveValues(max = manhat2$yrange[2])
      
      observeEvent(input$m_zoomin2, {
        m_ylim2$max <- pmax(m_ylim2$max * 0.88, 5)
        yr <- c(manhat2$yrange[1], m_ylim2$max)
        yr <- yr + diff(yr) * c(-0.05, 0.05)
        plotlyProxy("manhattan2", session) %>%
          plotlyProxyInvoke("relayout",
                            list(yaxis = list(range = yr,
                                              title = man_ylab[2],
                                              ticks = "outside",
                                              zeroline = FALSE, showline = TRUE)))
      })
      
      observeEvent(input$m_zoomout2, {
        m_ylim2$max <- pmin(m_ylim2$max / 0.88, manhat2$yrange[2])
        yr <- c(manhat2$yrange[1], m_ylim2$max)
        yr <- yr + diff(yr) * c(-0.05, 0.05)
        plotlyProxy("manhattan2", session) %>%
          plotlyProxyInvoke("relayout",
                            list(yaxis = list(range = yr,
                                              title = man_ylab[2],
                                              ticks = "outside",
                                              zeroline = FALSE, showline = TRUE)))
      })
      
      # zoom chrom2 y axis
      chr_y2 <- reactiveValues(max = 0, range = c(0, 0))
      
      observeEvent(input$chr_zoomin2, {
        chr_y2$max <- pmax(chr_y2$max * 0.88, 5)
        yr <- c(chr_y2$range[1], chr_y2$max)
        yr <- yr + diff(yr) * c(-0.05, 0.05)
        plotlyProxy("chrom2", session) %>%
          plotlyProxyInvoke("relayout",
                            list(yaxis = list(range = yr,
                                              title = man_ylab[2],
                                              ticks = "outside",
                                              zeroline = FALSE, showline = TRUE)))
      })
      
      observeEvent(input$chr_zoomout2, {
        chr_y2$max <- pmin(chr_y2$max / 0.88, chr_y2$range[2])
        yr <- c(chr_y2$range[1], chr_y2$max)
        yr <- yr + diff(yr) * c(-0.05, 0.05)
        plotlyProxy("chrom2", session) %>%
          plotlyProxyInvoke("relayout",
                            list(yaxis = list(range = yr,
                                              title = man_ylab[2],
                                              ticks = "outside",
                                              zeroline = FALSE, showline = TRUE)))
      })
      
      # chrom2 highlight
      observeEvent(coords$xrange, {
        req(man2, input$show_chrom, coords$chr)
        plotlyProxy("chrom2", session) %>%
          plotlyProxyInvoke("relayout",
                            list(shapes = list(
                              list(type = "rect",
                                   line = list(width = 1, color = "red"),
                                   x0 = coords$xrange[1] / 1e6,
                                   x1 = coords$xrange[2] / 1e6, y0 = 0, y1 = 1,
                                   xref = "x", yref = "paper", layer = "below")
                            )))
      })
      
    }  # end of 2nd manhattan section
    
    input_biotype <- reactive({input$biotype}) %>% debounce(2000)
    
    loc <- reactiveValues(i = NULL)
    locv2 <- reactiveValues(i = NULL)
    ntrace <- reactiveVal()
    genes <- reactiveValues(x = NULL)
    ld_snp <- reactiveVal(NULL)
    cur_index <- reactiveVal(NULL)
    
    output$locus <- renderPlotly({
      req(coords$chr, coords$xrange)
      # temporary fix for plotly minallowed not working 
      req(coords$xrange[1] >= 0)
      loc1 <- locus(data = data, xrange = coords$xrange,
                     seqname = coords$chr, ens_db = ens_db,
                     chrom = chrom[1], pos = pos[1], p = p[1], labs = labs[1])
      # validate(need(loc1$data, "Locus contains no SNPs/datapoints"))
      validate(need(nrow(loc1$data) < 1.5e5, "Too many datapoints. Zoom in."))
      loc1$TX$fullname <- expandGenes(loc1$TX, fullnames)
      if (!is.null(recomb) && input$recomb) {
        loc1 <- link_recomb(loc1, recomb = recomb)
      }
      
      if (man2) {
        loc2 <- locus(data = data2, xrange = coords$xrange,
                      seqname = coords$chr, ens_db = ens_db,
                      chrom = chrom[2], pos = pos[2], p = p[2], labs = labs[2],
                      tx = FALSE)
        if (!is.null(recomb) && input$recomb) {
          loc2 <- link_recomb(loc2, recomb = recomb)
        }
        locv2$i <- loc2
      }
      
      isolate(cur_index(loc1$index_snp))
      pin <- ld_snp()
      ld_msg <- NULL
      if (!is.null(pin) && pin %in% loc1$data[, labs[1]]) {
        loc1$index_snp <- pin
        loc1b <- withCallingHandlers(
          try(link_LD(loc1, token = ld_token, pop = ld_pop)),
          message = function(m) {
            txt <- conditionMessage(m)
            ld_msg <<- trimws(txt)
          })
        if (inherits(loc1b, "try-error")) {
          ld_msg <- attr(loc1b, "condition")$message
        } else loc1 <- loc1b
        removeNotification("ld_busy")
        if (!"ld" %in% colnames(loc1$data)) {
          showNotification(
            paste0("LD failed for ", pin,
                   if (is.null(ld_msg)) "" else paste0(" - ", ld_msg)),
            type = "error", duration = 10)
        }
      }
      loc$i <- loc1
      
      if (!is.null(eqtl_gene)) {
        ind <- loc1$data[, p[1]] < pcutoff
        eqtls <- loc1$data[ind, eqtl_gene]
        genes$x <- genes1 <- unique(eqtls)
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
        h <- if (!man2) c(364, 24 * needrow + 80) else c(230, 230, 20 * needrow + 40)
        maxrows <- NULL
      }
      hideFeedback("tex")
      p <- locus_plotly(loc1, h, filter_gene_biotype = biotype, pcutoff = pcutoff,
                   width = width, eqtl_gene = eqtl_gene, beta = eqtl_beta,
                   add_hover = add_hover, scheme = locscheme, maxrows = maxrows,
                   loc2 = if (man2) loc2 else NULL,
                   ylab = man_ylab)
      ntrace(length(p$x$data) -2)
      p
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
    
    # temporary fix for R/plotly layout.xaxis.minallowed not working
    observeEvent(coords$xrange, {
      if (coords$xrange[1] < 0) {
        coords$xrange <- signif(coords$xrange - coords$xrange[1], 3)
      }
    })
    
    output$pos <- renderText({
      req(coords$chr %in% chr_set[[1]], coords$xrange)
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
          SeqNameFilter(seq_filter)))
        if (length(loc) > 1) loc <- loc[1]
        chr <- names(seqlengths(loc))
        m <- mean(c(start(loc), end(loc)))
        xr <- as.integer(c(m - 5e5, m + 5e5))
      } else if (grepl("^rs", input$tex)) {
        w <- which(data[, labs[1]] == input$tex)
        if (length(w) > 0) {
          chr <- data[w[1], chrom[1]]
          xr <- data[w[1], pos[1]] + c(-5e5, 5e5)
        } else {
          showFeedback("tex", "SNP not found")
          return()
        }
      } else {
        if (nchar(tex) > 1) showFeedback("tex", "not found")
        return()
      }
      xr <- as.integer(pmax(xr, 0))
      
      if (chr %in% chr_set[[1]]) {
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
      req(coords$chr %in% chr_set[[1]], coords$xrange)
      s <- event_data("plotly_relayout", source = "plotly_locus")
      req(c("xaxis.range[0]", "xaxis.range[1]") %in% names(s))
      xr <- c(s$`xaxis.range[0]`, s$`xaxis.range[1]`)
      xd <- diff(xr) * (1 - 1/1.02) / 2
      xr <- xr + c(xd, -xd)
      coords$xrange <- as.integer(xr * 1e6)
    })
    
    loc_width <- reactiveVal(600)
    
    observe({
      loc_width(session$clientData$output_locus_width)
    })
    
    # redo gene tracks only
    observeEvent(c(loc_width(), input_biotype()), {
      req(loc$i)
      maxrows <- if (input$alltracks) NULL else 8
      gt <- genetrack_ly(loc$i, filter_gene_biotype = input_biotype(),
                         width = loc_width(), blanks = "show", plot = FALSE,
                         maxrows = maxrows)
      if (nrow(gt$TX) == 0) {
        # blank gene tracks
        p <- plotlyProxy("locus", session) %>%
          plotlyProxyInvoke("restyle",
                            list(x = list(NULL), y = list(NULL), text = list(NULL),
                                 hoverinfo = "text"),
                            list(ntrace())) %>%
          plotlyProxyInvoke("update",
                            list(x = list(NULL), y = list(NULL),
                                 text = list(NULL), hoverinfo = "none"),
                            list(shapes = NULL),
                            list(ntrace() + 1L)) %>%
          plotlyProxyInvoke("relayout",
                            list(annotations = list(list(x = 0, y = 0.008,
                                                    text = "No gene tracks",
                                                    xref = "paper", yref = "paper",
                                                    showarrow = FALSE))))
        return(p)
      }
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
      rec <- !is.null(recomb) && input$recomb
      yref <- paste0("y", rec + man2 + rec * man2 +2)
      
      y0 <- -EX$row - 0.15
      y1 <- -EX$row + 0.15
      shapes <- lapply(seq_len(nrow(EX)), function(i) {
        list(type = "rect", fillcolor = exon_col, line = list(color = exon_border,
                                                              width = 0.5),
             x0 = EX$start[i], x1 = EX$end[i], xref = "x",
             y0 = y0[i], y1 = y1[i], yref = yref)
      })
      ok <- !is.na(TX$gene_name2)
      if (sum(ok) > 0) {
        xtex <- TX$tx[ok]
        ytex <- TX$ty[ok]
        ttext <- TX$gene_name2[ok]
      } else {
        xtex <- ytex <- 0
        ttext <- ""
      }
      
      plotlyProxy("locus", session) %>%
        plotlyProxyInvoke("restyle",
                          list(x = list(lx), y = list(ly), text = list(ht),
                               hoverinfo = "text"),
                          list(ntrace())) %>%
        plotlyProxyInvoke("update",
                          list(x = list(xtex), y = list(ytex),
                               text = list(ttext), hoverinfo = "none"),
                          list(shapes = shapes),
                          list(ntrace() + 1L)) %>%
        plotlyProxyInvoke("relayout",
                          list(annotations = list(NULL)))
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
    
    get_ld <- reactiveVal(FALSE)
    
    # retrieve LD
    observeEvent(input$ld_get, {
      snp <- cur_index()
      if (is.null(snp) || is.na(snp)) {
        showNotification("No index SNP in view", type = "warning")
        return()
      }
      showNotification(paste0("Fetching LD for ", snp, " - click a point to re-base"),
                       id = "ld_busy", duration = NULL)
      ld_snp(snp)
    })
    
    observeEvent(input$ld_clear, {
      ld_snp(NULL)
      removeNotification("ld_busy")
    })
    
    observe({
      s <- event_data("plotly_click", source = "plotly_locus")
      req(s, !is.null(s$key))
      cur <- isolate(ld_snp())
      req(!is.null(cur))
      snp <- as.character(s$key)[1]
      req(!is.na(snp), nzchar(snp))
      if (identical(snp, cur)) return()
      showNotification(paste0("Fetching LD for ", snp),
                       id = "ld_busy", duration = NULL)
      ld_snp(snp)
    })
    
    output$ld_status <- renderText({
      req(ld_snp() %in% loc$i$data[, labs])
      paste0("LD: ", ld_snp(), " (", ld_pop, ")")
    })
    
    # inline conditional UI
    output$save_ui <- renderUI({
      req(coords$chr)
      downloadButton("save", NULL, icon = icon("floppy-disk"))
    })
    
    # save files
    output$save <- downloadHandler(filename = function() {
      paste0("loc_", dat_name, "_", coords$chr, "_", coords$xrange[1], ".pdf")
    },
    content = function(file) {
      req(loc$i)
      pdf(file)
      if (!man2) {
        locus_plot(loc$i)
      } else {
        oldpar <- set_layers(2)
        on.exit(par(oldpar))
        scatter_plot(loc$i, xticks = FALSE, bty = "u",
                     ylab = bquote(.(traits[1]) ~ -log[10] ~ P))
        scatter_plot(locv2$i, xticks = FALSE, bty = "u",
                     ylab = bquote(.(traits[2]) ~ -log[10] ~ P))
        genetracks(loc$i, blanks = "hide")
      }
      dev.off()
    })
    
  }
  
  runApp(list(ui = ui, server = server)) %>%
    suppress_warnings("please add `event_register\\(p")
}
