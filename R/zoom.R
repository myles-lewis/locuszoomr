
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
#' @param scheme Vector of 4 colours for main Manhattan plot: 1st, 2nd colours
#'   for alternating chromosomes, 3rd colour for significant points, 4th colour
#'   for highlighting chromosomes.
#' @param scheme2 Vector of colours for 2nd Manhattan plot.
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param eqtl_gene Determines which column in `data` contains eQTL genes.
#'   Currently `data2` cannot be used to show eQTL colour data.
#' @param beta Optional column name for beta coefficient to display upward
#'   triangles for positive beta and downward triangles for negative beta
#'   (significant SNPs only). If `data2` is supplied, then a vector can be used
#'   to specify different beta columns in `data` and `data2`; use `NA` to
#'   indicate no `beta`.
#' @param eqtl_scheme Colour scheme for eQTL genes.
#' @param add_hover Optional vector of column names in `data` to add to the
#'   plotly hover text for scatter points. Not available for `data2`.
#' @param mh_points Number of points to display in manhattan plot. Default is
#'   `1e5`.
#' @param recomb Optional `GRanges` class object of recombination data.
#' @param gwas_catalog Optional local GWAS Catalog, shown as an extra track at
#'   every locus: `TRUE` to load the copy saved by [download_gwas_catalog()],
#'   or an object from [prep_gwas_catalog()] or the path of an RDS file of it.
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
#' @param align Logical, whether to align the chromosomes in the main Manhattan
#'   plots and chromosome subplots when comparing 2 GWAS.
#' @returns No return value. Opens an interactive shiny window.
#' @importFrom plotly plotlyOutput renderPlotly event_data config plotlyProxy
#' @importFrom plotly plotlyProxyInvoke layout
#' @importFrom shiny fluidPage tabsetPanel tabPanel fluidRow column actionButton 
#' @importFrom shiny icon uiOutput checkboxInput textOutput splitLayout req
#' @importFrom shiny textInput conditionalPanel h5 runApp debounce isolate
#' @importFrom shiny renderUI reactiveValues reactive observe observeEvent radioButtons
#' @importFrom shiny reactiveVal validate need renderText updateTextInput outputOptions
#' @importFrom shiny showNotification removeNotification HTML downloadButton downloadHandler
#' @importFrom shiny showTab hideTab
#' @importFrom htmltools tags br
#' @importFrom gtools mixedsort
#' @importFrom stats as.formula setNames
#' @importFrom grDevices dev.off pdf
#' @export

zoom <- function(data, ens_db,
                 chrom = NULL, pos = NULL, p = NULL, labs = NULL,
                 data2 = NULL,
                 traits = NULL,
                 scheme = c('royalblue', 'skyblue', 'red', "#00CD00"),
                 scheme2 = c("#33a02c", "#b2df8a", "purple", "red"),
                 pcutoff = 5e-8,
                 eqtl_gene = NULL,
                 beta = NULL,
                 eqtl_scheme = c("#FF0000", "#00FFFF", "#FF9000", "#0080FF", "#FFFF00",
                                 "#0000FF", "#80DD00", "#8000FF", "#009900", "#FF00FF"),
                 add_hover = NULL,
                 mh_points = 1e5,
                 recomb = NULL,
                 gwas_catalog = NULL,
                 ld_token = Sys.getenv("LDLINK_TOKEN"),
                 ld_pop = "EUR",
                 seq_filter = c(1:22, 'X', 'Y'),
                 AnnotationDb = "org.Hs.eg.db",
                 align = TRUE) {
  if (!requireNamespace("shinycssloaders", quietly = TRUE))
    stop("Need shinycssloaders package", call. = FALSE)
  if (!requireNamespace("shinyWidgets", quietly = TRUE))
    stop("Need shinyWidgets package", call. = FALSE)
  if (!requireNamespace("shinyFeedback", quietly = TRUE))
    stop("Need shinyFeedback package", call. = FALSE)
  if (!requireNamespace("DT", quietly = TRUE))
    stop("Need DT package", call. = FALSE)
  
  # how zoom() was called, for the R code button
  data_sym <- arg_name(substitute(data))
  data2_sym <- arg_name(substitute(data2))
  recomb_sym <- arg_name(substitute(recomb))
  ens_source <- ens_db_source(substitute(ens_db), ens_db)
  gcat_source <- catalog_source(substitute(gwas_catalog), gwas_catalog)
  gwas_catalog <- as_gwas_catalog(gwas_catalog)
  dat_name <- deparse(substitute(data))
  raw_names <- colnames(data)
  data <- data.frame(data)
  # autodetect headings
  dc <- detect_cols(data, chrom[1], pos[1], p[1], labs[1])
  
  # autodetect headings gwas2
  dc2 <- NULL
  man2 <- FALSE
  if (!is.null(data2)) {
    dat2_name <- deparse(substitute(data2))
    raw_names2 <- colnames(data2)
    data2 <- data.frame(data2)
    if (!is.null(chrom)) chrom <- rep_len(chrom, 2)
    if (!is.null(pos)) pos <- rep_len(pos, 2)
    if (!is.null(p)) p <- rep_len(p, 2)
    if (!is.null(labs)) labs <- rep_len(labs, 2)
    if (!is.null(beta)) beta <- rep_len(beta, 2)
    dc2 <- detect_cols(data2, chrom[2], pos[2], p[2], labs[2])
    if (is.null(traits)) traits <- c(dat_name, dat2_name)
    man2 <- TRUE
  }
  
  chrom <- c(dc$chrom, dc2$chrom)
  pos <- c(dc$pos, dc2$pos)
  p <- c(dc$p, dc2$p)
  labs <- c(dc$labs, dc2$labs)
  
  message("Generating Manhattan plot", (if (man2) " 1"))
  data_prep <- prep_changes(data, raw_names, chrom[1], pos[1], p[1], labs[1])
  data <- complete_data(data, chrom[1], pos[1], p[1])
  
  if (is.null(eqtl_gene)) {
    data[, labs[1]] <- unique_snps(data, labs[1], chrom[1])
  } else {
    data[, labs[1]] <- unique_snps(data, labs[1], eqtl_gene)
  }
  # currently eqtl_gene can only apply to data1
  
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
    eqtl_colour <- eqtl_colours(data[data[, p[1]] < pcutoff, ], chrom[1], pos[1],
                                eqtl_gene, eqtl_scheme)
  }
  
  show_ld <- nzchar(ld_token)
  show_eqtl <- show_ld & is.null(eqtl_gene)
  # local Catalog button: always shows the Catalog, this just toggles its pickers
  gcat_btn <- actionButton("gcat_show", "GWAS Catalog", class = "btn-success btn-sm")
  ot_btn <- actionButton("ot_get",
                         shiny::tagList(tags$img(
                           src = "https://platform.opentargets.org/favicon.png",
                           alt = ""), "Open Targets"),
                         class = "btn-sm")
  
  # apply min_p_snp to data for manhat?
  # smallest floating point
  data[which(data[, p[1]] < 5e-324), p[1]] <- 5e-324
  manhat <- manhattan(data, chrom[1], pos[1], p[1], labs[1], pcutoff = pcutoff,
                      npoints = mh_points)
  chr_set <- list()
  chr_set[[1]] <- manhat$chrom_list
  man_ylab <- "-log<sub>10</sub> P"
  height <- c(300, 220, 624)
  chrom_lim <- manhat$chrom_lim
  man_lim <- NULL
  
  if (man2) {
    message("Generating Manhattan plot 2")
    data2_prep <- prep_changes(data2, raw_names2, chrom[2], pos[2], p[2], labs[2])
    data2 <- complete_data(data2, chrom[2], pos[2], p[2])
    data2[, labs[2]] <- unique_snps(data2, labs[2], chrom[2])
    data2[which(data2[, p[2]] < 5e-324), p[2]] <- 5e-324
    manhat2 <- manhattan(data2, chrom[2], pos[2], p[2], labs[2], pcutoff = pcutoff,
                         npoints = mh_points)
    chr_set[[2]] <- manhat2$chrom_list
    full_chr_set <- unique(unlist(chr_set))
    man_ylab <- paste(traits, man_ylab)
    height <- c(220, 180, 824)
    if (align) {
      message("Aligning Manhattans")
      aligned <- align_manhats(list(manhat, manhat2))
      manhat <- aligned$manhats[[1]]
      manhat2 <- aligned$manhats[[2]]
      chrom_lim <- aligned$chrom_lim
      man_lim <- aligned$man_lim
    } else {
      chrom_lim <- NULL
    }
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
      }")),
      tags$style(HTML("
      #ld_get, #eqtl_get, #trait_get, #gcat_show, #ot_get, #ld_clear { margin-bottom: 4px; }
      #ot_get { background-color: #fff; border-color: #3489ca; color: #3489ca; }
      #ot_get:hover, #ot_get:focus { border-color: #18405e; color: #18405e; }
      #ot_get img { height: 14px; margin-right: 4px; vertical-align: -2px; }")),
      # with a local Catalog its pickers can show alongside eQTL/QTL pickers
      (if (!is.null(gwas_catalog)) {
        tags$style(HTML("
      #ui_link_traits:not(:empty), #ui_link_qtl:not(:empty), #ui_link_genes:not(:empty) {
        border-top: 1px solid #ddd; margin-top: 8px; padding-top: 4px; }"))
      })
    ),
    shinyFeedback::useShinyFeedback(),
    tabsetPanel(id = "topTabs",
      tabPanel("Plot",
               fluidRow(
                 column(11,
                        shinycssloaders::withSpinner(
                          plotlyOutput("manhattan", width = "85vw", height = height[1]),
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
                          shinycssloaders::withSpinner(
                            plotlyOutput("manhattan2", width = "85vw", height = height[1]),
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
                                                  shinycssloaders::withSpinner(
                                                    plotlyOutput("chrom", width = "85vw", height = height[2]),
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
                                                    shinycssloaders::withSpinner(
                                                      plotlyOutput("chrom2", width = "85vw", height = height[2]),
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
                        uiOutput("save_ui", inline = T),
                        uiOutput("fig_code_ui", inline = T)
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
                        shinyWidgets::dropdown(
                          (if (!is.null(recomb)) {
                            checkboxInput("recomb", "show recombination rate", value = TRUE)
                          } else NULL),
                          (if (!is.null(gwas_catalog)) {
                            checkboxInput("gwas_catalog", "show GWAS Catalog", value = TRUE)
                          }),
                          checkboxInput("alltracks", "show all gene tracks"),
                          shinyWidgets::pickerInput("biotype", h5("Select gene biotypes"),
                                      choices = biotypes, selected = biotypes,
                                      multiple = TRUE,
                                      options = shinyWidgets::pickerOptions(actionsBox = TRUE,
                                                              selectedTextFormat = 'count > 1')),
                          (if (show_ld) {
                            fluidRow(
                              column(12,
                                     h5("LDlink tools"),
                                     actionButton("ld_get", "Get LD", icon = icon("circle-nodes"),
                                                  class = "btn-primary btn-sm"),
                                     (if (show_eqtl) {
                                       actionButton("eqtl_get", "eQTL", icon = icon("compass"),
                                                    class = "btn-info btn-sm")
                                     }),
                                     (if (is.null(gwas_catalog)) {
                                       actionButton("trait_get", "GWAS trait",
                                                    class = "btn-success btn-sm")
                                     } else gcat_btn),
                                     ot_btn,
                                     actionButton("ld_clear", "Clear",
                                                  class = "btn-default btn-sm")
                              ))
                          }),
                          # the local Catalog and Open Targets need no LDlink token
                          (if (!show_ld && !is.null(gwas_catalog)) {
                            fluidRow(
                              column(12,
                                     h5("GWAS Catalog"),
                                     gcat_btn, ot_btn,
                                     actionButton("ld_clear", "Clear",
                                                  class = "btn-default btn-sm")))
                          }),
                          (if (show_ld || !is.null(gwas_catalog)) {
                            conditionalPanel(
                              if (is.null(gwas_catalog)) "true" else "output.gcat_open",
                            fluidRow(
                              column(12,
                                     uiOutput("ui_link_traits"),
                                     conditionalPanel("output.has_trait",
                                                      checkboxInput("trait_collapse", "collapse duplicate associations",
                                                                    value = TRUE)
                                     ),
                                     (if (show_ld && !is.null(gwas_catalog)) {
                                       shiny::tagList(
                                         shiny::sliderInput("trait_r2_min", "minimum r\u00b2 with the index SNP",
                                                            min = 0, max = 1, value = 0, step = 0.05),
                                         checkboxInput("trait_ld_colour", "colour by LD", value = FALSE))
                                     })
                              )))
                          }),
                          (if (show_ld || !is.null(gwas_catalog)) {
                            uiOutput("ui_link_qtl")
                          }),
                          (if (show_eqtl) {
                            uiOutput("ui_link_genes")
                          }),
                          (if (!is.null(eqtl_gene)) {
                            uiOutput("ui_genes")
                          }),
                          radioButtons("export", h5("Export"),
                                       list(pdf = "pdf", plotly = "rds"), inline = TRUE),
                          right = TRUE, icon = icon("gear")
                        ))
                 ),
                 fluidRow(
                   column(12,
                          plotlyOutput("locus", width = "95vw", height = height[3])
                   )
                 )
      ),
      tabPanel("Table",
               fluidRow(
                 column(12, br(),
                        conditionalPanel('output.coords_ok',
                                         radioButtons("tab1_type", label = NULL,
                                                      choices = c("local", "genome"),
                                                      inline = TRUE)
                        ),
                        DT::dataTableOutput("table")
                 ))
               ),
      (if (man2) {
        tabPanel("Table 2",
                 fluidRow(
                   column(12, br(),
                          conditionalPanel('output.coords_ok',
                                           radioButtons("tab2_type", label = NULL,
                                                        choices = c("local", "genome"),
                                                        inline = TRUE)
                          ),
                          DT::dataTableOutput("table2")
                   ))
                 )
      }),
      (if (!is.null(gwas_catalog)) {
        tabPanel("GWAS Catalog",
                 fluidRow(column(12, br(), DT::dataTableOutput("catalog_table"))))
      }),
      tabPanel(title = uiOutput("extra_tab_nm"), value = "linktab",
               fluidRow(
                 column(12, br(), DT::dataTableOutput("link_table")
                 )))
    )
  )
  
  server <- function(input, output, session) {
    
    output$manhattan <- renderPlotly({
      plotly_manhattan(manhat, ylab = man_ylab[1], pcutline = NULL,
                       scheme = scheme, xlim = man_lim)
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
                       xlim = chrom_lim[coords$chr, ],
                       source = "plotly_chrom") %>%
        layout(margin = list(t = 5),
               shapes = list(
                 list(type = "rect",
                      line = list(width = 1, color = scheme[4]),
                      x0 = xr[1] / 1e6,
                      x1 = xr[2] / 1e6, y0 = 0, y1 = 1,
                      xref = "x", yref = "paper", layer = "above")))
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
    
    ###########
    # 2nd gwas
    if (man2) {
      output$manhattan2 <- renderPlotly({
        req(man2)
        plotly_manhattan(manhat2, ylab = man_ylab[2], pcutline = NULL,
                         scheme = scheme2, xlim = man_lim,
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
                         xlim = chrom_lim[coords$chr, ],
                         source = "plotly_chrom2") %>%
          layout(margin = list(t = 5),
                 shapes = list(
                   list(type = "rect",
                        line = list(width = 1, color = scheme2[4]),
                        x0 = xr[1] / 1e6,
                        x1 = xr[2] / 1e6, y0 = 0, y1 = 1,
                        xref = "x", yref = "paper", layer = "above")))
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
      
      # Table 2 tab
      output$table2 <- DT::renderDataTable({
        dat <- data2
        if (input$tab2_type == "local" && !is.null(coords$chr)) {
          dat <- locv2$i$data
          validate(need(!is.null(dat), "No data at locus"))
        }
        fcols <- colnames(dat)[sapply(dat, class) == "numeric"]
        fcols <- fcols[!fcols %in% c(chrom[2], pos[2])]
        DT::datatable(dat, caption = traits[2], rownames = FALSE) %>%
          DT::formatSignif(fcols, digits = 3)
      })
      
    }  # end of 2nd manhattan section
    
    #########
    
    input_biotype <- reactive({input$biotype}) %>% debounce(2000)
    
    loc <- reactiveValues(i = NULL)
    locv2 <- reactiveValues(i = NULL)
    ntrace <- reactiveVal()
    genes <- reactiveValues(x = NULL)
    link_snp <- reactiveVal(NULL)
    link_type <- reactiveVal(NULL)
    cur_index <- reactiveVal(NULL)
    link_eqtl <- reactiveValues(genes = NULL, tissues = NULL, traits = NULL,
                                categories = NULL, qtl_types = NULL,
                                qtl_genes = NULL)
    link_data <- reactiveValues(LDtrait = NULL, trait_category = NULL,
                                LDexp = NULL, OTqtl = NULL)
    ld_colour <- reactiveVal(FALSE)
    save_plotly <- reactiveValues(p = NULL)
    
    output$locus <- renderPlotly({
      req(coords$chr, coords$xrange)
      # temporary fix for plotly minallowed not working 
      req(coords$xrange[1] >= 0)
      loc1 <- locus(data = data, xrange = coords$xrange,
                     seqname = coords$chr, ens_db = ens_db,
                     chrom = chrom[1], pos = pos[1], p = p[1], labs = labs[1])
      validate(need(is.null(loc1$data) || nrow(loc1$data) < 1.5e5,
                    "Too many datapoints. Zoom in."))
      loc1$TX$fullname <- expandGenes(loc1$TX, fullnames)
      if (!is.null(recomb) && input$recomb) {
        loc1 <- link_recomb(loc1, recomb = recomb)
      }
      
      if (man2) {
        loc2 <- locus(data = data2, xrange = coords$xrange,
                      seqname = coords$chr, ens_db = ens_db,
                      chrom = chrom[2], pos = pos[2], p = p[2], labs = labs[2],
                      tx = FALSE)
        validate(need(is.null(loc2$data) || nrow(loc2$data) < 1e5,
                      "Too many datapoints. Zoom in."))
        if (!is.null(recomb) && input$recomb) {
          loc2 <- link_recomb(loc2, recomb = recomb)
        }
      }
      
      # retrieve LD
      isolate(cur_index(loc1$index_snp))
      pin <- link_snp()
      ld_msg <- NULL
      if (link_type() == "LD" && !is.null(pin) &&
          (pin %in% loc1$data[, labs[1]] || (man2 && pin %in% loc2$data[, labs[2]]))) {
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
        if (man2 && !is.null(pin) && "ld" %in% colnames(loc1$data)) {
          loc2$index_snp <- pin
          loc2b <- try(link_LD(loc2, token = ld_token, pop = ld_pop))
          if (!inherits(loc2b, "try-error")) loc2 <- loc2b
        }
      }
      
      # retrieve eQTL
      eqtl_scheme <- NULL
      if (link_type() == "eQTL" && !is.null(pin) &&
          (pin %in% loc1$data[, labs[1]] || (man2 && pin %in% loc2$data[, labs[2]]))) {
        loc1$index_snp <- pin
        loc1b <- try(link_eqtl(loc1, token = ld_token), silent = TRUE)
        if (inherits(loc1b, "try-error")) {
          ld_msg <- trimws(attr(loc1b, "condition")$message)
          message("\n", ld_msg)
        } else loc1 <- loc1b
        removeNotification("ld_busy")
        if (is.null(loc1$LDexp)) {
          showNotification(
            paste0("eQTL request failed for ", pin,
                   if (is.null(ld_msg)) "" else paste0(" - ", ld_msg)),
            type = "error", duration = 10)
        } else {
          LDX1 <- loc1$LDexp[loc1$LDexp$inData, ]
          g <- unique(LDX1$Gene_Symbol[order(LDX1$P_value)])
          tiss <- unique(LDX1$Tissue)
          if (man2 && !is.null(pin) && !is.null(loc1$LDexp)) {
            loc2$index_snp <- pin
            loc2b <- try(suppressMessages(link_eqtl(loc2, token = ld_token)))
            if (!inherits(loc2b, "try-error")) {
              loc2 <- loc2b
              LDX2 <- loc2$LDexp[loc2$LDexp$inData, ]
              g <- union(g, LDX2$Gene_Symbol[order(LDX2$P_value)])
              tiss <- union(tiss, LDX2$Tissue)
              eqtl_scheme <- setNames(rainbow(length(g)), g)  # synchronise colours
            }
          }
          link_eqtl$genes <- sort(g)
          link_eqtl$tissues <- sort(tiss)
        }
      } else {
        link_eqtl$genes <- NULL
      }
      
      loc$i <- loc1
      if (man2) locv2$i <- loc2
      
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
      
      # retrieve trait
      if (link_type() == "GWAS trait" && !is.null(pin) &&
          (pin %in% loc1$data[, labs[1]] || (man2 && pin %in% loc2$data[, labs[2]]))) {
        loc1$index_snp <- pin
        loc1b <- try(link_trait(loc1, token = ld_token, category = TRUE))
        if (inherits(loc1b, "try-error")) {
          ld_msg <- attr(loc1b, "condition")$message
        } else loc1 <- loc1b
        removeNotification("ld_busy")
        if (is.null(loc1$LDtrait)) {
          showNotification(
            paste0("GWAS trait request failed for ", pin,
                   if (is.null(ld_msg)) "" else paste0(" - ", ld_msg)),
            type = "error", duration = 10)
        }
        link_eqtl$traits <- sort(unique(loc1$LDtrait$GWAS_Trait))
        link_eqtl$categories <- sort(unique(loc1$trait_category$category))
      } else if (!is.null(gwas_catalog) && isTRUE(input$gwas_catalog)) {
        loc1b <- try(link_gwas_catalog(loc1, gwas_catalog))
        if (inherits(loc1b, "try-error")) {
          showNotification(attr(loc1b, "condition")$message,
                           id = "gcat_error", type = "error", duration = 10)
        } else {
          loc1 <- loc1b
          traits_new <- sort(unique(loc1$LDtrait$GWAS_Trait))
          # the pickers re-render with every trait selected, so drop filters
          # from the previous locus now rather than in a second render
          if (!identical(traits_new, isolate(link_eqtl$traits))) {
            trait_filter(NULL)
            trait_category_filter(NULL)
          }
          link_eqtl$traits <- traits_new
          link_eqtl$categories <- sort(unique(loc1$trait_category$category))
        }
      } else {
        link_eqtl$traits <- NULL
        link_eqtl$categories <- NULL
      }
      
      # retrieve Open Targets QTL credible sets
      if (identical(link_type(), "Open Targets") && !is.null(pin) &&
          (pin %in% loc1$data[, labs[1]] || (man2 && pin %in% loc2$data[, labs[2]]))) {
        loc1$index_snp <- pin
        ld_msg <- NULL
        loc1b <- try(link_ot(loc1))
        if (inherits(loc1b, "try-error")) {
          ld_msg <- attr(loc1b, "condition")$message
        } else loc1 <- loc1b
        removeNotification("ld_busy")
        if (is.null(loc1$OTqtl)) {
          showNotification(
            paste0("No Open Targets QTL data for ", pin,
                   if (is.null(ld_msg)) "" else paste0(" - ", ld_msg)),
            type = "warning", duration = 10)
        }
        link_eqtl$qtl_types <- sort(unique(loc1$OTqtl$type))
        link_eqtl$qtl_genes <- sort(unique(loc1$OTqtl$gene))
      } else {
        link_eqtl$qtl_types <- NULL
        link_eqtl$qtl_genes <- NULL
      }
      link_data$LDtrait <- loc1$LDtrait
      link_data$trait_category <- loc1$trait_category
      link_data$LDexp <- loc1$LDexp
      link_data$OTqtl <- loc1$OTqtl
      
      isolate(width <- loc_width())
      isolate(biotype <- input_biotype())
      h <- c(0.6, 0.4)
      maxrows <- 8
      if (input$alltracks) {
        cex.width <- 0.7 * par("pin")[1] * 80 / (width - 250)
        tryTX <- mapRow(loc1$TX, xlim = loc1$xrange, cex.text = cex.width,
                        blanks = "show")
        needrow <- pmax(max(tryTX$row, na.rm = TRUE), 8)
        h <- if (!man2) c(364, 24 * needrow + 80) else c(305, 305, 20 * needrow + 40)
        maxrows <- NULL
      }
      shinyFeedback::hideFeedback("tex")
      p <- locus_plotly(loc1, h, filter_gene_biotype = biotype, pcutoff = pcutoff,
                   width = width, eqtl_gene = eqtl_gene, beta = beta,
                   add_hover = add_hover, scheme = locscheme, maxrows = maxrows,
                   loc2 = if (man2) loc2 else NULL,
                   ylab = man_ylab,
                   gene_filter = eqtl_gene_filter(),
                   tissue_filter = eqtl_tissue_filter(),
                   eqtl_scheme = eqtl_scheme,
                   trait_filter = trait_filter(),
                   trait_category_filter = trait_category_filter(),
                   trait_r2_min = if ("ld" %in% colnames(loc1$data)) r2_min(),
                   trait_ld_scheme = if (!is.null(gwas_catalog) && ld_colour() &&
                                         "ld" %in% colnames(loc1$data))
                     eval(formals(scatter_plotly)$LD_scheme),
                   trait_size_by_r2 = is.null(gwas_catalog),
                   trait_collapse = !isFALSE(input$trait_collapse),
                   qtl_type_filter = qtl_type_filter(),
                   qtl_gene_filter = qtl_gene_filter())
      ntrace(length(p$x$data) -2)
      save_plotly$p <- p
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
    
    output$ui_link_genes <- renderUI({
      req(link_eqtl$genes)
      fluidRow(column(12,
        shinyWidgets::pickerInput("eqtl_gene_filter", h5("Select eQTL genes"),
                    choices = link_eqtl$genes, selected = link_eqtl$genes,
                    multiple = TRUE,
                    options = shinyWidgets::pickerOptions(actionsBox = TRUE,
                                            selectedTextFormat = 'count > 3')),
        shinyWidgets::pickerInput("eqtl_tissue_filter", h5("Select eQTL tissues"),
                    choices = link_eqtl$tissues, selected = link_eqtl$tissues,
                    multiple = TRUE,
                    options = shinyWidgets::pickerOptions(actionsBox = TRUE,
                                            selectedTextFormat = 'count > 1',
                                            liveSearch = TRUE))
      ))
    })
    
    output$ui_link_traits <- renderUI({
      req(link_eqtl$traits)
      cat_ok <- !(is.null(link_eqtl$categories) || 
          (length(link_eqtl$categories) == 1 && is.na(link_eqtl$categories)))
      fluidRow(
        column(12,
               (if (cat_ok) {
                 shinyWidgets::pickerInput(
                   "trait_category_filter", h5("Select trait categories"),
                   choices = link_eqtl$categories,
                   selected = link_eqtl$categories,
                   multiple = TRUE,
                   options = shinyWidgets::pickerOptions(
                     actionsBox = TRUE,
                     selectedTextFormat = 'count > 1',
                     dropdownAlignRight = 'auto'))
               }),
               shinyWidgets::pickerInput("trait_filter", h5("Select traits"),
                                         choices = link_eqtl$traits,
                                         selected = link_eqtl$traits,
                                         multiple = TRUE,
                                         options = shinyWidgets::pickerOptions(
                                           actionsBox = TRUE,
                                           selectedTextFormat = 'count > 1',
                                           dropdownAlignRight = 'auto',
                                           liveSearch = TRUE))
        ))
    })
    
    # conditionalPanel for trait_collapse checkbox
    output$has_trait <- reactive({!is.null(link_eqtl$traits)})
    outputOptions(output, "has_trait", suspendWhenHidden = FALSE)
    
    output$ui_link_qtl <- renderUI({
      req(link_eqtl$qtl_types)
      fluidRow(column(12,
        shinyWidgets::pickerInput("qtl_type_filter", h5("Select QTL types"),
                                  choices = link_eqtl$qtl_types,
                                  selected = link_eqtl$qtl_types,
                                  multiple = TRUE,
                                  options = shinyWidgets::pickerOptions(
                                    actionsBox = TRUE,
                                    selectedTextFormat = 'count > 2')),
        (if (length(link_eqtl$qtl_genes) > 1) {
          shinyWidgets::pickerInput("qtl_gene_filter", h5("Select QTL genes"),
                                    choices = link_eqtl$qtl_genes,
                                    selected = link_eqtl$qtl_genes,
                                    multiple = TRUE,
                                    options = shinyWidgets::pickerOptions(
                                      actionsBox = TRUE,
                                      selectedTextFormat = 'count > 3'))
        })
      ))
    })
    
    # local GWAS Catalog button: show/hide the trait pickers; if the Catalog
    # track is switched off or replaced by Open Targets, bring it back
    gcat_open <- reactiveVal(FALSE)
    output$gcat_open <- reactive(gcat_open())
    outputOptions(output, "gcat_open", suspendWhenHidden = FALSE)
    observeEvent(input$gcat_show, {
      ot <- identical(link_type(), "Open Targets")
      if (!isTRUE(input$gwas_catalog) || ot) {
        shiny::updateCheckboxInput(session, "gwas_catalog", value = TRUE)
        if (ot) {
          link_snp(NULL)
          link_type(NULL)
        }
        gcat_open(TRUE)
      } else gcat_open(!gcat_open())
    })
    
    eqtl_gene_filter <- reactiveVal(NULL)
    eqtl_tissue_filter <- reactiveVal(NULL)
    trait_filter <- reactiveVal(NULL)
    trait_category_filter <- reactiveVal(NULL)
    qtl_type_filter <- reactiveVal(NULL)
    qtl_gene_filter <- reactiveVal(NULL)
    
    deb_gene_filter <- reactive({input$eqtl_gene_filter}) %>% debounce(500)
    deb_tissue_filter <- reactive({input$eqtl_tissue_filter}) %>% debounce(500)
    deb_trait_filter <- reactive({input$trait_filter}) %>% debounce(500)
    deb_category_filter <- reactive({input$trait_category_filter}) %>%
      debounce(500)
    trait_r2_min <- reactive({input$trait_r2_min}) %>% debounce(500)
    # held here so Clear etc. reset it in the same render, not after the
    # slider's round trip
    r2_min <- reactiveVal(0)
    deb_qtl_filter <- reactive({input$qtl_type_filter}) %>% debounce(500)
    deb_qtl_gene_filter <- reactive({input$qtl_gene_filter}) %>% debounce(500)
    
    # prevents ui update double trigger output$locus
    # NULL special meaning = all categories ticked
    observeEvent(deb_gene_filter(), {
      i <- deb_gene_filter()
      geneset <- link_eqtl$genes
      req(i, geneset)
      if (!identical(i, geneset)) {
        eqtl_gene_filter(i)
      } else eqtl_gene_filter(NULL)
    })
    
    observeEvent(deb_tissue_filter(), {
      i <- deb_tissue_filter()
      tset <- link_eqtl$tissues
      req(i, tset)
      if (!identical(i, tset)) {
        eqtl_tissue_filter(i)
      } else eqtl_tissue_filter(NULL)
    })
    
    observeEvent(deb_trait_filter(), {
      i <- deb_trait_filter()
      tset <- link_eqtl$traits
      req(i, tset)
      if (!identical(i, tset)) {
        trait_filter(i)
      } else trait_filter(NULL)
    })
    
    observeEvent(deb_category_filter(), {
      i <- deb_category_filter()
      cset <- link_eqtl$categories
      req(i, cset)
      if (!identical(i, cset)) {
        trait_category_filter(i)
      } else trait_category_filter(NULL)
    })
    
    observeEvent(deb_qtl_filter(), {
      i <- deb_qtl_filter()
      tset <- link_eqtl$qtl_types
      req(tset)
      if (is.null(i)) {
        qtl_type_filter(character(0))
      } else if (!identical(i, tset)) {
        qtl_type_filter(i)
      } else qtl_type_filter(NULL)
    }, ignoreNULL = FALSE)
    
    observeEvent(deb_qtl_gene_filter(), {
      i <- deb_qtl_gene_filter()
      gset <- link_eqtl$qtl_genes
      req(length(gset) > 1)
      if (is.null(i)) {
        qtl_gene_filter(character(0))
      } else if (!identical(i, gset)) {
        qtl_gene_filter(i)
      } else qtl_gene_filter(NULL)
    }, ignoreNULL = FALSE)
    
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
      shinyFeedback::hideFeedback("tex")
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
          shinyFeedback::showFeedback("tex", "SNP not found")
          return()
        }
      } else {
        if (nchar(tex) > 1) shinyFeedback::showFeedback("tex", "not found")
        return()
      }
      xr <- as.integer(pmax(xr, 0))
      
      if (chr %in% chr_set[[1]]) {
        coords$chr <- chr
        if (any(is.na(xr))) {
          shinyFeedback::showFeedback("tex", "invalid entry")
          return()
        }
        coords$xrange <- xr
        shinyFeedback::hideFeedback("tex")
      } else {
        shinyFeedback::showFeedback("tex", "not present")
      }
    })
    
    # Table tab
    output$table <- DT::renderDataTable({
      dat <- data
      if (input$tab1_type == "local" && !is.null(coords$chr)) {
        dat <- loc$i$data
        validate(need(!is.null(dat), "No data at locus"))
      }
      fcols <- colnames(dat)[sapply(dat, class) == "numeric"]
      fcols <- fcols[!fcols %in% c(chrom[1], pos[1])]
      DT::datatable(dat, caption = traits[1], rownames = FALSE) %>%
        DT::formatSignif(fcols, digits = 3)
    })
    
    # Extra table tab
    output$link_table <- DT::renderDataTable({
      # a local GWAS Catalog has its own tab, so this one is for eQTL
      if (!is.null(link_data$LDtrait) && is.null(gwas_catalog)) {
        tab <- trait_table(shiny::reactiveValuesToList(link_data), trait_filter(),
                           trait_category_filter(),
                           !isFALSE(input$trait_collapse))
        tab$PMID <- pubmed_links(tab$PMID)
        num <- c("r2", "P")
      } else if (!is.null(link_data$OTqtl)) {
        tab <- qtl_table(shiny::reactiveValuesToList(link_data),
                         qtl_type_filter(), qtl_gene_filter())
        num <- c("P", "Beta", "PIP")
      } else {
        tab <- eqtl_table(shiny::reactiveValuesToList(link_data), eqtl_gene_filter(),
                          eqtl_tissue_filter())
        num <- c("r2", "Effect size", "P")
      }
      req(tab)
      DT::datatable(tab, rownames = FALSE,
                    escape = setdiff(colnames(tab), "PMID")) %>%
        DT::formatSignif(num, digits = 3)
    })
    
    # local GWAS Catalog tab
    output$catalog_table <- DT::renderDataTable({
      validate(need(!is.null(link_data$LDtrait), "GWAS Catalog track is switched off"))
      tab <- trait_table(shiny::reactiveValuesToList(link_data), trait_filter(),
                         trait_category_filter(), !isFALSE(input$trait_collapse))
      tab$PMID <- pubmed_links(tab$PMID)
      DT::datatable(tab, rownames = FALSE,
                    escape = setdiff(colnames(tab), "PMID"),
                    caption = if (!is.null(gwas_catalog$date))
                      paste("GWAS Catalog associations added up to", gwas_catalog$date)) %>%
        DT::formatSignif(c("r2", "P"), digits = 3)
    })
    
    # control extra table tab title
    output$extra_tab_nm <- renderText({
      if (!is.null(link_data$LDtrait) && is.null(gwas_catalog)) return("GWAS Catalog associations")
      if (!is.null(link_data$LDexp)) return("eQTL associations")
      if (!is.null(link_data$OTqtl)) return("QTL credible sets (Open Targets)")
      ""
    })
    
    # show/hide extra table tab
    observe({
      if ((!is.null(link_data$LDtrait) && is.null(gwas_catalog)) ||
          !is.null(link_data$LDexp) || !is.null(link_data$OTqtl)) {
        showTab("topTabs", "linktab")
      } else {
        hideTab("topTabs", "linktab")
      }
    })
    
    # detect change to x axis range
    observeEvent(event_data("plotly_relayout", source = "plotly_locus"), {
      req(coords$chr %in% chr_set[[1]], coords$xrange)
      s <- event_data("plotly_relayout", source = "plotly_locus")
      req(c("xaxis.range[0]", "xaxis.range[1]") %in% names(s))
      xr <- c(s$`xaxis.range[0]`, s$`xaxis.range[1]`)
      xd <- diff(xr) * 0.01 / 1.02
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
      yref <- save_plotly$p$x$data[[ntrace() +2]]$yaxis
      
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
    
    # main manhattan highlight
    observeEvent(coords$chr, {
      req(coords$chr)
      # manhattan1
      if (coords$chr %in% chr_set[[1]]) {
        chrom_range <- manhat$chrom_range[coords$chr, ]
        shapes <- list(
          list(type = "rect",
               line = list(width = 1, color = scheme[4]),
               x0 = chrom_range[1],
               x1 = chrom_range[2], y0 = 0, y1 = 1,
               xref = "x", yref = "paper", layer = "above")
        )
      } else {
        shapes <- list()
      }
      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout",
                          list(shapes = shapes))
      
      # manhattan2 highlight
      if (man2) {
        if (coords$chr %in% chr_set[[2]]) {
          chrom_range <- manhat2$chrom_range[coords$chr, ]
          shapes <- list(
            list(type = "rect",
                 line = list(width = 1, color = scheme2[4]),
                 x0 = chrom_range[1],
                 x1 = chrom_range[2], y0 = 0, y1 = 1,
                 xref = "x", yref = "paper", layer = "above")
          )
        } else {
          shapes <- list()
        }
        plotlyProxy("manhattan2", session) %>%
          plotlyProxyInvoke("relayout",
                            list(shapes = shapes))
      }
    })
    
    # chrom highlight
    observeEvent(coords$xrange, {
      req(input$show_chrom, coords$chr)
      plotlyProxy("chrom", session) %>%
        plotlyProxyInvoke("relayout",
                          list(shapes = list(
                            list(type = "rect",
                                 line = list(width = 1, color = scheme[4]),
                                 x0 = coords$xrange[1] / 1e6,
                                 x1 = coords$xrange[2] / 1e6, y0 = 0, y1 = 1,
                                 xref = "x", yref = "paper", layer = "above")
                          )))
      # chrom2 highlight
      if (man2) {
        plotlyProxy("chrom2", session) %>%
          plotlyProxyInvoke("relayout",
                            list(shapes = list(
                              list(type = "rect",
                                   line = list(width = 1, color = scheme2[4]),
                                   x0 = coords$xrange[1] / 1e6,
                                   x1 = coords$xrange[2] / 1e6, y0 = 0, y1 = 1,
                                   xref = "x", yref = "paper", layer = "above")
                            )))
      }
    })
    
    # retrieve LD
    observeEvent(input$ld_get, {
      snp <- cur_index()
      if (is.null(snp) || is.na(snp)) {
        showNotification("No index SNP in view", type = "warning")
        return()
      }
      showNotification(paste0("Fetching LD for ", snp, " - click a point to re-base"),
                       id = "ld_busy", duration = NULL)
      link_snp(NULL)  # so a repeat press re-renders
      link_snp(snp)
      link_type("LD")
    })
    
    # Clear, eQTL and Open Targets drop LD, so reset the LD-based Catalog options
    observeEvent(list(input$ld_clear, input$eqtl_get, input$ot_get), {
      shiny::updateSliderInput(session, "trait_r2_min", value = 0)
      shiny::updateCheckboxInput(session, "trait_ld_colour", value = FALSE)
      r2_min(0)
      ld_colour(FALSE)
    }, ignoreInit = TRUE)
    
    # fetch LD for the index SNP in view, as Get LD would
    fetch_index_ld <- function(why) {
      snp <- cur_index()
      if (is.null(snp) || is.na(snp)) {
        showNotification("No index SNP in view", type = "warning")
        return()
      }
      showNotification(paste0("Fetching LD for ", snp, " ", why),
                       id = "ld_busy", duration = NULL)
      # re-setting the same SNP would not re-render, e.g. after a failed fetch
      link_snp(NULL)
      link_snp(snp)
      link_type("LD")
    }
    
    # the r2 cut-off and LD colouring need LD, so fetch it if missing; run
    # before output$locus so it renders once, with LD
    has_ld <- function() "ld" %in% colnames(isolate(loc$i$data))
    
    observeEvent(trait_r2_min(), {
      r2_min(as.numeric(trait_r2_min()))
      req(trait_r2_min() > 0, !has_ld())
      ld_colour(TRUE)
      shiny::updateCheckboxInput(session, "trait_ld_colour", value = TRUE)
      fetch_index_ld("for the r\u00b2 cut-off")
    }, priority = 1)
    
    # with a local Catalog, Get LD also colours it by LD (in the same render)
    observeEvent(input$ld_get, {
      req(!is.null(gwas_catalog), !is.null(cur_index()), !is.na(cur_index()))
      ld_colour(TRUE)
      shiny::updateCheckboxInput(session, "trait_ld_colour", value = TRUE)
    }, priority = 1)
    
    observeEvent(input$trait_ld_colour, {
      ld_colour(isTRUE(input$trait_ld_colour))
      req(isTRUE(input$trait_ld_colour), !has_ld())
      fetch_index_ld("to colour by LD")
    }, priority = 1)
    
    observeEvent(input$ld_clear, {
      link_snp(NULL)
      link_type(NULL)
      link_eqtl$genes <- NULL
      link_eqtl$tissues <- NULL
      removeNotification("ld_busy")
    })
    
    # LDlink rebase click on locus
    observe({
      s <- event_data("plotly_click", source = "plotly_locus")
      req(s, !is.null(s$key))
      cur <- isolate(link_snp())
      ltype <- isolate(link_type())
      req(!is.null(cur))
      snp <- as.character(s$key)[1]
      req(!is.na(snp), nzchar(snp))
      # trait-track points can be SNPs absent from the data, which LD can't use
      req(snp %in% isolate(loc$i$data[, labs[1]]) ||
            (man2 && snp %in% isolate(locv2$i$data[, labs[2]])))
      if (identical(snp, cur)) return()
      showNotification(paste0("Fetching ", ltype, " for ", snp),
                       id = "ld_busy", duration = NULL)
      link_snp(snp)
    })
    
    output$ld_status <- renderText({
      req(link_snp())
      if (man2) {
        req(link_snp() %in% c(loc$i$data[, labs[1]], locv2$i$data[, labs[2]]))
      } else {
        req(link_snp() %in% loc$i$data[, labs])
      }
      if (link_type() == "LD") {
        paste0("LD: ", link_snp(), " (", ld_pop, ")")
      } else if (link_type() == "eQTL") {
        "eQTL overlay"
      }
    })
    
    # retrieve eQTL
    observeEvent(input$eqtl_get, {
      snp <- cur_index()
      if (is.null(snp) || is.na(snp)) {
        showNotification("No index SNP in view", type = "warning")
        return()
      }
      showNotification(paste0("Fetching eQTL data"),
                       id = "ld_busy", duration = NULL)
      link_snp(NULL)  # so a repeat press re-renders
      link_snp(snp)
      link_type("eQTL")
    })
    
    # retrieve trait
    observeEvent(input$trait_get, {
      snp <- cur_index()
      if (is.null(snp) || is.na(snp)) {
        showNotification("No index SNP in view", type = "warning")
        return()
      }
      showNotification(paste0("Fetching GWAS trait data"),
                       id = "ld_busy", duration = NULL)
      link_snp(NULL)  # so a repeat press re-renders
      link_snp(snp)
      link_type("GWAS trait")
    })
    
    # retrieve Open Targets QTL credible sets
    observeEvent(input$ot_get, {
      snp <- cur_index()
      if (is.null(snp) || is.na(snp)) {
        showNotification("No index SNP in view", type = "warning")
        return()
      }
      showNotification("Fetching Open Targets QTL data",
                       id = "ld_busy", duration = NULL)
      link_snp(NULL)  # so a repeat press re-renders
      link_snp(snp)
      link_type("Open Targets")
    })
    
    # inline conditional UI
    output$save_ui <- renderUI({
      req(coords$chr)
      downloadButton("save", NULL, icon = icon("floppy-disk"))
    })
    
    # save files
    output$save <- downloadHandler(filename = function() {
      paste0("loc_", dat_name, "_", coords$chr, "_", coords$xrange[1], ".",
             input$export)
    },
    content = function(file) {
      req(loc$i)
      if (input$export == "pdf") {
        pdf(file)
        if (!man2) {
          locus_plot(loc$i, eqtl_gene = eqtl_gene, beta = beta,
                     blanks = "hide")
        } else {
          oldpar <- set_layers(2)
          on.exit(par(oldpar))
          bty <- if (!is.null(recomb) && input$recomb) "u" else "l" 
          scatter_plot(loc$i, xticks = FALSE, bty = bty,
                       ylab = bquote(.(traits[1]) ~ -log[10] ~ P),
                       eqtl_gene = eqtl_gene, beta = beta[1])
          scatter_plot(locv2$i, xticks = FALSE, bty = bty, beta = beta[2],
                       ylab = bquote(.(traits[2]) ~ -log[10] ~ P))
          genetracks(loc$i, blanks = "hide")
        }
        dev.off()
      } else {
        saveRDS(save_plotly$p, file)
      }
    })
    
    # R code for a static figure of the current view
    current_figure_state <- function() isolate({
      pin <- link_snp()
      in_view <- !is.null(pin) &&
        (pin %in% loc$i$data[, labs[1]] || (man2 && pin %in% locv2$i$data[, labs[2]]))
      bt <- input_biotype()
      if (setequal(bt, biotypes)) bt <- NULL
      figure_state(data = data_sym, data2 = if (man2) data2_sym, man2 = man2,
                   data_prep = data_prep, data2_prep = if (man2) data2_prep,
                   qtl_found = !is.null(link_data$OTqtl),
                   ens_db = ens_source, recomb = recomb_sym,
                   recomb_on = !is.null(recomb) && isTRUE(input$recomb),
                   gwas_catalog = if (!is.null(gwas_catalog)) gcat_source,
                   catalog_on = isTRUE(input$gwas_catalog),
                   chrom = chrom, pos = pos, p = p, labs = labs,
                   chr = coords$chr, xrange = coords$xrange,
                   link_type = link_type(), link_snp = pin, link_snp_in_view = in_view,
                   ld_pop = ld_pop, pcutoff = pcutoff, traits = if (man2) traits,
                   eqtl_gene = eqtl_gene, beta = beta,
                   select_gene = input$select_gene, biotype = bt,
                   maxrows = if (isTRUE(input$alltracks)) NULL else 8,
                   trait_filter = trait_filter(),
                   category_filter = trait_category_filter(),
                   collapse = !isFALSE(input$trait_collapse),
                   r2_min = r2_min(), ld_colour = ld_colour(),
                   eqtl_gene_filter = eqtl_gene_filter(),
                   eqtl_tissue_filter = eqtl_tissue_filter(),
                   qtl_type_filter = qtl_type_filter(),
                   qtl_gene_filter = qtl_gene_filter())
    })
    
    output$fig_code_ui <- renderUI({
      req(coords$chr)
      actionButton("fig_code", NULL, icon = icon("code"),
                   title = "R code for this view")
    })
    
    observeEvent(input$fig_code, {
      code <- figure_code(current_figure_state())
      shiny::showModal(shiny::modalDialog(
        tags$pre(id = "fig_code_text", code),
        title = "R code for this view", size = "l", easyClose = TRUE,
        footer = shiny::tagList(
          tags$button("Copy", class = "btn btn-primary",
                      onclick = paste0("navigator.clipboard.writeText(",
                                       "document.getElementById('fig_code_text').innerText)")),
          shiny::modalButton("Close"))))
    })
    
  }
  
  runApp(list(ui = ui, server = server)) %>%
    suppress_warnings("please add `event_register\\(p")
}
