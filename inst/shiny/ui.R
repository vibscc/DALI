suppressPackageStartupMessages(library(DALI))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinycssloaders))

fillPage(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "css/dali.min.css"),
    ),
    tags$div(class = "container-fluid header mb-2 p-0",
             tags$img(src = "images/dali.png", class = "header-logo"),
             tags$div(class = "col-sm-4",
                      uiOutput("headerUI"),
                      htmlOutput("dataset.metrics", container = tags$div, class = "metrics")
             )
    ),
    fluidPage(
        tabsetPanel(id  = "VDJ",
                    type = "pills",
                    tabPanel("General view",
                             fluidRow(
                                 column(12, uiOutput("reduction.tabs.chain.usage") %>% withSpinner())
                             ),
                             fluidRow(
                                 sidebarPanel(width = 3,
                                              selectInput("chain.usage.chain", label = "Chain", choices = c("VDJ", "VJ")),
                                              selectInput("chain.usage.region", label = "Region", choices = c("V", "J", "C")),
                                              checkboxInput("chain.usage.add.missing.families", label = "Show missing families", value = F),
                                              checkboxInput("chain.usage.cluster.cols", label = "Cluster groups based on VDJ genes", value = F)
                                 ),
                                 column(9, plotOutput("chain.usage.heatmap") %>% withSpinner())
                             )
                    ),
                    tabPanel("Clone view",
                             fluidRow(
                                 column(12, uiOutput("reduction.tabs.expansion") %>% withSpinner())
                             ),
                             fluidRow(
                                 column(8,
                                        fluidRow(
                                            sidebarPanel(width = 2,
                                                         selectInput("clonotype.group.by", label = "Group data by", choices = NULL),
                                                         selectizeInput("clonotype.group", label = "Group", choices = NULL),
                                                         sliderInput("cdr3.frequency.threshold", value = 1, min = 0, max = 250, label = "Highlight threshold"),
                                                         checkboxInput("cdr3.frequency.show.missing", label = "Show cells without VDJ data")
                                            )
                                        ),
                                        fluidRow(
                                            column(3, plotOutput("cdr3.frequency")),
                                            column(5, tableOutput("top.clonotypes"))
                                        )),
                                 column(4,
                                        sidebarPanel(width = 12,
                                                     selectizeInput("featureplot.clonotype", label = "Clonotype location", choices = NULL),
                                                     selectizeInput("featureplot.reduction", label = "Reduction", choices = NULL),
                                        ),
                                        plotOutput("featureplot.clonotype") %>% withSpinner()
                                 )
                             )

                    ),
                    tabPanel("Population comparison",
                             fluidRow(
                                 column(4, uiOutput("reduction.tabs.comparison") %>% withSpinner()),
                                 column(8,
                                        fluidRow(
                                            column(8, plotOutput("barplot.comparison") %>% withSpinner()),
                                            sidebarPanel(width = 4,
                                                         selectInput("compare.group.by", label = "Group data by", choices = NULL),
                                                         selectizeInput("compare.ident.1", label = "Ident 1 (yellow)", choices = NULL, multiple = T),
                                                         selectizeInput("compare.ident.2", label = "Ident 2 (red)", choices = NULL, multiple = T),
                                                         selectInput("compare.chain", label = "Chain", choices = c("VDJ", "VJ"), multiple = F),
                                                         selectInput("compare.region", label = "Region", choices = c("V", "J", "C"), multiple = F)
                                            )
                                        ),
                                        plotOutput("spectratypeplot") %>% withSpinner()
                                 )
                             )
                    ),
                    tabPanel("Clonotypes",
                             tabsetPanel(
                                 tabPanel("Table & Lineage",
                                          DT::DTOutput("clonotypes.table"),
                                          fluidRow(
                                              column(4, uiOutput("clonotype.lineage.ui")),
                                              column(4, plotOutput("clonotype.lineage")),
                                              column(4, plotOutput("clonotype.smh") %>% withSpinner())
                                          ),
                                 ),
                                 tabPanel("Family & Chain useage",
                                          fluidRow(
                                              column(width = 6,
                                                     plotOutput("circosplot.genes", height = 700, width = 700) %>% withSpinner(),
                                              ),
                                              column(width = 6,
                                                     plotOutput("circosplot.chains", height = 700, width = 700) %>% withSpinner(),
                                              ),

                                          ),
                                          fluidRow(
                                              column(width = 6,
                                                     sidebarPanel( width = 10,
                                                                   h4("Family to V-gene"),
                                                                   fluidRow(
                                                                       selectInput("subsetby.circos.genes", label = "Group data by ", choices = NULL),
                                                                       selectInput("gene.subset.group", label = "Select group: ", choices = NULL)
                                                                   ),
                                                     ),
                                              ),
                                              column(width = 6,
                                                     sidebarPanel(width = 10,
                                                                  h4("V D J chain useage"),
                                                                  fluidRow(
                                                                      selectInput("subsetby.circos.chains", label = "Group data by ", choices = NULL),
                                                                      selectInput("chain.subset.group", label = "Select group: ", choices = NULL )
                                                                  ),
                                                     ),
                                              )
                                          )
                                 ),
                             )
                    ),
                    tabPanel("Transcriptomics",
                             fluidRow(
                                 column(12,
                                        sidebarPanel(width = 6,
                                                     fluidRow(
                                                         column(6, selectInput("transcriptomics.assay", label = "Assay", choices = NULL)),
                                                         column(6, selectInput("transcriptomics.reduction", label = "Reduction", choices = NULL))
                                                     ),
                                                     fluidRow(
                                                        column(6, selectizeInput("transcriptomics.feature", label = "Feature", choices = NULL)),
                                                        column(6, checkboxInput("order", label = "Set positive cells in foreground", value = F))
                                                     )
                                        ),
                                        sidebarPanel(width = 6,
                                                     selectizeInput("transcriptomics.clonotype", label = "Clonotype", choices = NULL)
                                        )
                                 )
                             ),
                             fluidRow(
                                 column(6,
                                        plotOutput("transcriptomics.featureplot") %>% withSpinner()
                                 ),
                                 column(6,
                                        plotOutput("transcriptomics.clonotype.featureplot") %>% withSpinner()
                                 )
                             )
                    ),
                    tabPanel("DEG",
                             tabsetPanel(
                                 tabPanel("Select groups",
                                          div(class = "d-flex justify-content-around",
                                              div(class = "col-sm-5 well",
                                                  h4("Specify group 1"),
                                                  fluidRow(
                                                      column(6, sliderInput(
                                                          "sig.P",
                                                          "significant P-value",
                                                          0.01,
                                                          0.1,
                                                          0.05)
                                                      ),
                                                      column(6, sliderInput(
                                                          "sig.logFC",
                                                          "significant Log2(FC)",
                                                          0,
                                                          5,
                                                          0.6,
                                                          step = 0.01)
                                                      ),
                                                  ),
                                                  fluidRow(
                                                      column(4, selectInput("deg.assay", label = "Assay for results", choices = NULL)),
                                                      column(8, selectizeInput("deg.group.by", label = "Metadata column", choices = NULL))
                                                  ),
                                                  fluidRow(
                                                      column(12, selectizeInput("deg.ident.1", label = "Values", multiple = T, choices = NULL))
                                                  )
                                              ),
                                              div(class = "col-sm-5 well",
                                                  h4("Specify group 2"),
                                                  fluidRow(
                                                      column(4, radioButtons("deg.ident.2.choice", "", c("All other cells" = 1, "Selected cells" = 2), inline = T)),
                                                      column(8, selectizeInput("deg.ident.2", label = "Values", multiple = T, choices = NULL))
                                                  ),
                                                  fluidRow(
                                                      column(4, actionButton("deg.calculate", "Calculate DEG"))
                                                  )
                                              )
                                          )
                                 ),
                                 tabPanel("Results",
                                          DT::DTOutput("deg.output"),
                                          tags$br(),
                                          plotOutput("deg.volcano.plot")
                                 )
                             )
                    )
        ),
        tabsetPanel(id = "Seurat",
                    type = "pills",
                    tabPanel("Clustering & Transcriptomics",
                             fluidRow(
                                 column(12,
                                        column(width = 6,
                                               plotOutput("transcriptomics.featureplot.novdj") %>% withSpinner()
                                        ),
                                        sidebarPanel(width = 6,
                                                     fluidRow(
                                                         column(6, selectInput("transcriptomics.assay.novdj", label = "Assay", choices = NULL)),
                                                         column(6, selectInput("transcriptomics.reduction.novdj", label = "Reduction", choices = NULL))
                                                     ),
                                                     selectizeInput("transcriptomics.feature.novdj", label = "Feature", choices = NULL)
                                        )
                                 )
                             ),
                             fluidRow(
                                 column(6,
                                        plotOutput("dim.reduction") %>% withSpinner()
                                 )
                             )
                    ),
                    tabPanel("DEG Selector",
                             div(class = "d-flex justify-content-around",
                                 div(class = "col-sm-5 well",
                                     h4("Specify group 1"),
                                     fluidRow(
                                         column(6, sliderInput(
                                             "sig.P.novdj",
                                             "significant P-value",
                                             0.01,
                                             0.1,
                                             0.05)
                                         ),
                                         column(6, sliderInput(
                                             "sig.logFC.novdj",
                                             "significant Log2(FC)",
                                             0,
                                             5,
                                             0.6,
                                             step = 0.01)
                                         )
                                     ),
                                     fluidRow(
                                         column(4, selectInput("deg.assay.novdj", label = "Assay for results", choices = NULL)),
                                         column(8, selectizeInput("deg.group.by.novdj", label = "Metadata column", choices = NULL)),
                                     ),
                                     fluidRow(
                                         column(12, selectizeInput("deg.ident.1.novdj", label = "Values", multiple = T, choices = NULL))
                                     )
                                 ),
                                 div(class = "col-sm-5 well",
                                     h4("Specify group 2"),
                                     fluidRow(
                                         column(4, radioButtons("deg.ident.2.choice.novdj", "", c("All other cells" = 1, "Selected cells" = 2), inline = T)),
                                         column(8, selectizeInput("deg.ident.2.novdj", label = "Values", multiple = T, choices = NULL))
                                     ),
                                     fluidRow(
                                         column(4, actionButton("deg.calculate.novdj", "Calculate DEG"))
                                     )
                                 )
                             )
                    ),
                    tabPanel("DEG Results",
                             DT::DTOutput("deg.output.novdj"),
                             tags$br(),
                             plotOutput("deg.volcano.plot.novdj")
                    )
        )
    )
)
