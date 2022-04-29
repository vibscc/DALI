suppressPackageStartupMessages(library(DALI))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinycssloaders))

fillPage(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "css/dali.min.css"),
    ),
    tags$div(
        class = "container-fluid header mb-2 p-0",
        tags$img(src = "images/dali.png", class = "header-logo"),
        tags$div(class = "col-sm-3",
            tags$div(class = "form-group row col-sm-12",
                tags$label("Assay", class = "col-sm-3 text-right col-form-label"),
                tags$div(class = "col-sm-9",
                    tags$select(name = "active.assay", id = "active.assay", class = "form-control rounded-all-90")
                ),
            ),
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
        		        selectInput("chain.usage.color", label = "Colorscheme", choices = c("coolwarm","viridis")),
		                checkboxInput("chain.usage.add.missing.families", label = "Show missing families", value = F),
                        checkboxInput("chain.usage.cluster.cols", label = "Cluster groups based on VDJ genes", value = F)
                        ),
                    column(9, plotOutput("chain.usage.heatmap") %>% withSpinner()))
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
                DT::DTOutput("clonotypes.table"),
                fluidRow(
                    column(4, uiOutput("clonotype.lineage.ui")),
                    column(8, plotOutput("clonotype.lineage"))
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
                            selectizeInput("transcriptomics.feature", label = "Feature", choices = NULL)
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
                                    column(4, selectizeInput("deg.group.by", label = "Metadata column", choices = NULL)),
                                    column(8, selectizeInput("deg.ident.1", label = "Values", multiple = T, choices = NULL))
                                ),
                                fluidRow(
                                    column(4, selectInput("deg.assay", label = "Assay for results", choices = NULL)),
                                    column(4, actionButton("deg.calculate", "Calculate DEG"))
                                )
                            ),
                            div(class = "col-sm-5 well",
                                h4("Specify group 2"),
                                fluidRow(
                                    column(4, radioButtons("deg.ident.2.choice", "", c("All other cells" = 1, "Selected cells" = 2), inline = T)),
                                    column(4, selectizeInput("deg.ident.2", label = "Values", multiple = T, choices = NULL))
                                )
                            )
                        )
                    ),
                    tabPanel("Results", DT::DTOutput("deg.output"))
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
                             uiOutput("dim.reduction.tabs") %>% withSpinner()
                      )
                  )
             ),
             tabPanel("DEG Selector",
                  div(class = "d-flex justify-content-around",
                      div(class = "col-sm-5 well",
                          h4("Specify group 1"),
                          fluidRow(
                              column(4, selectizeInput("deg.group.by.novdj", label = "Metadata column", choices = NULL)),
                              column(8, selectizeInput("deg.ident.1.novdj", label = "Values", multiple = T, choices = NULL))
                          ),
                          fluidRow(
                              column(4, selectInput("deg.assay.novdj", label = "Assay for results", choices = NULL)),
                              column(4, actionButton("deg.calculate.novdj", "Calculate DEG"))
                          )
                      ),
                      div(class = "col-sm-5 well",
                          h4("Specify group 2"),
                          fluidRow(
                              column(4, radioButtons("deg.ident.2.choice.novdj", "", c("All other cells" = 1, "Selected cells" = 2), inline = T)),
                              column(4, selectizeInput("deg.ident.2.novdj", label = "Values", multiple = T, choices = NULL))
                          )
                      )
                  )
             ),
             tabPanel("DEG Results",
                      DT::DTOutput("deg.output.novdj")
            )
        )
    )
)
