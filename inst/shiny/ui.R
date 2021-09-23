suppressPackageStartupMessages(library(Diversity))
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
        tabsetPanel(
            type = "pills",
            tabPanel("General view",
                fluidRow(
                    column(12, uiOutput("reduction.tabs.chain.usage") %>% withSpinner())
                ),
                fluidRow(
                    column(3,
                        selectInput("chain.usage.chain", label = "Chain", choices = c("VDJ", "VJ")),
                        selectInput("chain.usage.region", label = "Region", choices = c("V", "D", "J", "C")),
                        checkboxInput("chain.usage.add.missing.families", label = "Show missing families", value = F)
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
                            column(2,
                               selectInput("clonotype.group.by", label = "Group data by", choices = NULL),
                               selectizeInput("clonotype.group", label = "Group", choices = NULL),
                               sliderInput("cdr3.frequency.threshold", value = 1, min = 0, max = 250, label = "Highlight threshold"),
                               checkboxInput("cdr3.frequency.show.missing", label = "Show cells without VDJ data")
                            )
                        ),
                        fluidRow(
                            column(3, plotOutput("cdr3.frequency")),
                            column(5, tableOutput("top.clonotypes"))
                        )
                    ),
                    column(4,
                        selectizeInput("featureplot.clonotype", label = "Clonotype location", choices = NULL),
                        selectizeInput("featureplot.reduction", label = "Reduction", choices = NULL),
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
                            column(4,
                                selectInput("compare.group.by", label = "Group data by", choices = NULL),
                                selectizeInput("compare.ident.1", label = "Ident 1 (yellow)", choices = NULL, multiple = T),
                                selectizeInput("compare.ident.2", label = "Ident 2 (red)", choices = NULL, multiple = T)
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
                    column(6,
                        fluidRow(
                            column(4, selectInput("transcriptomics.assay", label = "Assay", choices = NULL)),
                            column(4, selectizeInput("transcriptomics.feature", label = "Feature", choices = NULL)),
                            column(4, selectInput("transcriptomics.reduction", label = "Reduction", choices = NULL)),
                        ),
                        plotOutput("transcriptomics.featureplot")
                    ),
                    column(6,
                        selectizeInput("transcriptomics.clonotype", label = "Clonotype", choices = NULL),
                        plotOutput("transcriptomics.clonotype.featureplot")
                    )
                )
            ),
            tabPanel("DEG",
                fluidRow(
                    column(6,
                        h3("Specify group 1"),
                        fluidRow(
                            column(4, selectizeInput("deg.group.by", label = "Metadata column", choices = NULL)),
                            column(8, selectizeInput("deg.ident.1", label = "Values", multiple = T, choices = NULL))
                        ),
                        fluidRow(
                            column(4, selectInput("deg.assay", label = "Assay for results", choices = NULL)),
                            column(4, actionButton("deg.calculate", "Calculate DEG"))
                        )
                    ),
                    column(6,
                       h3("Specify group 2"),
                       fluidRow(
                           column(4, radioButtons("deg.ident.2.choice", "", c("All other cells" = 1, "Selected cells" = 2), inline = T)),
                           column(4, selectizeInput("deg.ident.2", label = "Values", multiple = T, choices = NULL))
                       )
                    )
                ),
                DT::DTOutput("deg.output")
            )
        )
    )


)
