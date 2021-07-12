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
                tags$label("Assay", class = "col-sm-3 text-right"),
                tags$div(class = "col-sm-9",
                    tags$select(name = "active.assay", id = "active.assay", class = "form-control rounded-all-90")
                ),
            ),
            htmlOutput('dataset.metrics', container = tags$div, class = "metrics")
        )
    ),
    fluidPage(
        tabsetPanel(
            type = "pills",
            tabPanel("General view",
                fluidRow(
                    column(12, uiOutput('reduction.tabs.chain.usage') %>% withSpinner())
                ),
                fluidRow(
                    column(3,
                        selectInput('chain.usage.chain', label = "Chain", choices = NULL),
                        selectInput('chain.usage.region', label = "Region", choices = c("V", "D", "J", "C")),
                        checkboxInput('chain.usage.add.missing.families', label = "Show missing families", value = F)
                    ),
                    column(9, plotOutput('chain.usage.barplot') %>% withSpinner())
                )
            ),
            tabPanel("Clone view",
                fluidRow(
                    column(12, uiOutput('reduction.tabs.expansion') %>% withSpinner())
                ),
                fluidRow(
                    column(8,
                        fluidRow(
                            column(4,
                                   selectInput("clonotype.group.by", label = "Group data by", choices = list("seurat_clusters")),
                                   selectizeInput("clonotype.group", label = "Group", choices = NULL)
                            ),
                            column(4,
                                   sliderInput("cdr3.frequency.threshold", value = 1, min = 0, max = 250, label = "Highlight threshold"),
                                   checkboxInput("cdr3.frequency.show.missing", label = "Show cells without VDJ data")
                            )
                        ),
                        fluidRow(
                            column(3, plotOutput('cdr3.frequency')),
                            column(9, tableOutput('top.clonotypes'))
                        )
                    ),
                    column(4,
                        selectizeInput('featureplot.clonotype', label = "Clonotype location", choices = NULL),
                        selectizeInput('featureplot.reduction', label = "Reduction", choices = NULL),
                        plotOutput('featureplot.clonotype') %>% withSpinner()
                    )
                )
            ),
            tabPanel("Population comparison",
                fluidRow(
                    column(4, uiOutput('reduction.tabs.comparison') %>% withSpinner()),
                    column(8,
                        fluidRow(
                            column(8, plotOutput('barplot.comparison') %>% withSpinner()),
                            column(4,
                                selectInput("compare.group.by", label = "Group data by", choices = list("seurat_clusters")),
                                selectizeInput("compare.ident.1", label = "Ident 1 (red)", choices = NULL, multiple = T),
                                selectizeInput("compare.ident.2", label = "Ident 2 (blue)", choices = NULL, multiple = T)
                            )
                        ),
                        plotOutput('spectratypeplot') %>% withSpinner()
                    )
                )
            )
        )
    )


)
