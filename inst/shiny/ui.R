suppressPackageStartupMessages(library(Diversity))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))

fluidPage(
    navbarPage("Diversity",
        tabPanel("Exploration",
            sidebarLayout(
                sidebarPanel(width = 2,
                    h4("Scatterplot"),
                    selectInput('scatterplot.chain', label = "Chain", choices = list("Heavy" = "H", "Light" = "L")),
                    selectInput('scatterplot.region', label = "Region", choices = c("V", "D", "J", "C")),
                    checkboxInput('scatterplot.by.family', label = 'By family', value = T),
                    h4("VH barplot | CDR3 length"),
                    selectInput('group.highlight', label = "Group to highlight", choices = NULL),
                    uiOutput("dataset.metrics")
                ),
                mainPanel(
                    uiOutput("reduction.tabs"),
                    fluidRow(
                        column(6, plotOutput('barplot')),
                        column(6, plotOutput('cdr3.length'))
                    ),
                    fluidRow(
                        column(6, plotOutput('circosplot')),
                        column(6, plotOutput('cdr3.frequency'))
                    )
                )
            )
        ),
        tabPanel("Comparison",
            sidebarLayout(
                sidebarPanel(width = 2,
                    selectInput("compare.group.by", label = "Group data by", choices = list("seurat_clusters")),
                    # selectInput("compare.group.by", label = "Group data by", choices = NULL),
                    selectizeInput("compare.ident.1", label = "Ident 1", choices = NULL, multiple = T),
                    selectizeInput("compare.ident.2", label = "Ident 2", choices = NULL, multiple = T),
                    checkboxInput("compare.grid", label = "Show as grid", value = F),
                    checkboxInput("compare.legend", label = "Show legend", value = F)
                ),
                mainPanel(
                    plotOutput('barplot.comparison')
                )
            )
        )
    )
)
