suppressPackageStartupMessages(library(Diversity))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))

fluidPage(
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
)
