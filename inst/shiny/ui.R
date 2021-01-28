suppressPackageStartupMessages(library(Diversity))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))

fluidPage(
    sidebarLayout(
        sidebarPanel(width = 2,
            selectInput('group.highlight', label = "Group to highlight", choices = NULL),
            uiOutput("dataset.metrics")
        ),
        mainPanel(
            uiOutput("reduction.tabs"),
            fluidRow(
                column(6, plotOutput('barplot')),
                column(6, plotOutput('lineplot'))
            )
        )
    )
)
