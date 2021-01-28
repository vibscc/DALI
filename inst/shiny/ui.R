library(Diversity)
library(htmlwidgets)
library(plotly)
library(shiny)

fluidPage(
    sidebarLayout(
        sidebarPanel(width = 2,
            selectInput('group.highlight', label = "Group to highlight", choices = NULL)
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
