library(shiny)
library(plotly)
library(htmlwidgets)

fluidPage(
    sidebarLayout(
        sidebarPanel(width = 2,
            textInput('group.highlight', value = 0, label = "Group to highlight")
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
