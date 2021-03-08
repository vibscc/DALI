suppressPackageStartupMessages(library(Diversity))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinycssloaders))

fluidPage(
    navbarPage("Diversity",
        tabPanel("Explore",
            sidebarLayout(
                sidebarPanel(width = 2,
                    selectInput('active.assay', label = "Assay", choices = NULL),
                    selectInput('scatterplot.chain', label = "Chain", choices = list("Heavy" = "H", "Light" = "L","Alpha" = "A", "Beta" = "B")),
                    selectInput('scatterplot.region', label = "Region", choices = c("V", "D", "J", "C")),
                    checkboxInput('scatterplot.by.family', label = 'By family', value = T),
                    selectizeInput('featureplot.clonotype', label = "Clonotype", choices = NULL),
                    uiOutput("dataset.metrics")
                ),
                mainPanel(
                    fluidRow(
                        column(6, uiOutput("reduction.tabs")  %>% withSpinner()),
                        column(6, plotOutput("featureplot.clonotype") %>% withSpinner())
                    ),
                    fluidRow(
                        column(6, plotOutput('circosplot') %>% withSpinner()),
                        column(6, plotOutput('cdr3.frequency') %>% withSpinner())
                    )
                )
            )
        ),
        tabPanel("Explore subset",
            sidebarLayout(
                sidebarPanel(width = 2,
                    selectInput('group.by', label = "Group by", choices = NULL),
                    selectInput('group.highlight', label = "Group to highlight", choices = NULL),
                    selectInput('sequence.type', label = "Sequence", choices = c("AA", "NT"))
                ),
                mainPanel(
                    # uiOutput("reduction.tabs"),
                    fluidRow(
                        column(6, plotOutput('barplot') %>% withSpinner()),
                    ),
                    fluidRow(
                        column(6, plotOutput('cdr3.frequency.subset') %>% withSpinner()),
                        column(6, plotOutput('cdr3.length') %>% withSpinner())
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
                    plotOutput('barplot.comparison') %>% withSpinner()
                )
            )
        )
    )
)
