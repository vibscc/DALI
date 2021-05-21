suppressPackageStartupMessages(library(Diversity))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinycssloaders))

fluidPage(
    navbarPage("Diversity",
        tabPanel("Explore data",
            sidebarLayout(
                sidebarPanel(width = 2,
                    selectInput('active.assay', label = "Assay", choices = NULL),
                    uiOutput('dataset.metrics')
                ),
                mainPanel(
                    tabsetPanel(
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
                                    selectInput("clonotype.group.by", label = "Group data by", choices = list("seurat_clusters")),
                                    selectizeInput("clonotype.group", label = "Group", choices = NULL),
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
        )
    )
)
