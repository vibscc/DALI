library(promises)
library(future)
library(shinyFiles)
plan(multisession)

function(input, output, session) {

    # ======================================================================= #
    # Initialize app
    # ======================================================================= #
    volumes <- c("Home (~)" = fs::path_home(), "/" = "/")
    shinyFileChoose(input, "seurat_rds", session = session, roots = volumes)
    shinyDirChoose(input, "vdj_dir", session = session, roots = volumes)
    vals <- reactiveValues(data = .GlobalEnv$.data.object.VDJ)
    upload <- reactiveValues(
        seurat.rds = NULL,
        vdj.dir = NULL
    )

    app.initialize <- function() {
        renderReductionPlotsChainUsage(isolate(vals$data))
        renderReductionPlotsExpansion(isolate(vals$data))
        renderReductionPlotsComparison(isolate(vals$data))

        groups <- levels(isolate(vals$data@meta.data$seurat_clusters))
        updateSelectInput(session, "group.highlight", choices = groups)

        # updateSelectInput(session, "compare.group.by", choices = colnames(isolate(vals$data@meta.data)), selected = "seurat_clusters")
        updateSelectizeInput(session, "compare.ident.1", choices = groups)
        updateSelectizeInput(session, "compare.ident.2", choices = groups)

        assays <- names(isolate(vals$data@misc$VDJ))
        updateSelectInput(session, "active.assay", choices = assays, selected = DefaultAssayVDJ(isolate(vals$data)))

        metadata.columns <- colnames(isolate(vals$data@meta.data))
        updateSelectInput(session, "group.by", choices = metadata.columns, selected = "seurat_clusters")

        updateSelectInput(session, "chain.usage.chain", choices = Diversity:::AvailableChainsList(isolate(vals$data)))
        selected <- if ('umap' %in% names(isolate(vals$data@reductions))) 'umap' else if ('tsne' %in% names(isolate(vals$data@reductions))) 'tsne' else NULL
        updateSelectizeInput(session, "featureplot.reduction", choices = names(isolate(vals$data@reductions)), selected = selected)
    }

    # ======================================================================= #
    # Function definitions
    # ======================================================================= #

    # Render DimPlots for each reduction
    renderReductionPlotsChainUsage <- function(object) {
        for (reduction in names(object@reductions)) {
            # Make all data available in a local scope, since plots are not rendered instantly.
            # Without the local scope, each plot would be identical to the last plot
            local({
                r <- reduction

                output[[paste0('reduction.plot.', r)]] <- renderPlot({
                    Seurat::DimPlot(
                        object,
                        reduction = r,
                        group.by = "seurat_clusters",
                        cols = Diversity:::GetCategoricalColorPalette(object@meta.data$seurat_clusters)
                    ) + theme(
                        axis.line = element_blank(),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text = element_blank()
                    ) + ggtitle("Clustering")
                })

                output[[paste0('reduction.plot.vdj.', r)]] <- renderPlot({
                    DimplotChainRegion(
                        object,
                        grid = F,
                        reduction = r,
                        chain = input$chain.usage.chain,
                        region = input$chain.usage.region
                    ) + theme(
                        axis.line = element_blank(),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text = element_blank()
                    ) + ggtitle("Chain usage")
                })
            })
        }
    }

    renderReductionPlotsExpansion <- function(object) {
        for (reduction in names(object@reductions)) {
            local({
                r <- reduction

                output[[paste0('expansion.reduction.plot.', r)]] <- renderPlot({
                    Seurat::DimPlot(
                        object,
                        reduction = r,
                        group.by = "seurat_clusters"
                    ) + theme(
                        axis.line = element_blank(),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text = element_blank()
                    ) + ggtitle("Clustering")
                })

                output[[paste0('expansion.reduction.plot.', r, '.exp')]] <- renderPlot({
                    ExpansionPlot(
                        object,
                        reduction = r,
                        threshold = 2,
                        negative.alpha = 0.7
                    ) + theme(
                        axis.line = element_blank(),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text = element_blank(),
                        plot.title = element_text(hjust = 0.5, face = "bold", vjust = 1, size = 16, margin = margin(0,0,7,7))
                    ) + ggtitle("Expansion plot")
                })

                # output[[paste0('graph.', r)]] <- renderPlot({
                #     CloneConnGraph(
                #         object = object,
                #         reduction = r
                #     ) + theme(
                #         legend.position = "none"
                #     ) + ggtitle("Clonal connection plot")
                #
                # })
            })
        }
    }

    renderReductionPlotsComparison <- function(object) {
        for (reduction in names(object@reductions)) {
            local({
                r <- reduction

                output[[paste0('dimred.', r)]] <- renderPlot({
                    Seurat::DimPlot(
                        object,
                        reduction = r
                    ) + theme(
                        axis.line = element_blank(),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text = element_blank()
                    )
                })
            })
        }
    }

    # ======================================================================= #
    # Load data if necessary
    # ======================================================================= #

    dataUploadModal <- function(error = NULL) {
        modalDialog(
            div(
                strong("Select a Seurat Rds file:"),
                shinyFilesButton("seurat_rds", "Browse...", "Choose an Rds file to load",  multiple = F, filetype = list(data = c("Rds", "rds"))),
                textOutput("seurat.rds.path.text", inline = T),
            ),
            div(
                strong("(Optional) Cellranger output folder", title = "This is only required when the selected Seurat object does NOT contain VDJ data yet (loaded by Diversity::Read10X_vdj() or other Diversity method)"),
                shinyDirButton("vdj_dir", "Browse...", "Select cellranger output folder containing VDJ data", multiple = F),
                textOutput("vdj.dir.text", inline = T)
            ),
            div(
                uiOutput("vdj.type.select")
            ),
            if (!is.null(error)) {
                div(
                    h4("Error:", class = "text-danger"),
                    p(error)
                )
            },
            title = "Load data",
            footer = tagList(
                actionButton("load", "Load")
            ),
            easyClose = F,
            size = "l"
        )
    }

    loadingModal <- function() {
        modalDialog(
            "Data is loading. This can take a while for larger datasets!",
            footer = NULL,
            easyClose = F
        )
    }

    if (is.null(isolate(vals$data)) || !IsValidSeuratObject(isolate(vals$data))) {
        showModal(dataUploadModal())
    } else {
        app.initialize()
    }

    observeEvent(input$seurat_rds, {
        upload$seurat.rds <- shinyFiles::parseFilePaths(volumes, input$seurat_rds)
    })

    observeEvent(input$vdj_dir, {
        upload$vdj.dir <- shinyFiles::parseDirPath(volumes, input$vdj_dir)
    })

    observeEvent(input$load, {
        if (length(upload$seurat.rds) == 0 || is.null(upload$seurat.rds)) {
            showModal(dataUploadModal(error = "Missing Seurat Rds file!"))
            return()
        }

        if (!file.exists(upload$seurat.rds$datapath)) {
            showModal(dataUploadModal(error = paste0("Could not find file '", upload$seurat.rds$datapath, "'")))
        }

        removeModal()
        showModal(loadingModal())

        data <- readRDS(upload$seurat.rds$datapath)

        if (length(upload$vdj.dir) > 0) {
            data <- tryCatch({
                Read10X_vdj(data, upload$vdj.dir, type = input$vdj_type)
            }, error = function(e) {
                showModal(dataUploadModal(error = paste0(e, " Make sure the selected VDJ data matches the Seurat object + the correct VDJ type is selected")))
                return(NA)
            })
            if (is.na(data)) {
                return()
            }
        }


        if (IsValidSeuratObject(data)) {
            vals$data <- data
            removeModal()
            app.initialize()
            return()
        } else {
            showModal(
                dataUploadModal(
                    error = "Selected Seurat-object does not contain VDJ (correct) information. Please select a cellranger output folder for the given Seurat object AND make sure the VDJ type is set correctly."
                )
            )
        }
    })

    output$seurat.rds.path.text <- renderText({
        req(upload$seurat.rds)

        gsub("/+", "/", upload$seurat.rds$datapath)
    })

    output$vdj.dir.text <- renderText({
        req(upload$vdj.dir)

        gsub("/+", "/", upload$vdj.dir)
    })

    output$vdj.type.select <- renderUI({
        req(upload$vdj.dir)

        selectInput("vdj_type", label = "VDJ data type", choices = c("TCR", "BCR"))
    })

    # ======================================================================= #
    # Reduction plots UI tabs
    # ======================================================================= #

    output$dataset.metrics <- renderUI({
        req(vals$data)

        if ("h.v_gene" %in% colnames(vals$data@meta.data)) {
            cells.with.VDJ <- vals$data@meta.data %>% filter(!is.na(.data$h.v_gene) | !is.na(.data$l.v_gene) ) %>% nrow()
        }
        if ("a.v_gene" %in% colnames(vals$data@meta.data)) {
            cells.with.VDJ <- vals$data@meta.data %>% filter(!is.na(.data$a.v_gene) | !is.na(.data$b.v_gene) ) %>% nrow()
        }

        list(
            div("# cells: ", strong(ncol(vals$data))),
            div("# cells with VDJ info: ", strong(cells.with.VDJ))
        )
    })

    # ======================================================================= #
    # Reduction plots UI tabs
    # ======================================================================= #

    # Create tabsetPanel with tabPanel for each dimensionality reducion in the dataset
    output$reduction.tabs.chain.usage <- renderUI({
        req(vals$data)

        tabs <- lapply(names(vals$data@reductions), function(reduction) {
            plotname.dimred <- paste0('reduction.plot.', reduction)
            plotname.dimred.vdj <- paste0('reduction.plot.vdj.', reduction)

            tabPanel(
                Diversity:::FormatDimred(reduction),
                fluidRow(
                    column(6, plotOutput(plotname.dimred) %>% withSpinner()),
                    column(6, plotOutput(plotname.dimred.vdj) %>% withSpinner())
                )
            )
        })

        # TODO: select default tab in a more elegant way. PCA should be avoided as default tab, since this is the least informative
        tabs[['selected']] <- if ('umap' %in% names(vals$data@reductions)) 'UMAP' else if ('tsne' %in% names(vals$data@reductions)) 'tSNE' else NULL
        do.call(tabsetPanel, tabs)
    })

    output$reduction.tabs.expansion <- renderUI({
        req(vals$data)

        tabs <- lapply(names(vals$data@reductions), function(reduction) {
            plotname.dimred <- paste0('expansion.reduction.plot.', reduction)
            plotname.dimred.expansion <- paste0('expansion.reduction.plot.', reduction, '.exp')
            plotname.graph <- paste0('graph.', reduction)

            tabPanel(
                Diversity:::FormatDimred(reduction),
                fluidRow(
                    column(6, plotOutput(plotname.dimred) %>% withSpinner()),
                    column(6, plotOutput(plotname.dimred.expansion) %>% withSpinner()),
                    # column(4, plotOutput(plotname.graph) %>% withSpinner())
                )
            )
        })

        tabs[['selected']] <- if ('umap' %in% names(vals$data@reductions)) 'UMAP' else if ('tsne' %in% names(vals$data@reductions)) 'tSNE' else NULL
        do.call(tabsetPanel, tabs)
    })

    output$reduction.tabs.comparison <- renderUI({
        req(vals$data)

        tabs <- lapply(names(vals$data@reductions), function(reduction) {
            plotname <- paste0('dimred.', reduction)

            tabPanel(
                Diversity:::FormatDimred(reduction),
                plotOutput(plotname) %>% withSpinner()
            )
        })

        tabs[['selected']] <- if ('umap' %in% names(vals$data@reductions)) 'UMAP' else if ('tsne' %in% names(vals$data@reductions)) 'tSNE' else NULL
        do.call(tabsetPanel, tabs)
    })

    # ======================================================================= #
    # Chain usage
    # ======================================================================= #

    output$chain.usage.barplot <- renderPlot({
        req(vals$data, input$chain.usage.chain, input$chain.usage.region)

        BarplotChainRegion(
            vals$data,
            chain = input$chain.usage.chain,
            region = input$chain.usage.region,
            add.missing.families = input$chain.usage.add.missing.families,
            # group.by = input$chain.usage.group.by
        ) + ggtitle("Chain usage")
    })

    # ======================================================================= #
    # Ridge CDR3-length
    # ======================================================================= #

    output$spectratypeplot <- renderPlot({
        req(vals$data, input$compare.group.by)

        CDR3Plot(
            vals$data,
            group.by = input$compare.group.by,
            sequence.type = "AA",
            plot.type = "ridge"
        )
    })

    # ======================================================================= #
    # Frequency CDR3 AA sequences
    # ======================================================================= #

    output$cdr3.frequency <- renderPlot({
        req(vals$data, input$clonotype.group.by, input$clonotype.group, input$cdr3.frequency.threshold)

        ClonotypeFrequency(
            vals$data,
            chain = NULL,
            use.sequence = F,
            group.by = input$clonotype.group.by,
            subset = input$clonotype.group,
            threshold = input$cdr3.frequency.threshold,
            show.missing = input$cdr3.frequency.show.missing
        )
    })

    output$top.clonotypes <- renderTable({
        req(vals$data, input$clonotype.group.by, input$clonotype.group)

        n.cells <- sum(vals$data@meta.data[input$clonotype.group.by] == input$clonotype.group)

        top.clonotypes <- Diversity:::CalculateFrequency(vals$data, 'clonotype', input$clonotype.group.by, F) %>%
            filter(.data[[input$clonotype.group.by]] == input$clonotype.group) %>%
            arrange(desc(freq)) %>%
            select(c(clonotype, freq)) %>%
            head(n = 25) %>%
            mutate(perc = round((freq/n.cells) * 100, digits = 1))

        vals$top.clonotypes <- top.clonotypes

        h_seqs <- c()
        l_seqs <- c()
        for (clonotype in top.clonotypes$clonotype) {
            h_seqs <- c(h_seqs, Diversity:::ClonotypeToSequence(vals$data, clonotype, "H"))
            l_seqs <- c(l_seqs, Diversity:::ClonotypeToSequence(vals$data, clonotype, "L"))
        }

        top.clonotypes$h_seq <- h_seqs
        top.clonotypes$l_seq <- l_seqs

        colnames(top.clonotypes) <- c("Clonotype", "Cells", "Pct.group", "H CDR3 AA seq", "L CDR3 AA seq")
        top.clonotypes
    })

    # ======================================================================= #
    # Featureplot clonotype
    # ======================================================================= #

    output$featureplot.clonotype <- renderPlot({
        req(vals$data, input$featureplot.clonotype, input$featureplot.reduction)

        FeaturePlotChainRegion(vals$data, input$featureplot.reduction, input$featureplot.clonotype) + theme(
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            plot.title = element_blank()
        )
    })

    # ======================================================================= #
    # Observers
    # ======================================================================= #

    # Update ident choices on group.by change

    observeEvent(input$group.by, {
        req(vals$data, input$group.by)

        updateSelectizeInput(session, "group.highlight", choices = NULL, selected = NULL)
        groups <- isolate(vals$data@meta.data[, input$group.by]) %>% as.character() %>% unique() %>% gtools::mixedsort(x = .)
        updateSelectizeInput(session, "group.highlight", choices = groups)
    }, priority = 10)

    observeEvent(input$compare.group.by, {
        req(vals$data, input$compare.group.by)

        groups <- vals$data@meta.data[, input$compare.group.by] %>% as.character() %>% unique() %>% gtools::mixedsort(x = .)
        updateSelectizeInput(session, "compare.ident.1", choices = groups, selected = groups[[1]])
        updateSelectizeInput(session, "compare.ident.2", choices = groups, selected = groups[[2]])
    })

    observeEvent(input$clonotype.group.by, {
        req(vals$data, input$clonotype.group.by)

        groups <- vals$data@meta.data[, input$compare.group.by] %>% as.character() %>% unique() %>% gtools::mixedsort(x = .)
        updateSelectizeInput(session, "clonotype.group", choices = groups, selected = groups[[1]])
    })

    # Update available regions on chain change

    observeEvent(input$chain.usage.chain, {
        req(vals$data, input$chain.usage.chain)

        updateSelectInput(session, "chain.usage.region", choices = Diversity:::AvailableRegions(input$chain.usage.chain))
    })

    # Top clonotypes change

    observeEvent(vals$top.clonotypes, {
        updateSelectInput(session, "featureplot.clonotype", choices = vals$top.clonotypes$clonotype)
    })

    # ======================================================================= #
    # Barplot to compare groups
    # ======================================================================= #

    output$barplot.comparison <- renderPlot({
        req(vals$data, input$compare.ident.1, input$compare.ident.2)

        BarplotChainRegion(
            vals$data,
            group.by = input$compare.group.by,
            ident.1 = input$compare.ident.1,
            ident.2 = input$compare.ident.2,
            legend = F
        )
    })
}
