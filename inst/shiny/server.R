function(input, output, session) {

    # ======================================================================= #
    # Initialize app
    # ======================================================================= #

    vals <- reactiveValues(data = .GlobalEnv$.data.object.VDJ)

    app.initialize <- function() {
        renderReductionPlots(isolate(vals$data))

        groups <- levels(isolate(vals$data@meta.data$seurat_clusters))
        updateSelectInput(session, "group.highlight", choices = groups)

        # updateSelectInput(session, "compare.group.by", choices = colnames(isolate(vals$data@meta.data)), selected = "seurat_clusters")
        updateSelectizeInput(session, "compare.ident.1", choices = groups)
        updateSelectizeInput(session, "compare.ident.2", choices = groups)

        assays <- names(isolate(vals$data@misc$VDJ))
        updateSelectInput(session, "active.assay", choices = assays, selected = DefaultAssayVDJ(isolate(vals$data)))

        metadata.columns <- colnames(isolate(vals$data@meta.data))
        updateSelectInput(session, "group.by", choices = metadata.columns, selected = "seurat_clusters")

        clonotypes <- isolate(vals$data@meta.data$clonotype) %>% unique() %>% gtools::mixedsort(.)
        updateSelectizeInput(session, "featureplot.clonotype", choices = clonotypes, server = T)
    }

    # ======================================================================= #
    # Function definitions
    # ======================================================================= #

    # Render DimPlots for each reduction
    renderReductionPlots <- function(object) {
        for (reduction in names(object@reductions)) {
            # Make all data available in a local scope, since plots are not rendered instantly.
            # Without the local scope, each plot would be identical to the last plot
            local({
                r <- reduction

                plotname <- paste0('reduction.plot.', r)
                output[[plotname]] <- renderPlotly({
                    ggplotly(DimPlot_vh(
                        object,
                        grid = F,
                        reduction = r,
                        chain = input$scatterplot.chain,
                        region = input$scatterplot.region,
                        by.family = input$scatterplot.by.family)
                    ) %>% onRender("
                        function(el) {
                            el.on('plotly_legendclick', function(d) {

                                // const input = document.getElementById('group.highlight');

                                // Create fake keyup-event to trigger shiny
                                // const ev = document.createEvent('Event');
                                // ev.initEvent('keyup');
                                // ev.which = ev.keyCode = 13;

                                // input.value = d.curveNumber;
                                // input.dispatchEvent(ev);

                                // Prevent other handlers from firing
                                // return false;
                            })
                        }
                    ")
                })
            })
        }
    }

    # ======================================================================= #
    # Load data if necessary
    # ======================================================================= #

    dataUploadModal <- function(failed = F) {
        modalDialog(
            fileInput("file", "Choose an rds file to load", accept = ".rds"),
            if (failed) {
                div("Invalid file! Make sure the object already contains VDJ data. This can be done with the `Read10X_vdj` function of `Diversity`")
            },

            footer = tagList(
                actionButton("load", "Load")
            ),
            easyClose = F
        )
    }

    if (is.null(isolate(vals$data)) || !isValidSeuratObject(isolate(vals$data))) {
        showModal(dataUploadModal())
    } else {
        app.initialize()
    }

    observeEvent(input$load, {
        data <- readRDS(input$file$datapath)

        if (isValidSeuratObject(data)) {
            vals$data <- readRDS(input$file$datapath)
            app.initialize()
            removeModal()
        } else {
            showModal(dataUploadModal(failed = T))
        }
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
            h4("Dataset metrics"),
            div("# cells: ", ncol(vals$data)),
            div("# cells with VDJ info: ", cells.with.VDJ)
        )
    })

    # ======================================================================= #
    # Reduction plots UI tabs
    # ======================================================================= #

    # Create tabsetPanel with tabPanel for each dimensionality reducion in the dataset
    output$reduction.tabs <- renderUI({
        req(vals$data)

        tabs <- lapply(names(vals$data@reductions), function(reduction) {
            plotname <- paste0('reduction.plot.', reduction)
            tabPanel(
                Diversity:::formatDimred(reduction),
                plotlyOutput(plotname) %>% withSpinner()
            )
        })

        # TODO: select default tab in a more elegant way. PCA should be avoided as default tab, since this is the least informative
        tabs[['selected']] <- if ('umap' %in% names(vals$data@reductions)) 'UMAP' else if ('tsne' %in% names(vals$data@reductions)) 'tSNE' else NULL

        do.call(tabsetPanel, tabs)
    })

    # ======================================================================= #
    # Barplot
    # ======================================================================= #

    output$barplot <- renderPlot({
        req(input$group.by, input$group.highlight)

        barplot_vh(
            vals$data,
            group.by = input$group.by,
            ident.1 = input$group.highlight,
            chain = input$scatterplot.chain
        )
    })

    # ======================================================================= #
    # Lineplot CDR3-length
    # ======================================================================= #

    output$cdr3.length <- renderPlot({
        req(input$group.by, input$group.highlight, input$sequence.type)

        SpectratypePlot(
            vals$data,
            group.by = input$group.by,
            subset = input$group.highlight,
            sequence.type = input$sequence.type,
            plot.type = "line"
        )
    })

    # ======================================================================= #
    # Circosplot
    # ======================================================================= #

    output$circosplot <- renderPlot({
        req(vals$data)

        circosplot(vals$data)
    })

    output$circosplot.subset <- renderPlot({
        req(vals$data, input$group.by, input$group.highlight)

        circosplot(vals$data, group.by = input$group.by, subset = input$group.highlight)
    })

    # ======================================================================= #
    # Frequency CDR3 AA sequences
    # ======================================================================= #

    output$cdr3.frequency <- renderPlot({
        req(vals$data)

        CDR3freq(vals$data, chain = NULL)
    })

    output$cdr3.frequency.subset <- renderPlot({
        req(vals$data, input$group.by, input$group.highlight, input$sequence.type)

        CDR3freq(
            vals$data,
            chain = NULL,
            group.by = input$group.by,
            subset = input$group.highlight,
            sequence.type = input$sequence.type
        )
    })

    # ======================================================================= #
    # Featureplot clonotype
    # ======================================================================= #

    output$featureplot.clonotype <- renderPlot({
        req(vals$data, input$featureplot.clonotype)

        FeaturePlot_vdj(vals$data, input$featureplot.clonotype) + theme(
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            plot.title = element_blank()
        )
    })

    # ======================================================================= #
    # Update ident choices on group.by change
    # ======================================================================= #

    observeEvent(input$group.by, {
        req(vals$data, input$group.by)

        updateSelectizeInput(session, "group.highlight", choices = NULL, selected = NULL)
        groups <- isolate(vals$data@meta.data[, input$group.by]) %>% as.character() %>% unique() %>% gtools::mixedsort(x = .)
        updateSelectizeInput(session, "group.highlight", choices = groups)
    }, priority = 10)

    observeEvent(input$compare.group.by, {
        req(vals$data, input$compare.group.by)

        groups <- vals$data@meta.data[, input$compare.group.by] %>% as.character() %>% unique() %>% gtools::mixedsort(x = .)
        updateSelectizeInput(session, "compare.ident.1", choices = groups)
        updateSelectizeInput(session, "compare.ident.2", choices = groups)
    })

    # ======================================================================= #
    # Barplot to compare groups
    # ======================================================================= #

    output$barplot.comparison <- renderPlot({
        req(vals$data, input$compare.ident.1, input$compare.ident.2)

        barplot_vh(vals$data, group.by = input$compare.group.by, ident.1 = input$compare.ident.1, ident.2 = input$compare.ident.2, grid = input$compare.grid, legend = input$compare.legend)
    })
}
