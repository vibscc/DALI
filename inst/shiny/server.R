function(input, output, session) {

    # ======================================================================= #
    # Initialize app
    # ======================================================================= #

    vals <- reactiveValues(data = .GlobalEnv$.data.object.VDJ)

    app.initialize <- function() {
        updateSelectInput(session, "group.highlight", choices = levels(isolate(vals$data@meta.data$seurat_clusters)))
        renderReductionPlots(isolate(vals$data))
    }

    # ======================================================================= #
    # Function defenitions
    # ======================================================================= #

    # Render DimPlots for each reduction
    # TODO: replace with reduction plot, colored by VDJ data (v_gene, c_gene...)
    renderReductionPlots <- function(object) {
        for (reduction in names(object@reductions)) {
            local({
                plotname <- paste0('reduction.plot.', reduction)
                r <- reduction
                output[[plotname]] <- renderPlotly({
                    ggplotly(Seurat::DimPlot(object, reduction = r)) %>%
                        onRender("
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

        cells.with.VDJ <- vals$data@meta.data %>% filter(!is.na(.data$h.v_gene) | !is.na(.data$l.v_gene)) %>% nrow()

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
            tabPanel(Diversity:::formatDimred(reduction), plotlyOutput(plotname))
        })
        # TODO: select default tab in a more elegant way. PCA should be avoided as default tab, since this is the least informative
        tabs[['selected']] <- if('umap' %in% names(vals$data@reductions)) 'UMAP' else if('tsne' %in% names(vals$data@reductions)) 'tSNE' else NULL
        do.call(tabsetPanel, tabs)
    })

    # ======================================================================= #
    # Barplot
    # ======================================================================= #

    output$barplot <- renderPlot({
        req(input$group.highlight)
        barplot_vh(vals$data, groups.to.plot = input$group.highlight)
    })

    # ======================================================================= #
    # Lineplot CDR3-length
    # ======================================================================= #

    output$lineplot <- renderPlot({
        req(input$group.highlight)

        cdr3length(vals$data, subset = input$group.highlight)
    })
}
