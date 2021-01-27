function(input, output) {

    vals <- reactiveValues(data = .GlobalEnv$.data.object.VDJ)

    # ======================================================================= #
    # Load data if necessary
    # ======================================================================= #

    dataUploadModal <- function(failed = F) {
        modalDialog(
            fileInput("file", "Choose an rds file to load", accept = ".rds"),
            if (failed) {
                div("Invalid file!")
            },

            footer = tagList(
                actionButton("load", "Load")
            ),
            easyClose = F
        )
    }

    if (is.null(isolate(vals$data)) || !isValidSeuratObject(isolate(vals$data))) {
        showModal(dataUploadModal())
    }

    observeEvent(input$load, {
        data <- readRDS(input$file$datapath)

        if (isValidSeuratObject(data)) {
            vals$data <- readRDS(input$file$datapath)
            removeModal()
        } else {
            showModal(dataUploadModal(failed = T))
        }
    })

    # ======================================================================= #
    # Reduction plots
    # ======================================================================= #

    # Create tabsetPanel with tabPanel for each dimensionality reducion in the dataset
    output$reduction.tabs <- renderUI({
        tabs <- lapply(names(vals$data@reductions), function(reduction) {
            plotname <- paste0('reduction.plot.', reduction)
            tabPanel(formatDimred(reduction), plotlyOutput(plotname))
        })
        # TODO: select default tab in a more elegant way. PCA should be avoided as default tab, since this is the least informative
        tabs[['selected']] <- if('umap' %in% names(vals$data@reductions)) 'UMAP' else if('tsne' %in% names(vals$data@reductions)) 'tSNE' else NULL
        do.call(tabsetPanel, tabs)
    })

    # Render DimPlots for each reduction
    # TODO: replace with reduction plot, colored by VDJ data (v_gene, c_gene...)
    for (reduction in names(isolate(vals$data@reductions))) {
        local({
            plotname <- paste0('reduction.plot.', reduction)
            r <- reduction
            output[[plotname]] <- renderPlotly({
                ggplotly(Seurat::DimPlot(vals$data, reduction = r)) %>%
                onRender("
                    function(el) {
                        el.on('plotly_click', function(d) {
                            console.log('Click: ', d);
                        })
                        el.on('plotly_legendclick', function(d) {

                            console.log(d);
                            const input = document.getElementById('group.highlight');

                            // Create fake keyup-event to trigger shiny
                            const ev = document.createEvent('Event');
                            ev.initEvent('keyup');
                            ev.which = ev.keyCode = 13;

                            input.value = d.curveNumber;
                            input.dispatchEvent(ev);

                            // Prevent other handlers from firing
                            return false;
                        })
                    }
                ")
            })
        })
    }

    # ======================================================================= #
    # Barplot
    # ======================================================================= #

    output$barplot <- renderPlot(barplot_vh(vals$data, groups.to.plot = input$group.highlight))

    # ======================================================================= #
    # Lineplot CDR3-length
    # ======================================================================= #

    output$lineplot <- renderPlot(cdr3length(vals$data, subset = input$group.highlight))

    # ======================================================================= #
    # Event handling
    # ======================================================================= #

}
