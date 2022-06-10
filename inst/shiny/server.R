library(shinyFiles)
library(dplyr)
suppressPackageStartupMessages(library(circlize))

function(input, output, session) {

    # ======================================================================= #
    # Initialize app
    # ======================================================================= #
    volumes <- c("Home (~)" = fs::path_home(), "/" = "/")
    shinyFileChoose(input, "seurat_rds", session = session, roots = volumes)
    shinyDirChoose(input, "bcr_dir", session = session, roots = volumes)
    shinyDirChoose(input, "tcr_dir", session = session, roots = volumes)
    shinyFileChoose(input, "reference_fasta", session = session, roots = volumes)

    vals <- reactiveValues(
        data = .GlobalEnv$.data.object.VDJ,
        lineage.tab = NULL,
        has_vdj = FALSE,
        loaded_data = !is.null(.GlobalEnv$.data.object.VDJ)
    )

    upload <- reactiveValues(
        seurat.rds = NULL,
        bcr.dir = NULL,
        tcr.dir = NULL
    )

    app.initialize <- function() {
        vals$data <- Seurat::AddMetaData(isolate(vals$data), metadata = Seurat::Idents(isolate(vals$data)), col.name = "default.clustering")
        vals$has_vdj <- IsValidSeuratObject(isolate(vals$data))

        renderReductionPlotsChainUsage(isolate(vals$data))
        renderReductionPlotsExpansion(isolate(vals$data))
        renderReductionPlotsComparison(isolate(vals$data))

        groups <- levels(isolate(vals$data@meta.data$default.clustering))
        reductions <- names(isolate(vals$data@reductions))

        categorical.metadata <- c()
        categorical.metadata.comparable <- c()
        rows <- nrow(isolate(vals$data@meta.data))
        for (column in colnames(isolate(vals$data@meta.data))) {
            coldata <- isolate(vals$data@meta.data[[column]])
            unique.count <- coldata %>% unique() %>% length()

            if ((is.factor(coldata) || is.character(coldata) || is.logical(coldata)) & unique.count < 0.75 * rows) {
                categorical.metadata <- c(categorical.metadata, column)
                if (unique.count > 1) {
                    categorical.metadata.comparable <- c(categorical.metadata.comparable, column)
                }
            }
        }
        categorical.metadata <- gtools::mixedsort(categorical.metadata)
        metadata.default <- if ("default.clustering" %in% categorical.metadata) "default.clustering" else NULL

        assays <- names(isolate(vals$data@assays))
        assays.default <- Seurat::DefaultAssay(isolate(vals$data))
        assays.vdj <- names(isolate(vals$data@misc$VDJ))
        assay.selection <- DefaultAssayVDJ(isolate(vals$data))

        updateSelectInput(session, "group.highlight", choices = groups)

        updateSelectInput(session, "compare.group.by", choices = categorical.metadata.comparable, selected = metadata.default)
        categorical.metadata.no.clonotype <- categorical.metadata[categorical.metadata != "clonotype"]
        updateSelectInput(session, "clonotype.group.by", choices = categorical.metadata.no.clonotype, selected = metadata.default)

        updateSelectInput(session, "active.assay", choices = assays.vdj, selected = assay.selection)

        metadata.columns <- colnames(isolate(vals$data@meta.data))
        updateSelectInput(session, "group.by", choices = metadata.columns, selected = "default.clustering")

        updateSelectInput(session, "subsetby.circos.genes",choices = categorical.metadata, selected = metadata.default )
        updateSelectInput(session, "subsetby.circos.chains",choices = categorical.metadata, selected = metadata.default)

        selected <- if ("umap" %in% reductions) "umap" else if ("tsne" %in% reductions) "tsne" else NULL
        updateSelectizeInput(session, "featureplot.reduction", choices = reductions, selected = selected)

        updateSelectInput(session, "transcriptomics.assay", choices = assays, selected = assays.default)
        updateSelectInput(session, "transcriptomics.reduction", choices = reductions, selected = selected)

        clonotypes <- unique(isolate(vals$data@meta.data$clonotype)) %>% gtools::mixedsort(x = .)
        updateSelectizeInput(session, "transcriptomics.clonotype", choices = clonotypes, server = T)

        updateSelectizeInput(session, "deg.group.by", choices = categorical.metadata, selected = metadata.default, server = T)
        updateSelectInput(session, "deg.assay", choices = assays, selected = assays.default)

        updateSelectInput(session, "transcriptomics.assay.novdj", choices = assays, selected = assays.default)
        updateSelectInput(session, "transcriptomics.reduction.novdj", choices = reductions, selected = selected)

        updateSelectizeInput(session, "deg.group.by.novdj", choices = categorical.metadata, selected = metadata.default, server = T)
        updateSelectInput(session, "deg.assay.novdj", choices = assays, selected = assays.default)

        vals$categorical.metadata <- categorical.metadata
    }

    output$headerUI <- renderUI({
        if (!vals$has_vdj) {
            tags$div(class = "form-group row col-sm-12",
                     tags$div(class = "col-sm-9"),
                     div(class = "col-sm-3",
                         actionButton("more.files", label = "File Management", class = "files")
                     )
                )
        } else {
            tags$div(class = "form-group row col-sm-12",
              tags$label("Assay", class = "col-sm-3 text-right col-form-label"),
              tags$div(class = "col-sm-6",
                       tags$select(name = "active.assay", id = "active.assay", class = "form-control rounded-all-90"),
              ),
              div(class = "col-sm-3",
                  actionButton("more.files", label = "File Management", class = "files")
              )
            )
        }
    })

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
                        group.by = "default.clustering",
                        cols = DALI:::GetCategoricalColorPalette(object@meta.data$default.clustering)
                    ) + theme(
                        axis.line = element_blank(),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text = element_blank()
                    ) + ggtitle("Clustering")
                })


                    output[[paste0('reduction.plot.vdj.', r)]] <- renderPlot({
                        if (!vals$has_vdj) {
                            return()
                        }

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
                            axis.text = element_blank(),
                            panel.grid.major = element_blank()
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
                        group.by = "default.clustering",
                        cols = DALI:::GetCategoricalColorPalette(object@meta.data$default.clustering)
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
                    ) + ggtitle("Clonal expansion")
                })
                output[[paste0('graph.', r)]] <- renderPlot({
                    CloneConnGraph(
                        object = object,
                        reduction = r
                    ) + theme(
                        legend.position = "none",
                        plot.title = element_text(hjust = 0.5, face = "bold", vjust = 1, size = 16, margin = margin(0,0,7,7))
                    ) + ggtitle("Clonal connection plot")

                })
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
                        reduction = r,
                        cols = DALI:::GetCategoricalColorPalette(object@meta.data$default.clustering)
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
                tags$br(),
                shinyFilesButton("seurat_rds", "Browse...", "Choose an Rds file to load",  multiple = F, filetype = list(data = c("Rds", "rds"))),
                actionButton("seurat_rds_remove", "Remove", class = "upload-remove"),
                textOutput("seurat.rds.path.text", inline = T),
            ),
            div(
                h4("-- BCR DATA (optional) --"),
                shinyDirButton("bcr_dir", "Browse...", "Select cellranger output folder containing vdj_b (BCR) data", multiple = F),
                actionButton("bcr_dir_remove", "Remove", class = "upload-remove"),
                textOutput("bcr.dir.text", inline = T)
            ),
            div(
                h4("-- TCR DATA (optional) --"),
                shinyDirButton("tcr_dir", "Browse...", "Select cellranger output folder containing vdj_t (TCR) data", multiple = F),
                actionButton("tcr_dir_remove", "Remove", class = "upload-remove"),
                textOutput("tcr.dir.text", inline = T)
            ),
            if (!is.null(error)) {
                div(
                    h4("Error:", class = "text-danger"),
                    p(error)
                )
            },
            title = "Load data",
            footer = tagList(
                actionButton("load", "Load"),
                actionButton("close", "Close", class = "closer")
            ),
            easyClose = F,
            size = "l"
        )
    }

    loadingModal <- function() {
        modalDialog(
            "Data is loading. This can take a while for larger datasets!",
            title = "Loading...",
            footer = NULL,
            easyClose = F
        )
    }

    if (is.null(isolate(vals$data))) {
        showModal(dataUploadModal())
    } else {
        app.initialize()
    }

    observeEvent(input$seurat_rds_remove, {
        upload$seurat.rds <- NULL
        vals$loaded_data <- FALSE
    })

    observeEvent(input$tcr_dir_remove, {
        upload$tcr.dir <- NULL
    })

    observeEvent(input$bcr_dir_remove, {
        upload$bcr.dir <- NULL
    })

    observeEvent(input$seurat_rds, {
        upload$seurat.rds <- shinyFiles::parseFilePaths(volumes, input$seurat_rds)
    })

    observeEvent(input$bcr_dir, {
        upload$bcr.dir <- shinyFiles::parseDirPath(volumes, input$bcr_dir)
    })

    observeEvent(input$tcr_dir, {
        upload$tcr.dir <- shinyFiles::parseDirPath(volumes, input$tcr_dir)
    })

    observeEvent(input$close, {
        if (vals$loaded_data == TRUE) {
            removeModal()
        } else {
            showModal(dataUploadModal(error = "Load the data first!"))
            return()
        }
    })

    observeEvent(input$load, {
        has_new_seurat_rds <- length(upload$seurat.rds) > 0 && !is.null(upload$seurat.rds)

        if (is.null(vals$data) && !has_new_seurat_rds) {
            showModal(dataUploadModal(error = "Missing Seurat Rds file!"))
            return()
        }

        if (has_new_seurat_rds && !file.exists(upload$seurat.rds$datapath)) {
            showModal(dataUploadModal(error = paste0("Could not find file '", upload$seurat.rds$datapath, "'")))
        }


        removeModal()
        showModal(loadingModal())
        if (has_new_seurat_rds) {
            data <- readRDS(upload$seurat.rds$datapath)
        } else {
            data <- vals$data
        }

        if (length(upload$bcr.dir) > 0) {
            data <- tryCatch({
                Read10X_vdj(data, upload$bcr.dir, assay = "BCR")
            }, error = function(e) {
                showModal(dataUploadModal(error = paste0(e, " Make sure the selected VDJ data matches the Seurat object + is of the type vdj_b (BCR)")))
                return(NULL)
            })
            if (!isS4(data) & is.null(data)) {
                return()
            }
        }

        if (length(upload$tcr.dir) > 0) {
            data <- tryCatch({
                Read10X_vdj(data, upload$tcr.dir, assay = "TCR")
            }, error = function(e) {
                showModal(dataUploadModal(error = paste0(e, " Make sure the selected VDJ data matches the Seurat object + is of the type vdj_t (TCR)")))
                return(NULL)
            })
            if (!isS4(data) & is.null(data)) {
                return()
            }
        }

        vals$data <- data
        removeModal()
        vals$loaded_data <- TRUE
        app.initialize()
        return()
    })

    output$seurat.rds.path.text <- renderText({
        req(upload$seurat.rds)

        gsub("/+", "/", upload$seurat.rds$datapath)
    })

    output$bcr.dir.text <- renderText({
        req(upload$bcr.dir)

        gsub("/+", "/", upload$bcr.dir)
    })

    output$tcr.dir.text <- renderText({
        req(upload$tcr.dir)

        gsub("/+", "/", upload$tcr.dir)
    })

    # ======================================================================= #
    # Reduction plots UI tabs
    # ======================================================================= #
    output$dataset.metrics <- renderUI({
        req(vals$data)


        if (!vals$has_vdj) {
            cells.with.VDJ <- "Not applicable"
            .data$vdj.v_gene <- NULL
        } else {
            cells.with.VDJ <- vals$data@meta.data %>% filter(!is.na(.data$vdj.v_gene) | !is.na(.data$vj.v_gene) ) %>% nrow()
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
            if (vals$has_vdj) {
                plotname.dimred.vdj <- paste0('reduction.plot.vdj.', reduction)
            } else {
                plotname.dimred.vdj <- NULL
            }

            tabPanel(
                DALI:::FormatDimred(reduction),
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
                DALI:::FormatDimred(reduction),
                fluidRow(
                    column(4, plotOutput(plotname.dimred) %>% withSpinner()),
                    column(4, plotOutput(plotname.dimred.expansion) %>% withSpinner()),
                    column(4, plotOutput(plotname.graph) %>% withSpinner())
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
                DALI:::FormatDimred(reduction),
                plotOutput(plotname) %>% withSpinner()
            )
        })

        tabs[['selected']] <- if ('umap' %in% names(vals$data@reductions)) 'UMAP' else if ('tsne' %in% names(vals$data@reductions)) 'tSNE' else NULL
        do.call(tabsetPanel, tabs)
    })


    output$dim.reduction.tabs <- renderUI({
        req(vals$data)

        tabs <- lapply(names(vals$data@reductions), function(reduction) {
            plotname <- paste0('dimred.', reduction)

            tabPanel(
                DALI:::FormatDimred(reduction),
                plotOutput(plotname) %>% withSpinner()
            )
        })

        tabs[['selected']] <- if ('umap' %in% names(vals$data@reductions)) 'UMAP' else if ('tsne' %in% names(vals$data@reductions)) 'tSNE' else NULL
        do.call(tabsetPanel, tabs)
    })

    # ======================================================================= #
    # Chain usage
    # ======================================================================= #


    output$chain.usage.heatmap <- renderPlot({
        req(vals$data,vals$data@meta.data$vdj.v_fam, input$chain.usage.chain, input$chain.usage.region)
        HeatmapChainRegion(
            vals$data,
            color = input$chain.usage.color,
            chain = input$chain.usage.chain,
            region = input$chain.usage.region,
            add.missing.families = input$chain.usage.add.missing.families,
            show.missing.values = F,
            cluster.cols = input$chain.usage.cluster.cols
        )
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

        if (!input$clonotype.group %in% vals$data@meta.data[, input$clonotype.group.by]) {
            return()
        }

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

        top.clonotypes <- DALI:::CalculateFrequency(vals$data, 'clonotype', input$clonotype.group.by, F) %>%
            filter(.data[[input$clonotype.group.by]] == input$clonotype.group) %>%
            arrange(desc(freq)) %>%
            select(c(clonotype, freq)) %>%
            head(n = 25) %>%
            mutate(perc = round((freq/n.cells) * 100, digits = 1))

        vals$top.clonotypes <- top.clonotypes

        h_seqs <- c()
        l_seqs <- c()
        for (clonotype in top.clonotypes$clonotype) {
            h_seqs <- c(h_seqs, DALI:::ClonotypeToSequence(vals$data, clonotype, "VDJ"))
            l_seqs <- c(l_seqs, DALI:::ClonotypeToSequence(vals$data, clonotype, "VJ"))
        }

        top.clonotypes$h_seq <- h_seqs
        top.clonotypes$l_seq <- l_seqs

        colnames(top.clonotypes) <- c("Clonotype", "Cells", "pct.group", "VDJ CDR3 AA seq", "VJ CDR3 AA seq")
        top.clonotypes
    })


    # ======================================================================= #
    # Featureplot clonotype
    # ======================================================================= #

    output$featureplot.clonotype <- renderPlot({
        req(vals$data, input$featureplot.clonotype, input$featureplot.reduction)

        FeaturePlotClonotype(vals$data, input$featureplot.reduction, input$featureplot.clonotype, size = 0.8, missing.alpha = 0.4) + theme(
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

        updateSelectizeInput(session, "compare.ident.1", choices = groups, selected = groups[[1]], server = T)
        updateSelectizeInput(session, "compare.ident.2", choices = groups, selected = groups[[2]], server = T)
    })

    observeEvent(input$clonotype.group.by, {
        req(vals$data, input$clonotype.group.by)

        groups <- vals$data@meta.data[, input$clonotype.group.by] %>% as.character() %>% unique() %>% gtools::mixedsort(x = .)
        updateSelectizeInput(session, "clonotype.group", choices = groups, selected = groups[[1]], server = T)
    })

    observeEvent(input$subsetby.circos.genes, {
        req(vals$data, input$subsetby.circos.genes)

        gene.subsets <- vals$data@meta.data[, input$subsetby.circos.genes] %>% as.character() %>% unique() %>% gtools::mixedsort(x = .)
        gene.subsets <- c("Select group...", gene.subsets)
        updateSelectizeInput(session, "gene.subset.group", choices = gene.subsets, selected = "Select group...", server = T)
    })

    observeEvent(input$subsetby.circos.chains, {
        req(vals$data, input$subsetby.circos.chains)

        chain.subsets <- vals$data@meta.data[, input$subsetby.circos.chains] %>% as.character() %>% unique() %>% gtools::mixedsort(x = .)
        chain.subsets <- c("Select group...", chain.subsets)
        updateSelectizeInput(session, "chain.subset.group", choices = chain.subsets, selected = "Select group...", server = T)
    })

    # Top clonotypes change

    observeEvent(vals$top.clonotypes, {
        updateSelectInput(session, "featureplot.clonotype", choices = vals$top.clonotypes$clonotype)
    })

    # VDJ assay change

    observeEvent(input$active.assay, {
        req(vals$data)
        DefaultAssayVDJ(vals$data) <- input$active.assay
        app.initialize()
    })

    # Remove UI elements (when no optional files present)
    # add them again when we upload a new file
    observe({
        if (vals$has_vdj) {
            showTab("VDJ", "General view")
            showTab("VDJ", "Population comparison")
            showTab("VDJ", "Clonotypes")
            showTab("VDJ", "Transcriptomics")
            showTab("VDJ", "DEG")
            showTab("VDJ", "Clone view")
            hideTab("Seurat", "Clustering & Transcriptomics")
            hideTab("Seurat", "DEG Selector")
            hideTab("Seurat", "DEG Results")
        } else {
            hideTab("VDJ", "General view")
            hideTab("VDJ", "Population comparison")
            hideTab("VDJ", "Clonotypes")
            hideTab("VDJ", "Transcriptomics")
            hideTab("VDJ", "DEG")
            hideTab("VDJ", "Clone view")
            showTab("Seurat", "Clustering & Transcriptomics")
            showTab("Seurat", "DEG Selector")
            showTab("Seurat", "DEG Results")
        }
    })

    # upload new files

    observeEvent(input$more.files, {
        showModal(dataUploadModal())
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
            region = input$compare.region,
            chain = input$compare.chain,
            legend = F
        )
    })


    # ======================================================================= #
    # Clonotypes table
    # ======================================================================= #

    output$clonotypes.table <- DT::renderDT({
        req(vals$data)

        data <- vals$data@meta.data %>% filter(!is.na(clonotype))

        cols.for.nt.sequence <- c("fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt")
        cols.for.aa.sequence <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")

        vdj.prefix <- "vdj"
        vj.prefix <- "vj"

        vdj.seq.nt.cols <- paste0(vdj.prefix, ".", cols.for.nt.sequence)
        vdj.seq.aa.cols <- paste0(vdj.prefix, ".", cols.for.aa.sequence)
        vj.seq.nt.cols <- paste0(vj.prefix, ".", cols.for.nt.sequence)
        vj.seq.aa.cols <- paste0(vj.prefix, ".", cols.for.aa.sequence)

        vars <- list("vdj.seq.nt" = vdj.seq.nt.cols, "vdj.seq.aa" = vdj.seq.aa.cols, "vj.seq.nt" = vj.seq.nt.cols, "vj.seq.aa" = vj.seq.aa.cols)
        i <- 1
        for (cols in vars) {
            for (col in cols) {
                colname <- names(vars)[i]
                style.classname <- (strsplit((strsplit(col, "\\.") %>% unlist())[2], "_") %>% unlist())[1]

                if (!colname %in% colnames(data)) {
                    data[, colname] <- ""
                }
                if (!col %in% colnames(data)) {
                    data[, col] <- "-"
                }

                data[is.na(data[,col]), col] <- "-"

                data[,colname] <- paste0(data[,colname], paste0("<span class='", style.classname, "' title='", style.classname, "'>", data[,col], "</span>"))
            }
            i <- i + 1
        }

        vals$clonotype.table <- data %>%
            dplyr::group_by(clonotype) %>%
            summarize(
                n.cells = n(),
                vdj.seq.nt = vdj.seq.nt %>% unique() %>% paste(collapse = "<br>"),
                vdj.seq.aa = vdj.seq.aa %>% unique() %>% paste(collapse = "<br>"),
                vj.seq.nt = vj.seq.nt %>% unique() %>% paste(collapse = "<br>"),
                vj.seq.aa = vj.seq.aa %>% unique() %>% paste(collapse = "<br>"),
            )

        DT::datatable(
            vals$clonotype.table,
            escape = F,
            rownames = F,
            options = list(scrollX = T),
            selection = "single"
        ) %>% DT::formatStyle(names(vars), `font-family` = "monospace")
    })

    observeEvent(input$clonotypes.table_rows_selected, {
        req(vals$clonotype.table)

        vals$clonotype.table.selected <- vals$clonotype.table[input$clonotypes.table_rows_selected, "clonotype"] %>% pull("clonotype")
    })

    observeEvent(input$reference_fasta, {
        vals$reference <- shinyFiles::parseFilePaths(volumes, input$reference_fasta)
    })

    output$reference.fasta.path.text <- renderText({
        req(vals$reference)

        gsub("/+", "/", vals$reference$datapath)
    })

    output$clonotype.lineage.ui <- renderUI({
        if (!DALI::DefaultAssayVDJ(vals$data) == "BCR") {
            return()
        }

        assays <- names(isolate(vals$data@assays))

        div(class = "well",
            h3("B-cell lineage"),
            strong("Select VDJ reference fasta:"),
            shinyFilesButton("reference_fasta", "Browse...", "Choose the VDJ reference fasta for loaded dataset",  multiple = F, filetype = list(data = c("fa", "fasta", "fas"))),
            textOutput("reference.fasta.path.text", inline = T),
            checkboxInput("lineage.color.tips", label = "Color tips by metadata/feature", value = F),
            tabsetPanel(id = "lineage_tabs",
                tabPanel("Metadata",
                    selectizeInput("lineage.metadata", label = "Metadata", choices = vals$categorical.metadata, multiple = F)
                ),
                tabPanel("Feature",
                    selectInput("lineage.assay", label = "Assay", choices = assays, multiple = F),
                    selectizeInput("lineage.feature", label = "Feature", choices = NULL, multiple = F)
                )
            )
        )
    })

    output$clonotype.lineage <- renderPlot({
        req(vals$reference, vals$data, vals$clonotype.table.selected, vals$lineage.tab)

        if (length(vals$reference$datapath) == 0 | is.null(vals$reference)) {
            return()
        }

        if (DALI::DefaultAssayVDJ(vals$data) == "BCR") {
            color.tip.by <- NULL

            if (input$lineage.color.tips) {
                if (vals$lineage.tab == "metadata") {
                    color.tip.by <- input$lineage.metadata
                } else if (vals$lineage.tab == "feature" & nchar(input$lineage.feature) > 0) {
                    color.tip.by <- input$lineage.feature
                }
            }

            DALI::LineageTree(
                object = vals$data,
                clonotype = vals$clonotype.table.selected,
                reference = vals$reference$datapath,
                color.tip.by = color.tip.by
            )
        }
    })

    observeEvent(input$lineage.assay, {
        Seurat::DefaultAssay(vals$data) <- input$lineage.assay

        features <- rownames(vals$data) %>% sort()
        updateSelectizeInput(session, "lineage.feature", choices = features, server = T)
    })

    observe({
        if (!is.null(input$lineage_tabs)) {
            if (input$lineage_tabs == "Metadata") {
                vals$lineage.tab <- "metadata"
            } else if (input$lineage_tabs == "Feature") {
                vals$lineage.tab <- "feature"
            } else {
                vals$lineage.tab <- NULL
            }
        }
    })


    # ======================================================================= #
    # Transcriptomics
    # ======================================================================= #

    observeEvent(input$transcriptomics.assay, {
        req(vals$data, input$transcriptomics.assay)

        Seurat::DefaultAssay(vals$data) <- input$transcriptomics.assay

        features <- rownames(vals$data) %>% sort()
        updateSelectizeInput(session, "transcriptomics.feature", choices = features, selected = features[1], server = T)
    })

    output$transcriptomics.featureplot <- renderPlot({
        req(vals$data, input$transcriptomics.feature, input$transcriptomics.reduction, input$transcriptomics.assay)

        Seurat::FeaturePlot(vals$data, input$transcriptomics.feature, reduction = input$transcriptomics.reduction) + theme(
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
        ) + ggtitle(paste0("Featureplot (", input$transcriptomics.feature ,")"))
    })

    output$transcriptomics.clonotype.featureplot <- renderPlot({
        req(vals$data, input$transcriptomics.clonotype, input$transcriptomics.reduction)

        FeaturePlotClonotype(vals$data, input$transcriptomics.reduction, input$transcriptomics.clonotype) + theme(
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", vjust = 1, size = 16, margin = margin(0,0,7,7))
        ) + ggtitle(paste0("Featureplot clonotype (", input$transcriptomics.clonotype, ")"))
    })

    # ======================================================================= #
    # DEG
    # ======================================================================= #

    observeEvent(input$deg.group.by, {
        req(vals$data)

        choices <- vals$data@meta.data[[input$deg.group.by]] %>% unique() %>% gtools::mixedsort()

        updateSelectizeInput(session, "deg.ident.1", choices = choices, selected = NULL, server = T)
        updateSelectizeInput(session, "deg.ident.2", choices = choices, selected = NULL, server = T)
    })

    observeEvent(input$deg.calculate, {
        req(vals$data, input$deg.assay)

        vals$deg.results <- NULL

        ident.1 <- input$deg.ident.1

        if (input$deg.ident.2.choice == 1) {
            ident.2 <- NULL
        } else {
            ident.2 <- input$deg.ident.2
        }

        if (is.null(input$deg.ident.1)) {
            showNotification("Group 1 can't be empty", type = c("error"), session = session)
        } else if (intersect(ident.1, ident.2) %>% length() > 0) {
            showNotification("Group 1 and 2 should not overlap!", type = c("error"), session = session)
        } else {
            withProgress(message = "Calculating DEG", detail = "This may take a while", min = 0, max = 1, value = 1, {
                vals$deg.results <- Seurat::FindMarkers(vals$data, ident.1 = ident.1, ident.2 = ident.2, group.by = input$deg.group.by, assay = input$deg.assay)
            })
        }
    })

    output$deg.output <- DT::renderDT({
        req(vals$data, vals$deg.results)

        vals$deg.results
    })

    # ======================================================================= #
    # CircosPlot
    # ======================================================================= #

    output$circosplot.genes <- renderPlot({
        req(vals$data)
        circos.clear()
        CircosPlotGenes(object = vals$data,
                        group.by = input$subsetby.circos.genes,
                        subset = input$gene.subset.group,
                        seed = 0)
    })

    output$circosplot.chains <- renderPlot({
        req(vals$data)
        circos.clear()
        CircosPlotChains(object = vals$data,
                         group.by = input$subsetby.circos.chains,
                         subset = input$chain.subset.group)
    })

    # ####################################################################### #
    # Seurat Analysis
    # ####################################################################### #

    # Transcriptomics subtab

    output$dim.reduction <- renderPlot({
        req(vals$data, input$transcriptomics.reduction.novdj)
        Seurat::DimPlot(
            vals$data,
            group.by = "default.clustering",
            reduction = input$transcriptomics.reduction.novdj,
            cols = DALI:::GetCategoricalColorPalette(vals$data@meta.data$default.clustering)
        ) + theme(
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
        ) + ggtitle("Clustering")
    })


    observeEvent(input$transcriptomics.assay.novdj, {
        req(vals$data, input$transcriptomics.assay.novdj)

        Seurat::DefaultAssay(vals$data) <- input$transcriptomics.assay.novdj

        features <- rownames(vals$data) %>% sort()
        updateSelectizeInput(session, "transcriptomics.feature.novdj", choices = features, selected = features[1], server = T)
    })

    output$transcriptomics.featureplot.novdj <- renderPlot({
        req(vals$data, input$transcriptomics.feature.novdj, input$transcriptomics.reduction.novdj, input$transcriptomics.assay.novdj)

        Seurat::FeaturePlot(vals$data, input$transcriptomics.feature.novdj, reduction = input$transcriptomics.reduction.novdj) + theme(
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
        ) + ggtitle(paste0("Featureplot (", input$transcriptomics.feature.novdj ,")"))
    })

    # DEG analyses subtab

    observeEvent(input$deg.group.by.novdj, {
        req(vals$data)

        choices <- vals$data@meta.data[[input$deg.group.by.novdj]] %>% unique() %>% gtools::mixedsort()

        updateSelectizeInput(session, "deg.ident.1.novdj", choices = choices, selected = NULL, server = T)
        updateSelectizeInput(session, "deg.ident.2.novdj", choices = choices, selected = NULL, server = T)
    })

    observeEvent(input$deg.calculate.novdj, {
        req(vals$data, input$deg.assay.novdj)

        vals$deg.results.novdj <- NULL

        ident.1.novdj <- input$deg.ident.1.novdj

        if (input$deg.ident.2.choice.novdj == 1) {
            ident.2.novdj <- NULL
        } else {
            ident.2.novdj <- input$deg.ident.2.novdj
        }

        if (is.null(input$deg.ident.1.novdj)) {
            showNotification("Group 1 can't be empty", type = c("error"), session = session)
        } else if (intersect(ident.1.novdj, ident.2.novdj) %>% length() > 0) {
            showNotification("Group 1 and 2 should not overlap!", type = c("error"), session = session)
        } else {
            withProgress(message = "Calculating DEG", detail = "This may take a while", min = 0, max = 1, value = 1, {
                vals$deg.results.novdj <- Seurat::FindMarkers(vals$data, ident.1 = ident.1.novdj, ident.2 = ident.2.novdj, group.by = input$deg.group.by.novdj, assay = input$deg.assay.novdj)
            })
        }
    })

    output$deg.output.novdj <- DT::renderDT({
        req(vals$data, vals$deg.results.novdj)

        vals$deg.results.novdj
    })
}
