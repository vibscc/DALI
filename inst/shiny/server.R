library(shinyFiles)
library(dplyr)

function(input, output, session) {

    # ======================================================================= #
    # Initialize app
    # ======================================================================= #
    volumes <- c("Home (~)" = fs::path_home(), "/" = "/")
    shinyFileChoose(input, "seurat_rds", session = session, roots = volumes)
    shinyDirChoose(input, "bcr_dir", session = session, roots = volumes)
    shinyDirChoose(input, "tcr_dir", session = session, roots = volumes)
    shinyFileChoose(input, "reference_fasta", session = session, roots = volumes)

    vals <- reactiveValues(data = .GlobalEnv$.data.object.VDJ)
    upload <- reactiveValues(
        seurat.rds = NULL,
        bcr.dir = NULL,
        tcr.dir = NULL
    )

    app.initialize <- function() {
        vals$data <- Seurat::AddMetaData(isolate(vals$data), metadata = Seurat::Idents(isolate(vals$data)), col.name = "default.clustering")

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

        updateSelectInput(session, "group.highlight", choices = groups)

        updateSelectInput(session, "compare.group.by", choices = categorical.metadata.comparable, selected = metadata.default)
        updateSelectInput(session, "clonotype.group.by", choices = categorical.metadata, selected = metadata.default)

        assays.vdj <- names(isolate(vals$data@misc$VDJ))
        updateSelectInput(session, "active.assay", choices = assays.vdj, selected = DefaultAssayVDJ(isolate(vals$data)))

        metadata.columns <- colnames(isolate(vals$data@meta.data))
        updateSelectInput(session, "group.by", choices = metadata.columns, selected = "default.clustering")

        selected <- if ("umap" %in% reductions) "umap" else if ("tsne" %in% reductions) "tsne" else NULL
        updateSelectizeInput(session, "featureplot.reduction", choices = reductions, selected = selected)

        assays <- names(isolate(vals$data@assays))
        assays.default <- Seurat::DefaultAssay(isolate(vals$data))
        updateSelectInput(session, "transcriptomics.assay", choices = assays, selected = assays.default)
        updateSelectInput(session, "transcriptomics.reduction", choices = reductions, selected = selected)

        clonotypes <- unique(isolate(vals$data@meta.data$clonotype)) %>% gtools::mixedsort(x = .)
        updateSelectizeInput(session, "transcriptomics.clonotype", choices = clonotypes, server = T)

        updateSelectizeInput(session, "deg.group.by", choices = categorical.metadata, selected = metadata.default, server = T)
        updateSelectInput(session, "deg.assay", choices = assays, selected = assays.default)
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
                        group.by = "default.clustering",
                        cols = Diversity:::GetCategoricalColorPalette(object@meta.data$default.clustering)
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
                        cols = Diversity:::GetCategoricalColorPalette(object@meta.data$default.clustering)
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
                        cols = Diversity:::GetCategoricalColorPalette(object@meta.data$default.clustering)
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
                h4("-- BCR DATA (optional) --"),
                shinyDirButton("bcr_dir", "Browse...", "Select cellranger output folder containing vdj_b (BCR) data", multiple = F),
                textOutput("bcr.dir.text", inline = T)
            ),
            div(
                h4("-- TCR DATA (optional) --"),
                shinyDirButton("tcr_dir", "Browse...", "Select cellranger output folder containing vdj_t (TCR) data", multiple = F),
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
                actionButton("load", "Load")
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

    if (is.null(isolate(vals$data)) || !IsValidSeuratObject(isolate(vals$data))) {
        showModal(dataUploadModal())
    } else {
        app.initialize()
    }

    observeEvent(input$seurat_rds, {
        upload$seurat.rds <- shinyFiles::parseFilePaths(volumes, input$seurat_rds)
    })

    observeEvent(input$bcr_dir, {
        upload$bcr.dir <- shinyFiles::parseDirPath(volumes, input$bcr_dir)
    })

    observeEvent(input$tcr_dir, {
        upload$tcr.dir <- shinyFiles::parseDirPath(volumes, input$tcr_dir)
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

        cells.with.VDJ <- vals$data@meta.data %>% filter(!is.na(.data$vdj.v_gene) | !is.na(.data$vj.v_gene) ) %>% nrow()

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

    output$chain.usage.heatmap <- renderPlot({
        req(vals$data, input$chain.usage.chain, input$chain.usage.region)

        HeatmapChainRegion(
            vals$data,
            chain = input$chain.usage.chain,
            region = input$chain.usage.region,
            add.missing.families = input$chain.usage.add.missing.families,
            show.missing.values = F
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
            h_seqs <- c(h_seqs, Diversity:::ClonotypeToSequence(vals$data, clonotype, "VDJ"))
            l_seqs <- c(l_seqs, Diversity:::ClonotypeToSequence(vals$data, clonotype, "VJ"))
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

    # VDJ assay change

    observeEvent(input$active.assay, {
        DefaultAssayVDJ(vals$data) <- input$active.assay
        app.initialize()
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
        if (!Diversity::DefaultAssayVDJ(vals$data) == "BCR") {
            return()
        }

        div(
            strong("Select VDJ reference fasta:"),
            shinyFilesButton("reference_fasta", "Browse...", "Choose the VDJ reference fasta for loaded dataset",  multiple = F, filetype = list(data = c("fa", "fasta", "fas"))),
            textOutput("reference.fasta.path.text", inline = T)
        )
    })

    output$clonotype.lineage <- renderPlot({
        req(vals$reference, vals$data, vals$clonotype.table.selected)

        if (length(vals$reference) == 0 || is.null(vals$reference)) {
            return()
        }

        if (Diversity::DefaultAssayVDJ(vals$data) == "BCR") {
            Diversity::LineageTree(vals$data, vals$clonotype.table.selected, vals$reference$datapath)
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
}
