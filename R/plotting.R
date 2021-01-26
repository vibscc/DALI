#' Barplot for IGHV-family per group
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param chain Chain to plot. Can be any of c('l', 'light', 'h', 'heavy')
#'
#' @importFrom dplyr %>% arrange count rename
#' @importFrom ggplot2 ggplot geom_col labs aes scale_fill_manual geom_label theme element_rect element_blank element_text unit
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr spread
#'
#' @export

barplot_vh <- function(object, group.by = NULL, chain = "h") {
    chain <- tolower(chain)
    if (chain == 'heavy' || chain == 'h') {
        chain.column <- 'h.V.fam'
    } else if (chain == 'light' || chain == 'l') {
        chain.column <- 'l.V.fam'
    } else {
        chain.column <- 'l.c_gene'
    }

    if (is.null(group.by)) {
        group.by <- 'seurat_clusters'
    }

    if (!group.by %in% colnames(object@meta.data)) {
        stop("Invalid group.by column ", group.by)
    }


    families <- object@meta.data[, chain.column] %>% na.omit() %>% unique()
    families <- families %>% gtools::mixedsort(decreasing = sum(grepl('-', .)) > 0)

    data <- object@meta.data %>%
        count(.data[[chain.column]], .data[[group.by]]) %>%
        na.omit() %>%
        spread(.data[[group.by]], n) %>%
        replace(is.na(.), 0) %>%
        arrange(factor(.data[[chain.column]], levels = families)) %>%
        column_to_rownames(chain.column)

    plots <- list()

    for (group in colnames(data)) {
        plot.data <- data[, group] %>%
            as.data.frame() %>%
            rename(freq = .data[["."]])
        plot.data$fam <- factor(rownames(data), levels = families)

        plots[[group]] <- ggplot(plot.data, aes(x = fam, y = freq, fill = fam)) +
            geom_col() +
            labs(y = "Cell number", x = "Family", title = group) +
            scale_fill_manual(values = colorRampPalette(c("darkblue", "lightblue"))(nrow(data))) +
            geom_label(data = NULL, aes(label = freq), fill = "white", size = 2) +
            theme(
                panel.background = element_rect(fill = "white"), # bg of the panel
                plot.background = element_rect(fill = "white"), # bg of the plot
                panel.grid.major = element_blank(), # get rid of major grid
                panel.grid.minor = element_blank(), # get rid of minor grid
                legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
                legend.position="none",
                legend.key.size = unit(0.1, "cm"),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            )

    }

    gridExtra::grid.arrange(grobs = plots, ncol = 3)
}

#' Circosplot for family to gene distribution
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param subset Subset data to these groups
#'
#' @importFrom circlize chordDiagram circos.track circos.text CELL_META
#' @importFrom dplyr %>% select
#' @importFrom graphics strwidth
#' @importFrom stats na.omit
#'
#' @export
circosplot <- function(object, group.by = NULL, subset = NULL) {
    if (is.null(group.by)) {
        group.by <- "seurat_clusters"
    }

    if (!group.by %in% colnames(object@meta.data)) {
        stop("Invalid group.by column ", group.by)
    }

    if (!is.null(subset)) {
        cells <- rownames(object@meta.data)[object@meta.data[[group.by]] %in% subset]
        object <- subset(object, cells = cells)
    }

    plot.data <- object@meta.data %>%
                    select(h.V.fam, l.v_gene) %>%
                    na.omit() %>%
                    table() %>%
                    as.matrix()

    plot.data <- plot.data[gtools::mixedsort(rownames(plot.data), decreasing = T), ]

    chordDiagram(plot.data, annotationTrack = "grid",
                 preAllocateTracks = list(track.height = plot.data %>% dimnames() %>% unlist() %>% strwidth() %>% max()/3))
    circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
    }, bg.border = NA)
}

#' Lineplot with the CDR3 lenght per group
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param subset Subset data to these groups
#'
#' @importFrom dplyr %>% mutate select full_join count
#' @importFrom ggplot2 ggplot aes geom_line labs theme element_rect element_line element_blank unit element_text
#' @importFrom stats na.omit
#'
#' @export

cdr3length <- function(object, group.by = NULL, subset = NULL) {
    if (is.null(group.by)) {
        group.by <- "seurat_clusters"
    }

    if (!group.by %in% colnames(object@meta.data)) {
        stop("Invalid group.by column ", group.by)
    }

    if (!is.null(subset)) {
        cells <- rownames(object@meta.data)[object@meta.data[[group.by]] %in% subset]
        object <- subset(object, cells = cells)
    }

    plots <- list()
    groups <- unique(object@meta.data[[group.by]]) %>% gtools::mixedsort()

    for(group in groups) {
        cells <- rownames(object@meta.data)[object@meta.data[[group.by]] == group]
        subset <- subset(object, cells = cells)

        plot.data.h <- subset@meta.data %>%
            mutate(len = nchar(h.cdr3)) %>%
            count(len) %>%
            na.omit() %>%
            mutate(freq = n/sum(n) * 100) %>%
            select(len, freq)
        plot.data.l <- subset@meta.data %>%
            mutate(len = nchar(l.cdr3)) %>%
            count(len) %>%
            na.omit() %>%
            mutate(freq = n/sum(n) * 100) %>%
            select(len, freq)

        plot.data <- full_join(plot.data.h, plot.data.l, by = "len") %>% replace(is.na(.), 0)
        colnames(plot.data) <- c("cdr3.length", "heavy.chain", "light.chain")

        plots[[group]] <- ggplot(plot.data, aes(x = cdr3.length)) +
            geom_line(aes(y = heavy.chain), color = "black") +
            geom_line(aes(y = light.chain), color = "red") +
            labs(x = "CDR3 length (AA)", y = "Frequency of cells", title = paste0("CDR3 length - ", group)) +
            theme(
                panel.background = element_rect(fill = "white"), # bg of the panel
                plot.background = element_rect(fill = "white"), # bg of the plot
                axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(), # get rid of major grid
                panel.grid.minor = element_blank(), # get rid of minor grid
                legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
                legend.key.size = unit(0.1, "cm"),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            )
    }

    gridExtra::grid.arrange(grobs = plots, ncol = 3)
}
