#' Barplot for IGHV-family per group
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
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
        tidyr::spread(.data[[group.by]], n) %>%
        replace(is.na(.), 0) %>%
        arrange(factor(.data[[chain.column]], levels = families)) %>%
        column_to_rownames(chain.column)

    plots <- list()

    for (group in colnames(data)) {
        plot.data <- data[, group] %>%
            as.data.frame() %>%
            rename(freq = .data[['.']])
        plot.data$family <- factor(rownames(data), levels = families)

        plots[[group]] <- ggplot(plot.data, aes(x = family, y = freq, fill = family)) +
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

    grid.arrange(grobs = plots, ncol = 3)
}
