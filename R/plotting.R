#' Barplot for IGHV-family per group
#'
#' @param object Seurat object
#' @param ident.1 Identy class(es) to plot
#' @param ident.2 Second identity class(es). This class will be used as comparison. If NULL, ident.1 will just be used to subset the data
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param region Region to plot. Available options: 'V'(ariable) or 'C'(onstant)
#' @param chain Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR
#' @param by.family Group genes of 1 family together. Default = TRUE
#' @param legend Should the legend be included in the plot. Default = TRUE
#' @param grid Organize plots in grid. Each plot contains information about 1 group. Default = FALSE
#'
#' @importFrom dplyr arrange case_when count rename  %>%
#' @importFrom ggplot2 aes element_blank element_rect element_text facet_grid geom_bar geom_text ggplot labs position_dodge scale_fill_manual theme unit ylim
#' @importFrom grDevices colorRampPalette
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr spread
#'
#' @export

barplot_vh <- function(object, ident.1 = NULL, ident.2 = NULL, group.by = NULL, region = c("V", "C"), chain = availableChains(object), by.family = T, legend = T, grid = F) {
    region <- match.arg(region) %>% tolower()
    chain <- match.arg(chain) %>% tolower()

    data.column <- paste0(chain, '.', region, '_')

    if (by.family && region == 'v') {
        data.column <- paste0(data.column, 'fam')
    } else {
        data.column <- paste0(data.column, 'gene')
    }

    if (is.null(group.by)) {
        group.by <- 'seurat_clusters'
    }

    if (!group.by %in% colnames(object@meta.data)) {
        stop("Invalid group.by column ", group.by)
    }

    families <- object@meta.data[, data.column] %>% na.omit() %>% unique()

    # Add missing families
    if (by.family) {
        families.completed <- c()
        prefixes <- gsub("/","",gsub("[0-9-]", "", families) %>% unique())

        for (prefix in prefixes) {
            families.with.prefix <- families[grepl(prefix, families)]

            if(data.column=="b.v_fam" | data.column=="l.v_fam" | data.column=="h.v_fam"){
                family.numbers <- c()
                for (family in families.with.prefix) {
                    family.numbers <- c(family.numbers, gsub("[A-Za-z-]", "", family)) %>% as.numeric()
                }

                families.completed <- c(families.completed, paste0(prefix, '-', seq(1,max(family.numbers))))
            }
            if(data.column=="a.v_fam"){
                family.numbers <- c()
                for (family in families.with.prefix) {
                    family.numbers <- c(family.numbers, gsub("[A-Za-z-]", "", family))
                }

                # family.numbers<-strsplit(family.numbers, "/")[[1]]

                families.completed <- c(families.completed, paste0(prefix, '-', family.numbers))
            }


        }
        families <- families.completed
    }

    families <- families %>% gtools::mixedsort(x = ., decreasing = sum(grepl('-', .)) > 0)

    data.filtered <- object@meta.data

    if (!is.null(ident.1) || !is.null(ident.2)) {
        data.filtered <- data.filtered %>%
            filter(.data[[group.by]] %in% c(ident.1, ident.2))
    }

    if (!is.null(ident.2)) {
        data.filtered[[group.by]] <- data.filtered[[group.by]] %>% as.character()

        for (ident in list(ident.1, ident.2)) {
            for (group in ident) {
                data.filtered[[group.by]] <- gsub(paste0("^", group, "$"), paste0(group.by, " (", paste(ident, collapse = ","),")"), data.filtered[[group.by]])
            }
        }
    }

    data <- data.filtered %>%
        count(.data[[data.column]], .data[[group.by]]) %>%
        group_by(.data[[group.by]]) %>%
        mutate(freq = prop.table(.data$n) * 100) %>%
        mutate(freq = round(.data$freq, 2)) %>%
        select(.data[[data.column]], .data[[group.by]], .data$freq) %>%
        na.omit() %>%
        spread(.data[[group.by]], .data$freq) %>%
        replace(is.na(.), 0)

    missing.families <- setdiff(families, data[[data.column]])

    if (length(missing.families) > 0) {
        missing.data <- matrix(0, nrow = length(missing.families), ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
        missing.data <- cbind(missing.families, missing.data)
        colnames(missing.data) <- c(data.column, colnames(data)[-1])
        data <- rbind(data, missing.data)
    }

    plot.data <- data %>%
        melt(data = ., id.vars = c(data.column), variable.name = "group", value.name = "freq") %>%
        mutate(freq = as.numeric(.data$freq))

    plot.data[[data.column]] <- factor(plot.data[[data.column]], levels = families)

    plot <- ggplot(plot.data, aes(x = .data[[data.column]], y = .data$freq, fill = .data$group)) +
        geom_bar(position = "dodge", stat = "identity") +
        ylim(0, NA) +
        labs(y = "Cell number", x = "Family", title = .data$group) +
        # geom_text(data = NULL, aes(label = .data$freq), size = 2, position = position_dodge(width = 1), vjust = -0.5) +
        theme(
            panel.background = element_rect(fill = "white"), # bg of the panel
            plot.background = element_rect(fill = "white"), # bg of the plot
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            legend.background = element_rect(fill = "transparent"), # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
            legend.key.size = unit(0.1, "cm"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        )

    if (grid) {
        plot <- plot + facet_grid(~ .data$group)
    }

    if (!legend) {
        plot <- plot + theme(
            legend.position = "none"
        )
    }

    return(plot)
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
#' @importFrom rlang .data
#' @importFrom stats na.omit
#'
#' @export
circosplot <- function(object, group.by = NULL, subset = NULL) {
    type = object@misc$default.assay.VDJ

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

    if(type=="TCR"){    Names <- c("a.v_fam","b.v_gene")  }
    if(type=="BCR"){    Names <- c("h.v_fam","l.v_gene")  }

    plot.data <- object@meta.data %>%
        select(.data[[Names[1]]], .data[[Names[2]]]) %>%
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
#' @importFrom rlang .data
#' @importFrom stats na.omit
#'
#' @export

cdr3length <- function(object, group.by = NULL, subset = NULL) {
    type = object@misc$default.assay.VDJ

    if (is.null(group.by)) {
        group.by <- "seurat_clusters"
    }

    if (!group.by %in% colnames(object@meta.data)) {
        stop("Invalid group.by column ", group.by)
    }

    if (!is.null(subset) && subset != '') {
        cells <- rownames(object@meta.data)[object@meta.data[[group.by]] %in% subset]
        object <- subset(object, cells = cells)
    }

    plots <- list()
    groups <- unique(object@meta.data[[group.by]]) %>% gtools::mixedsort(x = .)

    if(type=="TCR"){    Names <- c("a.cdr3","b.cdr3")  }
    if(type=="BCR"){    Names <- c("h.cdr3","l.cdr3")  }

    for (group in groups) {
        cells <- rownames(object@meta.data)[object@meta.data[[group.by]] == group]
        subset <- subset(object, cells = cells)

        plot.data.h <- subset@meta.data %>%
            mutate(len = nchar(.data[[Names[1]]])) %>%
            count(.data$len) %>%
            na.omit() %>%
            mutate(freq = .data$n/sum(.data$n) * 100) %>%
            select(.data$len, .data$freq)
        plot.data.l <- subset@meta.data %>%
            mutate(len = nchar(.data[[Names[2]]])) %>%
            count(.data$len) %>%
            na.omit() %>%
            mutate(freq = .data$n/sum(.data$n) * 100) %>%
            select(.data$len, .data$freq)

        plot.data <- full_join(plot.data.h, plot.data.l, by = "len") %>% replace(is.na(.), 0)
        colnames(plot.data) <- c("cdr3.length", "heavy.chain", "light.chain")

        plots[[group]] <- ggplot(plot.data, aes(x = .data$cdr3.length)) +
            geom_line(aes(y = .data$heavy.chain), color = "black") +
            geom_line(aes(y = .data$light.chain), color = "red") +
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

    gridExtra::grid.arrange(grobs = plots, ncol = min(length(plots), 3))
}




#' Dimplot for IGXV-family
#'
#' @param object Seurat object
#' @param region Region to plot. Available options: 'V'(ariable), 'D', 'J'(unction), 'C'(onstant).
#' @param chain Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR
#' @param by.family Group genes of 1 family together. Only effective with the V-gene. Default = TRUE
#' @param grid If TRUE, show per gene type in grid. If FALSE, show all genes types together on plot. Default = TRUE
#' @param ... Extra parameters passed to Seurat::Dimplot
#'
#' @importFrom dplyr %>%
#' @importFrom Seurat DimPlot
#'
#' @export

DimPlot_vh <- function(object, region = c("V", "D", "J", "C"), chain = availableChains(object), by.family = T, grid = T, ...) {

  region <- match.arg(region) %>% tolower()
  chain <- match.arg(chain) %>% tolower()

  data.column <- paste0(chain, '.', region, '_')

  if (by.family && region == 'v') {
    data.column <- paste0(data.column, 'fam')
  } else {
    data.column <- paste0(data.column, 'gene')
  }

  split <- data.column

  if (!grid) {
    split <- NULL
  }

  families <- object@meta.data[, data.column] %>% na.omit() %>% unique()
  families <- families %>% gtools::mixedsort(x = ., decreasing = sum(grepl('-', .)) > 0)
  Seurat::DimPlot(object, group.by = data.column, split.by = split, order = rev(families), ...)
}


#' Dimplot for IGXV-family
#'
#' @param object Seurat object
#' @param chain Chain to plot, available options: "L", "H", NULL (= both)
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#'
#' @importFrom ggplot2 aes_string element_blank element_line element_rect geom_bar ggplot ggtitle theme xlab ylab
#' @importFrom dplyr group_by n summarise ungroup %>%
#'
#' @export

CDR3freq <- function(object, chain = availableChains(object), group.by = NULL) {

  if (!is.null(chain)) {
    chain <- match.arg(chain) %>% tolower()
  }

  if (is.null(group.by)) {
    group.by <- "seurat_clusters"
  }

  chains <- chain

  if (is.null(chains)) {
    chains <- availableChains(object) %>% tolower()
  }

  plots <- list()
  for (chain in chains) {
    plot.title <- paste0("CDR3 ", toupper(chain), "-chain AA sequence frequency")
    data.column <- paste0(chain, ".cdr3")

    plot.data <- object@meta.data %>%
      group_by(.data[[data.column]], .data[[group.by]]) %>%
      summarise(freq = n()) %>%
      ungroup()

    plots[[data.column]] <- ggplot(plot.data, aes(x = .data[[group.by]], y = .data$freq, fill = .data[[data.column]])) +
      geom_bar(position = "fill", stat = "identity") +
      ylab("Frequency") +
      xlab("Cluster") +
      ggtitle(plot.title) +
      theme(legend.position = "none",
            panel.background = element_rect("white"),
            axis.line = element_line(color = "black", size = 0.4),
            axis.text.y = element_blank()
      )
  }

  gridExtra::grid.arrange(grobs = plots, ncol = min(length(plots), 2))
}

#' Plot of clonotype expansion
#'
#' @param object Seurat object
#' @param reduction Which dimensionality reduction to use
#' @param clonotype.column Metadata column with clonotype information. Default = 'clonotype'
#' @param min.color Color for cells without clonotype. Default = grey
#' @param max.color Color for cells with the highest number of clonotypes. Default = red
#'
#' @importFrom dplyr %>% add_count all_of select
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom ggplot2 aes geom_point scale_color_gradient
#'
#' @export

plot_expansion <- function(object, reduction, clonotype.column = 'clonotype', min.color = NULL, max.color = NULL) {

  if (!clonotype.column %in% colnames(object@meta.data)) {
    stop("Invalid clonotype column ", clonotype.column, call. = F)
  }

  if (is.null(min.color)) {
    min.color <- 'grey'
  }

  if (is.null(max.color)) {
    max.color <- 'red'
  }

  if (!reduction %in% names(object@reductions)) {
    stop("Invalid reduction ", reduction, call. = F)
  }

  coordinates <- object@reductions[[reduction]]@cell.embeddings[,c(1,2)] %>% as.data.frame()
  x.name <- colnames(coordinates)[[1]]
  y.name <- colnames(coordinates)[[2]]

  data <- object@meta.data %>%
    rownames_to_column("barcode") %>%
    select(all_of(c("barcode", clonotype.column))) %>%
    na.omit() %>%
    add_count(.data$clonotype) %>%
    column_to_rownames("barcode")

  plot.data <- coordinates
  plot.data$clonotype_count <- 0
  plot.data[rownames(data), 'clonotype_count'] <- data$n %>% as.numeric()

  ggplot() +
    geom_point(data = subset(plot.data, clonotype_count == 0), aes(x = .data[[x.name]], y = .data[[y.name]], color = .data$clonotype_count)) +
    geom_point(data = subset(plot.data, clonotype_count > 0), aes(x = .data[[x.name]], y = .data[[y.name]], color = .data$clonotype_count)) +
    scale_color_gradient(low = min.color, high = max.color)
}
