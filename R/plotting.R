#' Barplot for IGHV-family per group
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param groups.to.plot Which groups of the group by column should be plotted. Default = NULL (all)
#' @param region Region to plot. Available options: 'V'(ariable) or 'C'(onstant)
#' @param chain Chain to plot. Options: 'H'(eavy) or 'L'(ight)
#' @param by.family Group genes of 1 family together. Default = TRUE
#'
#' @importFrom dplyr %>% arrange count case_when rename
#' @importFrom ggplot2 ggplot geom_col labs aes scale_fill_manual geom_label theme element_rect element_blank element_text unit ylim
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr spread
#'
#' @export

barplot_vh <- function(object, group.by = NULL, groups.to.plot = NULL, region = c("V", "C"), chain = c("H", "L"), by.family = T) {
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
        prefix <- gsub("[0-9-]", "", families[1])
        family.numbers <- c()
        for (family in families) {
            family.numbers <- c(family.numbers, gsub("[A-Za-z-]", "", family) %>% as.numeric())
        }
        families <- paste0(prefix, '-', seq(1,max(family.numbers)))
    }

    families <- families %>% gtools::mixedsort(x = ., decreasing = sum(grepl('-', .)) > 0)

    data <- object@meta.data %>%
        filter(case_when(!is.null(groups.to.plot) ~ .data[[group.by]] %in% groups.to.plot,
                         T ~ T)) %>%
        count(.data[[data.column]], .data[[group.by]]) %>%
        na.omit() %>%
        spread(.data[[group.by]], .data$n) %>%
        replace(is.na(.), 0)

    missing.families <- setdiff(families, data[[data.column]])

    if (length(missing.families) > 0) {
        missing.data <- matrix(0, nrow = length(missing.families), ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
        missing.data <- cbind(missing.families, missing.data)
        colnames(missing.data) <- c(data.column, colnames(data)[-1])
        data <- rbind(data, missing.data)
    }

    data <- data %>%
        arrange(factor(.data[[data.column]], levels = families))  %>%
        column_to_rownames(data.column)

    plots <- list()

    for (group in colnames(data)) {
        plot.data <- data[, group] %>%
            as.data.frame() %>%
            rename(freq = .data[["."]]) %>%
            mutate(freq = as.numeric(.data$freq))
        plot.data$fam <- factor(rownames(data), levels = families)

        plots[[group]] <- ggplot(plot.data, aes(x = .data$fam, y = .data$freq, fill = .data$fam)) +
            geom_col() +
            ylim(0, NA) +
            labs(y = "Cell number", x = "Family", title = group) +
            scale_fill_manual(values = colorRampPalette(c("darkblue", "lightblue"))(nrow(data))) +
            geom_label(data = NULL, aes(label = .data$freq), fill = "white", size = 2) +
            theme(
                panel.background = element_rect(fill = "white"), # bg of the panel
                plot.background = element_rect(fill = "white"), # bg of the plot
                panel.grid.major = element_blank(), # get rid of major grid
                panel.grid.minor = element_blank(), # get rid of minor grid
                legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
                legend.position = "none",
                legend.key.size = unit(0.1, "cm"),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            )

    }

    gridExtra::grid.arrange(grobs = plots, ncol = min(length(plots), 3))
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
                    select(.data$h.v_fam, .data$l.v_gene) %>%
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

    for (group in groups) {
        cells <- rownames(object@meta.data)[object@meta.data[[group.by]] == group]
        subset <- subset(object, cells = cells)

        plot.data.h <- subset@meta.data %>%
            mutate(len = nchar(.data$h.cdr3)) %>%
            count(.data$len) %>%
            na.omit() %>%
            mutate(freq = .data$n/sum(.data$n) * 100) %>%
            select(.data$len, .data$freq)
        plot.data.l <- subset@meta.data %>%
            mutate(len = nchar(.data$l.cdr3)) %>%
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
#' @param fam Heavy or light chain family. Options are heavy or light, default = NULL
#' @param chain Heavy or light chain. Options are heavy/h or light/l, default = NULL
#' Is used in combination with gene.
#' @param gene c/j/d gene. Options are c, j or d, default = NULL
#' Is used in combination with chain.
#' @param reduction Specify which reduction to use, default = "tsne"
#' @param grid If TRUE, show per gene type in grid. If FALSE, show all genes types together on plot. Default=T
#'
#' @importFrom dplyr %>%
#'
#' @export

DimPlot_vh <- function(object, fam = NULL, chain = NULL, gene = NULL, reduction = "tsne", grid = T) {

  if (is.null(fam)) {
    if (is.null(chain) | is.null(gene)) {
      stop("Either fam or either chain and gene should be given a value")
    } else{
      if (chain == 'heavy') {chain <- "h"}
      if (chain == 'light') {chain <- "l"}
      columnname <- paste0(chain,".",gene,"_gene")

      if (grid) {split <- columnname} else {split <- NULL}

      families <- object@meta.data[, columnname] %>% na.omit() %>% unique()
      families <- families %>% gtools::mixedsort(x = ., decreasing = sum(grepl('-', .)) > 0)
      Seurat::DimPlot(object, group.by = columnname, split.by = split, ncol = 3,order = rev(families))

    }
  } else {
    if (fam == 'heavy' || fam == 'h') { columnname <- "h.v_fam" }
    if (fam == 'light' || fam == 'l') { columnname <- "l.v_fam" }

    if (grid) {split <- columnname} else {split <- NULL}

    families <- object@meta.data[, columnname] %>% na.omit() %>% unique()
    families <- families %>% gtools::mixedsort(x = ., decreasing = sum(grepl('-', .)) > 0)
    Seurat::DimPlot(object, group.by = columnname, split.by = split, ncol = 4,order = rev(families))

  }
}


#' Dimplot for IGXV-family
#'
#' @param object Seurat object
#' @param columnname Chain to plot, available options:"l.cdr3","h.cdr3","cdr3" (=both)
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#'
#' @importFrom ggplot2 aes_string element_blank element_line element_rect geom_bar ggplot ggtitle theme xlab ylab
#' @importFrom dplyr count
#'
#' @export

CDR3freq <- function(object, columnname = c("l.cdr3","h.cdr3","cdr3"), group.by = 'seurat_clusters') {

  columnname <- match.arg(columnname)

  if (columnname == "l.cdr3") {TITLE <- "light chain"}
  if (columnname == "h.cdr3") {TITLE <- "heavy chain"}

  if (columnname == "h.cdr3" | columnname == "l.cdr3") {

    data <- count(object@meta.data, c(columnname, group.by) )

    p <- ggplot(data, aes_string(fill = columnname, y = .data$freq, x = group.by)) +
      geom_bar(position = "fill", stat = "identity") +
      ylab("Frequency") +
      xlab("Cluster") +
      ggtitle(paste0("CDR3 ", TITLE ," frequency" )) +
      theme(legend.position = "none",
            panel.background = element_rect("white"),
            axis.line = element_line(color = "black", size = 0.4),
            axis.text.y = element_blank())
  } else if (columnname == "cdr3") {
    data1 <- count(object@meta.data, c("l.cdr3", "seurat_clusters") )
    colnames(data1)[1] <- "cdr3"
    data2 <- count(object@meta.data, c("h.cdr3", "seurat_clusters") )
    colnames(data2)[1] <- "cdr3"


    p1 <- ggplot(data1, aes_string(fill = columnname, y = .data$freq, x = group.by)) +
      geom_bar(position = "fill", stat = "identity") +
      ylab("Frequency light chain") +
      # ggtitle(paste0("CDR3 ",TITLE ," frequency" )) +
      theme(legend.position = "none",
            panel.background = element_rect("White"),
            axis.line.y = element_line(color = "black", size = 0.4),
            axis.text = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank() )
    p2 <- ggplot(data2, aes_string(fill = columnname, y = .data$freq, x = group.by)) +
      geom_bar(position = "fill", stat = "identity") +
      ylab("Frequency heavy chain") +
      xlab("Cluster") +
      # ggtitle(paste0("CDR3 ",TITLE ," frequency" )) +
      theme(legend.position = "none",
            panel.background = element_rect("White"),
            axis.line = element_line(color = "black", size = 0.4),
            axis.text.y = element_blank())

    p <- gridExtra::grid.arrange(p1, p2, ncol = 1)

  }
  return(p)
}


