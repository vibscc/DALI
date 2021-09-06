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
#' @param add.missing.families Should missing families be added to the plot. Default = TRUE
#' @param percent.total Should the fraction of cells be calculated from the total number or cells in the group or just the cells with VDJ info. Default = TRUE (= from total)
#' @param show.missing.values Should missing values be shown in the plot. Default = FALSE
#' @param cols Colors to use
#'
#' @importFrom dplyr arrange case_when count rename  %>%
#' @importFrom ggplot2 aes element_blank element_rect element_text facet_grid geom_bar geom_text ggplot labs position_dodge scale_fill_manual theme unit ylim
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr spread
#'
#' @export

BarplotChainRegion <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  group.by = NULL,
  region = c("V", "C"),
  chain = AvailableChains(object),
  by.family = T,
  legend = T,
  grid = F,
  add.missing.families = T,
  percent.total = T,
  show.missing.values = F,
  cols = NULL
) {
  region <- match.arg(region) %>% tolower()
  chain <- match.arg(chain) %>% tolower()

  data.column <- paste0(chain, '.', region, '_')

  if (by.family && region == 'v') {
      data.column <- paste0(data.column, 'fam')
  } else {
      data.column <- paste0(data.column, 'gene')
      by.family <- F
  }

  if (is.null(group.by)) {
      group.by <- 'seurat_clusters'
  }

  if (!group.by %in% colnames(object@meta.data)) {
      stop("Invalid group.by column ", group.by)
  }

  if (!is.null(ident.1) && !ident.1 %in% object@meta.data[[group.by]]) {
    stop("Invalid ident.1")
  }

  if (!is.null(ident.2) && !ident.2 %in% object@meta.data[[group.by]]) {
    stop("Invalid ident.2")
  }

  if (is.null(ident.1) && !is.null(ident.2)) {
    stop("Can't specify ident.2 without ident.1")
  }

  families <- object@meta.data[, data.column] %>% na.omit() %>% unique()

  # Add missing families
  if (by.family && add.missing.families) {
      families.completed <- c()
      prefixes <- gsub("/.*$", "", gsub("[0-9-]", "", families)) %>% unique()

      for (prefix in prefixes) {
        families.with.prefix <- families[grepl(prefix, families)]

        # Ignore all families that contain a / in the name
        # These families will just be appended to the final families without attempting to complete the missing families
        families.ignored <- families.with.prefix[grepl("/", families.with.prefix)]

        # Only keep families without / in the name
        families.with.prefix <- setdiff(families.with.prefix, families.ignored)

        family.numbers <- c()
        for (family in families.with.prefix) {
            family.numbers <- c(family.numbers, gsub("[A-Za-z-]", "", family)) %>% as.numeric()
        }

        if (length(family.numbers) > 0) {
          families.completed <- c(families.completed, paste0(prefix, '-', seq(1,max(family.numbers))))
        }

        families.completed <- c(families.completed, families.ignored)
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
      filter(case_when((!percent.total && !show.missing.values) ~ !is.na(.data[[data.column]]),
                       T ~ T)) %>%
      group_by(.data[[group.by]]) %>%
      mutate(freq = prop.table(.data$n) * 100) %>%
      mutate(freq = round(.data$freq, 2)) %>%
      select(.data[[data.column]], .data[[group.by]], .data$freq) %>%
      filter(case_when(!show.missing.values ~ !is.na(.data[[data.column]]),
                        T ~ T)) %>%
      spread(.data[[group.by]], .data$freq) %>%
      mutate(family = replace(.data[[data.column]], is.na(.data[[data.column]]), "UNKNOWN")) %>%
      select(-all_of(data.column)) %>%
      replace(is.na(.), 0)

  if (nrow(data) == 0) {
    stop("Provided identities don't have any VDJ data. Can't plot without data!")
  }

  missing.families <- setdiff(families, data$family)

  if (length(missing.families) > 0) {
      missing.data <- matrix(0, nrow = length(missing.families), ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
      missing.data <- cbind(missing.data, missing.families)
      colnames(missing.data) <- colnames(data)
      data <- rbind(data, missing.data)
  }

  plot.data <- data %>%
      melt(data = ., id.vars = "family", variable.name = "group", value.name = "freq") %>%
      mutate(freq = as.numeric(.data$freq))

  plot.data[[data.column]] <- factor(plot.data$family, levels = families)

  if (is.null(cols)) {
    cols <- GetCategoricalColorPalette(plot.data$group)
  }

  plot <- ggplot(plot.data, aes(x = .data[[data.column]], y = .data$freq, fill = .data$group)) +
      geom_bar(position = "dodge", stat = "identity") +
      ylim(0, NA) +
      labs(y = "Percentage cells", x = "Family", title = .data$group) +
      theme(
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          legend.key.size = unit(0.1, "cm"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      scale_fill_manual(values = cols)

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

#' Barplot with clonotype distribution
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param subset Subset data to these groups
#' @param clonotypes Clonotypes to plot. Default = top 10
#' @param position Position of the bars in the plots. Options = stack or dodge
#' @param cols Colors to use
#'
#' @importFrom dplyr %>% filter group_by n select summarise
#' @importFrom ggplot2 aes element_text geom_bar ggplot theme
#' @importFrom rlang .data
#' @importFrom stats na.omit
#'
#' @export

BarplotClonotypes <- function(object, group.by = NULL, subset = NULL, clonotypes = NULL, position = c("stack", "dodge"), cols = NULL) {

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

  if (is.null(clonotypes)) {
    clonotypes.ordered <- table(object@meta.data$clonotype) %>% sort(decreasing = T) %>% names()
    clonotypes <- clonotypes.ordered[1:10]
  }

  position <- match.arg(position)

  clonotypes.missing <- setdiff(clonotypes, unique(object@meta.data$clonotype))

  if (length(clonotypes.missing) == length(clonotypes)) {
    stop("Could not find any of the provided clonotypes")
  }

  if (length(clonotypes.missing) > 0) {
    message(paste0("Could not find following clonotypes: ", paste(clonotypes.missing, collapse = ", ")))
  }

  plot.data <- object@meta.data %>%
    select(.data[[group.by]], .data$clonotype) %>%
    filter(.data$clonotype %in% clonotypes) %>%
    na.omit() %>%
    group_by(.data[[group.by]], .data$clonotype) %>%
    summarise(n = n())

  plot.data$clonotype <- factor(plot.data$clonotype, levels = clonotypes)

  if (is.null(cols)) {
    cols <- GetCategoricalColorPalette(plot.data[[group.by]])
  }

  ggplot(plot.data, aes(x = .data$clonotype, y = .data$n, fill = .data[[group.by]])) +
    geom_bar(position = position, stat = "identity") +
    theme(
      axis.text.x = element_text(angle = 90)
    ) +
    scale_fill_manual(values = cols)
}

#' Circosplot for family to gene distribution
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param subset Subset data to these groups
#' @param seed Random seed to use. Using the same seed ensures that colors in the plot will be identical. Default = NULL
#'
#' @importFrom circlize chordDiagram circos.track circos.text CELL_META
#' @importFrom dplyr %>% select
#' @importFrom graphics strwidth
#' @importFrom rlang .data
#' @importFrom stats na.omit
#'
#' @export

CircosPlot <- function(object, group.by = NULL, subset = NULL, seed = NULL) {

    if (!is.null(seed)) {
      set.seed(seed)
    }

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

    type <- object@misc$default.assay.VDJ

    family.column <- paste0(if (type == "TCR") "a" else "h", ".v_fam")
    gene.column <- paste0(if (type == "TCR") "b" else "l", ".v_gene")

    plot.data <- object@meta.data %>%
        select(.data[[family.column]], .data[[gene.column]]) %>%
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

#' Plot length of CDR3 AA/NT sequences
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param subset Subset data to these groups
#' @param plot.type Type of plot. Options = ridge, line
#' @param sequence.type Which CDR3 sequence should be used? Options: "AA" or "NT"
#' @param cols Colors to use. Only has an effect on the ridgeplots
#'
#' @export

CDR3Plot <- function(object, group.by = NULL, subset = NULL, plot.type = c("ridge", "line"), sequence.type = c("AA", "NT"), cols = NULL) {

    plot.type <- match.arg(plot.type)
    sequence.type <- match.arg(sequence.type)

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

    type <- DefaultAssayVDJ(object)
    cdr3.sequence <- if (sequence.type == "AA") ".cdr3" else ".cdr3_nt"
    heavy.cdr3.column <- paste0(if (type == "TCR") "a" else "h", cdr3.sequence)
    light.cdr3.column <- paste0(if (type == "TCR") "b" else "l", cdr3.sequence)

    if (plot.type == "line") {
      plots <- CDR3Plot.line(object, group.by, heavy.cdr3.column, light.cdr3.column, sequence.type)
    } else if (plot.type == "ridge") {
      plots <- CDR3Plot.ridge(object, group.by, heavy.cdr3.column, light.cdr3.column, cols)
    }

    gridExtra::grid.arrange(grobs = plots, ncol = min(length(plots), 3))
}

#' Plot length of CDR3 aa sequences as a line plot
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param heavy.cdr3.column Column name for heavy cdr3 data
#' @param light.cdr3.column Column name for light cdr3 data
#' @param sequence.type AA or NT
#'
#' @importFrom dplyr %>% mutate select full_join count
#' @importFrom ggplot2 ggplot aes geom_line labs theme element_rect element_line element_blank unit element_text
#' @importFrom rlang .data
#' @importFrom stats na.omit

CDR3Plot.line <- function(object, group.by, heavy.cdr3.column, light.cdr3.column, sequence.type) {
  plots <- list()
  groups <- unique(object@meta.data[[group.by]]) %>% gtools::mixedsort(x = .)

  for (group in groups) {
    cells <- rownames(object@meta.data)[object@meta.data[[group.by]] == group]
    subset <- subset(object, cells = cells)

    plot.data.h <- subset@meta.data %>%
      mutate(len = nchar(.data[[heavy.cdr3.column]])) %>%
      count(.data$len) %>%
      na.omit() %>%
      mutate(freq = .data$n/sum(.data$n) * 100) %>%
      select(.data$len, .data$freq)
    plot.data.l <- subset@meta.data %>%
      mutate(len = nchar(.data[[light.cdr3.column]])) %>%
      count(.data$len) %>%
      na.omit() %>%
      mutate(freq = .data$n/sum(.data$n) * 100) %>%
      select(.data$len, .data$freq)

    plot.data <- full_join(plot.data.h, plot.data.l, by = "len") %>% replace(is.na(.), 0)
    colnames(plot.data) <- c("cdr3.length", "heavy.chain", "light.chain")

    plots[[group]] <- ggplot(plot.data, aes(x = .data$cdr3.length)) +
      geom_line(aes(y = .data$heavy.chain), color = "black") +
      geom_line(aes(y = .data$light.chain), color = "red") +
      labs(x = paste0("CDR3 length (", sequence.type,")"), y = "Frequency of cells", title = paste0("CDR3 length - ", group)) +
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

  return(plots)
}

#' Plot length of CDR3 aa sequences as a ridgeplot
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param heavy.cdr3.column Column name for heavy cdr3 data
#' @param light.cdr3.column Column name for light cdr3 data
#' @param cols Colors to use
#'
#' @importFrom dplyr %>% mutate
#' @importFrom ggplot2 ggplot aes coord_cartesian scale_y_discrete scale_x_continuous
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom rlang .data

CDR3Plot.ridge <- function(object, group.by, heavy.cdr3.column, light.cdr3.column, cols) {
  plots <- list()

  for (column in c(heavy.cdr3.column, light.cdr3.column)) {

    plot.data <- object@meta.data %>%
      mutate(len = nchar(.data[[column]]))

    if (column == heavy.cdr3.column) {
      if (DefaultAssayVDJ(object) == 'TCR') {
        title <- "TCR-a"
      } else {
        title <- "Heavy"
      }
    } else {
      if (DefaultAssayVDJ(object) == 'TCR') {
        title <- "TCR-b"
      } else {
        title <- "Light"
      }
    }

    title <- paste0(title, " chain")

    if (is.null(cols)) {
      cols <- GetCategoricalColorPalette(plot.data[[group.by]])
    }

    plots[[column]] <- ggplot(plot.data, aes(x = .data$len, y = .data[[group.by]], fill = .data[[group.by]])) +
      geom_density_ridges() +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = cols) +
      coord_cartesian(clip = "off") +
      theme_ridges(grid = F, center_axis_labels = T) +
      ggtitle(title)
  }

  return(plots)
}

#' Dimplot for IGXV-family
#'
#' @param object Seurat object
#' @param region Region to plot. Available options: 'V'(ariable), 'D', 'J'(unction), 'C'(onstant).
#' @param chain Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR
#' @param by.family Group genes of 1 family together. Only effective with the V-gene. Default = TRUE
#' @param grid If TRUE, show per gene type in grid. If FALSE, show all genes types together on plot. Default = TRUE
#' @param highlight Family or gene to highlight. Default = NULL
#' @param ... Extra parameters passed to Seurat::Dimplot
#'
#' @importFrom dplyr %>%
#'
#' @export

DimplotChainRegion <- function(object, region = c("V", "D", "J", "C"), chain = AvailableChains(object), by.family = T, grid = T, highlight = NULL, cols = NULL, ...) {

  region <- match.arg(region) %>% tolower()
  chain <- match.arg(chain) %>% tolower()

  if (length(highlight) > 1) {
    stop("Can only select 1 family/gene to highlight", call. = F)
  }

  data.column <- paste0(chain, '.', region, '_')

  if (by.family && region == 'v') {
    data.column <- paste0(data.column, 'fam')
  } else {
    data.column <- paste0(data.column, 'gene')
  }

  if (is.na(object@meta.data[, data.column]) %>% sum() == nrow(object@meta.data)) {
    stop("No data for the combination ", toupper(region), "-region + ", toupper(chain), "-chain", call. = F)
  }

  families <- object@meta.data[, data.column] %>% na.omit() %>% unique()
  families <- families[families %>% gsub(pattern = "-", replacement = "_") %>% gtools::mixedorder(x = .)]

  if (!is.null(highlight) && !highlight %in% families) {
    stop("Invalid highlight for selected region/chain combination", call. = F)
  }

  split <- data.column
  if (!grid || !is.null(highlight)) {
    split <- NULL
  }

  cells.highlight = NULL
  if (!is.null(highlight)) {
    cells.highlight = rownames(object@meta.data)[object@meta.data[[data.column]] %in% highlight]
  }

  if (is.null(cols)) {
    cols <- GetCategoricalColorPalette(object@meta.data[[data.column]])
  }

  Seurat::DimPlot(object, group.by = data.column, split.by = split, cells.highlight = cells.highlight, order = rev(families), cols = cols, ...) +
    theme(
      panel.grid.major = element_line(color = "grey"),
      axis.line = element_blank()
    )
}


#' Plot the frequency of clonotypes in (subset of) clusters
#'
#' @param object Seurat object
#' @param chain Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR. NULL = both
#' @param group.by Metadata column to group the family data by. Default = seurat_clusters
#' @param subset Subset data to these groups
#' @param use.sequence Use AA/NT sequence instead of clonotype. Default = FALSE
#' @param sequence.type What sequences to use, available options: "AA" or "NT". Only functional when use.sequence = T
#' @param clonotype.column Metadata column with clonotype information. Default = 'clonotype'
#' @param bulk Group all cells together and handle as bulk. Default = FALSE
#' @param show.missing Show missing values in plot. Default = FALSE
#' @param plot.type Type of plot, options = "bar" or "violin"
#' @param threshold Highlight clonotypes with more or equal cells than threshold. Only works with barplot. Default = 1
#'
#' @importFrom dplyr %>%
#'
#' @export

ClonotypeFrequency <- function(
  object,
  chain = AvailableChains(object),
  group.by = NULL,
  subset = NULL,
  use.sequence = F,
  sequence.type = c("AA", "NT"),
  clonotype.column = NULL,
  bulk = F,
  show.missing = F,
  plot.type = c("bar", "violin"),
  threshold = 1
) {

  if (!is.null(chain)) {
    chain <- match.arg(chain) %>% tolower()
  }

  if (is.null(group.by)) {
    group.by <- "seurat_clusters"
  }

  if (is.null(clonotype.column)) {
    clonotype.column <- "clonotype"
  }

  plot.type <- match.arg(plot.type)

  if (!use.sequence && !clonotype.column %in% colnames(object@meta.data)) {
    stop("Invalid clonotype column ", clonotype.column, call. = F)
  }

  if (!is.null(subset) && subset != '') {
    cells <- rownames(object@meta.data)[object@meta.data[[group.by]] %in% subset]
    object <- subset(object, cells = cells)
  }

  sequence.type <- match.arg(sequence.type)

  chains <- chain

  if (is.null(chains)) {
    chains <- AvailableChains(object) %>% tolower()
  }

  if (plot.type == "bar") {
    plots <- ClonotypeFrequency.bar(object, chains, group.by, use.sequence, sequence.type, clonotype.column, show.missing, bulk, threshold)
  } else if (plot.type == "violin") {
    plots <- ClonotypeFrequency.violin(object, chains, group.by, use.sequence, sequence.type, clonotype.column, show.missing, bulk, cols)
  }

  gridExtra::grid.arrange(grobs = plots, ncol = min(length(plots), 2))
}

#' Plot the frequency of clonotypes in (subset of) clusters in a barplot
#'
#' @param object Seurat object
#' @param chains Vector of chain(s) to plot
#' @param group.by Metadata column to group the family data by
#' @param use.sequence Use AA/NT sequence instead of clonotype
#' @param sequence.type What sequences to use, available options: "AA" or "NT"
#' @param clonotype.column Metadata column with clonotype information
#' @param show.missing Show missing values in plot
#' @param bulk Group all cells together and handle as bulk
#' @param threshold Highlight clonotypes with more or equal cells than threshold
#'
#' @importFrom ggplot2 aes element_line element_rect geom_bar ggplot ggtitle scale_color_manual theme xlab ylab
#' @importFrom dplyr %>% arrange case_when group_by n summarise ungroup

ClonotypeFrequency.bar <- function(
  object,
  chains,
  group.by,
  use.sequence,
  sequence.type,
  clonotype.column,
  show.missing,
  bulk,
  threshold
) {
  plots <- list()

  for (chain in chains) {
    data.column <- clonotype.column
    plot.title <- paste0("Clonotype frequency")

    if (use.sequence) {
      plot.title <- paste0("CDR3 ", toupper(chain), "-chain ", sequence.type, " sequence frequency")
      data.column <- paste0(chain, ".cdr3")

      if (sequence.type == "NT") {
        data.column <- paste0(data.column, "_nt")
      }
    }

    plot.data <- CalculateFrequency(object, data.column, group.by, show.missing) %>% arrange(.data$freq)

    plot.data[[data.column]] <- factor(plot.data[[data.column]], levels = unique(plot.data[[data.column]]))

    if (bulk) {
      plot.data[[group.by]] <- "dataset"
    }

    y.lab <- paste0("Frequency (", if (show.missing) "all cells" else "cells with VDJ data", ")")
    plots[[data.column]] <- ggplot(plot.data, aes(x = .data[[group.by]], y = .data$freq, color = (!is.na(.data$clonotype) & .data$freq >= threshold))) +
    # plots[[data.column]] <- ggplot(plot.data, aes(x = .data[[group.by]], y = .data$freq, color = .data$freq > threshold)) +
      geom_bar(position = "fill", stat = "identity", alpha = 0.1) +
      ylab(y.lab) +
      xlab("Group") +
      ggtitle(plot.title) +
      theme(
        legend.position = "none",
        panel.background = element_rect("white"),
        axis.line = element_line(color = "black", size = 0.4),
      ) +
      scale_color_manual(values = c(
        "TRUE" = "red",
        "FALSE" = "lightgrey"
      ))
  }

  return(plots)
}

#' Plot the frequency of clonotypes in (subset of) clusters in a violinplot
#'
#' @param object Seurat object
#' @param chains Vector of chain(s) to plot
#' @param group.by Metadata column to group the family data by
#' @param use.sequence Use AA/NT sequence instead of clonotype
#' @param sequence.type What sequences to use, available options: "AA" or "NT"
#' @param clonotype.column Metadata column with clonotype information
#' @param show.missing Show missing values in plot
#' @param bulk Group all cells together and handle as bulk
#'
#' @importFrom ggplot2 aes element_line geom_jitter geom_violin ggplot ggtitle theme xlab ylab
#' @importFrom dplyr %>% arrange case_when group_by n summarise ungroup

ClonotypeFrequency.violin <- function(
  object,
  chains,
  group.by,
  use.sequence,
  sequence.type,
  clonotype.column,
  show.missing,
  bulk
) {
  plots <- list()

  for (chain in chains) {

    data.column <- clonotype.column
    plot.title <- paste0("Clonotype frequency")

    if (use.sequence) {
      plot.title <- paste0("CDR3 ", toupper(chain), "-chain ", sequence.type, " sequence frequency")
      data.column <- paste0(chain, ".cdr3")

      if (sequence.type == "NT") {
        data.column <- paste0(data.column, "_nt")
      }
    }

    plot.data <- CalculateFrequency(object, data.column, group.by, show.missing) %>%
      group_by(.data[[group.by]]) %>%
      mutate(freq = prop.table(.data$n) * 100) %>%
      mutate(freq = round(.data$freq, 2))


    plots[[data.column]] <- ggplot(plot.data, aes(x = .data[[group.by]], y = .data$freq)) +
      geom_violin(scale = "width") +
      geom_jitter(aes(color = .data[[data.column]])) +
      ylab("Frequency") +
      xlab("Cluster") +
      ggtitle(plot.title) +
      theme(
        legend.position = if (legend) "right" else "none"
      )
  }

  return(plots)
}

#' Plot of clonotype expansion
#'
#' @param object Seurat object
#' @param reduction Which dimensionality reduction to use
#' @param clonotype.column Metadata column with clonotype information. Default = 'clonotype'
#' @param color.low Color for cells without clonotype. Default = lightgrey
#' @param color.mid Color for cells with the medium number of clonotypes. Default = red
#' @param color.high Color for cells with the highest number of clonotypes. Default = darkred
#' @param threshold Cells with clonotype count < threshold are colored with the min.color. Default = 1
#' @param positive.size Size of the dots with clonotype count > threshold. Default = 0.5
#' @param negative.size Size of the dots with clonotype count <= threshold. Default = 0.5
#' @param positive.alpha Alpha of the dots with clonotype count > threshold. Default = 1
#' @param negative.alpha Alpha of the dots with clonotype count < threshold. Default = 1
#'
#' @importFrom dplyr %>% add_count all_of select
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom ggplot2 aes geom_point scale_color_gradient2 scale_color_gradientn theme_classic
#'
#' @export

ExpansionPlot <- function(
  object,
  reduction,
  clonotype.column = 'clonotype',
  color.low = 'lightgrey',
  color.mid = 'red',
  color.high = 'darkred',
  threshold = 1,
  positive.size = 0.5,
  negative.size = 0.5,
  positive.alpha = 1,
  negative.alpha = 1
) {

  if (!clonotype.column %in% colnames(object@meta.data)) {
    stop("Invalid clonotype column ", clonotype.column, call. = F)
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

  max_value <- max(plot.data$clonotype_count)

  midpoint <- (max_value + threshold) / 2

  plot <- ggplot() +
    geom_point(
      data = subset(plot.data, clonotype_count <= threshold),
      aes(
        x = .data[[x.name]], y = .data[[y.name]],
        color = .data$clonotype_count
      ),
      size = negative.size,
      alpha = negative.alpha
    ) +
    geom_point(
      data = subset(plot.data, clonotype_count > threshold),
      aes(
        x = .data[[x.name]],
        y = .data[[y.name]],
        color = .data$clonotype_count
      ),
      size = positive.size,
      alpha = positive.alpha
    ) +
    theme_classic()

  if (threshold > 0) {
    plot <- plot + scale_color_gradientn(colours = c(color.low, color.low, color.mid, color.high), values = c(0, threshold / max_value, midpoint / max_value, 1))
  } else {
    plot <- plot + scale_color_gradient2(low = color.low, mid = color.mid, high = color.high, midpoint = midpoint)
  }

  plot
}

#' Featureplot of clonotypes
#'
#' @param object Seurat object
#' @param reduction Which dimensionality reduction to use
#' @param clonotypes Clonotypes to plot
#' @param clonotype.column Metadata column with clonotype information. Default = 'clonotype'
#' @param size Size of the dots with highlighted clonotype. Default = 0.5
#' @param missing.size Size of the dots without clonotype. Default = 0.5
#' @param alpha Alpha of the dots with highlighted clonotype. Default = 1
#' @param missing.alpha Alpha of the dots without clonotype. Default = 1
#' @param missing.color Color of dots without clonotype. Default = lightgrey
#'
#' @importFrom dplyr %>% mutate
#' @importFrom ggplot2 aes geom_point ggplot scale_color_discrete theme_classic
#'
#' @export

FeaturePlotChainRegion <- function(
  object,
  reduction,
  clonotypes,
  clonotype.column = NULL,
  size = 0.5,
  missing.size = 0.5,
  alpha = 1,
  missing.alpha = 1,
  missing.color = "lightgrey"
  ) {

  if (is.null(clonotype.column)) {
    clonotype.column <- "clonotype"
  }

  if (!clonotype.column %in% colnames(object@meta.data)) {
    stop("Invalid clonotype column ", clonotype.column, call. = F)
  }

  if (!reduction %in% names(object@reductions)) {
    stop("Invalid reduction ", reduction, call. = F)
  }

  invalid.clonotypes <- setdiff(clonotypes, unique(object@meta.data[[clonotype.column]]))
  if (length(invalid.clonotypes) > 0) {
    stop("Invalid clonotypes: ", paste(invalid.clonotypes, collapse = ", "), call. = F)
  }

  coordinates <- object@reductions[[reduction]]@cell.embeddings[,c(1,2)] %>% as.data.frame()
  x.name <- colnames(coordinates)[[1]]
  y.name <- colnames(coordinates)[[2]]

  object@meta.data <- object@meta.data %>%
    mutate(clonotypes = ifelse(.data[[clonotype.column]] %in% clonotypes, .data[[clonotype.column]], NA))

  plot.data <- coordinates
  plot.data$clonotypes <- object@meta.data[rownames(plot.data), 'clonotypes']

  ggplot() +
    geom_point(
      data = subset(plot.data, is.na(clonotypes)),
      aes(
        x = .data[[x.name]], y = .data[[y.name]],
        color = .data$clonotypes
      ),
      size = missing.size,
      alpha = missing.alpha
    ) +
    geom_point(
      data = subset(plot.data, !is.na(clonotypes)),
      aes(
        x = .data[[x.name]],
        y = .data[[y.name]],
        color = .data$clonotypes
      ),
      size = size,
      alpha = alpha
    ) +
    theme_classic() +
    scale_color_discrete(na.value = missing.color)
}

#' Display connections between clusters
#'
#' @param object Seurat object
#' @param reduction Dimensionality reduction
#' @param group.by  Metadata column to group the family data by. Default = seurat_clusters
#' @param groups.highlight Groups for which to highlight the edges in the graph
#' @param clonotype.column Metadata column with clonotype information. Default = clonotype
#' @param cols Colors to use
#'
#' @importFrom dplyr %>% bind_rows distinct filter group_by group_map n
#' @importFrom ggplot2 aes element_blank geom_point scale_color_manual theme
#' @importFrom ggraph geom_edge_arc geom_node_point ggraph scale_edge_alpha scale_edge_color_manual scale_edge_width
#' @importFrom rlang .data
#' @importFrom tidygraph as_tbl_graph
#'
#' @export

CloneConnGraph <- function(object, reduction, group.by = NULL, groups.highlight = NULL, clonotype.column = NULL, cols = NULL) {

  if (is.null(group.by)) {
    group.by <- "seurat_clusters"
  }

  if (!reduction %in% names(object@reductions)) {
    stop("Invalid redution ", reduction, call. = F)
  }

  if (!group.by %in% colnames(object@meta.data)) {
    stop("Invalid group.by column ", group.by)
  }

  if (is.null(clonotype.column)) {
    clonotype.column <- "clonotype"
  }

  if (!clonotype.column %in% colnames(object@meta.data)) {
    stop("Invalid clonotype.column ", group.by)
  }

  object@meta.data[, group.by] <- as.character(object@meta.data[, group.by])

  edges <- object@meta.data[, c(clonotype.column, group.by)] %>%
    na.omit() %>%
    group_by(clonotype) %>%
    distinct() %>%
    filter(n() > 1) %>%
    group_map(~ gtools::combinations(nrow(.x), 2, as.numeric(.x[[1]])) %>% as.data.frame()) %>%
    bind_rows()

  colnames(edges) <- c("from", "to")
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)

  edges <- edges %>% group_by(.data$from, .data$to) %>% count()

  if (!is.null(groups.highlight)) {
    edges$highlight <- F
    edges[edges$from %in% groups.highlight, 'highlight'] <- T
    edges[edges$to %in% groups.highlight, 'highlight'] <- T
  }

  getCenters <- function(groups, column) {
    centers <- c()
    for (group in groups) {
      cells <- rownames(object@meta.data)[object@meta.data[, group.by] == group]
      test <- mean(object@reductions[[reduction]]@cell.embeddings[cells, column])

      centers <- c(centers, test)
    }
    return(centers)
  }

  graph <- as_tbl_graph(edges, directed = F) %>%
    mutate(x = getCenters(.data$name, 1), y = getCenters(.data$name, 2))

  dimred <- object@reductions[[reduction]]@cell.embeddings %>% as.data.frame()
  dimred$group <- factor(object@meta.data[[group.by]], levels = gtools::mixedsort(unique(object@meta.data[[group.by]])))

  label.x.axis <- colnames(dimred)[[1]]
  label.y.axis <- colnames(dimred)[[2]]

  if (is.null(cols)) {
    cols <- GetCategoricalColorPalette(dimred$group)
  }

  plot <- ggraph(graph, layout = "manual", x = .data$x, y = .data$y) +
    geom_point(data = dimred, aes(x = .data[[label.x.axis]], y = .data[[label.y.axis]], color = .data$group)) +
    scale_color_manual(values = cols) +
    scale_edge_alpha(range = c(0.1, 1)) +
    scale_edge_width(range = c(0.5, 2)) +
    theme(
      panel.background = element_blank()
    )

  if (!is.null(groups.highlight)) {
    plot <- plot +
      geom_edge_arc(aes(color = .data$highlight, alpha = .data$n, width = .data$n), strength = 0.2) +
      scale_edge_color_manual(values = c("grey", "red"))
  } else {
    plot <- plot +
      geom_edge_arc(aes(alpha = .data$n, width = .data$n), strength = 0.2)
  }

  plot + geom_node_point()
}
