#' Barplot for IGHV-family per group
#'
#' @param object Seurat object
#' @param ident.1 Identy class(es) to plot
#' @param ident.2 Second identity class(es). This class will be used as comparison. If NULL, ident.1 will just be used to subset the data
#' @param group.by Metadata column to group the family data by.
#' @param region Region to plot. Available options: 'V'(ariable) or 'C'(onstant)
#' @param chain Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR
#' @param by.family Group genes of 1 family together. Default = TRUE
#' @param legend Should the legend be included in the plot. Default = TRUE
#' @param grid Organize plots in grid. Each plot contains information about 1 group. Default = FALSE
#' @param add.missing.families Should missing families be added to the plot. Default = TRUE
#' @param percent.total Should the fraction of cells be calculated from the total number or cells in the group or just the cells with VDJ info. Default = TRUE (= from total)
#' @param show.missing.values Should missing values be shown in the plot. Default = FALSE
#' @param color.theme Color theme to use. Default = DALI
#' @param colors Colors to use. This overwrites the selected color theme
#'
#' @importFrom dplyr case_when count %>%
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
  region = c("V", "J", "C"),
  chain = c("VDJ", "VJ"),
  by.family = T,
  legend = T,
  grid = F,
  add.missing.families = T,
  percent.total = T,
  show.missing.values = F,
  color.theme = ColorThemes(),
  colors = NULL
) {
  region <- match.arg(region) %>% tolower()
  chain <- match.arg(chain) %>% tolower()
  color.theme <- match.arg(color.theme)

  data.column <- GetDataColumn(chain, region, by.family)

  if (!grepl("fam", data.column)) {
    by.family <- F
  }

  if (is.null(group.by)) {
    object <- Seurat::AddMetaData(object, Seurat::Idents(object), "default.clustering")
    group.by <- "default.clustering"
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

  if (by.family && add.missing.families) {
    families <- AddMissingVDJFamilies(families)
  }

  families <- families %>% gtools::mixedsort(x = ., decreasing = sum(grepl("-", .)) > 0)

  data.filtered <- object@meta.data

  if (!is.null(ident.1) || !is.null(ident.2)) {
      data.filtered <- data.filtered %>%
          filter(.data[[group.by]] %in% c(ident.1, ident.2))
  }

  if (!is.null(ident.2)) {
      data.filtered[[group.by]] <- data.filtered[[group.by]] %>% as.character()

      for (ident in list(ident.1, ident.2)) {
          for (group in ident) {
              data.filtered[[group.by]] <- gsub(paste0("^", group, "$"), paste0(group.by, " (", paste(ident, collapse = ","), ")"), data.filtered[[group.by]])
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

  if (is.null(colors)) {
      if (color.theme == "Colorblind") {
        colors <- c("#FDE725", "#423C81")
      } else {
          colors <- GetCategoricalColorPalette(plot.data$group, color.theme)
      }
  }

  plot <- ggplot(plot.data, aes(x = .data[[data.column]], y = .data$freq, fill = .data$group)) +
      geom_bar(position = "dodge", stat = "identity") +
      ylim(0, NA) +
      labs(y = "Percentage cells", x = "Family") +
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
      scale_fill_manual(values = colors)

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

#' Heatmap showing to percentage of cells with a specifiek chain/region combination
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by.
#' @param chain Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR
#' @param region Region to plot. Available options: 'V'(ariable) or 'C'(onstant)
#' @param by.family Group genes of 1 family together. Default = TRUE
#' @param add.missing.families Should missing families be added to the plot. Default = TRUE
#' @param percent.total Should the fraction of cells be calculated from the total number or cells in the group or just the cells with VDJ info. Default = TRUE (= from total)
#' @param show.missing.values Should missing values be shown in the plot. Default = FALSE
#' @param cluster.rows Should rows (genes) be clustered in the heatmap. Default = FALSE
#' @param cluster.cols Should columns (groups) be clustered in the heatmap. Default = FALSE
#' @param color.scheme Colorscheme to use for the heatmap. Options: "coolwarm", "viridis". Default = "coolwarm"
#' @param ... parameters to pass to pheatmap::pheatmap()
#'
#' @importFrom dplyr case_when count filter group_by mutate select %>%
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr spread
#'
#' @export

HeatmapChainRegion <- function(
  object,
  group.by = NULL,
  chain = c("VDJ", "VJ"),
  region = c("V", "J", "C"),
  by.family = T,
  add.missing.families = T,
  percent.total = T,
  show.missing.values = F,
  cluster.rows = F,
  cluster.cols = F,
  color.scheme = c("coolwarm", "viridis", "gray to blue", "turning red"),
  ...
) {

  region <- match.arg(region) %>% tolower()
  chain <- match.arg(chain) %>% tolower()
  color.scheme <- match.arg(color.scheme)

  if (is.null(group.by)) {
    object <- Seurat::AddMetaData(object, Seurat::Idents(object), "default.clustering")
    group.by <- "default.clustering"
  }

  data.column <- GetDataColumn(chain, region, by.family)

  if (!grepl("fam", data.column)) {
    by.family <- F
  }

  families <- object@meta.data[, data.column] %>% na.omit() %>% unique()

  if (by.family && add.missing.families) {
    families <- AddMissingVDJFamilies(families)
  }

  families <- families %>% gtools::mixedsort(x = ., decreasing = sum(grepl("-", .)) > 0)

  data <- object@meta.data %>%
    count(.data[[data.column]], .data[[group.by]]) %>%
    filter(case_when((!percent.total && !show.missing.values) ~ !is.na(.data[[data.column]]),
                     T ~ T)) %>%
    group_by(.data[[group.by]]) %>%
    mutate(freq = prop.table(.data$n) * 100) %>%
    mutate(freq = round(.data$freq, 2)) %>%
    select(.data[[data.column]], .data[[group.by]], .data$freq) %>%
    spread(.data[[group.by]], .data$freq) %>%
    mutate(family = replace(.data[[data.column]], is.na(.data[[data.column]]), "NA")) %>%
    select(-all_of(data.column)) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames("family")

  if (show.missing.values) {
    families <- c(families, "NA")
  }

  plot.data <- data[families, ]
  rownames(plot.data) <- families

  return(pheatmap::pheatmap(plot.data, color = ColorScale(color.scheme), cluster_rows = cluster.rows, cluster_cols = cluster.cols, angle_col = 0))
}

#' Barplot with clonotype distribution
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by.
#' @param subset Subset data to these groups
#' @param clonotypes Clonotypes to plot. Default = top 10
#' @param position Position of the bars in the plots. Options = stack or dodge
#' @param color.theme Color theme to use. Default = DALI
#' @param colors Colors to use. This overwrites the selected color theme
#'
#' @importFrom dplyr %>% filter group_by n select summarise
#' @importFrom ggplot2 aes element_text geom_bar ggplot theme
#' @importFrom rlang .data
#' @importFrom stats na.omit
#'
#' @export

BarplotClonotypes <- function(
    object,
    group.by = NULL,
    subset = NULL,
    clonotypes = NULL,
    position = c("stack", "dodge"),
    color.theme = ColorThemes(),
    colors = NULL
) {
  color.theme <- match.arg(color.theme)

  if (is.null(group.by)) {
    object <- Seurat::AddMetaData(object, Seurat::Idents(object), "default.clustering")
    group.by <- "default.clustering"
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

  if (is.null(colors)) {
    colors <- GetCategoricalColorPalette(plot.data[[group.by]], color.theme)
  }

  ggplot(plot.data, aes(x = .data$clonotype, y = .data$n, fill = .data[[group.by]])) +
    geom_bar(position = position, stat = "identity") +
    theme(
      axis.text.x = element_text(angle = 90)
    ) +
    scale_fill_manual(values = colors)
}

#' Circosplot for family to gene distribution
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by.
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

CircosPlotGenes <- function(object, group.by = NULL, subset = NULL, seed = NULL) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (is.null(group.by)) {
        object <- Seurat::AddMetaData(object, Seurat::Idents(object), "default.clustering")
        group.by <- "default.clustering"
    }

    if (!group.by %in% colnames(object@meta.data)) {
        stop("Invalid group.by column ", group.by)
    }

    if (!is.null(subset)) {
        if (!(subset == "Select group...")) {
            cells <- rownames(object@meta.data)[object@meta.data[[group.by]] %in% subset]
            object <- subset(object, cells = cells)
        }
    }

    type <- object@misc$default.assay.VDJ

    family.column <- "vdj.v_fam"
    gene.column <- "vj.v_gene"

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


#' Circosplot for VDJ chain useage distribution
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by.
#' @param subset Subset data to these groups
#' @param color.theme Color theme to use. Default = DALI
#'
#' @importFrom circlize chordDiagram circos.track circos.text CELL_META
#' @importFrom dplyr %>% select
#' @importFrom graphics strwidth
#' @importFrom rlang .data
#' @importFrom stats na.omit
#'
#' @export

CircosPlotChains <- function(object, group.by = NULL, subset = NULL, color.theme = ColorThemes()) {
    color.theme <- match.arg(color.theme)

    if (is.null(group.by)) {
        object <- Seurat::AddMetaData(object, Seurat::Idents(object), "default.clustering")
        group.by <- "default.clustering"
    }

    if (!group.by %in% colnames(object@meta.data)) {
        stop("Invalid group.by column ", group.by)
    }

    if (!is.null(subset)) {
        if (!(subset == "Select group...")) {
            cells <- rownames(object@meta.data)[object@meta.data[[group.by]] %in% subset]
            object <- subset(object, cells = cells)
        }
    }

    v.gene.column <- "vdj.v_gene"
    d.gene.column <- "vdj.d_gene"
    j.gene.column <- "vdj.j_gene"

    plot.data_dj <- object@meta.data %>%
        select(.data[[d.gene.column]], .data[[j.gene.column]]) %>%
        na.omit() %>%
        table()
    plot.data_dv <- object@meta.data %>%
        select(.data[[d.gene.column]], .data[[v.gene.column]]) %>%
        na.omit() %>%
        table()
    plot.data_jv <- object@meta.data %>%
        select(.data[[j.gene.column]], .data[[v.gene.column]]) %>%
        na.omit() %>%
        table()

    plot.data_dj <- EqualiseTableRows(plot_data_x = plot.data_dj,plot_data_y = plot.data_dv)
    plot.data_dv <- EqualiseTableRows(plot_data_x = plot.data_dv,plot_data_y = plot.data_dj)
    plot.data_jv <- EqualiseTableRows(t(plot.data_jv),t(plot.data_dv)) %>% t()
    plot.data_dv <- EqualiseTableRows(t(plot.data_dv),t(plot.data_jv)) %>% t()

    plot.data_jj <- matrix(rep(0, nrow(plot.data_jv)*ncol(plot.data_dj)), nrow = nrow(plot.data_jv))
    colnames(plot.data_jj) <- colnames(plot.data_dj)
    rownames(plot.data_jj) <- rownames(plot.data_jv)

    right.data <- rbind(plot.data_dv,plot.data_jv)
    left.data <- rbind(plot.data_dj, plot.data_jj)
    plot.data <- cbind(left.data,right.data)

    nm <- unique(unlist(dimnames(plot.data)))
    group <- structure(gsub("[^VDJ]","" , nm), names = nm)

    col.vector <- unname(group) %>% as.factor()
    levels(col.vector) <- GetCategoricalColorPalette(levels(col.vector), color.theme)

    grid.col <- structure(col.vector,
                          names = names(group))

    plot.data <- plot.data[gtools::mixedsort(rownames(plot.data), decreasing = T), ]
    if (is.null(ncol(plot.data))) {
        stop("No enough data. Genes for at least two chains out of V/D/J is required")
    }

    chordDiagram(plot.data, group = group, grid.col = grid.col, annotationTrack = "grid",
                 preAllocateTracks = list(track.height = plot.data %>% dimnames() %>% unlist() %>% strwidth() %>% max()/3))
    circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
    }, bg.border = NA)
}

#' Plot length of CDR3 AA/NT sequences
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by.
#' @param subset Subset data to these groups
#' @param plot.type Type of plot. Options = ridge, line
#' @param sequence.type Which CDR3 sequence should be used? Options: "AA" or "NT"
#' @param color.theme Color theme to use. Default = DALI
#' @param colors Colors to use. This overwrites the selected color theme
#'
#' @export

CDR3Plot <- function(
    object,
    group.by = NULL,
    subset = NULL,
    plot.type = c("ridge", "line"),
    sequence.type = c("AA", "NT"),
    color.theme = ColorThemes(),
    colors = NULL
) {
    plot.type <- match.arg(plot.type)
    sequence.type <- match.arg(sequence.type)
    color.theme <- match.arg(color.theme)

    if (is.null(group.by)) {
        object <- Seurat::AddMetaData(object, Seurat::Idents(object), "default.clustering")
        group.by <- "default.clustering"
    }

    if (!group.by %in% colnames(object@meta.data)) {
        stop("Invalid group.by column ", group.by)
    }

    if (!is.null(subset) && subset != "") {
        cells <- rownames(object@meta.data)[object@meta.data[[group.by]] %in% subset]
        object <- subset(object, cells = cells)
    }

    cdr3.sequence <- if (sequence.type == "AA") ".cdr3" else ".cdr3_nt"
    vdj.cdr3.column <- paste0("vdj", cdr3.sequence)
    vj.cdr3.column <- paste0("vj", cdr3.sequence)

    if (plot.type == "line") {
      plots <- CDR3Plot.line(object, group.by, vdj.cdr3.column, vj.cdr3.column, sequence.type)
    } else if (plot.type == "ridge") {
      plots <- CDR3Plot.ridge(object, group.by, vdj.cdr3.column, vj.cdr3.column, color.theme = color.theme, colors = colors)
    }

    gridExtra::grid.arrange(grobs = plots, ncol = min(length(plots), 3))
}

#' Plot length of CDR3 aa sequences as a line plot
#'
#' @param object Seurat object
#' @param group.by Metadata column to group the family data by.
#' @param vdj.cdr3.column Column name for heavy/alpha cdr3 data
#' @param vj.cdr3.column Column name for light/beta cdr3 data
#' @param sequence.type AA or NT
#'
#' @importFrom dplyr %>% mutate select full_join count
#' @importFrom ggplot2 ggplot aes geom_line labs theme element_rect element_line element_blank unit element_text
#' @importFrom rlang .data
#' @importFrom stats na.omit

CDR3Plot.line <- function(object, group.by, vdj.cdr3.column, vj.cdr3.column, sequence.type) {
  plots <- list()
  groups <- unique(object@meta.data[[group.by]]) %>% gtools::mixedsort(x = .)

  for (group in groups) {
    cells <- rownames(object@meta.data)[object@meta.data[[group.by]] == group]
    subset <- subset(object, cells = cells)

    plot.data.h <- subset@meta.data %>%
      mutate(len = nchar(.data[[vdj.cdr3.column]])) %>%
      count(.data$len) %>%
      na.omit() %>%
      mutate(freq = .data$n / sum(.data$n) * 100) %>%
      select(.data$len, .data$freq)
    plot.data.l <- subset@meta.data %>%
      mutate(len = nchar(.data[[vj.cdr3.column]])) %>%
      count(.data$len) %>%
      na.omit() %>%
      mutate(freq = .data$n / sum(.data$n) * 100) %>%
      select(.data$len, .data$freq)

    plot.data <- full_join(plot.data.h, plot.data.l, by = "len") %>% replace(is.na(.), 0)
    colnames(plot.data) <- c("cdr3.length", "vdj.chain", "vj.chain")

    plots[[group]] <- ggplot(plot.data, aes(x = .data$cdr3.length)) +
      geom_line(aes(y = .data$vdj.chain), color = "black") +
      geom_line(aes(y = .data$vj.chain), color = "red") +
      labs(x = paste0("CDR3 length (", sequence.type, ")"), y = "Frequency of cells", title = paste0("CDR3 length - ", group)) +
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
#' @param group.by Metadata column to group the family data by.
#' @param vdj.cdr3.column Column name for heavy/alpha cdr3 data
#' @param vj.cdr3.column Column name for light/beta cdr3 data
#' @param color.theme Color theme to use. Default = DALI
#' @param colors Colors to use. This overwrites the selected color theme
#'
#' @importFrom dplyr %>% mutate
#' @importFrom ggplot2 ggplot aes coord_cartesian scale_y_discrete scale_x_continuous
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom rlang .data

CDR3Plot.ridge <- function(object, group.by, vdj.cdr3.column, vj.cdr3.column, color.theme, colors = NULL) {
  plots <- list()

  for (column in c(vdj.cdr3.column, vj.cdr3.column)) {

    plot.data <- object@meta.data %>%
      mutate(len = nchar(.data[[column]]))

    if (column == vdj.cdr3.column) {
      if (DefaultAssayVDJ(object) == "TCR") {
        title <- "TCR-a"
      } else {
        title <- "Heavy"
      }
    } else {
      if (DefaultAssayVDJ(object) == "TCR") {
        title <- "TCR-b"
      } else {
        title <- "Light"
      }
    }

    title <- paste0(title, " chain")

    if (is.null(colors)) {
        colors <- GetCategoricalColorPalette(plot.data[[group.by]], color.theme)
    }

    plots[[column]] <- ggplot(plot.data, aes(x = .data$len, y = .data[[group.by]], fill = .data[[group.by]])) +
      geom_density_ridges() +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = colors) +
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
#' @param color.theme Color theme to use. Default = DALI
#' @param colors Colors to use. This overwrites the selected color theme
#' @param ... Extra parameters passed to Seurat::Dimplot
#'
#' @importFrom dplyr %>%
#'
#' @export

DimplotChainRegion <- function(
  object,
  region = c("V", "D", "J", "C"),
  chain = c("VDJ", "VJ"),
  by.family = T,
  grid = T,
  highlight = NULL,
  color.theme = ColorThemes(),
  colors = NULL,
  ...
) {

  region <- match.arg(region) %>% tolower()
  chain <- match.arg(chain) %>% tolower()
  color.theme <- match.arg(color.theme)

  if (length(highlight) > 1) {
    stop("Can only select 1 family/gene to highlight", call. = F)
  }

  data.column <- paste0(chain, ".", region, "_")

  if (by.family && region == "v") {
    data.column <- paste0(data.column, "fam")
  } else {
    data.column <- paste0(data.column, "gene")
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

  cells.highlight <- NULL
  if (!is.null(highlight)) {
    cells.highlight <- rownames(object@meta.data)[object@meta.data[[data.column]] %in% highlight]
  }

  if (is.null(colors)) {
    colors <- GetCategoricalColorPalette(object@meta.data[[data.column]], color.theme)
  }

  Seurat::DimPlot(object, group.by = data.column, split.by = split, cells.highlight = cells.highlight, order = rev(families), cols = colors, ...) +
    theme(
      panel.grid.major = element_line(color = "grey"),
      axis.line = element_blank()
    )
}


#' Plot the frequency of clonotypes in (subset of) clusters
#'
#' @param object Seurat object
#' @param chain Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR. NULL = both
#' @param group.by Metadata column to group the family data by.
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
  chain = c("VDJ", "VJ"),
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
    object <- Seurat::AddMetaData(object, Seurat::Idents(object), "default.clustering")
    group.by <- "default.clustering"
  }

  if (is.null(clonotype.column)) {
    clonotype.column <- "clonotype"
  }

  if (clonotype.column == group.by) {
      stop("Cannot group the data by the clonotype column")
  }

  plot.type <- match.arg(plot.type)

  if (!use.sequence && !clonotype.column %in% colnames(object@meta.data)) {
    stop("Invalid clonotype column ", clonotype.column, call. = F)
  }

  if (!is.null(subset) && subset != "") {
    cells <- rownames(object@meta.data)[object@meta.data[[group.by]] %in% subset]
    object <- subset(object, cells = cells)
  }

  sequence.type <- match.arg(sequence.type)

  chains <- chain

  if (is.null(chains)) {
    chains <- c("vdj", "vj")
  }

  if (plot.type == "bar") {
    plots <- ClonotypeFrequency.bar(object, chains, group.by, use.sequence, sequence.type, clonotype.column, show.missing, bulk, threshold)
  } else if (plot.type == "violin") {
    plots <- ClonotypeFrequency.violin(object, chains, group.by, use.sequence, sequence.type, clonotype.column, show.missing, bulk)
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
        legend.position = "none"
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
  clonotype.column = "clonotype",
  color.low = "lightgrey",
  color.mid = "red",
  color.high = "darkred",
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

  coordinates <- object@reductions[[reduction]]@cell.embeddings[, c(1, 2)] %>% as.data.frame()
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
  plot.data[rownames(data), "clonotype_count"] <- data$n %>% as.numeric()

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

FeaturePlotClonotype <- function(
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

  coordinates <- object@reductions[[reduction]]@cell.embeddings[, c(1, 2)] %>% as.data.frame()
  x.name <- colnames(coordinates)[[1]]
  y.name <- colnames(coordinates)[[2]]

  object@meta.data <- object@meta.data %>%
    mutate(clonotypes = ifelse(.data[[clonotype.column]] %in% clonotypes, .data[[clonotype.column]], NA))

  plot.data <- coordinates
  plot.data$clonotypes <- object@meta.data[rownames(plot.data), "clonotypes"]

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
#' @param group.by  Metadata column to group the family data by.
#' @param groups.highlight Groups for which to highlight the edges in the graph
#' @param clonotype.column Metadata column with clonotype information. Default = clonotype
#' @param color.theme Color theme to use. Default = DALI
#' @param colors Colors to use. This overwrites the selected color theme
#'
#' @importFrom dplyr %>% bind_rows distinct filter group_by group_map n
#' @importFrom ggplot2 aes element_blank geom_point scale_color_manual theme
#' @importFrom ggraph geom_edge_arc geom_node_point ggraph scale_edge_alpha scale_edge_color_manual scale_edge_width
#' @importFrom rlang .data
#' @importFrom tidygraph as_tbl_graph
#'
#' @export

CloneConnGraph <- function(
    object,
    reduction,
    group.by = NULL,
    groups.highlight = NULL,
    clonotype.column = NULL,
    color.theme = ColorThemes(),
    colors = NULL
) {
  color.theme <- match.arg(color.theme)

  if (is.null(group.by)) {
    object <- Seurat::AddMetaData(object, Seurat::Idents(object), "default.clustering")
    group.by <- "default.clustering"
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
    group_by(.data[[clonotype.column]]) %>%
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
    edges[edges$from %in% groups.highlight, "highlight"] <- T
    edges[edges$to %in% groups.highlight, "highlight"] <- T
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

  if (is.null(colors)) {
      colors <- GetCategoricalColorPalette(dimred$group, color.theme)
  }

  plot <- ggraph(graph, layout = "manual", x = .data$x, y = .data$y) +
    geom_point(data = dimred, aes(x = .data[[label.x.axis]], y = .data[[label.y.axis]], color = .data$group), size = 1) +
    scale_color_manual(values = colors) +
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

#' Use Seurat's DEG Table to plot a volcano plot
#'
#' @param deg Table with DEG data from Seurat::Findmarkers()
#' @param sig.P Significant P value to extract over/under expressed genes from
#' @param sig.logFC significant Log2(FC) value to extract over/underexpressed genes from
#' @param color.scheme Color scheme to use. Options: "coolwarm", "viridis". Default = "coolwarm"
#'
#' @importFrom ggplot2 geom_point ggplot labs scale_color_manual theme_minimal ylim
#' @importFrom ggrepel geom_text_repel
#'
#' @export

VolcanoPlotDEG <- function(deg, sig.P = 0.05, sig.logFC = 0.6, color.scheme = c("coolwarm", "viridis", "gray to blue", "turning red")) {
    color.scheme <- match.arg(color.scheme)

    # For coloring: adding a column to indicate diff expression
    deg$expression <- "zero"
    deg$expression[deg$avg_log2FC > sig.logFC & deg$p_val_adj < sig.P] <- "up"
    deg$expression[deg$avg_log2FC < -sig.logFC & deg$p_val_adj < sig.P] <- "down"

    if (color.scheme == "coolwarm") {
        colors <- c("red", "blue", "black")
    } else if (color.scheme == "viridis") {
        colors <-  c("#FDE725", "#423C81", "#249F87")
    } else if (color.scheme == "gray to blue") {
        colors <- c("#c0c0c0", "#c0c0c0", "#4575B4")
    } else if (color.scheme == "turning red") {
        colors <- c("#c0c0c0", "#c0c0c0", "#D73027")
    }

    names(colors) <- c("up", "down", "zero")

    #get names of diff expressed genes
    deg$genes <- NA
    deg$genes[deg$expression != "zero"] <- rownames(deg[deg$expression != "zero",])

    ymax <- -log10(deg$p_val_adj[deg$p_val_adj > 0] %>% min()) + 20

    ggplot(data = deg, aes(x = .data$avg_log2FC, y = -log10(.data$p_val_adj), col = expression, label = .data$genes)) +
       geom_point(show.legend = F) +
       geom_text_repel(show.legend = F) +
       theme_minimal() +
       ylim(-5,ymax) +
       scale_color_manual(values = colors) +
       labs(x = expression(Log[2](FC)), y = expression(-log[10](Adj-P-value)))
}

#' Plots a Dimensionality reduction plot with cells colored by mutation rate
#' @param object Seurat object
#' @param  clonotype clonotype to plot its mutation rate. Default plots every clonotype
#' @param regions VDJ region to use sequence from.
#' @param reference path to reference fasta sequence file
#' @param chain Which chain to use
#' @param  clonotype.column name of the metadata column that contains the clonotypes
#' @param order Plot positive cells in the foreground. Default = True
#'
#' @importFrom Seurat AddMetaData FeaturePlot
#'
#' @export

PlotSHM <- function(
        object,
        clonotype = NULL,
        regions = "V",
        reference,
        chain = c("VDJ", "VJ"),
        clonotype.column = "clonotype",
        order = T
) {
    chain <- match.arg(chain)

    distances <- GetMutationRate(
        object = object,
        clonotype = clonotype,
        regions = regions,
        reference = reference,
        chain = chain,
        clonotype.column = clonotype.column
    )

    object <- AddMetaData(object, distances, "mutationrate")

    return(FeaturePlot(object = object, features = "mutationrate", cols = ColorScale("turning red", n = 2), order = order, pt.size = 1))
}

#' Function to create a Trajectory analyses from from BCR or VDJ cell data using the dynoplot package
#'
#' @param object Seurat Object
#' @param method Method to calculate trajectory analyses. See https://dynverse.org/reference/dynmethods/ for options
#' @param reduction Dimensionality reduction for the plot
#' @param assay BCR or TCR assay. If not specified it will take all the data.
#' @param color.by Meta.data column to color the datapoints by. If empty a gradient will be used.
#' @param min.occurence Minimum occurence of a clonotype for it to be indicated on the graph. Default values = 5.
#' @param color.scheme Colorscheme of the colorgradient used to color by pseudotime when no color.by argument is provided.
#' @param size_cells Size of the plotted dots. Default size = 3.
#' @param alpha_cells Transparancy of the plotted dots. Default alpha = 0.5.
#'
#' @importFrom ggplot2 scale_colour_gradient2
#' @importFrom utils installed.packages
#'
#' @export

TrajectoryPlot <- function(
        object,
        method,
        reduction,
        assay = NULL,
        min.occurence = 5,
        color.by = NULL,
        color.scheme = c("coolwarm", "viridis", "gray to blue", "turning red"),
        size_cells = 3,
        alpha_cells = 0.5)
{
    # Check if dynplot & dynwrap are installed
    if (sum(c("dynverse") %in% installed.packages()) > 0) {
        stop("Package dynverse needs to be installed to use this function")
    }

    color.scheme <- match.arg(color.scheme)

    # Keep only specified assay data if assay is specified
    if (!is.null(assay) && assay %in% c('TCR', 'BCR')) {
        # get cells that are TCR/BCR cells if a subset is asked
        object <- SubsetAssay(object = object, assay = assay)
    } else if (!is.null(assay)){
        stop("Unrecognized assay: ", assay)
    }

    # Check if the dimensionality reduction i spresent
    if (is.null(object@reductions[[reduction]])) {
        stop("Dimensionality reduction ", reduction," was not performed for the Seurat object")
    }
    dimred <- object@reductions[[reduction]]@cell.embeddings

    if (!is.null(color.by)) {
        # don't color clonotypes with hardly any representation (change to NA)
        grouping <- object[[color.by]]
        table_groups <- table(grouping)

        for (i in 1:nrow(grouping)) {
            if (is.na(grouping[i, 1])) {
                next
            }
            if (table_groups[grouping[i, 1]] < min.occurence ) {
                grouping[i, 1] <- NA
            }
        }
        grouping <- as.character(grouping[,1])
        names(grouping) <- object[[color.by]] %>% rownames()
    }

    dataset <- dynwrap::wrap_expression(
        counts = t(as.matrix(object@assays$RNA@counts)),
        expression = t(as.matrix(object@assays$RNA@data))
    )

    # Check if method is an existing dyno method
    model <- tryCatch(
        dynwrap::infer_trajectory(dataset, method)
    , error = function(e) {return(NA)})
    if (is.na(model)[1]) {
        stop("Method is not recognized. Choose a method from: \n", paste0("ti_", dynwrap::get_ti_methods()$id, collapse = "(), "), "." )
    }

    if (is.null(color.by)) {
        # colorscheme when no color.by specified; the default option
        spectrum <- ColorScale(color.scheme, n = 3)
        pseudotime <-  dynwrap::calculate_pseudotime(model)
        midpoint <- max(pseudotime)/2

        trajectoryplot <- dynplot::plot_dimred(model, "pseudotime",
                                      pseudotime = pseudotime,
                                      dimred = dimred, alpha_cells = alpha_cells,
                                      size_cells = size_cells) +
            scale_colour_gradient2(low = spectrum[1], mid = spectrum[2], high = spectrum[3], midpoint = midpoint)

    } else {
        trajectoryplot <- dynplot::plot_dimred(model, grouping = grouping, dimred = dimred, alpha_cells = alpha_cells, size_cells = size_cells)
    }

    return(trajectoryplot)
}

#' Plotfunction to plot a simple Trajectory Analysis plot for the DALI shiny app
#'
#' @param object Seurat object
#' @param reduction Dimensionality reduction as basis for the analysis and plot.
#' @param gene Gene to base trajectory analysis on. Default = the first in the expression matrix
#' @param start_cluster_id ID of the cluster to root the trajectory tree on
#' @param color.theme Colorscheme for the plot.
#'
#' @importFrom slingshot slingshot
#' @importFrom dplyr first
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom graphics legend lines

ShinyTrajectoryPlot <- function(object, reduction, gene, start_cluster_id, color.theme = "DALI") {
    # Check if the reduction is actually present
    if (is.null(object@reductions[[reduction]])) {
        stop(reduction, " reduction was not performed on the object")
    }

    set.seed(1)

    dimred <- object@reductions[[reduction]]@cell.embeddings
    clustering <- object@meta.data$seurat_clusters

    if (!is.null(gene)) {
        # Check if the gene is present in the expression matrix
        if (is.null(object@assays$RNA@counts[gene, ])) {
            stop("Gene ", gene, " was not detected in the expression matrix")
        }
        start_cluster_id <- tibble(
            expression = object@assays$RNA@counts[gene,],
            cluster_id = clustering
        ) %>%
            group_by(cluster_id) %>%
            summarise(expression = mean(expression)) %>%
            arrange(desc(expression)) %>%
            pull(cluster_id) %>%
            first() %>%
            as.character()
    } else if (is.null(start_cluster_id)) {
        stop("Please Specify a Start Cluster ID or Gene to base the Trajectory tree")
    } else if (!(start_cluster_id %in% levels(clustering))) {
        stop( start_cluster_id," is not a valid cluster ID")
    }

    lineage <- slingshot(dimred, clustering, start.clus = start_cluster_id, allow.breaks = T)

    colors <- GetCategoricalColorPalette(levels(object@meta.data$seurat_clusters))

    plot(dimred, col = colors[clustering], asp = 1, pch = 16)
    for (lineage in lineage@metadata$curves) {
        lines(lineage, lwd = 3, col = 'black')
    }
}
