#' Calculate frequencies of column and group
#'
#' @param object Seurat object
#' @param data.column Metadata column to use
#' @param group.by Metadata column with group info
#' @param show.missing Should missing values be included or not

CalculateFrequency <- function(object, data.column, group.by, show.missing) {
    out <- object@meta.data %>%
        group_by(.data[[data.column]], .data[[group.by]]) %>%
        summarise(freq = n()) %>%
        ungroup() %>%
        filter(case_when(!show.missing ~ !is.na(.data[[data.column]]),
                         T ~ T))
    return(out)
}

#' Convert clonotype to sequences
#'
#' @param object Seurat object
#' @param ct clonotype to convert
#' @param chain Which chain to get sequence for
#' @param clonotype.column Metadata column with clonotype information. Default = clonotype
#' @param sequence.type One of AA or NT. Default = AA
#' @param collapse Character to use to paste multiple sequence together

ClonotypeToSequence <- function(object, ct, chain, clonotype.column = 'clonotype', sequence.type = 'AA', collapse = "\n") {
    sequence.column <- paste0(tolower(chain), ".cdr3", if (tolower(sequence.type) == 'nt') '_nt' else '')

    if (!sequence.column %in% colnames(object@meta.data)) {
        stop("Invalid metadata column ", sequence.column, call. = F)
    }

    data.filtered <- object@meta.data %>%
        filter(.data[[clonotype.column]] == ct)

    sequences <- data.filtered[, sequence.column] %>% unique() %>% na.omit()

    return(paste(sequences, collapse = collapse))
}

#' Get field name for given column name
#'
#' @param column Column for which you want the field
#' @param fields List of available fields
#' @param columns List of available columns

FieldForColumn <- function(column, fields, columns) {
    index = which(columns == column[[1]])

    return(fields[[index]])
}

#' Generate ggplot default colors
#'
#' @param n Number of colors needed
#' @param hue Minimum and maximum hue
#'
#' @importFrom grDevices hcl

ggplotColors <- function(n, h = c(15, 375)) {
    if ((diff(h) %% 360) < 1) {
        h[2] <- h[2] - 360/n
    }
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

#' Get a colorpalette for given categorical data
#'
#' @param data Categorical data
#'
#' @importFrom dplyr %>%
#' @importFrom Polychrome sky.colors

GetCategoricalColorPalette <- function(data) {
    n.categories <- unique(data) %>% length()

    if (n.categories <= 24) {
        return(sky.colors(max(c(3, n.categories))) %>% unname())
    }

    return(ggplotColors(n = n.categories))
}
