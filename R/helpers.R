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

ClonotypeToSequence <- function(object, ct, chain, clonotype.column = 'clonotype', sequence.type = 'AA') {
    sequence.column <- paste0(tolower(chain), ".cdr3", if (tolower(sequence.type) == 'nt') '_nt' else '')

    if (!sequence.column %in% colnames(object@meta.data)) {
        stop("Invalid metadata column ", sequence.column, call. = F)
    }

    data.filtered <- object@meta.data %>%
        filter(.data[[clonotype.column]] == ct)

    sequences <- data.filtered[, sequence.column] %>% unique() %>% na.omit()

    return(paste(sequences, collapse = " | "))
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
