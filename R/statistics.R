#' Calculate diversity in VDJ data
#'
#' @param object Seurat object
#' @param algorithm Algorithm used for diversity calculation. Options = "gini" or "shannon"
#' @param use.sequence Use AA/NT sequence instead of clonotype. Default = FALSE
#' @param sequence.type What sequences to use, available options: "AA" or "NT". Only functional when use.sequence = T
#' @param chain  Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR
#' @param clonotype.column Metadata column with clonotype information. Default = 'clonotype'
#'
#' @importFrom dplyr %>% group_by n summarise
#'
#' @export

ClonotypeDiversity <- function(
    object,
    algorithm = c("gini", "shannon"),
    use.sequence = F,
    sequence.type = c("AA", "NT"),
    chain = AvailableChains(object),
    clonotype.column = NULL
) {
    algorithm <- match.arg(algorithm)
    chain <- match.arg(chain) %>% tolower()
    sequence.type <- match.arg(sequence.type)

    if (is.null(clonotype.column)) {
        clonotype.column <- "clonotype"
    }

    if (!use.sequence && !clonotype.column %in% colnames(object@meta.data)) {
        stop("Invalid clonotype column ", clonotype.column, call. = F)
    }

    data.column <- clonotype.column

    if (use.sequence) {
        data.column <- paste0(chain, ".cdr3")

        if (sequence.type == "NT") {
            data.column <- paste0(data.column, "_nt")
        }
    }

    frequency.table <- object@meta.data %>%
        group_by(.data[[data.column]]) %>%
        summarise(freq = n()) %>%
        na.omit()

    if (algorithm == "gini") {
        return(CalculateDiversity.gini(frequency.table$freq))
    } else if (algorithm == "shannon") {
        return(CalculateDiversity.shannon(frequency.table$freq))
    } else {
        stop("Invalid algorithm ", algorithm)
    }
}

#' Calculate diversity using the gini index
#'
#' @param frequencies Vector of frequencies

CalculateDiversity.gini <- function(frequencies) {
    return(reldist::gini(frequencies, weights = rep(1, length = length(frequencies))))
}

#' Calculate diversity using the shannon index
#'
#' @param frequencies Vector of frequencies

CalculateDiversity.shannon <- function(frequencies) {
    return(vegan::diversity(frequencies))
}
