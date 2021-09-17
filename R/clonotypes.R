#' Define clonotypes
#'
#' @param object Seurat object
#' @param sequence Use the amino acis or nucleotide sequence. Default = aa
#' @param metric Distance metric to use. Default = identity
#' @param column Column to add to metadata. Default = clonotype.<metric>
#'
#' @importFrom dplyr %>% mutate na_if pull
#' @export

DefineClonotypes <- function(object, sequence = c("aa", "nt"), metric = c("identity"), column = NULL) {

    if (!IsValidSeuratObject(object)) {
        stop("Invalid Seurat object", call. = F)
    }

    sequence <- match.arg(sequence)
    metric <- match.arg(metric)

    if (is.null(column)) {
        column <- paste0("clonotype.", metric)
    }

    cdr3.column <- if (sequence == "aa") "cdr3" else "cdr3_nt"
    heavy.column <- paste0("h.", cdr3.column)
    light.column <- paste0("l.", cdr3.column)

    object@meta.data <- object@meta.data %>% mutate(tmp_cdr3_concat = paste0(.data[[heavy.column]], .data[[light.column]]))
    object@meta.data$tmp_cdr3_concat <- gsub("NA", "_", object@meta.data$tmp_cdr3_concat) %>% na_if(y = "__")

    sequences <- object@meta.data %>%
        pull(tmp_cdr3_concat) %>%
        na.omit() %>%
        unique()

    dist.mat <- CalculateDistances(sequences, metric)
    adj.mat <- DistToAdjacency(dist.mat)

    g <- igraph::graph_from_adjacency_matrix(adj.mat, "undirected")
    communities <- igraph::cluster_louvain(g)

    object@meta.data[[column]] <- NA

    i <- 1
    for (membership in communities$membership) {
        sequence <- sequences[i]

        object@meta.data[!is.na(object@meta.data$tmp_cdr3_concat) & object@meta.data$tmp_cdr3_concat == sequence, column] <- paste0("clonotype", membership)

        i <- i + 1
    }

    object@meta.data$tmp_cdr3_concat <- NULL

    return(object)
}

#' Convert a sparse distance matrix to an adjacency matrix
#'
#' @param dist.mat dgRMatrix with distance info

DistToAdjacency <- function(dist.mat) {
    max.val <- max(slot(dist.mat, "x"))
    slot(dist.mat, "x") <- slot(dist.mat, "x") - 1

    slot(dist.mat, "x") <- (max.val - slot(dist.mat, "x")) / max.val

    return(dist.mat)
}

#' Calculate distances between sequences with the given distance metric
#'
#' @param sequences Vector of sequences
#' @param metric Distance metric to use

CalculateDistances <- function(sequences, metric) {
    if (metric == "identity") {
        return(CalculateDistances.identity(sequences))
    } else {
        stop("Invalid distance metric ", metric, call. = F)
    }
}

#' Create distance matrix where identical sequences get following distances:
#' identical: 0
#' non-identical: 1
CalculateDistances.identity <- function(sequences) {
    # Distances are stored in a sparse matrix. therefor 1 should be added to each distance so we also keep 0 distances.
    # For the identity metric, this results in an identity matrix
    dist <- diag(length(sequences)) %>% as("dgRMatrix")

    return(dist)
}
