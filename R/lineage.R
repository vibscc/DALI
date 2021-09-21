#' Create a phylogenetic tree of the v-gene sequence of either the heavy or light chain
#'
#' @param object Seurat object
#' @param clonotype.name Clonotype to plot
#' @param clonotype.column Metadata column with clonotype information. Default = clonotype
#' @param airr Path to airr_rearrangement.tsv. This file can be found in cellranger multi output (> v6.1.1)
#' @param reference Path to reference fasta. This file can either be found in the output of cellranger multi or can be the input reference file on which cellranger multi was ran.
#' @param chain Which chain to use?
#'
#' @importFrom dplyr %>% filter pull
#'
#' @export

LineageTreeVGene <- function(object, clonotype.name, clonotype.column = NULL, airr, reference, chain = Diversity:::AvailableChains(object)) {
    chain <- match.arg(chain) %>% tolower()

    if (is.null(clonotype.column)) {
        clonotype.column <- "clonotype"
    }

    data <- object@meta.data %>% filter(!is.na(.data[[clonotype.column]]) & .data[[clonotype.column]] == clonotype.name)
    cells <- rownames(data)

    if (nrow(data) < 2) {
        stop("At least 2 cells need to be defined as given clonotype!", call. = F)
    }

    v_call <- data %>% pull(paste0(chain, ".v_gene")) %>% unique()

    if (length(v_call) > 1) {
        stop("Found more than 1 v_gene for given clonotype.", call. = F)
    }

    c_call.regex <- if (chain == "h") "IGH" else "IG[KL]"
    data.airr <- read.csv(airr, sep = "\t") %>% filter(cell_id %in% cells & grepl(c_call.regex, c_call))

    sequences <- GetGermline(reference, v_call = v_call)
    for (i in seq(1, nrow(data.airr))) {
        v_start <- data.airr[i, "v_sequence_start"]
        v_end <- data.airr[i, "v_sequence_end"]
        sequences <- c(sequences, substr(data.airr[i, "sequence_alignment"], v_start, v_end))
    }

    germline.name <- paste0("germline (", v_call, ")")
    distance.mat <- stringdist::stringdistmatrix(sequences, sequences, method = "lv")
    colnames(distance.mat) <- c(germline.name, data.airr$cell_id)
    rownames(distance.mat) <- c(germline.name, data.airr$cell_id)

    tree <- ape::nj(distance.mat) %>% ape::root(germline.name)
    plot(tree)
}

#' Extract the germline sequence from a reference fasta with given v, d and/or j call
#'
#' @param reference Path to reference fasta.
#' @param v_call Call for the v-gene. Default = NULL
#' @param d_call Call for the d-gene. Default = NULL
#' @param j_call Call for the j-gene. Default = NULL

GetGermline <- function(reference, v_call = NULL, d_call = NULL, j_call = NULL) {
    calls <- c(v_call, d_call, j_call)

    if (is.null(calls)) {
        stop("No genes specified! Specify either v_call, d_call, j_call or a combination.", call. = F)
    }

    reference.sequences <- Read10XVDJReference(reference)

    germline <- ""

    for (call in calls) {
        if (!call %in% names(reference.sequences)) {
            stop("Could not find ", call, " in given reference", call. = F)
        }

        germline <- paste0(germline, reference.sequences[[call]])
    }

    return(germline)
}

#' Read 10x VDJ reference FASTA and extract sequences for each gene
#'
#' @param fasta Path to fasta reference

Read10XVDJReference <- function(fasta) {
    con <- file(fasta, "r")

    sequences <- list()

    while (T) {
        line <- readLines(con, n = 1)
        if (length(line) == 0) {
            break
        }

        if (grepl("^>", line)) {
            header <- Parse10XVDJReferenceHeader(line)
        } else {
            if (!is.null(description)) {
                if (!header[["gene_name"]] %in% names(sequences)) {
                    sequences[[header[["gene_name"]]]] <- line
                } else {
                    sequences[[header[["gene_name"]]]] <- paste0(sequences[[header[["gene_name"]]]], line)
                }
            }
        }
    }
    close(con)

    return(sequences)
}

#' Extract information from FASTA header of 10x VDJ reference
#'
#' @param header Fasta header string

Parse10XVDJReferenceHeader <- function(header) {
    data <- list()

    parts <- strsplit(substr(description, 2, nchar(description)), "\\|") %>% unlist()

    if (!grepl("REGION", parts[[4]], ignore.case = T)) {
        return(NULL)
    }

    field2 <- strsplit(parts[[2]], " ") %>% unlist()
    data[["display_name"]] <- field2[[1]]
    data[["gene_name"]] <- parts[[3]]
    data[["region_type"]] <- parts[[4]]
    data[["chain_type"]] <- parts[[5]]
    data[["chain"]] <- parts[[6]]

    return(data)
}
