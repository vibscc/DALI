#' Build a phylogenetic tree based on v,d,j and/or c gene
#'
#' @param object Seurat object
#' @param clonotype Clonotype to plot
#' @param regions VDJ regions to use sequence for. Default = V
#' @param reference Path to reference fasta. This file can either be found in the output of cellranger multi or can be the input reference file on which cellranger multi was ran.
#' @param chain Which chain to use?
#' @param clonotype.column Metadata column with clonotype information. Default = clonotype
#'
#' @importFrom dplyr %>% filter pull
#'
#' @export

LineageTree <- function(object, clonotype, regions = "V", reference, chain = c("VDJ", "VJ"), clonotype.column = NULL) {
    if (DefaultAssayVDJ(object) != "BCR") {
        stop("Lineages can only be computed on BCR data!", call. = F)
    }

    chain <- match.arg(chain) %>% tolower()
    clonotype.name <- clonotype

    if (is.null(clonotype.column)) {
        clonotype.column <- "clonotype"
    }

    data <- object@meta.data %>% filter(!is.na(.data[[clonotype.column]]) & .data[[clonotype.column]] == clonotype.name)
    cells <- rownames(data)

    if (nrow(data) < 2) {
        stop("Clonotype needs at least 2 cells to build a tree!", call. = F)
    }

    calls <- list("v_call" = NULL, "d_call" = NULL, "j_call" = NULL, "c_call" = NULL)

    for (region in regions) {
        if (!region %in% AvailableRegions(chain)) {
            stop("Invalid region ", region, call. = F)
        }
        region <- tolower(region)
        region.call <- data %>% pull(paste0(chain, ".", region, "_gene")) %>% unique()

        if (is.na(region.call)) {
            next
        }

        if (length(region.call) > 1) {
            stop("Found more than 1 ", region, "_gene for given clonotype.", call. = F)
        }

        calls[[paste0(region, "_call")]] <- region.call
    }

    sequences <- GetGermline(reference, v_call = calls[["v_call"]], d_call = calls[["d_call"]], j_call = calls[["j_call"]], c_call = calls[["c_call"]])

    for (cell in cells) {
        sequences <- c(sequences, GetSequence(object, cell, "BCR", chain = chain %>% toupper(), regions = regions))
    }

    germline.name <- paste0("germline (", unname(unlist(calls)) %>% paste(collapse = ";"), ")")
    names(sequences) <- c(germline.name, cells)

    tree <- BuildLineageTree(sequences, outgroup = germline.name, root = T)
    plot(tree)
}

#' Build a lineage tree for the given sequences
#'
#' @param sequences Named vector of sequences
#' @param distance.method Method to use to compute distances between the sequences. Default = lv (Levenshtein)
#' @param root Should the tree be rooted. Default = FALSE
#' @param outgroup Outgroup to root the tree

BuildLineageTree <- function(sequences, distance.method = "lv", outgroup = NULL, root = F) {
    if (root & is.null(outgroup)) {
        stop("Missing outgroup to root the tree.", call. = F)
    }

    if (length(sequences) <= 2) {
        stop("Need at least 3 sequences to build a tree.", call. = F)
    }

    if (is.null(names(sequences))) {
        stop("Sequence vector should be named", call. = F)
    }

    distance.mat <- stringdist::stringdistmatrix(sequences, sequences, method = distance.method)
    colnames(distance.mat) <- names(sequences)
    rownames(distance.mat) <- names(sequences)

    tree <- ape::nj(distance.mat)

    if (root) {
        tree <- ape::root(tree, outgroup = outgroup)
    }

    return(tree)
}

#' Extract the germline sequence from a reference fasta with given v, d and/or j call
#'
#' @param reference Path to reference fasta.
#' @param v_call Call for the v-gene. Default = NULL
#' @param d_call Call for the d-gene. Default = NULL
#' @param j_call Call for the j-gene. Default = NULL
#' @param c_call Call for the c-gene. Default = NULL

GetGermline <- function(reference, v_call = NULL, d_call = NULL, j_call = NULL, c_call = NULL) {
    calls <- c(v_call, d_call, j_call, c_call)

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
            if (!is.null(header)) {
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

    parts <- strsplit(substr(header, 2, nchar(header)), "\\|") %>% unlist()

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
