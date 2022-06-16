#' Build a phylogenetic tree based on v,d,j and/or c gene
#'
#' @param object Seurat object
#' @param clonotype Clonotype to plot
#' @param regions VDJ regions to use sequence for. Default = V
#' @param reference Path to reference fasta. This file can either be found in the output of cellranger multi or can be the input reference file on which cellranger multi was ran.
#' @param chain Which chain to use?
#' @param color.tip.by Feature or metadata column to color tips of the tree by
#' @param clonotype.column Metadata column with clonotype information. Default = clonotype
#'
#' @importFrom dplyr %>% filter pull
#' @importFrom ggplot2 aes scale_color_gradient2
#' @importFrom ggtree geom_tiplab geom_tippoint ggtree
#' @importFrom scales muted
#'
#' @export

LineageTree <- function(object, clonotype, regions = "V", reference, chain = c("VDJ", "VJ"), color.tip.by = NULL, clonotype.column = "clonotype") {
    if (DefaultAssayVDJ(object) != "BCR") {
        print("Lineages can only be computed on BCR data!", call. = F)
    }
    clonotype.name <- clonotype

    chain <- match.arg(chain)

    metadata <- object@meta.data %>% filter(!is.na(.data[[clonotype.column]]) && .data[[clonotype.column]] == clonotype.name)
    sequences <- GetSequences(object = object, metadata = metadata, regions = regions, reference = reference, chain = chain )
    if (length(sequences) < 2) {
        stop("Clonotype needs at least 2 cells to build a tree!", call. = F)
    }

    tree <- BuildLineageTree(sequences, outgroup = names(sequences)[1], root = T)

    color.tip <- F
    continues.scale <- F

    if (!is.null(color.tip.by)) {
        if (color.tip.by %in% rownames(object)) {
            feature <- color.tip.by
            metadata <- NULL
            continues.scale <- T
        } else if (color.tip.by %in% colnames(object@meta.data)) {
            feature <- NULL
            metadata <- color.tip.by
        } else {
            stop("Could not find ", color.tip.by, " in either the features or the metadata", call. = F)
        }

        tree <- AddTreeMetadata(tree, object, metadata = metadata, features = feature)
        color.tip <- TRUE
    }

    tree <- ggtree(tree) + geom_tiplab(offset = 0.03)

    if (color.tip) {
        tree <- tree +
            geom_tippoint(aes(color = .data[[color.tip.by]]), size = 4)

        if (continues.scale) {
            vals <- tree$data %>% na.omit() %>% pull(color.tip.by)
            midpoint <- (max(vals) + min(vals)) / 2
            tree <- tree + scale_color_gradient2(low = muted("blue"), mid = "#ffffa8", high = muted("red"), midpoint = midpoint, na.value = "#000000")
        }
    }

    return(tree)
}

#' Build a lineage tree for the given sequences
#'
#' @param sequences Named vector of sequences
#' @param root Should the tree be rooted. Default = FALSE
#' @param outgroup Outgroup to root the tree

BuildLineageTree <- function(sequences, outgroup = NULL, root = F) {
    if (root & is.null(outgroup)) {
        stop("Missing outgroup to root the tree.", call. = F)
    }

    if (length(sequences) <= 2) {
        stop("Need at least 3 sequences to build a tree.", call. = F)
    }

    if (is.null(names(sequences))) {
        stop("Sequence vector should be named", call. = F)
    }

    distance.mat <- GetDistance(sequences = sequences,  distance.method = "lv" )

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
    header <- NULL
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

#' Add metadata to phylogenetic tree
#'
#' @param tree Tree of class phylo
#' @param object Seurat object
#' @param metadata Metadata column(s) to add to the tree
#' @param features Feature(s) to add to the tree
#'
#' @importFrom tibble tibble
#' @importFrom treeio full_join

AddTreeMetadata <- function(tree, object, metadata = NULL, features = NULL) {
    if (is.null(metadata) & is.null(features)) {
        stop("Specify either group.by or feature", call. = F)
    }

    tree.tibble <- tibble(label = tree$tip.label)

    if (sum(tree.tibble$label %in% colnames(object)) != nrow(tree.tibble) - 1) {
        stop("Cells from tree not all found in object", call. = F)
    }

    if (!is.null(metadata)) {
        for (column in metadata) {
            tree.tibble[[column]] <- object@meta.data[tree.tibble$label, column]
        }
    }

    if (!is.null(features)) {
        for (feature in features) {
            tree.tibble[[feature]] <- object@assays[[Seurat::DefaultAssay(object)]]@data[feature, ][tree.tibble$label]
        }
    }

    return(full_join(tree, tree.tibble, by = "label"))
}

#' Returns the mutation rate of each cell's VDJ chain to "mutation.rate" metadata column
#'
#' @param object Seurat object#'
#' @param clonotype Clonotype to get mutation rate from. Default uses every clonotype available.
#' @param regions VDJ region to use sequence for. DEfault = V.
#' @param reference Path to fasta reference sequence file
#' @param chain Which chain to use
#' @param clonotype.column name of the metadata column that contains the clonotypes
#'
#' @importFrom dplyr %>% filter

GetMutationRate <- function(
        object,
        clonotype = NULL,
        regions = "V",
        reference,
        chain = c("VDJ", "VJ"),
        clonotype.column = "clonotype"
        ) {

    if (DefaultAssayVDJ(object) != "BCR") {
        print("Mutationrate can only be computed on BCR data!", call. = F)
    }

    data <- object@meta.data %>% filter(!is.na(.data[[clonotype.column]]))
    if (is.null(clonotype)) {
        clonotypes <- data[[clonotype.column]] %>% unique()
    } else {
        clonotypes <- clonotype
    }

    chain <- match.arg(chain)

    distances <- c()

    for (clonotype.name in clonotypes) {
        clonotypedata <- data %>% filter(.data[[clonotype.column]] == clonotype.name)
        sequences <- GetSequences(object = object, metadata = clonotypedata, regions = regions, reference = reference, chain = chain)

        distance.mat <- GetDistance(sequences = sequences,  distance.method = "lv" )
        distance <- distance.mat[1,2:length(distance.mat[1,])]
        names(distance) <- colnames(distance.mat)[2:length(distance.mat[1,])]
        distances <- c(distances, distance)
    }
    return(distances)
}

#' Get a distance matrix with distances calculated from a reference
#'
#' @param sequences Named vector of sequences
#' @param distance.method Method to use to calculate distance between sequences. Defaul = "lv" (Levenstein)

GetDistance <- function(sequences, distance.method = "lv") {

    distance.mat <- stringdist::stringdistmatrix(sequences, sequences, method = distance.method)
    colnames(distance.mat) <- names(sequences)
    rownames(distance.mat) <- names(sequences)

    return(distance.mat)
}

#' Get sequences for reference fasta & clonotype
#'
#' @param object Seurat object
#' @param metadata Metadata of the cells of one clonotype
#' @param regions VDJ regions to use sequence for. Default = V
#' @param reference path to the reference fasta file
#' @param chain Which chain to use?

GetSequences <- function(object, metadata, regions = "V", reference, chain = c("VDJ", "VJ")) {
    cells <- rownames(metadata)
    chain <- match.arg(chain) %>% tolower()

    calls <- list("v_call" = NULL, "d_call" = NULL, "j_call" = NULL, "c_call" = NULL)

    for (region in regions) {
        if (!region %in% AvailableRegions(chain)) {
            stop("Invalid region ", region, call. = F)
        }
        region <- tolower(region)
        region.call <- metadata %>% pull(paste0(chain, ".", region, "_gene")) %>% unique()

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

    return(sequences)

}
