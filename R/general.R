#' Load 10x VDJ data in a seurat object
#'
#' @param object Seurat object
#' @param data.dir Cellranger output directory
#' @param type VDJ assay type for loaded data. This is automatically detected from input, but can be overwritten when something goes wrong.
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#' @param sort.by Column to sort the data to determine if chain is primary or secondary. Options = umis, reads
#' @param use.filtered Load filtered contig annotation. Default = TRUE
#'
#' @importFrom dplyr %>% add_count all_of arrange desc filter mutate rename_all select
#' @importFrom tibble column_to_rownames
#' @importFrom rlang .data
#' @importFrom utils read.csv
#'
#' @export

Read10X_vdj <- function(object, data.dir, type = NULL, force = F, sort.by = c('umis', 'reads'), use.filtered = T) {

    location.annotation.contig <- file.path(data.dir, paste0(if (use.filtered) "filtered" else "all", "_contig_annotations.csv"))
    location.metrics <- file.path(data.dir, "metrics_summary.csv")

    sort.by <- match.arg(sort.by)

    if (!file.exists(location.annotation.contig)) {
        stop("Contig annotation file (", location.annotation.contig, ") is missing!", call. = F)
    }

    if (!file.exists(location.metrics)) {
        stop("Metrics summary file (", location.metrics, ") is missing!", call. = F)
    }

    annotation.contig <- read.csv(location.annotation.contig, stringsAsFactors = F) %>% filter(.data$productive %in% c('True', 'true'))

    if (is.null(type)) {
        type <- getAssayForData(location.metrics)

        if (is.null(type)) {
            stop("Unable to determine if the data is TCR or BCR. Please check the input files or specify the type via the `type` parameter", call. = F)
        }
    }

    columns <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "reads", "umis", "dual_IR", "raw_clonotype_id")

    heavy.name <- if (type == "TCR") "a" else "h"
    light.name <- if (type == "TCR") "b" else "l"

    heavy <- annotation.contig %>%
        filter(grepl("^IGH|^TRA", .data$c_gene)) %>%
        add_count(.data$barcode) %>%
        mutate(dual_IR = .data$n == 2) %>%
        mutate(multichain = .data$n > 2) %>%
        select(all_of(columns)) %>%
        mutate(v_fam = get_v_families(.data$v_gene, type))

    heavy.primary <- heavy %>%
        arrange(desc(.data[[sort.by]]), .data$v_gene) %>%
        filter(!duplicated(.data$barcode)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(heavy.name, ".", .))

    heavy.secondary <- heavy %>%
        filter(.data$dual_IR) %>%
        arrange(.data[[sort.by]], desc(.data$v_gene)) %>%
        filter(!duplicated(.data$barcode)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(heavy.name, ".", .))

    colnames(heavy.primary) <- gsub(".\\.raw_clonotype_id", "clonotype", colnames(heavy.primary))
    colnames(heavy.secondary) <- gsub(".\\.raw_clonotype_id", "clonotype", colnames(heavy.secondary))

    light <- annotation.contig %>%
        filter(grepl("^IG[KL]|^TRB", .data$c_gene)) %>%
        add_count(.data$barcode) %>%
        mutate(dual_IR = .data$n == 2) %>%
        mutate(multichain = .data$n > 2) %>%
        select(all_of(columns)) %>%
        mutate(v_fam = get_v_families(.data$v_gene, type))

    light.primary <- light %>%
        arrange(desc(.data[[sort.by]]), .data$v_gene) %>%
        filter(!duplicated(.data$barcode)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(light.name, ".", .))

    light.secondary <- light %>%
        filter(.data$dual_IR) %>%
        arrange(.data[[sort.by]], desc(.data$v_gene)) %>%
        filter(!duplicated(.data$barcode)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(light.name, ".", .))

    colnames(light.primary) <- gsub(".\\.raw_clonotype_id", "clonotype", colnames(light.primary))
    colnames(light.secondary) <- gsub(".\\.raw_clonotype_id", "clonotype", colnames(light.secondary))

    object <- AddVDJDataForType(type, object, heavy.primary, heavy.secondary, light.primary, light.secondary, force)
    DefaultAssayVDJ(object) <- type

    return(object)
}

#' @method DefaultAssayVDJ Seurat
#'
#' @importFrom methods slot
#'
#' @export

DefaultAssayVDJ.Seurat <- function(object, ...) {

    return(slot(object, 'misc')[['default.assay.VDJ']])
}

#' @method DefaultAssayVDJ<- Seurat
#'
#' @importFrom methods slot slot<-
#'
#' @export

"DefaultAssayVDJ<-.Seurat" <- function(object, ..., value) {

    if (!value %in% names(x = slot(object, 'misc')[['VDJ']])) {
        stop("Cannot find assay ", value)
    }

    chain <- DefaultChainVDJ(object)

    object <- Seurat::AddMetaData(object, slot(object, 'misc')[['VDJ']][[value]][[paste0('heavy.', chain)]])
    object <- Seurat::AddMetaData(object, slot(object, 'misc')[['VDJ']][[value]][[paste0('light.', chain)]])

    slot(object, 'misc')[['default.assay.VDJ']] <- value

    return(object)
}

#' @method DefaultChainVDJ Seurat
#'
#' @importFrom methods slot
#'
#' @export

DefaultChainVDJ.Seurat <- function(object, ...) {

    return(slot(object, 'misc')[['default.chain.VDJ']])
}

#' @method DefaultChainVDJ<- Seurat
#'
#' @importFrom methods slot slot<-
#'
#' @export

"DefaultChainVDJ<-.Seurat" <- function(object, ..., value) {

    if (!value %in% c('primary', 'secondary')) {
        stop("Chain must either be primary or secondary")
    }

    assay <- DefaultAssayVDJ(object)

    object <- Seurat::AddMetaData(object, slot(object, 'misc')[['VDJ']][[assay]][[paste0('heavy.', value)]])
    object <- Seurat::AddMetaData(object, slot(object, 'misc')[['VDJ']][[assay]][[paste0('light.', value)]])

    slot(object, 'misc')[['default.chain.VDJ']] <- value

    return(object)
}

#' Determine assay type from metrics_summary.csv
#'
#' @param csv.path Path to metrics_summary.csv

getAssayForData <- function(csv.path) {
    metrics <- read.csv(csv.path, header = T)

    if (sum(grepl("IGK", colnames(metrics))) > 0 && sum(grepl("IGH", colnames(metrics))) > 0) {
        return("BCR")
    }

    if (sum(grepl("TRA", colnames(metrics))) > 0 && sum(grepl("TRB", colnames(metrics))) > 0) {
        return("TCR")
    }

    return(NULL)
}

#' Extract V-family from a vector of v-genes
#'
#' @param v_genes Vector of genes
#' @param type TCR/BCR
#'
#' @importFrom stringr str_replace_all


get_v_families <- function(v_genes, type) {
    v_families <- c()

    if (type == "TCR") { pattern <- "TR[ABD]" }
    if (type == "BCR" || is.null(type) ) {  pattern <- "^IG[KLH]V" }

    for (v_gene in v_genes) {
        if (!grepl(pattern, v_gene)) {
            v_families <- c(v_families, NA)
            next
        }

        family.full <- strsplit(v_gene, '-')[[1]][1]
        gene <- gsub("[0-9]", "", family.full)
        number <- gsub("[A-Za-z]", "", family.full)
        v_families <- c(v_families, paste0(gene, '-', number))
    }

    if (type == "TCR") {
        v_families <- str_replace_all(v_families, "TRAV/DV", "TRADV")
    }

    return(v_families)
}

#' Add VDJ metadata to misc slot in object
#'
#' @param type VDJ assay type
#' @param object Seurat object
#' @param heavy.primary Data frame with the metadata columns for the primary heavy chains
#' @param heavy.secondary Data frame with the metadata columns for the secondary heavy chains
#' @param light.primary Data frame with metadata columns for the primary light chains
#' @param light.secondary Data frame with the metadata columns for the secondary light chains
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE

AddVDJDataForType <- function(type, object, heavy.primary, heavy.secondary, light.primary, light.secondary, force = F) {
    if (!force) {
        overlap <- min(length(intersect(colnames(object), rownames(heavy.primary))) / nrow(heavy.primary), length(intersect(colnames(object), rownames(light.primary))) / nrow(light.primary))

        if (overlap < 0.50) {
            stop("Overlap in cell-barcodes is low. Please check if the barcodes in the Seurat object match the barcodes in the VDJ data.\n",
                 "To ignore this check, add the parameter `force = T`", call. = F)
        }
    }

    if (!'VDJ' %in% names(slot(object, 'misc'))) {
        slot(object, 'misc')[['VDJ']] <- list()
    }

    slot(object, 'misc')[['default.chain.VDJ']] <- 'primary'

    slot(object, 'misc')[['VDJ']][[type]] <- list(
        heavy.primary = heavy.primary,
        heavy.secondary = heavy.secondary,
        light.primary = light.primary,
        light.secondary = light.secondary
    )

    return(object)
}

#' Check validity of an object
#'
#' @param object Seurat object
#'
#' @importFrom methods is slot
#' @export

isValidSeuratObject <- function(object) {
    if (!is(object, "Seurat")) {
        return(F)
    }

    if (is.null(slot(object, 'misc')[['VDJ']])) {
        return(F)
    }

    return(T)
}

#' Get available chains from object
#'
#' @param object Seurat object

availableChains <- function(object) {
    assay <- DefaultAssayVDJ(object)

    if (grepl("BCR", assay, ignore.case = T)) {
        return(c("H", "L"))
    }

    if (grepl("TCR", assay, ignore.case = T)) {
        return(c("A", "B"))
    }

    return(c("H", "L", "A", "B"))
}
