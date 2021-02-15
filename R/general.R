#' Load 10x VDJ data in a seurat object
#'
#' @param object Seurat object
#' @param data.dir Cellranger output directory
#' @param type VDJ assay type for loaded data. This is automatically detected from input, but can be overwritten when something goes wrong.
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#'
#' @importFrom dplyr %>% all_of mutate rename_all select filter
#' @importFrom tibble column_to_rownames
#' @importFrom rlang .data
#' @importFrom utils read.csv
#'
#' @export

Read10X_vdj <- function(object, data.dir, type = NULL, force = F) {

    location.annotation.contig <- file.path(data.dir, "filtered_contig_annotations.csv")
    location.metrics <- file.path(data.dir, "metrics_summary.csv")

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

    columns <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt")

    heavy.name <- if (type == "TCR") "a" else "h"
    light.name <- if (type == "TCR") "b" else "l"

    heavy <- annotation.contig %>%
        filter(grepl("^IGH|^TRA", .data$c_gene)) %>%
        filter(!duplicated(.data$barcode)) %>%
        select(all_of(columns)) %>%
        mutate(v_fam = get_v_families(.data$v_gene, type)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(heavy.name, ".", .))

    light <- annotation.contig %>%
        filter(grepl("^IG[KL]|^TRB", .data$c_gene)) %>%
        filter(!duplicated(.data$barcode)) %>%
        select(all_of(columns)) %>%
        mutate(v_fam = get_v_families(.data$v_gene, type)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(light.name, ".", .))

    object <- AddVDJDataForType(type, object, heavy, light, force)
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

    object <- Seurat::AddMetaData(object, slot(object, 'misc')[['VDJ']][[value]][['heavy']])
    object <- Seurat::AddMetaData(object, slot(object, 'misc')[['VDJ']][[value]][['light']])

    slot(object, 'misc')[['default.assay.VDJ']] <- value

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
#' @param heavy Data frame with the metadata columns for the heavy chains
#' @param light Data frame with the metadata columns for the light chains
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE

AddVDJDataForType <- function(type, object, heavy, light, force = F) {
    if (!force) {
        overlap <- min(length(intersect(colnames(object), rownames(heavy))) / nrow(heavy), length(intersect(colnames(object), rownames(light))) / nrow(light))

        if (overlap < 0.50) {
            stop("Overlap in cell-barcodes is low. Please check if the barcodes in the Seurat object match the barcodes in the VDJ data.\n",
                 "To ignore this check, add the parameter `force = T`", call. = F)
        }
    }

    if (!'VDJ' %in% names(slot(object, 'misc'))) {
        slot(object, 'misc')[['VDJ']] <- list()
    }

    slot(object, 'misc')[['VDJ']][[type]] <- list(heavy = heavy, light = light)

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
