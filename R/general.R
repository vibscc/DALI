#' Load 10x VDJ data in a seurat object
#'
#' @param object Seurat object
#' @param data.dir Cellranger output directory
#' @param type VDJ assay type for loaded data. This is automatically detected from input, but can be overwritten when something goes wrong.
#'
#' @importFrom dplyr %>% all_of mutate rename_all select filter
#' @importFrom tibble column_to_rownames
#'
#' @export

Read10X_vdj <- function(object, data.dir, type = NULL) {

    location.annotation.contig <- file.path(data.dir, 'filtered_contig_annotations.csv')

    if (!file.exists(location.annotation.contig)) {
        stop("Contig annotation file (", location.annotation.contig, ") is missing!")
    }

    annotation.contig <- read.csv(location.annotation.contig, stringsAsFactors = F) %>% subset(productive == 'True')

    if (is.null(type)) {
        type <- getDefaultVDJAssay(data.dir)
    }

    columns <- c("barcode", "v_gene","d_gene","j_gene","c_gene", "cdr3","cdr3_nt")

    heavy <- annotation.contig %>%
        subset(grepl("^IGH", c_gene)) %>%
        filter(!duplicated(barcode)) %>%
        select(all_of(columns)) %>%
        mutate(V.fam = get_v_families(v_gene)) %>%
        column_to_rownames('barcode') %>%
        rename_all(~ paste0("h.", .))


    light <- annotation.contig %>%
        subset(grepl("^IG[KL]", c_gene)) %>%
        filter(!duplicated(barcode)) %>%
        select(all_of(columns)) %>%
        mutate(V.fam = get_v_families(v_gene)) %>%
        column_to_rownames('barcode') %>%
        rename_all(~ paste0("l.", .))

    object <- AddVDJDataForType(type, object, heavy, light)
    DefaultAssayVDJ(object) <- type

    return(object)
}

#' @method DefaultAssayVDJ Seurat
#'
#' @export

DefaultAssayVDJ.Seurat <- function(object, ...) {

    return(slot(object, 'misc')[['default.assay.VDJ']])
}

#' @method DefaultAssayVDJ<- Seurat
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

#' Determine assay type from cellranger output
#'
#' @param data.dir Cellranger output directory

getDefaultVDJAssay <- function(data.dir) {
    return('BCR')
}

#' Extract V-family from a vector of v-genes
#'
#' @param v_genes Vector of genes

get_v_families <- function(v_genes) {
    v_families <- c()

    for(v_gene in v_genes) {
        if (!grepl("^IG[KLH]V", v_gene)) {
            v_families <- c(v_families, NA)
            next
        }

        family.full <- strsplit(v_gene, '-')[[1]][1]
        gene <- gsub("[0-9]", "", family.full)
        number <- gsub("[A-Za-z]", "", family.full)
        v_families <- c(v_families, paste0(gene, '-', number))
    }

    return(v_families)
}

#' Add VDJ metadata to misc slot in object
#'
#' @param type VDJ assay type
#' @param object Seurat object
#' @param heavy Data frame with the metadata columns for the heavy chains
#' @param light Data frame with the metadata columns for the light chains

AddVDJDataForType <- function(type, object, heavy, light) {
    if (!'VDJ' %in% names(slot(object, 'misc'))) {
        slot(object, 'misc')[['VDJ']] <- list()
    }

    slot(object, 'misc')[['VDJ']][[type]] <- list(heavy = heavy, light = light)

    return(object)
}
