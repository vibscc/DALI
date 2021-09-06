#' Load 10x VDJ data in a seurat object
#'
#' @param object Seurat object
#' @param data.dir Cellranger output directory
#' @param type VDJ assay type for loaded data. This is automatically detected from input, but can be overwritten when something goes wrong.
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#' @param sort.by Column to sort the data to determine if chain is primary or secondary. Options = umis, reads
#' @param use.filtered Load filtered contig annotation. Default = TRUE
#'
#' @importFrom dplyr %>% filter
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

    annotation.contig <- read.csv(location.annotation.contig, stringsAsFactors = F) %>% filter(grepl('true', .data$productive, ignore.case = T))

    if (is.null(type)) {
        if (!file.exists(location.metrics)) {
            stop("Metrics summary file (", location.metrics, ") is missing!", call. = F)
        }

        type <- GetAssayForData(location.metrics)

        if (is.null(type)) {
            stop("Unable to determine if the data is TCR or BCR. Please check the input files or specify the type via the `type` parameter", call. = F)
        }
    } else {
        type <- match.arg(type, choices = c("BCR", "TCR"))
    }

    fields <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id")
    fields.extra <- c("fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt", "cdr2", "cdr2_nt", "fwr3", "fwr3_nt", "fwr4", "fwr4_nt")

    for (field in fields.extra) {
        if (field %in% colnames(annotation.contig)) {
            fields <- c(fields, field)
        }
    }
    columns <- gsub("raw_clonotype_id", "clonotype", fields)

    return(ReadData(object, type = type, data = annotation.contig, fields = fields, columns = columns, force = force, sort.by = sort.by))
}

#' Load 10x VDJ data in a seurat object from 10X AIRR rearrangement tsv
#'
#' @param object Seurat object
#' @param file 10X AIRR rearrangement tsv
#' @param type VDJ assay type for loaded data
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#'
#' @importFrom dplyr %>% filter
#' @importFrom rlang .data
#' @importFrom utils read.csv
#'
#' @export

Read10X_AIRR <- function(object, file, type, force = F) {
    data <- read.csv(file, sep = "\t") %>% filter(grepl('true', .data$productive, ignore.case = T))

    fields <- c("cell_id", "v_call", "d_call", "j_call", "c_call", "junction_aa", "junction", "consensus_count", "clone_id")
    columns <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "umis", "clonotype")

    return(ReadData(object, type, data, fields, columns, force, sort.by = 'umis'))
}

#' Load VDJ data from an AIRR rearrangement file
#'
#' @param object Seurat object
#' @param files AIRR rearrangement tsv. Can be multiple
#' @param type VDJ assay type for loaded data
#' @param fields Fields to keep from the AIRR file
#' @param columns Column names to map the field names to
#' @param only.productive Keep only productive rearrangements. Default = TRUE
#' @param productive.field Field containing productive information. Default = productive
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#'
#' @importFrom dplyr bind_rows
#' @importFrom rlang .data
#' @importFrom utils read.csv
#'
#' @export

Read_AIRR <- function(object, files, type, fields, columns, only.productive = T, productive.field = "productive", force = F) {
    data <- data.frame(matrix(ncol = 0, nrow = 0))
    for (f in files) {
        d <- read.csv(f, sep = "\t")
        data <- bind_rows(data, d)
    }

    if (only.productive) {
        data <- data %>% filter(grepl('true', .data[[productive.field]], ignore.case = T))
    }

    required <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "umis", "clonotype")

    for (required.field in required) {
        if (!required.field %in% columns) {
            stop("Missing required field ", required, call. = F)
        }
    }

    return(ReadData(object, type, data, fields, columns, force = force, sort.by = "umis"))
}

#' Load data in the Seurat object
#'
#' @param object Seurat object
#' @param type VDJ assay type for loaded data
#' @param data data frame containg all data
#' @param fields Fields to select from the data
#' @param columns Rename fields to these names. If not specified, just use the field names
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#' @param sort.by Column to sort the data to determine if chain is primary or secondary. Options = umis, reads
#'
#' @importFrom dplyr %>% add_count all_of arrange desc filter mutate mutate_all na_if rename_all select
#' @importFrom tibble column_to_rownames
#' @importFrom rlang .data

ReadData <- function(object, type, data, fields, columns = NULL, force = F, sort.by = c('umis', 'reads')) {

    if (!type %in% c("BCR", "TCR")) {
        stop("Invalid type '", type, "', must be one of TCR or BCR")
    }

    sort.by <- match.arg(sort.by)

    if (is.null(columns)) {
        columns <- fields
    } else if (length(fields) != length(columns)) {
        stop("Fields and columns should be equal length. Got ", length(fields), " and ", length(columns), " instead", call. = F)
    }

    for (field in fields) {
        if (!field %in% colnames(data)) {
            stop("Missing field ", field, " in data", call. = F)
        }
    }

    fields <- c(fields, "dual_IR")
    columns <- c(columns, "dual_IR")

    heavy.name <- if (type == "TCR") "a" else "h"
    light.name <- if (type == "TCR") "b" else "l"

    heavy.regex <- if (type == "TCR") "^TRA" else "^IGH"
    light.regex <- if (type == "TCR") "^TRB" else "^IG[KL]"

    data <- data %>% mutate_all(~ na_if(.x, ""))

    heavy <- data %>%
        filter(grepl(heavy.regex, .data[[FieldForColumn("c_gene", fields, columns)]])) %>%
        add_count(.data[[FieldForColumn("barcode", fields, columns)]]) %>%
        mutate(dual_IR = .data$n == 2) %>%
        mutate(multichain = .data$n > 2) %>%
        select(all_of(fields)) %>%
        `colnames<-`(columns) %>%
        mutate(v_fam = GetVFamilies(.data$v_gene, type))

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

    colnames(heavy.primary) <- gsub(".\\.clonotype", "clonotype", colnames(heavy.primary))
    colnames(heavy.secondary) <- gsub(".\\.clonotype", "clonotype", colnames(heavy.secondary))

    light <- data %>%
        filter(grepl(light.regex, .data[[FieldForColumn("c_gene", fields, columns)]])) %>%
        add_count(.data[[FieldForColumn("barcode", fields, columns)]]) %>%
        mutate(dual_IR = .data$n == 2) %>%
        mutate(multichain = .data$n > 2) %>%
        select(all_of(fields)) %>%
        `colnames<-`(columns) %>%
        mutate(v_fam = GetVFamilies(.data$v_gene, type))

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

    colnames(light.primary) <- gsub(".\\.clonotype", "clonotype", colnames(light.primary))
    colnames(light.secondary) <- gsub(".\\.clonotype", "clonotype", colnames(light.secondary))

    if (nrow(heavy.primary) == 0) {
        stop("Could not find heavy primary chains in the data!", call. = F)
    }

    if (nrow(light.primary) == 0) {
        stop("Could not find light primary chains in the data!", call. = F)
    }

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

GetAssayForData <- function(csv.path) {
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


GetVFamilies <- function(v_genes, type) {
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

IsValidSeuratObject <- function(object) {
    if (!is(object, "Seurat")) {
        return(F)
    }

    if (is.null(slot(object, "misc")[["VDJ"]])) {
        return(F)
    }

    assays <- c("TCR", "BCR")
    tables <- c("heavy.primary", "heavy.secondary", "light.primary", "light.secondary")

    missing <- 0
    for (assay in assays) {
        assay.slot <- slot(object, "misc")[["VDJ"]][[assay]]

        if (is.null(assay.slot)) {
            missing <- missing + 1
            next
        }

        for (table in tables) {
            if (nrow(assay.slot[[table]]) == 0) {
                return(F)
            }
        }
    }

    if (missing == length(assays)) {
        return(F)
    }

    return(T)
}

#' Get available chains from object
#'
#' @param object Seurat object

AvailableChains <- function(object) {
    assay <- DefaultAssayVDJ(object)

    if (grepl("BCR", assay, ignore.case = T)) {
        return(c("H", "L"))
    }

    if (grepl("TCR", assay, ignore.case = T)) {
        return(c("A", "B"))
    }

    return(c("H", "L", "A", "B"))
}
