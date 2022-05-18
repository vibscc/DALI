#' Load 10x VDJ data in a seurat object
#'
#' @param object Seurat object
#' @param data.dir Cellranger output directory
#' @param assay VDJ assay type for loaded data. This is automatically detected from input, but can be overwritten when something goes wrong.
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#' @param sort.by Column to sort the data to determine if chain is primary or secondary. Options = umis, reads
#' @param use.filtered Load filtered contig annotation. Default = TRUE
#' @param quiet Ignore warnings. Default = FALSE
#'
#' @importFrom dplyr %>% filter left_join
#' @importFrom rlang .data
#' @importFrom utils read.csv
#'
#' @export

Read10X_vdj <- function(object, data.dir, assay = NULL, force = F, sort.by = c("umis", "reads"), use.filtered = T, quiet = F) {

    location.annotation.contig <- file.path(data.dir, paste0(if (use.filtered) "filtered" else "all", "_contig_annotations.csv"))
    location.airr.rearrangement <- file.path(data.dir, "airr_rearrangement.tsv")

    sort.by <- match.arg(sort.by)

    if (!file.exists(location.annotation.contig)) {
        stop("Contig annotation file (", location.annotation.contig, ") is missing!", call. = F)
    }

    data <- read.csv(location.annotation.contig, stringsAsFactors = F) %>%
                            filter(grepl("true", .data$productive, ignore.case = T))

    fields <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id")
    fields.extra <- c("fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt", "cdr2", "cdr2_nt", "fwr3", "fwr3_nt", "fwr4", "fwr4_nt")

    if (file.exists(location.airr.rearrangement)) {
        airr.data <- read.csv(location.airr.rearrangement, sep = "\t")
        colnames(airr.data) <- gsub("sequence_id", "contig_id", colnames(airr.data))
        sequence.columns <- grep("sequence", colnames(airr.data), value = T)
        data <- left_join(data, airr.data[, c("contig_id", sequence.columns)], by = "contig_id")
        fields.extra <- c(fields.extra, sequence.columns)
    } else {
        if (!quiet) {
            warning("Could not find airr_rearrangement.tsv. Sequence information will not be loaded and some functionality for BCR lineage tracing will not be available", call. = F)
        }
    }

    for (field in fields.extra) {
        if (field %in% colnames(data)) {
            fields <- c(fields, field)
        }
    }
    columns <- gsub("raw_clonotype_id", "clonotype", fields)

    return(ReadData(object, assay = assay, data = data, fields = fields, columns = columns, force = force, sort.by = sort.by))
}

#' Load 10x VDJ data from multiple directories into a Seurat object containing multiple samples
#'
#' @param object Seurat object
#' @param id_column (Optional:) Metadata column that contains the information to distinguish samples.
#' @param data.dir directory containing Cellranger output directories for every sample. If column is specified this will be a named list
#' @param assay VDJ assay type for loaded data. This is automatically detected from input, but can be overwritten when something goes wrong.
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#' @param sort.by Column to sort the data to determine if chain is primary or secondary. Options = umis, reads
#' @param use.filtered Load filtered contig annotation. Default = TRUE
#' @param quiet Ignore warnings. Default = FALSE
#'
#' @importFrom dplyr %>% filter left_join
#' @importFrom rlang .data
#' @importFrom utils read.csv
#'
#' @export

Read10X_MultiVDJ <- function(object, data.dir, id_column = NULL , assay = NULL, force = F, sort.by = c("umis","reads"), use.filtered = T, quiet = F) {
    if (!is.null(id_column)) {
        # check if the sample names are unique
        if (length(unique(names(object@meta.data[[id_column]]))) != length(names(object@meta.data[[id_column]]))) {
            stop("Sample names are not unique:", unique(names(object@meta.data[[id_column]])))
        }

        samples <- as.factor(object@meta.data[[id_column]])
        directories <- data.dir
        overlap_matrix <- matrix(data = NA, nrow = length(levels(samples)), ncol = length(directories))
    } else {
        # Define samples using suffix & split barcodes into seperate elements in a list
        if (is.null(names(data.dir))) {
            stop("Please use a named vector for your directories")
        }
        samples <- strsplit(colnames(object), "_") %>% sapply("[", 2) %>% as.factor()
        directories <- list.dirs(data.dir, recursive = F)
        cellbarcodes <- strsplit(colnames(object), "_") %>% sapply("[", 1) %>% split(samples)

        # Check if there aren't too many VDJ directories
        if (length(cellbarcodes) < length(directories)) {
            stop("Too many VDJ data files for ",length(cellbarcodes)," samples!")
        }

        overlap_matrix <- matrix(data = NA, nrow = length(cellbarcodes), ncol = length(directories))
    }

    data <- list()
    fields <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id")
    fields.extra <- c("fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt", "cdr2", "cdr2_nt", "fwr3", "fwr3_nt", "fwr4", "fwr4_nt")

    ndata <- 1
    for  (vdjdir in directories) {
        # Convert vdjfiles to df
        vdjdir <- directories[1]
        location.annotation.contig <- file.path(vdjdir, paste0(if (use.filtered) "filtered" else "all", "_contig_annotations.csv"))
        location.airr.rearrangement <- file.path(vdjdir, "airr_rearrangement.tsv")

        if (!file.exists(location.annotation.contig)) {
            stop("Contig annotation file (", location.annotation.contig, ") is missing!", call. = F)
        }

        vdj_df <- read.csv(location.annotation.contig, stringsAsFactors = F) %>%
            filter(grepl("true", .data$productive, ignore.case = T))

        if (file.exists(location.airr.rearrangement)) {
            airr.data <- read.csv(location.airr.rearrangement, sep = "\t")
            colnames(airr.data) <- gsub("sequence_id", "contig_id", colnames(airr.data))
            sequence.columns <- grep("sequence", colnames(airr.data), value = T)
            vdj_df <- left_join(vdj_df, airr.data[, c("contig_id", sequence.columns)], by = "contig_id")
            fields.extra <- c(fields.extra, sequence.columns) %>% unique()
        } else {
            if (!quiet) {
                warning("Could not find airr_rearrangement.tsv. Sequence information will not be loaded and some functionality for BCR lineage tracing will not be available", call. = F)
            }
        }

        if (!is.null(id_column)) {
            # Change cellbarcode according to specified metadata column
            sample_id <- names(directories)[ndata]
            cellbarcodes <- grepl(sample_id,object@meta.data[[id_column]])
            cellbarcodes <- colnames(object)[cellbarcodes]
            suffix <- strsplit(colnames(object), "_") %>% sapply("[", 2)
            vdj_df$barcodes <- paste0(vdj_df$barcodes,"_",suffix)
        } else {
            vdj_df$ID <- ndata
        }

        data <- rbind(data, vdj_df)
        ndata <- ndata + 1
    }

    if (is.null(id_column)) {
        for (row in 1:length(cellbarcodes)) {
            for (col in 1:length(directories)) {
                # Check overlap of barcodes and store in overlapmatrix
                overlap_matrix[row,col] <- length(intersect(cellbarcodes[[row]], data[data$ID == col,]$barcode)) / nrow(data[data$ID == col,])
            }
        }
        # Now to check which vdj data belongs to which sample & Change the cellbarcode accordingly
        for (i in 1:ncol(overlap_matrix)) {
            df_id <- which.max(overlap_matrix[,i])
            data[data$ID == i,]$barcode <- paste0(data[data$ID == i,]$barcode,names(cellbarcodes[df_id]))
        }
    }

    for (field in fields.extra) {
        if (field %in% colnames(data)) {
            fields <- c(fields, field)
        }
    }
    columns <- gsub("raw_clonotype_id", "clonotype", fields)

    return(ReadData(object = object, assay = assay, data = data, fields = fields, columns = columns, force = force, sort.by = sort.by))
}

#' Load 10x VDJ data in a seurat object from 10X AIRR rearrangement tsv
#'
#' @param object Seurat object
#' @param file 10X AIRR rearrangement tsv
#' @param assay VDJ assay for loaded data
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#'
#' @importFrom dplyr %>% filter
#' @importFrom rlang .data
#' @importFrom utils read.csv
#'
#' @export

Read10X_AIRR <- function(object, file, assay, force = F) {
    data <- read.csv(file, sep = "\t") %>% filter(grepl("true", .data$productive, ignore.case = T))

    fields <- c("cell_id", "v_call", "d_call", "j_call", "c_call", "junction_aa", "junction", "consensus_count", "clone_id")
    columns <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "umis", "clonotype")

    sequence.columns <- grep("sequence", colnames(data), value = T)

    if (length(sequence.columns) > 0) {
        fields <- c(fields, sequence.columns)
        columns <- c(columns, sequence.columns)
    }

    return(ReadData(object, assay, data, fields, columns, force, sort.by = "umis"))
}

#' Load VDJ data from an AIRR rearrangement file
#'
#' @param object Seurat object
#' @param files AIRR rearrangement tsv. Can be multiple
#' @param assay VDJ assay for loaded data
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

Read_AIRR <- function(object, files, assay, fields, columns, only.productive = T, productive.field = "productive", force = F) {
    data <- data.frame(matrix(ncol = 0, nrow = 0))
    for (f in files) {
        d <- read.csv(f, sep = "\t")
        data <- bind_rows(data, d)
    }

    if (only.productive) {
        data <- data %>% filter(grepl("true", .data[[productive.field]], ignore.case = T))
    }

    required <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "umis", "clonotype")

    for (required.field in required) {
        if (!required.field %in% columns) {
            stop("Missing required field ", required.field, call. = F)
        }
    }

    return(ReadData(object, assay, data, fields, columns, force = force, sort.by = "umis"))
}

#' Load data in the Seurat object
#'
#' @param object Seurat object
#' @param assay VDJ assay for loaded data
#' @param data data frame containg all data
#' @param fields Fields to select from the data
#' @param columns Rename fields to these names. If not specified, just use the field names
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#' @param sort.by Column to sort the data to determine if chain is primary or secondary. Options = umis, reads
#'
#' @importFrom dplyr %>% add_count all_of arrange desc filter mutate mutate_all na_if rename_all select
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom rlang .data

ReadData <- function(object, assay, data, fields, columns = NULL, force = F, sort.by = c("umis", "reads")) {

    if (is.null(assay)) {
        if (sum(grepl("^TR[ABDG]", data$c_gene)) > 0) {
            assay <- "TCR"
        } else if (sum(grepl("^IG[HKL]", data$c_gene)) > 0) {
            assay <- "BCR"
        } else {
            stop("Could not determine if the data is TCR or BCR, please provide the data assay manually with the `assay` parameter", call. = F)
        }
    }

    if (!assay %in% c("BCR", "TCR")) {
        stop("Invalid assay '", assay, "', must be one of TCR or BCR")
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

    vdj.prefix <- "vdj"
    vj.prefix <- "vj"

    vdj.regex <- if (assay == "TCR") "^TR[BD]" else "^IGH"
    vj.regex <- if (assay == "TCR") "^TR[AG]" else "^IG[KL]"

    data <- data %>% mutate_all(~ na_if(.x, ""))

    vdj <- data %>%
        filter(grepl(vdj.regex, .data[[FieldForColumn("c_gene", fields, columns)]])) %>%
        add_count(.data[[FieldForColumn("barcode", fields, columns)]]) %>%
        mutate(dual_IR = .data$n == 2) %>%
        mutate(multichain = .data$n > 2) %>%
        select(all_of(fields)) %>%
        `colnames<-`(columns) %>%
        mutate(v_fam = GetVFamilies(.data$v_gene, assay))

    vdj.primary <- vdj %>%
        arrange(desc(.data[[sort.by]]), .data$v_gene) %>%
        filter(!duplicated(.data$barcode)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(vdj.prefix, ".", .)) %>%
        rownames_to_column(var = "barcode")

    vdj.secondary <- vdj %>%
        filter(.data$dual_IR) %>%
        arrange(.data[[sort.by]], desc(.data$v_gene)) %>%
        filter(!duplicated(.data$barcode)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(vdj.prefix, ".", .)) %>%
        rownames_to_column(var = "barcode")

    colnames(vdj.primary) <- gsub(".*\\.clonotype", "clonotype", colnames(vdj.primary))
    colnames(vdj.secondary) <- gsub(".*\\.clonotype", "clonotype", colnames(vdj.secondary))

    vj <- data %>%
        filter(grepl(vj.regex, .data[[FieldForColumn("c_gene", fields, columns)]])) %>%
        add_count(.data[[FieldForColumn("barcode", fields, columns)]]) %>%
        mutate(dual_IR = .data$n == 2) %>%
        mutate(multichain = .data$n > 2) %>%
        select(all_of(fields)) %>%
        `colnames<-`(columns) %>%
        mutate(v_fam = GetVFamilies(.data$v_gene, assay))

    vj.primary <- vj %>%
        arrange(desc(.data[[sort.by]]), .data$v_gene) %>%
        filter(!duplicated(.data$barcode)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(vj.prefix, ".", .)) %>%
        rownames_to_column(var = "barcode")

    vj.secondary <- vj %>%
        filter(.data$dual_IR) %>%
        arrange(.data[[sort.by]], desc(.data$v_gene)) %>%
        filter(!duplicated(.data$barcode)) %>%
        column_to_rownames("barcode") %>%
        rename_all(~ paste0(vj.prefix, ".", .)) %>%
        rownames_to_column(var = "barcode")

    colnames(vj.primary) <- gsub(".*\\.clonotype", "clonotype", colnames(vj.primary))
    colnames(vj.secondary) <- gsub(".*\\.clonotype", "clonotype", colnames(vj.secondary))

    if (nrow(vdj.primary) == 0 & !force) {
        stop("Could not find VDJ (heavy/alpha) primary chains in the data!", call. = F)
    }

    if (nrow(vj.primary) == 0 & !force) {
        stop("Could not find VJ (light/beta) primary chains in the data!", call. = F)
    }

    object <- AddVDJDataForAssay(assay, object, vdj.primary, vdj.secondary, vj.primary, vj.secondary, force)
    DefaultAssayVDJ(object) <- assay

    return(object)
}

#' @method DefaultAssayVDJ Seurat
#'
#' @importFrom methods slot
#'
#' @export

DefaultAssayVDJ.Seurat <- function(object, ...) {

    return(slot(object, "misc")[["default.assay.VDJ"]])
}

#' @method DefaultAssayVDJ<- Seurat
#'
#' @importFrom methods slot slot<-
#'
#' @export

"DefaultAssayVDJ<-.Seurat" <- function(object, ..., value) {

    if (!value %in% names(x = slot(object, "misc")[["VDJ"]])) {
        stop("Cannot find assay ", value)
    }

    chain <- DefaultChainVDJ(object)

    object <- Seurat::AddMetaData(object, GetInfoForMetadata(object, value, paste0("vj.", chain)))
    object <- Seurat::AddMetaData(object, GetInfoForMetadata(object, value, paste0("vdj.", chain)))

    slot(object, "misc")[["default.assay.VDJ"]] <- value

    return(object)
}

#' @method DefaultChainVDJ Seurat
#'
#' @importFrom methods slot
#'
#' @export

DefaultChainVDJ.Seurat <- function(object, ...) {

    return(slot(object, "misc")[["default.chain.VDJ"]])
}

#' @method DefaultChainVDJ<- Seurat
#'
#' @importFrom methods slot slot<-
#'
#' @export

"DefaultChainVDJ<-.Seurat" <- function(object, ..., value) {

    if (!value %in% c("primary", "secondary")) {
        stop("Chain must either be primary or secondary")
    }

    assay <- DefaultAssayVDJ(object)

    object <- Seurat::AddMetaData(object, GetInfoForMetadata(object, assay, paste0("vj.", value)))
    object <- Seurat::AddMetaData(object, GetInfoForMetadata(object, assay, paste0("vdj.", value)))

    slot(object, "misc")[["default.chain.VDJ"]] <- value

    return(object)
}

#' Extract V-family from a vector of v-genes
#'
#' @param v_genes Vector of genes
#' @param assay VDJ assay type. Options = TCR; BCR
#'
#' @importFrom stringr str_replace_all

GetVFamilies <- function(v_genes, assay) {
    v_families <- c()

    if (assay == "TCR") { pattern <- "TR[ABD]" }
    if (assay == "BCR" || is.null(assay) ) {  pattern <- "^IG[KLH]V" }

    for (v_gene in v_genes) {
        if (!grepl(pattern, v_gene)) {
            v_families <- c(v_families, NA)
            next
        }

        family.full <- strsplit(v_gene, "-")[[1]][1]
        gene <- gsub("[0-9]", "", family.full)
        number <- gsub("[A-Za-z]", "", family.full)
        v_families <- c(v_families, paste0(gene, "-", number))
    }

    if (assay == "TCR") {
        v_families <- str_replace_all(v_families, "TRAV/DV", "TRADV")
    }

    return(v_families)
}

#' Add VDJ metadata to misc slot in object
#'
#' @param assay VDJ assay type
#' @param object Seurat object
#' @param vdj.primary Data frame with the metadata columns for the primary heavy/alpha chains
#' @param vdj.secondary Data frame with the metadata columns for the secondary heavy/alpha chains
#' @param vj.primary Data frame with metadata columns for the primary light/beta chains
#' @param vj.secondary Data frame with the metadata columns for the secondary light/beta chains
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#'
#' @importFrom dplyr %>%

AddVDJDataForAssay <- function(assay, object, vdj.primary, vdj.secondary, vj.primary, vj.secondary, force = F) {
    if (!force) {
        overlap <- min(length(intersect(colnames(object), vdj.primary$barcode)) / nrow(vdj.primary), length(intersect(colnames(object), vj.primary$barcode)) / nrow(vj.primary))

        if (overlap < 0.50) {
            stop("Overlap in cell-barcodes is low. Please check if the barcodes in the Seurat object match the barcodes in the VDJ data.\n", call. = F)
        }
    }

    if (!"VDJ" %in% names(slot(object, "misc"))) {
        slot(object, "misc")[["VDJ"]] <- list()
    }

    if (assay %in% names(slot(object, "misc")[["VDJ"]])) {
        data <- list(
            vdj.primary = CreateVDJData(vdj.primary, slot(object, "misc")[["VDJ"]][[assay]][["vdj.primary"]]),
            vdj.secondary = CreateVDJData(vdj.secondary, slot(object, "misc")[["VDJ"]][[assay]][["vdj.secondary"]]),
            vj.primary = CreateVDJData(vj.primary, slot(object, "misc")[["VDJ"]][[assay]][["vj.primary"]]),
            vj.secondary = CreateVDJData(vj.secondary, slot(object, "misc")[["VDJ"]][[assay]][["vj.secondary"]])
        )

        intersect1 <- intersect(slot(object, "misc")[["VDJ"]][[assay]][["vdj.primary"]]$barcode, vdj.primary$barcode)
        intersect2 <- intersect(slot(object, "misc")[["VDJ"]][[assay]][["vdj.secondary"]]$barcode, vdj.secondary$barcode)
        intersect3 <- intersect(slot(object, "misc")[["VDJ"]][[assay]][["vj.primary"]]$barcode, vj.primary$barcode)
        intersect4 <- intersect(slot(object, "misc")[["VDJ"]][[assay]][["vj.secondary"]]$barcode, vj.secondary$barcode)

        overwritten <- unique(c(intersect1, intersect2, intersect3, intersect4)) %>% length()

        if (overwritten > 0) {
            warning("A total of ", overwritten, " cells have lost some original data for the VDJ assay '", assay, "'.", call. = F)
        }
    } else {
        data <- list(
            vdj.primary = CreateVDJData(vdj.primary),
            vdj.secondary = CreateVDJData(vdj.secondary),
            vj.primary = CreateVDJData(vj.primary),
            vj.secondary = CreateVDJData(vj.secondary)
        )
    }

    slot(object, "misc")[["default.chain.VDJ"]] <- "primary"

    slot(object, "misc")[["VDJ"]][[assay]] <- data

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
    tables <- c("vdj.primary", "vj.primary")

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

#' Get metadata for the given assay/chain combination
#'
#' @param object Seurat object
#' @param assay VDJ assay
#' @param chain vdj or vj chain
#'
#' @importFrom tibble column_to_rownames

GetInfoForMetadata <- function(object, assay, chain) {
    data <- slot(object, "misc")[["VDJ"]][[assay]][[chain]] %>% column_to_rownames("barcode")
    columns.to.ignore <- grep("sequence", colnames(data))

    if (length(columns.to.ignore) > 0) {
        return(data[, -columns.to.ignore])
    }

    return(data)
}

#' Merge seuratobjects and include their VDJ data
#'
#' @param ... seurat objects separated by comma
#'
#' @export

MergeVDJ <- function(...) {
    objects <- list(...)
    mergedobj <- merge(...)

    mergedobj <- AddVDJForAssay(mergedobj = mergedobj, objectlist = objects, assay = "BCR")
    mergedobj <- AddVDJForAssay(mergedobj = mergedobj, objectlist = objects, assay = "TCR")

    return(mergedobj)
}

#' Splits a Seurat object into a list of objects with their respective VDJ data in the misc slot
#'
#' @param object Seurat object
#' @param column metadata column to split the data on
#'
#' @export

SplitObject_VDJ <- function(object, column) {
    # Get barcodes belonging to each samples based on specified metadata column
    objects <- SplitObject(object, split.by = column)
    vdj_absence <- 0
    for (obj in objects) {
        barcodes <- colnames(obj)

        BCR <- subsetVDJ(object = object,barcodes = barcodes, assay = "BCR")
        TCR <- subsetVDJ(object = object,barcodes = barcodes, assay = "TCR")

        # Add misc data
        if (!is.null(BCR) & !is.null(TCR)) {
            Misc(object = obj, slot = "VDJ") <- list("TCR" = TCR, "BCR" = BCR)
            Misc(object = obj, slot = "default.assay.VDJ") <- "TCR"
            Misc(object = obj, slot = "default.chain.VDJ") <- DefaultChainVDJ(object)
            DefaultAssayVDJ(obj) <- "TCR"
        } else if (!is.null(BCR)) {
            Misc(object = obj, slot = "VDJ") <- list("BCR" = BCR)
            Misc(object = obj, slot = "default.assay.VDJ") <- "BCR"
            Misc(object = obj, slot = "default.chain.VDJ") <- DefaultChainVDJ(object)
            DefaultAssayVDJ(obj) <- "BCR"
        } else if (!is.null(TCR)) {
            Misc(object = obj, slot = "VDJ") <- list("TCR" = BCR)
            Misc(object = obj, slot = "default.assay.VDJ") <- "TCR"
            Misc(object = obj, slot = "default.chain.VDJ") <- DefaultChainVDJ(object)
            DefaultAssayVDJ(obj) <- "TCR"
        } else {
            vdj_absence <- vdj_absence + 1
            if (length(objects) == vdj_absence) {
                stop("No TCR or BCR assay present, use SplitObject() function from the Seurat Package")
            }
        }
    }
    #return a list of objects
    return(objects)
}
