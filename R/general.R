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

    sort.by <- match.arg(sort.by)

    fields <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id")
    fields.extra <- c("fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt", "cdr2", "cdr2_nt", "fwr3", "fwr3_nt", "fwr4", "fwr4_nt")

    sequence.columns <- GetAIRRSequenceColumns(data.dir, quiet)

    if (!is.null(sequence.columns)) {
        fields.extra <- c(fields.extra, sequence.columns) %>% unique()
    }

    data <- GetVDJ_Dataframe(data.dir, sequence.columns, use.filtered)

    for (field in fields.extra) {
        if (field %in% colnames(data)) {
            fields <- c(fields, field)
        }
    }
    columns <- gsub("raw_clonotype_id", "clonotype", fields)

    return(ReadData(object, assay = assay, data = data, fields = fields, columns = columns, force = force, sort.by = sort.by))
}

#' Load 10x VDJ data into a Seurat object containing multiple samples merged together
#'
#' @param object Seurat object
#' @param data.dir Cellranger output directory for every sample. Can be 1 directory or a named vector/list of multiple directories. The names must correspond to the sample ids in the metdata column specified by column.id
#' @param id.column Metadata column that contains the information to distinguish samples. Used to map the keys of the named vector to the sample in the object
#' @param assay VDJ assay type for loaded data. This is automatically detected from input, but can be overwritten when something goes wrong.
#' @param force Add VDJ data without checking overlap in cell-barcodes. Default = FALSE
#' @param sort.by Column to sort the data to determine if chain is primary or secondary. Options = umis, reads
#' @param use.filtered Load filtered contig annotation. Default = TRUE
#' @param quiet Ignore warnings. Default = FALSE
#'
#' @importFrom dplyr %>% bind_rows filter left_join
#' @importFrom rlang .data
#' @importFrom utils read.csv
#'
#' @export

Read10X_multi_vdj <- function(object, data.dir, id.column = NULL , assay = NULL, force = F, sort.by = c("umis","reads"), use.filtered = T, quiet = F) {

    if (!is.null(names(data.dir))) {
        # We got a named vector/list of directories

        if (is.null(id.column)) {
            stop("id.column cannot be null when providing multiple data directories as input")
        }

        if (!id.column %in% colnames(object@meta.data)) {
            stop("Invalid column '", id.column, "' provided as id.column. The column must be available in the metadata of the object")
        }

        if (length(unique(names(data.dir))) != length(names(data.dir))) {
            stop("Sample names are not unique for the provided data.dirs", call. = F)
        }
    } else if (length(data.dir) != 1) {
        stop("data.dir must either be a named list/vector of directories OR 1 directory", call. = F)
    }

    if (!is.null(id.column)) {
        sample.labels <- unique(object@meta.data[[id.column]])
        sample.suffixes <- c()

        for (sample.label in sample.labels) {
            cell.barcodes <- rownames(object@meta.data[object@meta.data[[id.column]] == sample.label, ])
            sample.suffix <- strsplit(cell.barcodes, "_") %>% sapply("[", 2) %>% unique()

            if (length(sample.suffix) != 1) {
                stop("Found multiple sample suffixes for label '", sample.label, "' from metadata column '", id.column, "'", call. = F)
            }
            sample.suffixes <- c(sample.suffixes, sample.suffix)
        }
        names(sample.suffixes) <- sample.labels
    } else {
        sample.suffixes <- strsplit(colnames(object), "_") %>% sapply("[", 2) %>% unique()
    }

    cell.barcodes.per.suffix <- list()
    for (suffix in sample.suffixes) {
        cell.barcodes.per.suffix[[suffix]] <- colnames(object)[grepl(paste0("_", suffix, "$"), colnames(object))]
        cell.barcodes.per.suffix[[suffix]] <- gsub("[0-9_-]", "", cell.barcodes.per.suffix[[suffix]])
    }

    data.all <- data.frame()
    fields <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id")
    fields.extra <- c("fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt", "cdr2", "cdr2_nt", "fwr3", "fwr3_nt", "fwr4", "fwr4_nt")

    for (i in 1:length(data.dir)) {
        directory <- data.dir[i]
        sample <- names(data.dir)[i]

        if (is.null(sample)) {
            sample <- i
        }

        sequence.columns <- GetAIRRSequenceColumns(directory)
        if (!is.null(sequence.columns)) {
            fields.extra <- c(fields.extra, sequence.columns) %>% unique()
        }

        data.sample <- GetVDJ_Dataframe(directory, sequence.columns, use.filtered)

        if (force && !is.null(id.column)) { # Force the data to the corresponding sample
            data.sample$barcode <- paste0(data.sample$barcode, "_", sample.suffixes[[sample]])
        } else { # Get the suffix with the highest overlap
            highest.overlap <- 0
            matching.suffix <- NULL
            for (suffix in names(cell.barcodes.per.suffix)) {
                overlap <- length(intersect(gsub("[0-9_-]", "", data.sample$barcode), cell.barcodes.per.suffix[[suffix]])) / nrow(data.sample)
                if (overlap > highest.overlap) {
                    highest.overlap <- overlap
                    matching.suffix <- suffix
                }
            }

            if (is.null(matching.suffix)) {
                message(names(cell.barcodes.per.suffix))
                stop("Could not find any suffix in the data that had any overlap with the cellbarcodes of the VDJ data!", call. = F)
            }

            if (!force && highest.overlap <= 0.1) {
                stop("Overlap in cell barcodes is not large enough between the provided VDJ data and the cells in the Seurat object", call. = F)
            }

            data.sample$barcode <- paste0(data.sample$barcode, "_", matching.suffix)
        }

        data.all <- dplyr::bind_rows(data.all, data.sample)
    }

    for (field in fields.extra) {
        if (field %in% colnames(data)) {
            fields <- c(fields, field)
        }
    }
    columns <- gsub("raw_clonotype_id", "clonotype", fields)

    return(ReadData(object = object, assay = assay, data = data.all, fields = fields, columns = columns, force = force, sort.by = sort.by))
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
        has.TCR <- FALSE
        has.BCR <- FALSE

        if (sum(grepl("^TR[ABDG]", data$c_gene)) > 0) {
            has.TCR <- TRUE
        }
        if (sum(grepl("^IG[HKL]", data$c_gene)) > 0) {
            has.BCR <- TRUE
        }

        if (has.TCR && has.BCR) {
            stop("Found both TCR and BCR data. VDJ data should be split by assay to use with DALI", call. = F)
        } else if (has.TCR) {
            assay <- "TCR"
        } else if (has.BCR) {
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

    if (!is.null(chain)) {
        object <- Seurat::AddMetaData(object, GetInfoForMetadata(object, value, paste0("vj.", chain)))
        object <- Seurat::AddMetaData(object, GetInfoForMetadata(object, value, paste0("vdj.", chain)))
    }

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

    if (!is.null(assay)) {
        object <- Seurat::AddMetaData(object, GetInfoForMetadata(object, assay, paste0("vj.", value)))
        object <- Seurat::AddMetaData(object, GetInfoForMetadata(object, assay, paste0("vdj.", value)))
    }

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
    if (!assay %in% names(slot(object, "misc")[["VDJ"]])) {
        stop("Invalid assay: ", assay, call. = F)
    }

    if (!chain %in% names(slot(object, "misc")[["VDJ"]][[assay]])) {
        stop("Invalid chain: ", chain, call. = F)
    }

    data <- slot(object, "misc")[["VDJ"]][[assay]][[chain]] %>% column_to_rownames("barcode")
    columns.to.ignore <- grep("sequence", colnames(data))

    if (length(columns.to.ignore) > 0) {
        return(data[, -columns.to.ignore])
    }

    return(data)
}

#' Merge multiple seurat objects and include their VDJ data
#'
#' @param ... Seurat objects separated by comma
#'
#' @export

MergeVDJ <- function(...) {
    objects <- list(...)
    merged.object <- merge(...) # Merge all data, except for the VDJ data in the misc

    i <- 1
    default.vdj.assay <- NULL
    default.vdj.chain <- NULL

    # Loop over all objects and append all uniqified vdj data into the new object
    for (object in objects) {
        for (assay in names(object@misc$VDJ)) {
            merged.object <- AppendVDJData(merged.object, Uniqify_VDJ(object, i, assay), assay)
        }
        i <- i + 1

        default.vdj.assay <- DefaultAssayVDJ(object)
        default.vdj.chain <- DefaultChainVDJ(object)
    }

    if (!is.null(default.vdj.assay) && !is.null(default.vdj.chain)) {
        DefaultChainVDJ(merged.object) <- default.vdj.chain
        DefaultAssayVDJ(merged.object) <- default.vdj.assay
    }

    return(merged.object)
}

#' Splits a Seurat object into a list of objects with their respective VDJ data in the misc slot
#'
#' @param object Seurat object
#' @param split.by Attribute for splitting.
#' @param quiet Supress warnings. Default = FALSE
#'
#' @importFrom Seurat SplitObject
#'
#' @export

SplitObject_VDJ <- function(object, split.by, quiet = F) {
    if (!split.by %in% colnames(object@meta.data)) {
        stop("Invalid column: ", split.by)
    }

    objects.split <- SplitObject(object, split.by = split.by)

    for (group in names(objects.split)) {
        barcodes <- colnames(objects.split[[group]])

        objects.split[[group]] <- ClearVDJ(objects.split[[group]])

        bcr.data <- SubsetVDJData(object, barcodes, assay = "BCR")
        tcr.data <- SubsetVDJData(object, barcodes, assay = "TCR")

        if (is.null(bcr.data) && is.null(tcr.data)) {
            if (!quiet) {
                warning("No VDJ data found for group '", group, "'")
            }
            next
        }

        if (!is.null(bcr.data)) {
            objects.split[[group]] <- AddVDJDataForAssay("BCR", objects.split[[group]], bcr.data[["vdj.primary"]], bcr.data[["vdj.secondary"]], bcr.data[["vj.primary"]], bcr.data[["vj.secondary"]])
        }

        if (!is.null(tcr.data)) {
            objects.split[[group]] <- AddVDJDataForAssay("TCR", objects.split[[group]], tcr.data[["vdj.primary"]], tcr.data[["vdj.secondary"]], tcr.data[["vj.primary"]], tcr.data[["vj.secondary"]])
        }

        # Set default assay to the default assay of the original object.
        # If this assay is not present for the subset, set the available assay as default
        if (DefaultAssayVDJ(object) %in% names(slot(objects.split[[group]], "misc")[["VDJ"]])) {
            DefaultAssayVDJ(objects.split[[group]]) <- DefaultAssayVDJ(object)
        } else {
            DefaultAssayVDJ(objects.split[[group]]) <- names(slot(objects.split[[group]], "misc")[["VDJ"]])[1]
        }
    }

    return(objects.split)
}

#' Subsets part of a Seurat Object specified by sample_id from a metadata column
#'
#' @param object Seurat object
#' @param subset.by Attribute for subsetting.
#' @param group.id Name of the group to subset for
#' @param assay VDJ assay (TCR or BCR) you want to keep. Default keeps both in the subset
#'
#' @export

SubsetObject_VDJ <- function(object, subset.by, group.id, assay = NULL) {

    if (is.null(object@misc$VDJ)) {
        warning("The object does not contain VDJ data")
    }

    if (!subset.by %in% colnames(object@meta.data)) {
        stop("Invalid subset.by column: ", subset.by, call. = F)
    }

    if (!group.id %in% unique(object@meta.data[[subset.by]])) {
        stop("The requested group (", group.id, ") is not present in metadata column ", subset.by, call. = F)
    }

    subset <- SplitObject(object, split.by = subset.by)[[group.id]]
    barcodes <- colnames(subset)

    bcr.data <- NULL
    tcr.data <- NULL
    default.vdj.assay <- assay

    if (is.null(assay)) {
        bcr.data <- SubsetVDJData(subset, barcodes, assay = "BCR")
        tcr.data <- SubsetVDJData(subset, barcodes, assay = "TCR")
        default.vdj.assay <- DefaultAssayVDJ(object)
    } else if (assay == "BCR") {
        bcr.data <- SubsetVDJData(subset, barcodes, assay = assay)
    } else if (assay == "TCR") {
        tcr.data <- SubsetVDJData(subset, barcodes, assay = assay)
    } else {
        stop("Invalid assay: ", assay, call. = F)
    }

    subset <- ClearVDJ(subset)

    if (!is.null(bcr.data)) {
        subset <- AddVDJDataForAssay("BCR", subset, bcr.data[["vdj.primary"]], bcr.data[["vdj.secondary"]], bcr.data[["vj.primary"]], bcr.data[["vj.secondary"]])
    }

    if (!is.null(tcr.data)) {
        subset <- AddVDJDataForAssay("TCR", subset, tcr.data[["vdj.primary"]], tcr.data[["vdj.secondary"]], tcr.data[["vj.primary"]], tcr.data[["vj.secondary"]])
    }

    DefaultAssayVDJ(subset) <- default.vdj.assay

    return(subset)
}
