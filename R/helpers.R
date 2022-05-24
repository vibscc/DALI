#' Calculate frequencies of column and group
#'
#' @param object Seurat object
#' @param data.column Metadata column to use
#' @param group.by Metadata column with group info
#' @param show.missing Should missing values be included or not

CalculateFrequency <- function(object, data.column, group.by, show.missing) {
    if (data.column == group.by) {
        stop("Can not calculate frequence when data column and group by are identical")
    }

    out <- object@meta.data %>%
        group_by(.data[[data.column]], .data[[group.by]]) %>%
        summarise(freq = n()) %>%
        ungroup() %>%
        filter(case_when(!show.missing ~ !is.na(.data[[data.column]]),
                         T ~ T))
    return(out)
}

#' Convert clonotype to sequences
#'
#' @param object Seurat object
#' @param ct clonotype to convert
#' @param chain Which chain to get sequence for
#' @param clonotype.column Metadata column with clonotype information. Default = clonotype
#' @param sequence.type One of AA or NT. Default = AA
#' @param collapse Character to use to paste multiple sequence together

ClonotypeToSequence <- function(object, ct, chain = c("VDJ", "VJ"), clonotype.column = "clonotype", sequence.type = "AA", collapse = "\n") {
    chain <- match.arg(chain) %>% tolower()
    sequence.column <- paste0(chain, ".cdr3", if (tolower(sequence.type) == "nt") "_nt" else "")

    if (!sequence.column %in% colnames(object@meta.data)) {
        stop("Invalid metadata column ", sequence.column, call. = F)
    }

    data.filtered <- object@meta.data %>%
        filter(.data[[clonotype.column]] == ct)

    sequences <- data.filtered[, sequence.column] %>% unique() %>% na.omit()

    return(paste(sequences, collapse = collapse))
}

#' Get field name for given column name
#'
#' @param column Column for which you want the field
#' @param fields List of available fields
#' @param columns List of available columns

FieldForColumn <- function(column, fields, columns) {
    index <- which(columns == column[[1]])

    return(fields[[index]])
}

#' Generate ggplot default colors
#'
#' @param n Number of colors needed
#' @param h Minimum and maximum hue
#'
#' @importFrom grDevices hcl

ggplotColors <- function(n, h = c(15, 375)) {
    if ((diff(h) %% 360) < 1) {
        h[2] <- h[2] - 360 / n
    }
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

#' Get a colorpalette for given categorical data
#'
#' @param data Categorical data
#'
#' @importFrom dplyr %>%
#' @importFrom Polychrome sky.colors

GetCategoricalColorPalette <- function(data) {
    n.categories <- unique(data) %>% length()

    if (n.categories <= 24) {
        return(sky.colors(max(c(3, n.categories))) %>% unname())
    }

    return(ggplotColors(n = n.categories))
}

#' Get metadata column for given chain and region
#'
#' @param chain VDJ chain
#' @param region VDJ region
#' @param by.family Data grouped by family. Default = F

GetDataColumn <- function(chain, region, by.family = F) {
    chain <- tolower(chain)
    region <- tolower(region)

    data.column <- paste0(chain, ".", region, "_")

    if (by.family && region == "v") {
        data.column <- paste0(data.column, "fam")
    } else {
        data.column <- paste0(data.column, "gene")
    }
}

#' Add missing VDJ families to a list of families.
#'
#' @param families Families to complete
#'
#' @importFrom dplyr %>%

AddMissingVDJFamilies <- function(families) {
    families.completed <- c()
    prefixes <- gsub("/.*$", "", gsub("[0-9-]", "", families)) %>% unique()

    for (prefix in prefixes) {
        families.with.prefix <- families[grepl(prefix, families)]

        # Ignore all families that contain a / in the name
        # These families will just be appended to the final families without attempting to complete the missing families
        families.ignored <- families.with.prefix[grepl("/", families.with.prefix)]

        # Only keep families without / in the name
        families.with.prefix <- setdiff(families.with.prefix, families.ignored)

        family.numbers <- c()
        for (family in families.with.prefix) {
            family.numbers <- c(family.numbers, gsub("[A-Za-z-]", "", family)) %>% as.numeric()
        }

        if (length(family.numbers) > 0) {
            families.completed <- c(families.completed, paste0(prefix, "-", seq(1, max(family.numbers))))
        }

        families.completed <- c(families.completed, families.ignored)
    }

    return(families.completed)
}

#' Get the sequence of a v,d,j or c gene for a cell
#'
#' @param object Seurat object
#' @param cell Cell id
#' @param assay VDJ assay to use. Default = DefaultAssayVDJ(object)
#' @param chain Chain to use. Options = VDJ (=heavy/alpha); VJ (=light/beta)
#' @param regions Regions to get sequence from
#'
#' @importFrom dplyr %>%
#'
#' @export

GetSequence <- function(object, cell, assay = NULL, chain = c("VDJ", "VJ"), regions = NULL) {
    chain <- match.arg(chain) %>% tolower()

    if (is.null(assay)) {
        assay <- DefaultAssayVDJ(object)
    }

    if (!assay %in% names(slot(object, "misc")[["VDJ"]])) {
        stop("Invalid assay ", assay, call. = F)
    }

    if (is.null(regions)) {
        regions <- AvailableRegions(chain)
    } else {
        for (region in regions) {
            if (!region %in% AvailableRegions(chain)) {
                stop("Invalid region ", region, " for VDJ assay ", DefaultAssayVDJ(object), call. = F)
            }
        }
    }

    regions <- tolower(regions)

    df.name <- paste0(chain, ".", DefaultChainVDJ(object))
    column.name.sequence <- paste0(chain, ".sequence")
    data.df <- slot(object, "misc")[["VDJ"]][[assay]][[df.name]] %>% column_to_rownames("barcode")

    if (!cell %in% rownames(data.df)) {
        stop("Invalid cell ", cell, call. = F)
    }

    sequence <- ""
    for (region in regions) {
        column.name.start <- paste0(chain, ".", region, "_sequence_start")
        column.name.end <- paste0(chain, ".", region, "_sequence_end")

        start <- data.df[cell, column.name.start]
        end <- data.df[cell, column.name.end]
        sequence.full <- data.df[cell, column.name.sequence]

        if (is.na(start) | is.na(end) | is.na(sequence.full)) {
            stop("Could not extract requested sequence for given cell. Missing information for region ", region, ".", call. = F)
        }

        if (end < start) {
            stop("Failed to extract sequence. End of sequence is before start for region ", region, ".", call. = F)
        }

        if (end > nchar(sequence.full)) {
            stop("Failed to extract sequence. End is outside the sequence for region ", region, ".", call. = F)
        }
        sequence <- paste0(sequence, substr(sequence.full, start, end))
    }

    return(sequence)
}

#' Merge 2 data frames and overwrite common columns with values from the second df
#'
#' @param x first data frame
#' @param y second data frame
#' @param by column to join both data frames on

MergeAndOverwrite <- function(x, y, by) {

    x[, setdiff(colnames(y), colnames(x))] <- NA
    y[, setdiff(colnames(x), colnames(y))] <- NA

    for (i in 1:nrow(x)) {
        if (!x[i, by] %in% y[, by]) {
            y <- rbind(y, x[i,])
        }
    }

    rownames(y) <- NULL

    return(y)
}

#' Create new VDJ data from by merging original and new content together
#'
#' @param new New VDJ data
#' @param orig Original VDJ data
#'
#' @importFrom dplyr %>% na_if
#' @importFrom stringr str_match

CreateVDJData <- function(new, orig = NULL) {
    new$clonotype <- new$clonotype %>% na_if("None")

    if (is.null(orig)) {
        return(new)
    }

    clonotype.suffixes <- str_match(orig$clonotype, "_S[0-9]$")[,1] %>% na.omit() %>% unique()

    if (length(clonotype.suffixes) > 0) {
        clonotype.suffix <- paste0("_", clonotype.suffixes %>% gsub(pattern = "_S", replacement = "") %>% as.numeric() %>% max() + 1)
    } else {
        clonotype.suffix <- "_2"
        orig$clonotype <- vapply(orig$clonotype, function(x) { ifelse(is.na(x), x, paste0(x, "_1"))}, character(1))
    }

    new$clonotype <- vapply(new$clonotype, function(x) { ifelse(is.na(x), x, paste0(x, clonotype.suffix))}, character(1))

    return(MergeAndOverwrite(orig, new, by = "barcode"))
}

#' Return colorscales for heatmaps from name
#'
#' @param name Colorscheme name
#' @param n Colors in spectrum. Default = 100

ColorScale <- function(name = c("coolwarm", "viridis"), n = 100) {
    if (is.null(name) || name == "coolwarm") {
        return(colorRampPalette(c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"))(n))
    } else if (name == "viridis") {
        return(colorRampPalette(c("#440154", "#443A83", "#31688E", "#21908C", "#35B779", "#8FD744", "#FDE725"))(n))
    } else {
        stop("invalid colorscheme ", name)
    }
}

#' Give cellbarcodes and clonotypes from VDJ data an identifier
#'
#' @param object Seurat object
#' @param data list of vdj dataframes
#' @param assay specify assay
#' @param index object name/id

Uniqify_VDJ <- function(object, data, assay = c("BCR","TCR"), index) {

    counter <- 1
    new_data <- list()
    assay <- match.arg(assay)
    assay_data <- object@misc$VDJ[[assay]]

    for (vdj in assay_data) {
        if (length(vdj$barcode) > 0) {

            vdj$clonotype <- paste0(vdj$clonotype, "_", as.character(index))
            vdj$barcode <- paste0(vdj$barcode, "_", as.character(index))
        }
        if (counter == 1) {
            converted_data <- list( "vdjp" = vdj)
        } else if  (counter == 2) {
            converted_data <- list( "vdjs" = vdj)
        } else if (counter == 3) {
            converted_data <- list( "vjp" = vdj)
        } else {
            converted_data <- list( "vjs" = vdj)
        }
        new_data <- append(new_data, converted_data)
        counter <- counter + 1
    }

    data$vdj.primary <- rbind(new_data$vdjp, data$vdj.primary)
    data$vdj.secondary <- rbind(new_data$vdjs, data$vdj.secondary)
    data$vj.primary <- rbind(new_data$vjp, data$vj.primary)
    data$vj.secondary <- rbind(new_data$vjs, data$vj.secondary)

    return(data)
}

#' Create a subset of VDJ data using given cellbarcodes
#'
#' @param object Seurat object to be split
#' @param barcodes barcodes of data to be retained
#' @param assay BCR or TCR data

SubsetVDJData <- function(object, barcodes, assay = c("TCR","BCR")) {
    if (is.null(object@misc$VDJ[[assay]])) {
        return(NULL)
    }
    assay <- match.arg(assay)
    data <- object@misc$VDJ[[assay]]
    VDJ <- list("vdj.primary" = data.frame(),
                "vdj.secondary" = data.frame(),
                "vj.primary" = data.frame(),
                "vj.secondary" = data.frame())

    index <- 1
    for (df in data) {
        transfers <- intersect(df$barcode,barcodes)

        if (index == 1) {
            VDJ$vdj.primary <- df[df$barcode %in% transfers,]
            rownames(VDJ$vdj.primary) <- NULL
        } else if (index == 2) {
            VDJ$vdj.secondary <-  df[df$barcode %in% transfers,]
            rownames(VDJ$vdj.secondary) <- NULL
        } else if (index == 3) {
            VDJ$vj.primary <-  df[df$barcode %in% transfers,]
            rownames(VDJ$vj.primary) <- NULL
        } else {
            VDJ$vj.secondary <-  df[df$barcode %in% transfers,]
            rownames(VDJ$vj.secondary) <- NULL
        }
        index <- index + 1
    }
    return(VDJ)
}

AddVDJForAssay <- function(mergedobj, objectlist, assay = c("BCR","TCR")) {
    VDJ_data <- list("vdj.primary" = data.frame(),
                     "vdj.secondary" = data.frame(),
                     "vj.primary" = data.frame(),
                     "vj.secondary" = data.frame())
    obj_idx <- 1
    exist <- F
    assay <- match.arg(assay)
    for (object in objectlist) {
        if (!is.null(object@misc$VDJ[[assay]])) {
            exist <- T
            VDJ_data <- Uniqify_VDJ(object, VDJ_data, assay, obj_idx)
        }
        obj_idx <- obj_idx + 1
    }
    mergedobj@misc$VDJ[[assay]] <- VDJ_data

    if (exist) {
        Misc(object = mergedobj, slot = "default.assay.VDJ") <- assay
        Misc(object = mergedobj, slot = "default.chain.VDJ") <- DefaultChainVDJ(object)
        DefaultAssayVDJ(mergedobj) <- assay
    }
    return(mergedobj)
}

#' Adds VDJ data from lists to the misc slot
#'
#' @param object Seurat object
#' @param TCR list with TCR data
#' @param BCR list with BCR data

AddMiscVDJData <- function(object,TCR = NULL ,BCR = NULL) {
    if (!is.null(BCR) & !is.null(TCR)) {
        Misc(object = object, slot = "VDJ") <- list("TCR" = TCR, "BCR" = BCR)
        Misc(object = object, slot = "default.assay.VDJ") <- "TCR"
        Misc(object = object, slot = "default.chain.VDJ") <- DefaultChainVDJ(object)
        DefaultAssayVDJ(object) <- "TCR"
    } else if (!is.null(BCR)) {
        Misc(object = object, slot = "VDJ") <- list("BCR" = BCR)
        Misc(object = object, slot = "default.assay.VDJ") <- "BCR"
        Misc(object = object, slot = "default.chain.VDJ") <- DefaultChainVDJ(object)
        DefaultAssayVDJ(object) <- "BCR"
    } else if (!is.null(TCR)) {
        Misc(object = object, slot = "VDJ") <- list("TCR" = BCR)
        Misc(object = object, slot = "default.assay.VDJ") <- "TCR"
        Misc(object = object, slot = "default.chain.VDJ") <- DefaultChainVDJ(object)
        DefaultAssayVDJ(object) <- "TCR"
    }
    return(object)
}

#' Get The sequence columns from an airr_rearrangement.tsv file
#'
#' @param data.dir path to directory containing the airr_rearrangement.tsv file
#' @param quiet Ignore warnings. Default = FALSE

GetAIRRSequenceColumns <- function(data.dir, quiet = F) {
    location.airr.rearrangement <- file.path(data.dir, "airr_rearrangement.tsv")

    if (file.exists(location.airr.rearrangement)) {
        airr.data <- read.csv(location.airr.rearrangement, sep = "\t")
        colnames(airr.data) <- gsub("sequence_id", "contig_id", colnames(airr.data))
        sequence.columns <- grep("sequence", colnames(airr.data), value = T)
    } else {
        if (!quiet) {
            warning("Could not find airr_rearrangement.tsv. Sequence information will not be loaded and some functionality for BCR lineage tracing will not be available", call. = F)
        }
    }
    return(sequence.columns)
}

#' Construct a dataframe with vdj data for a seurat object
#'
#' @param data.dir directory containing the vdj data
#' @param sequence.columns sequence columns from an AIRR file
#' @param use.filtered Load filtered contig annotation. Default = TRUE

GetVDJ_Dataframe <- function(data.dir, sequence.columns, use.filtered) {
    location.annotation.contig <- file.path(data.dir, paste0(if (use.filtered) "filtered" else "all", "_contig_annotations.csv"))
    location.airr.rearrangement <- file.path(data.dir, "airr_rearrangement.tsv")

    if (!file.exists(location.annotation.contig)) {
        stop("Contig annotation file (", location.annotation.contig, ") is missing!", call. = F)
    }
    vdj_df <- read.csv(location.annotation.contig, stringsAsFactors = F) %>%
        filter(grepl("true", .data$productive, ignore.case = T))
    airr.data <- read.csv(location.airr.rearrangement, sep = "\t")
    vdj_df <- left_join(vdj_df, airr.data[, c("contig_id", sequence.columns)], by = "contig_id")

    return(vdj_df)
}

