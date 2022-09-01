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

ColorScale <- function(name = c("coolwarm", "viridis", "gray to blue", "turning red"), n = 100) {
    name <- match.arg(name)

    if (name == "coolwarm") {
        return(colorRampPalette(c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"))(n))
    } else if (name == "viridis") {
        return(colorRampPalette(c("#440154", "#443A83", "#31688E", "#21908C", "#35B779", "#8FD744", "#FDE725"))(n))
    } else if (name == "gray to blue") {
        return(colorRampPalette(c("#c0c0c0", "#c0c0c0", "#c0c0c0", "#91BFDB", "#4575B4"))(n))
    } else if (name == "turning red") {
        return(colorRampPalette(c("#c0c0c0", "#c0c0c0", "#c0c0c0", "#FC8D59", "#D73027"))(n))
    }
}

#' Give cell barcodes and clonotypes from VDJ data an identifier to allow merging of VDJ data from multiple samples
#'
#' @param object Seurat object
#' @param sample.suffix Suffix to add to the barcodes/clonotypes to make them unique
#' @param assay VDJ assay to make unique.
#' @param columns.to.update List of columns that should be suffixed with sample.suffix. Default = c("clonotype", "barcode")

Uniqify_VDJ <- function(object, sample.suffix, assay = c("BCR","TCR"), columns.to.update = c("clonotype", "barcode")) {
    assay <- match.arg(assay)
    vdj.data <- object@misc$VDJ[[assay]]

    for (data.type in names(vdj.data)) {
        for (column in columns.to.update) {
            if (column %in% colnames(vdj.data[[data.type]])) {
                vdj.data[[data.type]][[column]] <- paste0(vdj.data[[data.type]][[column]], "_", sample.suffix)
            }
        }
    }

    return(vdj.data)
}

#' Create a subset of VDJ data using given cellbarcodes
#'
#' @param object Seurat object to be split
#' @param barcodes barcodes of data to be retained
#' @param assay BCR or TCR data
#'
#' @importFrom dplyr filter

SubsetVDJData <- function(object, barcodes, assay = c("TCR","BCR")) {
    if (is.null(object@misc$VDJ[[assay]])) {
        return(NULL)
    }

    assay <- match.arg(assay)
    data <- slot(object, "misc")[["VDJ"]][[assay]]

    vdj.data <- list(
        "vdj.primary" = data.frame(),
        "vdj.secondary" = data.frame(),
        "vj.primary" = data.frame(),
        "vj.secondary" = data.frame()
    )

    for (data.type in names(vdj.data)) {
        barcodes.to.transfer <- intersect(data[[data.type]]$barcode, barcodes)
        vdj.data[[data.type]] <- data[[data.type]] %>% filter(.data$barcode %in% barcodes.to.transfer)
    }

    if (sum(unlist(lapply(vdj.data, nrow))) == 0) {
        return(NULL)
    }

    return(vdj.data)
}

#' Append VDJ data to an object
#'
#' @param object Seurat object to append the VDJ data to
#' @param vdj.data List with 4 dataframes (vdj.primary, vdj.secondary, vj.primary, vj.secondary) to append
#' @param assay VDJ assay to add this new data to

AppendVDJData <- function(object, vdj.data, assay = c("BCR", "TCR")) {
    assay <- match.arg(assay)

    if (!assay %in% names(object@misc$VDJ)) {
        object@misc$VDJ[[assay]] <- list(
            "vdj.primary" = data.frame(),
            "vdj.secondary" = data.frame(),
            "vj.primary" = data.frame(),
            "vj.secondary" = data.frame()
        )
    }

    object@misc$VDJ[[assay]][["vdj.primary"]] <- rbind(object@misc$VDJ[[assay]][["vdj.primary"]], vdj.data[["vdj.primary"]])
    object@misc$VDJ[[assay]][["vdj.secondary"]] <- rbind(object@misc$VDJ[[assay]][["vdj.secondary"]], vdj.data[["vdj.secondary"]])
    object@misc$VDJ[[assay]][["vj.primary"]] <- rbind(object@misc$VDJ[[assay]][["vj.primary"]], vdj.data[["vj.primary"]])
    object@misc$VDJ[[assay]][["vj.secondary"]] <- rbind(object@misc$VDJ[[assay]][["vj.secondary"]], vdj.data[["vj.secondary"]])

    return(object)
}

#' Remove the VDJ data from a seurat object
#'
#' @param object Seurat object
#' @param assay VDJ assay to remove. Default = NULL (=> removes everything)

ClearVDJ <- function(object, assay = NULL) {
    if (!is.null(assay) && !assay %in% c("BCR", "TCR")) {
        stop("Invalid assay ", assay, call. = F)
    }

    # Assay not in object, so we return the original object
    if (!is.null(assay) && !assay %in% names(slot(object, "misc"))) {
        return(object)
    }

    # TODO: remove the metadata columns generated by DALI
    if (is.null(assay) || length(names(slot(object, "misc"))) == 1) { # We remove everything
        slot(object, "misc")[["VDJ"]] <- NULL
        slot(object, "misc")[["default.chain.VDJ"]] <- NULL
        slot(object, "misc")[["default.assay.VDJ"]] <- NULL
    } else {
        slot(object, "misc")[["VDJ"]][[assay]] <- NULL
        DefaultAssayVDJ(object) <- names(slot(object, "misc"))[1]
    }

    return(object)
}

#' Get The sequence columns from an airr_rearrangement.tsv file
#'
#' @param data.dir path to directory containing the airr_rearrangement.tsv file
#' @param quiet Ignore warnings. Default = FALSE

GetAIRRSequenceColumns <- function(data.dir, quiet = F) {
    if (!dir.exists(data.dir)) {
        stop("Invalid data directory: ", data.dir, call. = F)
    }

    location.airr.rearrangement <- file.path(data.dir, "airr_rearrangement.tsv")

    if (!file.exists(location.airr.rearrangement)) {
        if (!quiet) {
            warning("Could not find airr_rearrangement.tsv. Sequence information will not be loaded and some functionality for BCR lineage tracing will not be available", call. = F)
        }
        return(NULL)
    }

    airr.data <- read.csv(location.airr.rearrangement, sep = "\t")
    colnames(airr.data) <- gsub("sequence_id", "contig_id", colnames(airr.data))

    return(grep("sequence", colnames(airr.data), value = T))
}

#' Construct a dataframe with vdj data for a seurat object
#'
#' @param data.dir directory containing the vdj data
#' @param sequence.columns sequence columns from an AIRR file
#' @param use.filtered Load filtered contig annotation. Default = TRUE

GetVDJ_Dataframe <- function(data.dir, sequence.columns, use.filtered = T) {
    if (!dir.exists(data.dir)) {
        stop("Invalid data directory: ", data.dir, call. = F)
    }

    location.annotation.contig <- file.path(data.dir, paste0(if (use.filtered) "filtered" else "all", "_contig_annotations.csv"))
    location.airr.rearrangement <- file.path(data.dir, "airr_rearrangement.tsv")

    if (!file.exists(location.annotation.contig)) {
        stop("Contig annotation file (", location.annotation.contig, ") is missing!", call. = F)
    }

    vdj_df <- read.csv(location.annotation.contig, stringsAsFactors = F) %>%
        filter(grepl("true", .data$productive, ignore.case = T))

    if (file.exists(location.airr.rearrangement)) {
        airr.data <- read.csv(location.airr.rearrangement, sep = "\t")
        colnames(airr.data) <- gsub("sequence_id", "contig_id", colnames(airr.data))
        vdj_df <- left_join(vdj_df, airr.data[, c("contig_id", sequence.columns)], by = "contig_id")
    }

    return(vdj_df)
}

#' Returns vector of colors based on the selected color theme
#'
#' @param data Data to get categories from
#' @param theme Colorscheme to use. Default = DALI
#'
#' @importFrom Polychrome sky.colors
#'
#' @export

GetCategoricalColorPalette <- function(data, theme = ColorThemes()) {
    theme <- match.arg(theme)

    n <- unique(data) %>% length()

    if (theme == "DALI") {
        if (n <= 24) {
            return(sky.colors(max(c(3, n))) %>% rev() %>% unname())
        }

        return(ggplotColors(n = n))
    } else if (theme == "DALII") {
        if (n <= 15) {
            cols <- c("#FF932A" ,"#3CB0B5", "#5D3484",
                      "#54DAFF", "#FFCE54", "#E484ED",
                       "#0AFA1F", "#FF0000", "#967ADC",
                       "#0014CF", "#FF0CCD", "#0072B2",
                       "#2E8F26", "#02FFB9", "#B85100")
            return(cols[1:max(3, n)])
        }
        n.5 <- n/2
        cols <- ggplotColors(n.5, c(290, 180))
        cols <- c(cols, ggplotColors(n.5, c(170,45)))

        return(cols)
    } else if (theme == "Pastel") {
        if (n <= 10) {
            cols <- c("#78FFF1", "#FEB6DB", "#79CAFF",
                      "#F6F7A3", "#DAB6FE", "#C1DEC8",
                      "#F7D3BB", "#76D7D6",  "#D4F9E5",
                      "#F7E368")
            return(cols[1:max(3, n)])
        }

        return(hcl(h = (seq(0, 360, length = n)), c = 55, l = 70))
    } else if (theme == "Colorblind") {
        if (n <= 12) {
            cols <- c("#AA44AA", "#D55E00", "#44AAAA",
                      "#24FF24", "#FFB6DB", "#0072B2",
                      "#E69F00", "#FFFF99", "56B4E9",
                      "#000000", "#FFF600", "#FF007A")

            return(cols[1:max(3, n)])
        }

        return(ggplotColors(n = n))
    } else if (theme == "Spectrum") {
        return(ggplotColors(n = n) %>% rev())
    } else {
        stop("Invalid theme: ", theme)
    }
}

#' Get a list of the available color themes
#'
#' @export

ColorThemes <- function() {
    return(c("DALI", "DALII", "Pastel", "Colorblind", "Spectrum"))
}

#' Return table X, with the same rows as table Y
#'
#' @param plot_data_x table to be transformed
#' @param plot_data_y reference table

EqualiseTableRows <- function(plot_data_x, plot_data_y) {
    diff <- nrow(plot_data_y) - nrow(plot_data_x)
    if (diff > 0) {
        new.data <- matrix(rep(0, diff * ncol(plot_data_x)), nrow = diff)
        new.rownames <- c()
        for (row in rownames(plot_data_y)) {
            if (!(row %in% rownames(plot_data_x))) {
                new.rownames <- c(new.rownames, row)
            }
        }
        colnames(new.data) <- colnames(plot_data_x)
        rownames(new.data) <- new.rownames
        plot_data_x <- rbind(plot_data_x, new.data)
    }
    return(plot_data_x)
}

#' Helperfunction to subset BCR or TCR data in a Seurat object
#'
#' @param object Seurat object to be subsetted
#' @param assay BCR or TCR assay
#'
#' @importFrom dplyr %>%

SubsetAssay <- function(object, assay = c('TCR', 'BCR')) {
    assay <- match.arg(assay)

    barcodes <- c()
    for (df in object@misc$VDJ[[assay]]) {
        barcodes <- c(barcodes, df$barcode) %>% unique()
    }
    return(object[,barcodes])
}
