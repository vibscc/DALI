#' Explore VDJ data in an interactive way
#'
#' @param object Seurat object
#'
#' @export

Interactive_VDJ <- function(object = NULL) {
    app.directory <- system.file("shiny", package = "Diversity")
    if (app.directory == "") {
        stop("Could not shiny application directory. Try re-installing `Diversity`.", call. = FALSE)
    }

    # Load object in global environment so shiny app can use it
    if (!is.null(object)) {
        if (IsValidSeuratObject(object)) {
            .GlobalEnv$.data.object.VDJ <- object
        } else {
            stop("The given object is not valid. Please make sure it contains VDJ information, loaded by Diversity::Read10X_vdj()", call. = F)
        }
    } else {
        .GlobalEnv$.data.object.VDJ <- NULL
    }
    on.exit(rm('.data.object.VDJ', envir = .GlobalEnv))

    shiny::runApp(app.directory, display.mode = "normal")
}

#' Nicer formatting for the most common dimensionality reductions
#'
#' @param reduction Dimensionality reduction

FormatDimred <- function(reduction) {
    for (pattern in c('pca', 'tsne', 'umap')) {
        if (grepl(pattern, reduction, ignore.case = T)) {
            replacement = switch(pattern,
                'pca' = 'PCA',
                'tsne' = 'tSNE',
                'umap' = 'UMAP',
                pattern
            )

            return(sub(pattern, replacement, reduction, ignore.case = T))
        }
    }

    return(reduction)
}

#' Create options list of chains for dropdown
#'
#' @param object Seurat object

AvailableChainsList <- function(object) {
    chains <- list()

    for (chain in AvailableChains(object)) {
        key <- switch(chain,
            "H" = "Heavy",
            "L" = "Light",
            "A" = "Alpha",
            "B" = "Beta"
        )
        chains[[key]] = chain
    }

    return(chains)
}

#' Get regions for given chain
#'
#' @param chain One of H, L, A or B

AvailableRegions <- function(chain) {
    heavy <- c("H", "A")
    light <- c("L", "B")

    if (chain %in% heavy) {
        return(c("V", "D", "J", "C"))
    } else if (chain %in% light) {
        return(c("V", "J", "C"))
    } else {
        return(c())
    }
}
