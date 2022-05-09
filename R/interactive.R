#' Explore VDJ data in an interactive way
#'
#' @param object Seurat object
#' @param vdj Boolean to specify if object contains VDJ data. Default = True
#'
#' @export

Interactive_VDJ <- function(object = NULL, vdj = T) {
    app.directory <- system.file("shiny", package = "DALI")

    if (app.directory == "") {
        stop("Could not shiny application directory. Try re-installing `DALI`.", call. = FALSE)
    }

    # Load object in global environment so shiny app can use it
    if (!is.null(object)) {
        if (IsValidSeuratObject(object) | vdj) {
            .GlobalEnv$.data.object.VDJ <- object
        } else {
            stop("The given object is not valid. Please make sure it contains VDJ information, loaded by DALI::Read10X_vdj()", call. = F)
        }
    } else {
        .GlobalEnv$.data.object.VDJ <- NULL
    }
    .GlobalEnv$.no.vdj <- vdj
    .GlobalEnv$.loaded_data = TRUE

    on.exit(rm('.data.object.VDJ', envir = .GlobalEnv))
    on.exit(rm('.no.vdj', envir = .GlobalEnv))
    on.exit(rm('.loaded_data', envir = .GlobalEnv))
    shiny::runApp(app.directory, display.mode = "normal")
}

#' Nicer formatting for the most common dimensionality reductions
#'
#' @param reduction Dimensionality reduction

FormatDimred <- function(reduction) {
    for (pattern in c("pca", "tsne", "umap")) {
        if (grepl(pattern, reduction, ignore.case = T)) {
            replacement <- switch(pattern,
                "pca" = "PCA",
                "tsne" = "tSNE",
                "umap" = "UMAP",
                pattern
            )

            return(sub(pattern, replacement, reduction, ignore.case = T))
        }
    }

    return(reduction)
}

#' Get regions for given chain
#'
#' @param chain VDJ (heavy/alpha) OR VJ (light/beta)

AvailableRegions <- function(chain) {
    chain <- tolower(chain)
    if (chain == "vdj") {
        return(c("V", "D", "J", "C"))
    } else if (chain == "vj") {
        return(c("V", "J", "C"))
    } else {
        return(c())
    }
}
