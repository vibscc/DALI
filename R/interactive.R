#' Explore VDJ data in an interactive way
#'
#' @param object Seurat object
#' @param max.upload.size Maximum upload size for shiny application (in MB). Default = 1000
#'
#' @export

interactive_VDJ <- function(object = NULL, max.upload.size = 1000) {
    app.directory <- system.file("shiny", package = "Diversity")
    if (app.directory == "") {
        stop("Could not shiny application directory. Try re-installing `Diversity`.", call. = FALSE)
    }

    # Load object in global environment so shiny app can use it
    .GlobalEnv$.data.object.VDJ <- Seurat::DietSeurat(object, counts = F, data = T, scale.data = F, assays = NULL, dimreducs = names(object@reductions), graphs = NULL)
    on.exit(rm('.data.object.VDJ', envir = .GlobalEnv))

    options(shiny.maxRequestSize = max.upload.size*1024^2)
    shiny::runApp(app.directory, display.mode = "normal")
}

#' Nicer formatting for the most common dimensionality reductions
#'
#' @param reduction Dimensionality reduction

formatDimred <- function(reduction) {
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
