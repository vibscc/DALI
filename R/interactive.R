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
    .GlobalEnv$.data.object.VDJ <- object
    on.exit(rm('.data.object.VDJ', envir = .GlobalEnv))

    options(shiny.maxRequestSize = max.upload.size*1024^2)
    shiny::runApp(app.directory, display.mode = "normal")
}
