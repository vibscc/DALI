#' Explore VDJ data in an interactive way
#' @export
interactive <- function() {
    app.directory <- system.file("shiny", package = "Diversity")
    if (app.directory == "") {
        stop("Could not shiny application directory. Try re-installing `Diversity`.", call. = FALSE)
    }

    shiny::runApp(app.directory, display.mode = "normal")
}
