#' Get and set the default VDJ assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return The name of the default VDJ assay
#'
#' @rdname DefaultAssayVDJ
#' @export DefaultAssayVDJ
#'
DefaultAssayVDJ <- function(object, ...) {
    UseMethod(generic = 'DefaultAssayVDJ', object = object)
}

#' @inheritParams DefaultAssayVDJ
#' @param value Name of assay to set as default
#'
#' @return An object with the new default VDJ assay
#'
#' @rdname DefaultAssayVDJ
#' @export DefaultAssayVDJ<-
#'
"DefaultAssayVDJ<-" <- function(object, ..., value) {
    UseMethod(generic = 'DefaultAssayVDJ<-', object = object)
}
