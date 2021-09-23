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
    UseMethod(generic = "DefaultAssayVDJ", object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value New VDJ assay
#'
#' @return An object with the new default VDJ assay
#'
#' @rdname DefaultAssayVDJ
#' @export DefaultAssayVDJ<-
#'
"DefaultAssayVDJ<-" <- function(object, ..., value) {
    UseMethod(generic = "DefaultAssayVDJ<-", object = object)
}

#' Get and set the default VDJ assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return The name of the default VDJ assay
#'
#' @rdname DefaultChainVDJ
#' @export DefaultChainVDJ
#'
DefaultChainVDJ <- function(object, ...) {
    UseMethod(generic = "DefaultChainVDJ", object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value New VDJ assay
#'
#' @return An object with the new default VDJ assay
#'
#' @rdname DefaultChainVDJ
#' @export DefaultChainVDJ<-
#'
"DefaultChainVDJ<-" <- function(object, ..., value) {
    UseMethod(generic = "DefaultChainVDJ<-", object = object)
}
