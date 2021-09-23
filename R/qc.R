#' Detect cells which have both TCR and BCR annotations and are therefore likely doublets
#'
#' @param object Seurat object
#' @param reduction Dimensionality reduction to use
#' @param do.plot Show a plot indicating the doublets
#' @param return.object Return the Seurat object with a column added with doublet annotation
#' @param metadata.column Metadata column to add with doublet information. Default = vdj.doublets
#'
#' @importFrom dplyr %>%
#' @export

VDJDoublets <- function(object, reduction = NULL, do.plot = T, return.object = T, metadata.column = "vdj.doublets") {
    if (!"TCR" %in% names(object@misc$VDJ) || !"BCR" %in% names(object@misc$VDJ)) {
        stop("Object should contain both TCR and BCR to detect cells with both TCR and BCR information", call. = F)
    }

    cells.BCR <- rownames(object@misc$VDJ$BCR$vdj.primary)[!is.na(object@misc$VDJ$BCR$vdj.primary$clonotype)]
    cells.TCR <- rownames(object@misc$VDJ$TCR$vdj.primary)[!is.na(object@misc$VDJ$TCR$vdj.primary$clonotype)]

    doublets <- intersect(cells.BCR, cells.TCR)

    object@meta.data$vdj.doublets <- rownames(object@meta.data) %in% doublets

    if (do.plot) {
        print(Seurat::DimPlot(object, group.by = metadata.column, reduction = reduction))
    }

    if (return.object) {
        return(object)
    }
}
