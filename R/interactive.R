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
