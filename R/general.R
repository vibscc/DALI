#' Load 10x VDJ data in a seurat object
#'
#' @param object Seurat object
#' @param data.dir Cellranger output directory
#'
#' @export

Read10X_vdj <- function(object, data.dir) {

    location.annotation.contig <- file.path(data.dir, 'filtered_contig_annotations.csv')

    if (!file.exists(location.annotation.contig)) {
        stop("Contig annotation file (", location.annotation.contig, ") is missing!")
    }

    annotation.contig <- read.csv(location.annotation.contig, stringsAsFactors = F) %>% subset(productive == 'True')

    columns <- c("barcode", "v_gene","d_gene","j_gene","c_gene", "cdr3","cdr3_nt")

    heavy <- annotation.contig %>%
        subset(grepl("^IGH", c_gene)) %>%
        dplyr::filter(!duplicated(barcode)) %>%
        select(all_of(columns)) %>%
        mutate(V.fam = get_v_families(v_gene)) %>%
        column_to_rownames('barcode') %>%
        rename_all(~ paste0("h.", .))


    light <- annotation.contig %>%
        subset(grepl("^IG[KL]", c_gene)) %>%
        dplyr::filter(!duplicated(barcode)) %>%
        select(all_of(columns)) %>%
        mutate(V.fam = get_v_families(v_gene)) %>%
        column_to_rownames('barcode') %>%
        rename_all(~ paste0("l.", .))


    object <- AddMetaData(object, heavy)
    object <- AddMetaData(object, light)

    return(object)
}

#' Extract V-family from a vector of v-genes
#'
#' @param v_genes Vector of genes

get_v_families <- function(v_genes) {
    v_families <- c()

    for(v_gene in v_genes) {
        if (!grepl("^IG[KLH]V", v_gene)) {
            v_families <- c(v_families, NA)
            next
        }

        family.full <- strsplit(v_gene, '-')[[1]][1]
        gene <- gsub("[0-9]", "", family.full)
        number <- gsub("[A-Za-z]", "", family.full)
        v_families <- c(v_families, paste0(gene, '-', number))
    }

    return(v_families)
}
