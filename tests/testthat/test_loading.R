test_that("can load TCR data", {
    seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
    seuratObj <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC")

    metadata.columns <- colnames(seuratObj@meta.data)[7:ncol(seuratObj@meta.data)]

    expect_vector(metadata.columns, ptype = character(), size = 14)
    expect_setequal(metadata.columns, c("h.v_gene", "h.d_gene", "h.j_gene", "h.c_gene", "h.cdr3", "h.cdr3_nt", "h.V.fam", "l.v_gene", "l.d_gene", "l.j_gene", "l.c_gene", "l.cdr3", "l.cdr3_nt", "l.V.fam"))
})

