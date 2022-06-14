seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
seuratObj_BCR <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", quiet = T)

test_that("can split Seurat object without VDJ data", {
    expect_identical(SplitObject_VDJ(seuratObj, "seurat_clusters", quiet = T), Seurat::SplitObject(seuratObj, "seurat_clusters"))
})

test_that("can split Seurat object that include VDJ data", {
    split <- SplitObject_VDJ(seuratObj_BCR, "seurat_clusters", quiet = T)

    expect_length(split, 8)
    expect_identical(SplitObject_VDJ(seuratObj_BCR, "orig.ident")[[1]], seuratObj_BCR)
})


test_that("fails on wrong input", {
    expect_error(SplitObject_VDJ(seuratObj, "wrong.column"))
    expect_error(SubsetObject_VDJ(seuratObj_BCR, "wrong.column", "missing_id"))
    expect_error(SubsetObject_VDJ(seuratObj_BCR, "seurat_clusters", "missing_id"))
})
