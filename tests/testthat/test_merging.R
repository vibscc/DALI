seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
seuratObj_BCR <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", quiet = T)

test_that("can merge objects without VDJ data", {
    # Suppress warning 'Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.'
    suppressWarnings(
        expect_identical(MergeVDJ(seuratObj, seuratObj), merge(seuratObj, seuratObj))
    )
})

test_that("can merge seurat objects containing VDJ data", {
    options(warn = -1)
    merged <- MergeVDJ(seuratObj_BCR, seuratObj_BCR)
    options(warn = 0)

    expect_s4_class(merged, "Seurat")
    expect_equal(DefaultAssayVDJ(merged), "BCR")
    expect_equal(DefaultChainVDJ(merged), "primary")

    expect_equal(nrow(slot(merged, "misc")[["VDJ"]][["BCR"]][["vdj.primary"]]), nrow(slot(seuratObj_BCR, "misc")[["VDJ"]][["BCR"]][["vdj.primary"]]) * 2)
    expect_equal(ncol(slot(merged, "misc")[["VDJ"]][["BCR"]][["vdj.primary"]]), ncol(slot(seuratObj_BCR, "misc")[["VDJ"]][["BCR"]][["vdj.primary"]]))
    expect_equal(nrow(slot(merged, "misc")[["VDJ"]][["BCR"]][["vdj.secondary"]]), nrow(slot(seuratObj_BCR, "misc")[["VDJ"]][["BCR"]][["vdj.secondary"]]) * 2)
    expect_equal(ncol(slot(merged, "misc")[["VDJ"]][["BCR"]][["vdj.secondary"]]), ncol(slot(seuratObj_BCR, "misc")[["VDJ"]][["BCR"]][["vdj.secondary"]]))
    expect_equal(nrow(slot(merged, "misc")[["VDJ"]][["BCR"]][["vj.primary"]]), nrow(slot(seuratObj_BCR, "misc")[["VDJ"]][["BCR"]][["vj.primary"]]) * 2)
    expect_equal(ncol(slot(merged, "misc")[["VDJ"]][["BCR"]][["vj.primary"]]), ncol(slot(seuratObj_BCR, "misc")[["VDJ"]][["BCR"]][["vj.primary"]]))
    expect_equal(nrow(slot(merged, "misc")[["VDJ"]][["BCR"]][["vj.secondary"]]), nrow(slot(seuratObj_BCR, "misc")[["VDJ"]][["BCR"]][["vj.secondary"]]) * 2)
    expect_equal(ncol(slot(merged, "misc")[["VDJ"]][["BCR"]][["vj.secondary"]]), ncol(slot(seuratObj_BCR, "misc")[["VDJ"]][["BCR"]][["vj.secondary"]]))
})
