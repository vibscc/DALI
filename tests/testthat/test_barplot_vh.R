seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
seuratObj <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC")

test_that("barplot_vh catches invalid parameters", {
    expect_error(barplot_vh(seuratObj, ident.1 = "foo"), regexp = "Invalid ident.1")
    expect_error(barplot_vh(seuratObj, ident.2 = "foo"), regexp = "Invalid ident.2")
    expect_error(barplot_vh(seuratObj, ident.2 = 1), regexp = "Can't specify ident.2")
    expect_error(barplot_vh(seuratObj, region = "foo"))
    expect_error(barplot_vh(seuratObj, group.by = "foobar"), regexp = "Invalid group.by")
    expect_error(barplot_vh(seuratObj, chain = "foo"))
})

test_that("barplot_vh runs without errors", {
    expect_silent(barplot_vh(seuratObj))
    expect_silent(barplot_vh(seuratObj, ident.1 = 3))
    expect_error(barplot_vh(seuratObj, ident.1 = 1), regexp = "Provided identities")
    expect_silent(barplot_vh(seuratObj, by.family = F))
})

test_that("barplot_vh can hide legend", {
    plot <- barplot_vh(seuratObj, legend = F)

    expect_equal(plot$theme$legend.position, "none")
})

test_that("barplot_vh can show grid", {
    plot <- barplot_vh(seuratObj, grid = T)

    expect_gt(length(plot$facet$params), 0)
})
