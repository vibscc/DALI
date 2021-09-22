seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
seuratObj <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", quiet = T)

test_that("BarplotChainRegion catches invalid parameters", {
    expect_error(BarplotChainRegion(seuratObj, ident.1 = "foo"), regexp = "Invalid ident.1")
    expect_error(BarplotChainRegion(seuratObj, ident.2 = "foo"), regexp = "Invalid ident.2")
    expect_error(BarplotChainRegion(seuratObj, ident.2 = 1), regexp = "Can't specify ident.2")
    expect_error(BarplotChainRegion(seuratObj, region = "foo"))
    expect_error(BarplotChainRegion(seuratObj, group.by = "foobar"), regexp = "Invalid group.by")
    expect_error(BarplotChainRegion(seuratObj, chain = "foo"))
})

test_that("BarplotChainRegion runs without errors", {
    expect_silent(BarplotChainRegion(seuratObj))
    expect_silent(BarplotChainRegion(seuratObj, ident.1 = 3))
    expect_error(BarplotChainRegion(seuratObj, ident.1 = 1), regexp = "Provided identities")
    expect_silent(BarplotChainRegion(seuratObj, by.family = F))
})

test_that("BarplotChainRegion can hide legend", {
    plot <- BarplotChainRegion(seuratObj, legend = F)

    expect_equal(plot$theme$legend.position, "none")
})

test_that("BarplotChainRegion can show grid", {
    plot <- BarplotChainRegion(seuratObj, grid = T)

    expect_gt(length(plot$facet$params), 0)
})
