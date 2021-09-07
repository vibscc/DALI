seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
seuratObj_BCR <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC")

test_that("can load data", {

    expect_s4_class(seuratObj_BCR, "Seurat")

    metadata.columns <- colnames(seuratObj_BCR@meta.data)[7:ncol(seuratObj_BCR@meta.data)]

    expect_vector(metadata.columns, ptype = character(), size = 21)
    expect_setequal(metadata.columns, c("h.v_gene", "h.d_gene", "h.j_gene", "h.c_gene", "h.cdr3", "h.cdr3_nt", "h.v_fam", "h.reads", "h.umis", "h.dual_IR", "l.v_gene", "l.d_gene", "l.j_gene", "l.c_gene", "l.cdr3", "l.cdr3_nt", "l.v_fam", "l.reads", "l.umis", "l.dual_IR", "clonotype"))

    expect_named(seuratObj_BCR@misc, c("VDJ", "default.chain.VDJ", "default.assay.VDJ"), ignore.order = T)
    expect_equal(seuratObj_BCR@misc$default.chain.VDJ, "primary")

    expect_true(nrow(seuratObj_BCR@misc$VDJ$BCR$heavy.primary) > 0)
    expect_true(nrow(seuratObj_BCR@misc$VDJ$BCR$heavy.secondary) > 0)
    expect_true(nrow(seuratObj_BCR@misc$VDJ$BCR$light.primary) > 0)
    expect_true(nrow(seuratObj_BCR@misc$VDJ$BCR$light.secondary) > 0)
})


test_that("fails on invalid input", {
    expect_error(Read10X_vdj(seuratObj, "../testdata"))
    expect_error(Read10X_vdj(seuratobj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "foo"))
    expect_error(Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "TCR"))
})

test_that("can force load data", {
    obj <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "TCR", force = T)

    expect_equal(obj@misc$default.assay.VDJ, "TCR")
    expect_equal(nrow(obj@misc$VDJ$TCR$heavy.primary), 0)
    expect_equal(nrow(obj@misc$VDJ$TCR$heavy.secondary), 0)
    expect_equal(nrow(obj@misc$VDJ$TCR$light.primary), 0)
    expect_equal(nrow(obj@misc$VDJ$TCR$light.secondary), 0)
})

seuratObj_TCR <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "TCR", force = T)

test_that("sets type correctly", {
    expect_equal(seuratObj_BCR@misc$default.assay.VDJ, "BCR")
    expect_named(seuratObj_BCR@misc$VDJ, c("BCR"))

    expect_equal(seuratObj_TCR@misc$default.assay.VDJ, "TCR")
    expect_named(seuratObj_TCR@misc$VDJ, c("TCR"))
})

seuratObj_double <- Read10X_vdj(seuratObj_BCR, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "TCR", force = T)

test_that("can load multiple assays", {
    expect_named(seuratObj_double@misc$VDJ, c("BCR", "TCR"), ignore.order = T)
})
