seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
seuratObj_BCR <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", quiet = T)

test_that("can load data", {

    expect_s4_class(seuratObj_BCR, "Seurat")

    metadata.columns <- colnames(seuratObj_BCR@meta.data)[7:ncol(seuratObj_BCR@meta.data)]

    expect_vector(metadata.columns, ptype = character(), size = 21)
    expect_setequal(metadata.columns, c("vdj.v_gene", "vdj.d_gene", "vdj.j_gene", "vdj.c_gene", "vdj.cdr3", "vdj.cdr3_nt", "vdj.v_fam", "vdj.reads", "vdj.umis", "vdj.dual_IR", "vj.v_gene", "vj.d_gene", "vj.j_gene", "vj.c_gene", "vj.cdr3", "vj.cdr3_nt", "vj.v_fam", "vj.reads", "vj.umis", "vj.dual_IR", "clonotype"))

    expect_named(seuratObj_BCR@misc, c("VDJ", "default.chain.VDJ", "default.assay.VDJ"), ignore.order = T)
    expect_equal(seuratObj_BCR@misc$default.chain.VDJ, "primary")

    expect_true(nrow(seuratObj_BCR@misc$VDJ$BCR$vdj.primary) > 0)
    expect_true(nrow(seuratObj_BCR@misc$VDJ$BCR$vdj.secondary) > 0)
    expect_true(nrow(seuratObj_BCR@misc$VDJ$BCR$vj.primary) > 0)
    expect_true(nrow(seuratObj_BCR@misc$VDJ$BCR$vj.secondary) > 0)
})

test_that("shows warning for missing airr file", {
    expect_warning(Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "BCR"), regexp = ".*airr_rearrangement.tsv.*")
})


test_that("fails on invalid input", {
    expect_error(Read10X_vdj(seuratObj, "../testdata", quiet = T))
    expect_error(Read10X_vdj(seuratobj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "foo", quiet = T))
    expect_error(Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "TCR", quiet = T))
})

test_that("can force load data", {
    obj <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "TCR", force = T, quiet = T)

    expect_equal(obj@misc$default.assay.VDJ, "TCR")
    expect_equal(nrow(obj@misc$VDJ$TCR$vdj.primary), 0)
    expect_equal(nrow(obj@misc$VDJ$TCR$vdj.secondary), 0)
    expect_equal(nrow(obj@misc$VDJ$TCR$vj.primary), 0)
    expect_equal(nrow(obj@misc$VDJ$TCR$vj.secondary), 0)
})

seuratObj_TCR <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "TCR", force = T, quiet = T)

test_that("sets type correctly", {
    expect_equal(seuratObj_BCR@misc$default.assay.VDJ, "BCR")
    expect_named(seuratObj_BCR@misc$VDJ, c("BCR"))

    expect_equal(seuratObj_TCR@misc$default.assay.VDJ, "TCR")
    expect_named(seuratObj_TCR@misc$VDJ, c("TCR"))
})

seuratObj_double <- Read10X_vdj(seuratObj_BCR, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", type = "TCR", force = T, quiet = T)

test_that("can load multiple assays", {
    expect_named(seuratObj_double@misc$VDJ, c("BCR", "TCR"), ignore.order = T)
})
