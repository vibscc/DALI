seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
VDJseuratObj <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", quiet = T)
split <- SplitObject_VDJ(VDJseuratObj,"seurat_clusters")
split2 <- SplitObject_VDJ(VDJseuratObj,"orig.ident")
subset1 <- SubsetObject_VDJ(VDJseuratObj, "seurat_clusters", "1", "BCR")
subset2 <- SubsetObject_VDJ(VDJseuratObj, "seurat_clusters", "1", "TCR")

#can split objects by meta.data
test_that("can split Seurat objects that include VDJ data", {
    expect_equal(length(split) , 8)
    expect_equal(length(split2), 1)
    }
)

#can detect when meta.data column /sample_id does not exist
test_that("can detect wrong meta-data column", {
    expect_error(SplitObject_VDJ(seuratObj, "wrong.column"))
    expect_error(SubsetObject_VDJ(VDJseuratObj, "wrong.column", "missing_id"))
    expect_error(SubsetObject_VDJ(VDJseuratObj, "seurat_clusters", "missing_id"))
    }
)

#can detect when no vdj data is present
test_that("Can detect splitting an object without vdj data",{
    expect_error(SplitObject(seuratObj, "Seurat_clusters"))
    }
)

#can subset object by id column and sample id
test_that("Can subset object by metadata column & sample_id",{
    expect_warning(SubsetObject_VDJ(seuratObj, "seurat_clusters", "1"))
    expect_equal(length(subset1), 1)
    }
)

#can chooce to keep just one assay
test_that("can choose to subset and keep 1 assay", {
    expect_equal(length(subset1@misc$VDJ), 1)
    expect_equal(length(subset2@misc$VDJ), 0)
    }
)
