seuratObj <- readRDS("../testdata/seurat_objects/seuratObj_10x_sc5p_v2_hs_PBMC.rds")
VDJseuratObj <- Read10X_vdj(seuratObj, "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC", quiet = T)

# can merge data in objects
test_that("can merge seurat objects containing VDJ data", {
    expect_silent(suppressWarnings(MergeVDJ(VDJseuratObj,VDJseuratObj)))
    }
)

# can link directories and merge to object
test_that("can add VDJ data from directories to a merged seurat object", {
    metadata <- rep("ident1",times = length(colnames(seuratObj)))
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- colnames(seuratObj)
    seuratObj <- AddMetaData(seuratObj, metadata = metadata, col.name = "id")

    metadata <- rep("ident2",times = length(colnames(seuratObj)))
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- colnames(seuratObj)
    seuratObj2 <- AddMetaData(seuratObj, metadata = metadata, col.name = "id")

    expect_warning(mergedobj <- merge(seuratObj,seuratObj2))

    expect_silent(Read10X_MultiVDJ(mergedobj, id_column = "id", c("ident2" = "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC_Copy", "ident1" = "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC"), quiet =  T, force = T))
    }
)

# can detect too many dirs for the amount of available samples
test_that("gives an error when the #dirs is more than the #samples", {
    expect_error(Read10X_MultiVDJ(mergedobj, data.dir = "../testdata/empty_files", quiet = T))
    }
)

# gives warning when no airr
test_that("warns when there is no airr file", {
    metadata <- rep("ident1",times = length(colnames(seuratObj)))
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- colnames(seuratObj)
    seuratObj <- AddMetaData(seuratObj, metadata = metadata, col.name = "id")

    metadata <- rep("ident2",times = length(colnames(seuratObj)))
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- colnames(seuratObj)
    seuratObj2 <- AddMetaData(seuratObj, metadata = metadata, col.name = "id")

    expect_warning(mergedobj <- merge(seuratObj,seuratObj2))

    expect_warning(expect_warning(Read10X_MultiVDJ(mergedobj, id_column = "id", c("ident2" = "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC_Copy", "ident1" = "../testdata/cellranger_4.0.0/10x_sc5p_v2_hs_PBMC"), force  = T)))
    }
)

# can detect too similar samples
test_that("gives an error when sampels are too similar", {
    metadata <- rep("ident1",times = length(colnames(seuratObj)))
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- colnames(seuratObj)
    seuratObj <- AddMetaData(seuratObj, metadata = metadata, col.name = "id")

    metadata <- rep("ident2",times = length(colnames(seuratObj)))
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- colnames(seuratObj)
    seuratObj2 <- AddMetaData(seuratObj, metadata = metadata, col.name = "id")

    expect_warning(mergedobj <- merge(seuratObj,seuratObj2))

    expect_error(Read10X_MultiVDJ(mergedobj, data.dir = "../testdata/cellranger_4.0.0/", quiet = T))
    }
)

# can detect wrong id'd smaples (2 dirs for 1 sample)
test_that("gives an error when samples get linked to more than 1 directory", {
    expect_error(Read10X_MultiVDJ(mergedobj, data.dir = "../testdata/cellranger_4.0.0", quiet =  T, force = T))
    }
)

