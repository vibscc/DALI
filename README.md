# DALI

## Installation
To install the latest version of DALI, open R an install using `devtools`:

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("vibscc/DALI")
```

## Usage

### Starting the interactive shiny app
```
library("DALI")

seuratObj <- readRDS("<path/to/seurat_object.Rds>")

Interactive_VDJ(seuratObj)
```
If this seurat object does not have the VDJ data loaded yet using `Read10X_vdj()`, `Read10X_AIRR()` or `Read_AIRR()`, the application will prompt you to load in the data. Select the 10X cellranger output for either the BCR and/or TCR data linked to the same gene-expression data present in your seurat object. The app will then load in the data and start up.

### Loading 10X VDJ data in an existing seurat object
```
library("DALI")

seuratObj <- readRDS("<path/to/seurat_object.Rds>")
seuratObj <- Read10X_vdj(seuratObj, "<path/to/cellranger/bcr_or_tcr_out>", assay = "<assay>")

# <assay> can either be BCR or TCR
```

# Example data

Example data to use with this tool can be downloaded here: https://cloud.irc.ugent.be/public/index.php/s/9ys5czsaNtNQtSd

## FAQ

**Q**: What cellranger folder do I need to provide to load the vdj data?

**A**: If you used `cellranger vdj`, use the folder `outs` from the output. \
For `cellranger multi`, use the folder `vdj_b` (BCR) or `vdj_t` (TCR). This folder can be located in `outs/per_sample_outs/<sample>`.
