# DALI

## How to cite

Please cite as: Verstaen K, Lammens I, Roels J, Saeys Y, Lambrecht BN, Vandamme N, Vanhee S. DALI (Diversity AnaLysis Interface): a novel tool for the integrated analysis of multimodal single cell RNAseq data and immune receptor profiling. bioRxiv 2021.12.07.471549; doi: https://doi.org/10.1101/2021.12.07.471549

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

## Without a seurat object
## The application will allow you to upload all necessary data via the web browser

Interactive_VDJ()

## With a seurat object (1)

seuratObj <- readRDS("<path/to/seurat_object.Rds>")

Interactive_VDJ(seuratObj)
```
(1) If this seurat object does not have the VDJ data loaded yet using `Read10X_vdj()`, `Read10X_AIRR()` or `Read_AIRR()`, the application will prompt you to load in the data. Select the 10X cellranger output for either the BCR and/or TCR data linked to the same gene-expression data present in your seurat object. The app will then load in the data and start up.

### Loading 10X VDJ data in an existing seurat object
```
library("DALI")

seuratObj <- readRDS("<path/to/seurat_object.Rds>")
seuratObj <- Read10X_vdj(seuratObj, "<path/to/cellranger/bcr_or_tcr_out>", assay = "<assay>")

# <assay> can either be BCR or TCR
```

# Example data

Example data to use with this tool can be downloaded here: \
https://cloud.irc.ugent.be/public/index.php/s/9ys5czsaNtNQtSd

## FAQ

**Q**: What cellranger folder do I need to provide to load the vdj data?

**A**: If you used `cellranger vdj`, use the folder `outs` from the output. \
For `cellranger multi`, use the folder `vdj_b` (BCR) or `vdj_t` (TCR). This folder can be located in `outs/per_sample_outs/<sample>`.
