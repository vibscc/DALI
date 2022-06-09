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
To perform Trajectory analyses, installation of the Dynverse R packages and an installation of Docker is required.
Dynverse is installed using `devtools`:

```
devtools::install_github("dynverse/dyno")
```

## Usage

### Starting the interactive shiny app
```
library("DALI")

## Without a seurat object
## The application will allow you to upload all necessary data via the web browser

Interactive_DALI()

## With a seurat object (1)

seuratObj <- readRDS("<path/to/seurat_object.Rds>")
Interactive_DALI(seuratObj)
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

**Q**: How can I use DALI on merged samples?

**A**: Merging of the VDJ data (TCR and/or BCR) from multiple samples is currently not supported in DALI, but will come shortly! \
For now, this will require you to follow any of the following workflows:

1. If you used `cellranger multi` for the counts of the individual samples, you could aggregate the data into 1 datset using `cellranger aggr`. This will result in 1 count matrix and 1 set of VDJ related files which you can load into Seurat and DALI.
2. Manually concatenate the `filtered_contig_annotation.csv` or `all_contig_annotation.csv`. Make sure the clonotypes for each sample have a unique name! This can be done by simply appending `_<SAMPLE>` (replace `<SAMPLE>` with the actual name of your sample) to each clonotype definition (column `raw_clonotype_id`). This is necessary because clonotype 1 of sample 1 has no relation to clonotype 1 of sample 2. You could also define the clonotypes again on the combined dataset to find some common clonotypes between different samples. This could be done with tools like [immcantation](https://immcantation.readthedocs.io/en/stable/index.html)

**Q**: How can I process gamma/delta TCR data (gdTCR) using cellranger

**A**: Cellranger doesn't officially support this type of data (yet). There are however some instructions [here](https://kb.10xgenomics.com/hc/en-us/articles/360015793931-Can-I-detect-T-cells-with-delta-gamma-chains-in-my-V-D-J-data-) on how to work around the cellranger limitations.\
**NOTE**: Be careful if you want to load data into DALI processed thsi way if you already have 'normal' TCR (alpha + beta) data in your object. Make sure to change the naming of the V,D,J and C region back to its original name in the output files used by DALI: `airr_rearrangement.tsv` and `filtered/all_contig_annotations.csv`
