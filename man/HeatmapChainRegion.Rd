% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{HeatmapChainRegion}
\alias{HeatmapChainRegion}
\title{Heatmap showing to percentage of cells with a specifiek chain/region combination}
\usage{
HeatmapChainRegion(
  object,
  group.by = NULL,
  chain = c("VDJ", "VJ"),
  region = c("V", "J", "C"),
  by.family = T,
  add.missing.families = T,
  percent.total = T,
  show.missing.values = F,
  cluster.rows = F,
  cluster.cols = F,
  color.scheme = c("coolwarm", "viridis", "gray to blue", "turning red"),
  ...
)
}
\arguments{
\item{object}{Seurat object}

\item{group.by}{Metadata column to group the family data by.}

\item{chain}{Chain to plot. Options: 'H'(eavy), 'L'(ight) for BCR; 'A'(lpha), 'B'(eta) for TCR}

\item{region}{Region to plot. Available options: 'V'(ariable) or 'C'(onstant)}

\item{by.family}{Group genes of 1 family together. Default = TRUE}

\item{add.missing.families}{Should missing families be added to the plot. Default = TRUE}

\item{percent.total}{Should the fraction of cells be calculated from the total number or cells in the group or just the cells with VDJ info. Default = TRUE (= from total)}

\item{show.missing.values}{Should missing values be shown in the plot. Default = FALSE}

\item{cluster.rows}{Should rows (genes) be clustered in the heatmap. Default = FALSE}

\item{cluster.cols}{Should columns (groups) be clustered in the heatmap. Default = FALSE}

\item{color.scheme}{Colorscheme to use for the heatmap. Options: "coolwarm", "viridis". Default = "coolwarm"}

\item{...}{parameters to pass to pheatmap::pheatmap()}
}
\description{
Heatmap showing to percentage of cells with a specifiek chain/region combination
}
