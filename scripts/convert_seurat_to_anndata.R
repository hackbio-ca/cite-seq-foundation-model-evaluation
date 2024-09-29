library(remotes)

# Install the required packages
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("pmbio/MuDataSeurat")

# Convert the h5seurat data object to MuData
library(MuDataSeurat)
library(SeuratDisk)

# Load the h5seurat object
multi <- SeuratDisk::LoadH5Seurat("data/multi.h5seurat")
WriteH5MU(multi, "data/multi.h5mu")