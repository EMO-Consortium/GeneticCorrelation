## The data <br>

## Setup <br>
* install hdf5, download from website and then make install
* install hdf5r R package <br>

Shell
```bash
export LD_LIBRARY_PATH="/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib:$LD_LIBRARY_PATH"
```
R
```R
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/software/Rpackage/4.0/")
dyn.load("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib/libhdf5_hl.so.310", lib.loc = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib")
install.packages("hdf5r",configure.args="--with-hdf5=/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/bin/h5cc")
```

## Data preprocessing
```R
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/software/Rpackage/4.0/")
library(Seurat)
library(tidyverse)
library(here)
library(data.table)
library(RLinuxModules)
library(loomR)
library(SeuratDisk)

setwd("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/CD_Ramik/CO_EPI")
### Read the data
cds <- readRDS("CO_EPI_ctrl.RDS")

### Super important step!!! Tranforms seurat_obj$RNA@meta.features to a table. 
### We need seurat_obj$RNA@meta.features to have columns in order for the .loom file to be generated! 
seurat_obj <- FindVariableFeatures(object = cds)   
pfile <- as.loom(x = seurat_obj, filename = "CO_EPI.loom", verbose = TRUE, overwrite = TRUE)
### Close the .loom file after all
pfile$close_all()

```
