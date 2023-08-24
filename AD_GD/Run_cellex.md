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
remotes::install_github("mojaveazure/seurat-disk")

```

## Data preprocessing
```R
### R/4.2.2-foss-2022a-bare
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/software/Rpackage/4.0/")
dyn.load("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib/libhdf5_hl.so.310", lib.loc = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib")
library(Seurat)
library(tidyverse)
library(here)
library(data.table)
library(RLinuxModules)
library(loomR)
library(SeuratDisk)
library(R.utils)

setwd("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/CD_Ramik/CO_EPI")
### generate seurat object

## load metadata
metadata<-read.table("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/CD_Ramik/scp_metadata_combined.v2.txt",header = T)
metadata<-metadata[-1,]

## load data to seurat object
## make sure in the right folder
bdata <- Read10X(data.dir = filepath,gene.column = 1)
pbmc <- CreateSeuratObject(counts = bdata, project = "10x", min.cells = 3, min.features = 200)

## add metadata
mt<-metadata[match(rownames(pbmc@meta.data),metadata$NAME),]
pbmc@meta.data <-cbind(pbmc@meta.data,mt[,c("Celltype","disease__ontology_label")])

## subset and extract data
ctrl<-subset(pbmc, subset=disease__ontology_label=="normal")

### Read the data
# ctrl <- readRDS("CO_EPI_ctrl.RDS")

### Super important step!!! Tranforms seurat_obj$RNA@meta.features to a table. 
### We need seurat_obj$RNA@meta.features to have columns in order for the .loom file to be generated! 
ctrl <- FindVariableFeatures(object = ctrl)   
pfile <- as.loom(x = ctrl, filename = "CO_EPI.loom", verbose = TRUE, overwrite = TRUE)
### Close the .loom file after all
pfile$close_all()

```

* Additional preprocessing for Lung cells
```R
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/software/Rpackage/4.0/")
dyn.load("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib/libhdf5_hl.so.310", lib.loc = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib")
library(Seurat)
library(tidyverse)
library(here)
library(data.table)
library(RLinuxModules)
library(loomR)
library(SeuratDisk)
library(R.utils)

setwd("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/Lung_SciAdv")
ctrl <- readRDS("Ctrl.seurat.RDS")

## remove outlier cells
ctrl@meta.data$outlier<-F
ctrl@meta.data$outlier[grep("Outlier",ctrl@meta.data$Subclass_Cell_Identity)]<-T
ctrl.clean<-subset(ctrl, subset=outlier==F)

## Remove multiplet
ctrl.clean<-subset(ctrl.clean, subset=CellType_Category!="Multiplet")

ctrl.clean <- FindVariableFeatures(object = ctrl.clean)
pfile <- as.loom(x = ctrl.clean, filename = "Lung.loom", verbose = TRUE, overwrite = TRUE)
pfile$close_all()

## myeloid
sub<-subset(ctrl.clean, subset=CellType_Category=="Myeloid")
sub <- FindVariableFeatures(object = sub)
pfile <- as.loom(x = sub, filename = "Myeloid.loom", verbose = TRUE, overwrite = TRUE)
pfile$close_all()

other<-subset(ctrl.clean, subset=CellType_Category!="Myeloid")
other <- FindVariableFeatures(object = other)
pfile <- as.loom(x = other, filename = "Lung_other.loom", verbose = TRUE, overwrite = TRUE)
pfile$close_all()

```

* Additional preprocessing for blood cells
```R
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/software/Rpackage/4.0/")
dyn.load("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib/libhdf5_hl.so.310", lib.loc = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib")

library(Seurat)
library(tidyverse)
library(here)
library(data.table)
library(RLinuxModules)
library(loomR)
library(SeuratDisk)
library(R.utils)

setwd("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/PBMC")

## dataset1: covid_blish: https://www.covid19cellatlas.org/index.patient.html (this one:Peripheral Blood Mononuclear Cells (PBMCs), Blish lab)
## paper:https://www.nature.com/articles/s41591-020-0944-y
## can also check data from this paper: https://rupress.org/jem/article/218/8/e20210582/212379/Multi-omic-profiling-reveals-widespread
pbmc<-readRDS("./covid_blish/blish_covid.seu.rds")
ctrl<-subset(pbmc, subset=Status=="Healthy")

ctrl <- FindVariableFeatures(object = ctrl)
pfile <- as.loom(x = ctrl, filename = "PBMC_blish.loom", verbose = TRUE, overwrite = TRUE)
pfile$close_all()

## dataset2: onek1k
## only keep the young-mid adult, still to large
ctrl<-subset(pbmc, subset=age_group<65)

## another bigger dataset: The covid blood scRNAseq data was download from:https://data.humancellatlas.org/explore/projects/b963bd4b-4bc1-4404-8425-69d74bc636b8
## paper: https://www.nature.com/articles/s41591-021-01329-2
Convert("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/PBMC/covid/covid_portal.h5ad", ".h5seurat")
# You should make quit R and then load again so that you can read the .h5seurat
# This .d5seurat object can then be read manually
seuratObject <- LoadH5Seurat("/covid/.h5seurat") ## failed in converting the data

```

*Additional reprocessing for blood cells one1k
```R
### R/4.2.2-foss-2022a-bare
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/software/Rpackage/4.0/")
dyn.load("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib/libhdf5_hl.so.310", lib.loc = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/hdf5/lib")
library(Seurat)
library(tidyverse)
library(here)
library(data.table)
library(RLinuxModules)
library(loomR)
library(SeuratDisk)
library(R.utils)

setwd("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/PBMC")

dat<-readRDS("local.rds")
dat1<-subset(dat,subset=age<80)

### random 100 samples
metadata<-dat1@meta.data
set.seed(11111)
sample.id<-sample(as.character(metadata$donor_id), size=100, replace =F)

dat1@meta.data$sel<-0
dat1@meta.data$sel[which(dat1@meta.data$donor_id %in% sample.id)]<-1

### extract data
dat2<-subset(dat1,subset=sel==1)

df.meta<-dat2@meta.data
df.meta<-df.meta[,c("orig.ident","nCount_RNA","nFeature_RNA","donor_id","cell_type")]
dat2@meta.data<-df.meta

dat2 <- FindVariableFeatures(object = dat2)
pfile <- as.loom(x = dat2, filename = "PBMC_onek1k_100.loom", verbose = TRUE, overwrite = TRUE)
### Close the .loom file after all
pfile$close_all()
```

## Run cellex
```python
### ml Python 3.10.4
cd /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/CELLEX/cellex
python

### import packages
import sys
sys.path.append("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/CELLEX/")

para1 = int(sys.argv[1])

import loompy         # needed for importing data
import numpy as np    # needed for formatting data
import pandas as pd   # needed for formatting data
import cellex

### Set constants

dirOut = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/CD_Ramik/CO_EPI/out"                    # output directory for results and plots
prefixData = "CO_EPI"                      # prefix to prepend to files
pathData = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/scRNAseq/data/CD_Ramik/CO_EPI/CO_EPI.loom"

nameAnno = "Celltype"                  # metadata annotation column attribute name
nameId = "CellID"                             # metadata cell id column attribute name
nameClass = "Celltype"                  # possible classes to group by: trajectory_main, trajectory_sub

### re-format data 
with loompy.connect(pathData) as ds:
    rows = (ds.row_attrs["Gene"])
    cols = (ds.col_attrs[nameId])
    #our data
    data = pd.DataFrame(ds[:, :], index=rows, columns=cols)
    # the type-annotation for individual cells
    metadata = pd.DataFrame(data={"cell_type" : ds.col_attrs[nameAnno]}, index=ds.col_attrs[nameId])
    metadata_class = pd.DataFrame(data={"cell_class" : ds.col_attrs[nameClass]}, index=ds.col_attrs[nameAnno])

### Create ESObject and compute Expression Specificity
eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
eso.compute(verbose=True)

## Save results to disk
eso.save_as_csv(file_prefix=prefixData, path=dirOut, verbose=True)
    
```
