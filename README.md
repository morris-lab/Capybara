# R Package - Capybara
Capybara is a tool to measure cell identity and fate transitions. This approach is designed to measure cell identity as a continuum, at a single-cell resolution. Capybara enables classification of discrete entities as well as cells with multiple identities. This package have a dependency on R version (R >= 3.5.0). For details regarding the methods, usage and application, please refer to the following papaer - *\<Fill this in\>*.

## Installation
### Dependencies
Most dependency package can be installed along with Capybara through CRAN. The following dependency may need to be installed manually through BioConductor (Instructions can also be found here: https://bioconductor.org/).

Install BiocManager
```r
install.packages("BiocManager")
```
Install dependency packages
```r
BiocManager::install("limma")
```

### Install the package
Install devtools
```r
install.packages("devtools")
```
Install the package from GitHub.
```r
library("devtools")
devtools::install_github("morris-lab/Capybara")
```
Load the package
```r
library("Capybara")
```

## Step 1: Tissue-Level Classification
### Application of quadratic programming on reference and sample single-cell dataset using a bulk reference
Bulk transcriptome profiles of all tissues are mined from ARCHS4, a platform that contains most published RNA-seq and ChiP-seq datasets (Lachmann et al.,  2018). ARCHS4 obtains raw datasets from the Gene Expression Omnibus (GEO), realigned and processed through a uniform pipeline. We filtered to contain only poly-A and total RNA-seq data from C57BL/6 mice. With further filtering and preprocessing (details can be found in the method section of the paper), we landed with a reference of a total of 30 tissues. We provide our mined bulk references, including a matrix in raw counts and a matrix in reads per kilobase per million (RPKM), as a part of the Capybara package. Selection of your preferred normalization method can be applied to raw counts. Here, we will demonstrate the usage of the bulk raw counts in the pipeline.

**1. Load the bulk reference**
```r
# File path
bulk.raw.path <- system.file("extdata", "Bulk Reference Raw.Rds", package = "Capybara")
bulk.rpkm.path <- system.file("extdata", "Bulk Reference RPKM.Rds", package = "Capybara")
# Read the matrices
bulk.raw <- readRDS(bulk.raw.path)
bulk.rpkm <- readRDS(bulk.rpkm.path)
```

With the bulk reference, we next load the single-cell reference, such as a cell atlas, and the single-cell sample to be used. The datasets to be used should be in a matrix form with each row representing a gene and each column representing a cell. Here, we use the Mouse Cell Atlas (MCA) as background and single-cell RNA-seq data of mouse pancreatic cells (Baron et al., 2016) as examples for demonstration. MCA can be obtained from https://figshare.com/articles/MCA_DGE_Data/5435866. We included the mouse pancreatic dataset along with the package.

**2.Load the single-cell sample dataset and the corresponding meta data**

*Note: The meta data of this file contains 2 columns, where the first column represents cell.type and the second column represents barcode.*

```r
# Read in the pancreatic data file that come with the package
fpath <- system.file("extdata", "baron_dataset.zip", package = "Capybara")
extract.dir <- "."
# Extract the dataset
unzip(fpath, overwrite = FALSE, exdir = ".")
# Identify the full path
full.fpath.meta <- paste0(extract.dir, "/", "baron_et_al_pancreatic_meta.csv")
full.fpath.raw <- paste0(extract.dir, "/", "baron_et_al_pancreatic.csv")
# Load the count matrix and the meta data
baron.expr <- read.csv(full.fpath.raw, header = T, row.names = 1, stringsAsFactors = F)
baron.meta <- read.csv(full.fpath.meta, header = T, row.names = 1, stringsAsFactors = F)
```

**3. Application of QP on the sample single-cell data**

```r
single.round.QP.analysis(bulk.raw, baron.expr, scale.bulk.sc = "scale", unix.par = TRUE, 
                         force.eq = 1, n.cores = 4, save.to.path = "./", 
                         save.to.filename = "baron_bulk_classification_qp")
```

**4. Load the single-cell reference meta data**

*Note: The meta data of Mouse Cell Atlas contains 6 columns, including Cell.name, ClusterID, Tissue, Batch, Cell.Barcode, and Annotation. The annotation is what we used for high-resolution reference construction. We've included the version of meta data we used along with the package.*

```r
# Read the meta data
mca.meta.fpath <- system.file("extdata", "MCA_CellAssignments.csv", package = "Capybara")
mca <- read.csv(mca.meta.fpath, row.names = 1, header = T, stringsAsFactors = F)
# Clean up the meta data
mca.meta <- data.frame(row.names = mca$Cell.name, 
                       tissue = mca$Tissue,
                       cell.type = mca$Annotation,
                       stringsAsFactors = F)
```

**5. Load the single-cell reference atlas and apply QP tissue-by-tissue**

*Due to the large size of MCA count data, we did* ***NOT*** *include the counts along with the package. We further separated MCA into fetal/neonatal/embryonic and adult categories. The counts data were organized in the following manner.*

> Folder: MCA Counts

>> Tissue_1 Folder

>>> count.csv

>> Tissue_2 Folder

>>> count.csv

>> ...

```r
# List all possible files and tissues in the Mouse Cell Atlas
file.ls <- list.files("./MCA_Counts/", full.names = T)
base.nms <- basename(file.ls)

# Identify the tissues
unq.tissue <- unique(base.nms)

# Set a path to save all QP files for all tissues
general.path.to.save <- "./MCA_All_Tissue_QP/"
for (k in 1:length(unq.tissue)) {
  curr.tissue <- unq.tissue[k]
  curr.filename <- paste0("0", k, "_", curr.tissue, "_Bulk_ARCHS4")
  
  file.base.name <- base.nms[which(startsWith(base.nms, curr.tissue))][1]
  file.full <- file.ls[which(startsWith(base.nms, curr.tissue))][1]
  
  print(curr.tissue)
  
  sc.data <- read.csv(paste0(file.full, "/count.csv"), header = T, row.names = 1, stringsAsFactors = F)
  
  if (all(is.na(sc.data))) {
    print("There is no data in this counting matrix!")
  } else {
    single.round.QP.analysis(bulk.raw, sc.data, scale.bulk.sc = "scale", unix.par = TRUE, 
                             force.eq = 1, n.cores = 4, save.to.path = general.path.to.save, 
                             save.to.filename = curr.filename)
  }
}
```

**6. Selection of 90 cells from each tissue to construct a QP background*


### Identification of tissue correlate in the reference to the sample single-cell dataset

## Step 2: Generation of a High-Resolution Custom Reference, and Continuous Identity Measurement
After tissue-level classification, relevant cell types are selected from cell atlas and built as a single cell reference dataset. As an alternative, users could also use their own single-cell reference dataset to benchmark their samples.
### Systematic construction of a high-resolution reference

To alleviate the effect of technical variations, we construct pseudo-bulk references for each reference cell type. By default, 90 cells of each cell type would be used to build the reference. The construct.high.res.reference function returns a list containing expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference.
``` r
# Construction of a high-resolution reference
high.res.ref<-construct.high.res.reference(ref.mtx = capy.44h.mtx,coldata.df = capy.44h.meta)

```
### Application of quadratic programming on the self-established reference with the sample
## Step 3: Discrete Cell Type Classification and Multiple Identity Scoring

### Empirical p-value calculation
With the constructed single-cell reference, we apply QP to both the sample and reference single-cell datasets to generate continuous measurements of cell identity. The result of this step includes two lists of p-value matrices: one for the reference and the other for the sample. For each cell, each column of the p-value matrix denotes a cell type, while each row describes each round of 50 (default).
``` r
# Get expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference
ref.sc<-high.res.ref[[1]]
ref.meta<-high.res.ref[[2]]
ref.df<-high.res.ref[[3]]
# Read in query single-cell count data
sc.data<-capy.44h.mtx
# Measure cell identity in the reference dataset as a background 
single.round.QP.analysis(ref.df, ref.sc, n.cores = 1, save.to.path = "c:/Users/Yuheng Fu/Desktop/", save.to.filename = "03_MCA_Based_scClassifier_reference_mix90_normalize_select")
# Measure cell identity in the query dataset 
single.round.QP.analysis(ref.df, sc.data, n.cores = 1, save.to.path = "c:/Users/Yuheng Fu/Desktop/", save.to.filename = "04_MCA_Based_scClassifier_reference_mix90_test_normalize_select")

```
### Binarization with Mann-Whitney
A randomized test is performed using the background distributions as null to compute the occurrence probability or empirical p-values of each identity score. This test shapes the likelihood identity score occurrence as a continuous distribution, in which the cell type with the lowest likelihood rank is the classified identity. Capybara is also able to identify cells that harbor multiple identities, potentially representing cells transitioning between defined cell identities. To capture multiple cell identities, we use a Mann-Whitney (MW) test to compare the occurrence probabilities of the cell type with the lowest likelihood rank to that of other cell types, following the order from the second-lowest to the highest rank-sum. From this test, we calculate a p-value to determine whether two identities are equally likely to represent the identity of a specific cell. We stop our comparison when we identify the first cell type that is significantly (p-value < 0.05) less likely to represent one of the cell identities. A binarized matrix will be returned with each row representing a query cell and each column representing a possible cell type. 1 means inferred cell type in the matrix. 
``` r
# Read in background and testing identity scores
background.mtx <- read.csv("c:/Users/Yuheng Fu/Desktop/03_MCA_Based_scClassifier_reference_mix90_normalize_select_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("c:/Users/Yuheng Fu/Desktop/04_MCA_Based_scClassifier_reference_mix90_test_normalize_select_scale.csv", header = T, row.names = 1, stringsAsFactors = F)

# Conduct reference randomization to get empirical p-value matrix
col.sub <- ncol(background.mtx) - 2
system.time(ref.perc.list <- percentage.calc(background.mtx[,c(1:col.sub)], background.mtx[,c(1:col.sub)]))

# Conduct test randomization to get empirical p-value matrix
perc.list <- percentage.calc(as.matrix(mtx.test[,c(1:col.sub)]), as.matrix(background.mtx[,c(1:col.sub)]))
saveRDS(list(ref.perc.list, perc.list), "c:/Users/Yuheng Fu/Desktop/permutation_list.RDS")

# Binarization of inference results
bin.count <- binarization.mann.whitney(mtx = mtx.test[,c(1:col.sub)], ref.perc.ls = ref.perc.list, ref.meta = ref.meta, perc.ls = perc.list)

```
### Classification
Finally, we return a classification table of each query cell and its inferred cell type. Cells with multiple inferred identities are marked as "Multi_ID". Cells with no significant inferred identity are marked as "unassigned"
``` r
classification <- binary.to.classification(bin.count[,c(1:col.sub)])
```
## Analysis of Cells with Multiple Identities
``` r

```
