# R Package - Capybara <img src="/examples/Monocle_hat_colin.png" height="30" width="25">


Capybara is a tool to measure cell identity and fate transitions. This approach is designed to measure cell identity as a continuum, at a single-cell resolution. Capybara enables classification of discrete identities as well as cells with multiple identities. This package has a dependency on R version (R >= 3.5.0). For details regarding the methods, usage and application, please refer to the following paper: Kong et al., BioRxiv 2020 (https://www.biorxiv.org/content/10.1101/2020.02.17.947390v1)

## Installation
### Dependencies

Most dependencies can be installed along with Capybara through CRAN. The following dependencies may need to be installed manually through BioConductor (Instructions can also be found here: https://bioconductor.org/).

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

With the bulk reference, we next load the single-cell reference, such as a cell atlas, and the single-cell sample to be used. The datasets to be used should be in a matrix form with each row representing a gene and each column representing a cell. Here, we use the Mouse Cell Atlas (MCA) as background and single-cell RNA-seq data of mouse pancreatic cells (Baron et al., 2016) as examples for demonstration. MCA can be obtained from https://figshare.com/articles/MCA_DGE_Data/5435866. We included the mouse pancreatic dataset in the package.

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

Notice: For Windows users, please set unix.par=F and n.cores=1

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

**6. Selection of 90 cells from each tissue to construct a QP background**

With all QP scores calculated on the bulk transcriptome profiles of all tissues, we select 90 most relevant cells of each tissue in the MCA (90 highest scored cells in the MCA to each bulk tissue) as a QP background. We use this QP background further to map our sample single-cell data. Assuming that each cell in each cell type of the MCA takes a unique combination of QP scores to each tissue in ARCHS4, cells in the sample that share similar combination to those in MCA are marked to relate to the corresponding tissue in the MCA. Here, we demonstrate how we constructed the backgrounds as following. We included the background matrices along with the packages such that it can be directly used for convenience.

**(a) Get QP scores for all annotated cells**
```r
# Read the QP files from the directory
qp.files.to.read.clean <- list.files("./MCA_All_Tissue_QP/", full.names = T)

full.qp.mtx.known.annotation <- data.frame()
full.qp.mtx.unknown.annotation <- data.frame()
for (i in 1:length(qp.files.to.read.clean)) {
  curr.file <- qp.files.to.read.clean[i]
  curr.qp.rslt <- read.csv(curr.file, header = T, row.names = 1, stringsAsFactors = F)
  
  cells.to.keep <- intersect(rownames(mca.meta), rownames(curr.qp.rslt))
  cells.unlabel <- setdiff(rownames(curr.qp.rslt), cells.to.keep)
  
  curr.sub.mtx.to.keep <- curr.qp.rslt[cells.to.keep, ]
  curr.sub.mtx.unlabel <- curr.qp.rslt[cells.unlabel, ]
  
  if (nrow(full.qp.mtx.known.annotation) <= 0) {
    full.qp.mtx.known.annotation <- curr.sub.mtx.to.keep
    full.qp.mtx.unknown.annotation <- curr.sub.mtx.unlabel
  } else {
    full.qp.mtx.known.annotation <- rbind(full.qp.mtx.known.annotation, curr.sub.mtx.to.keep)
    full.qp.mtx.unknown.annotation <- rbind(full.qp.mtx.unknown.annotation, curr.sub.mtx.unlabel)
  }
}

full.qp.mtx.known.annotation.qp.score.only <- full.qp.mtx.known.annotation[,c(1:(ncol(full.qp.mtx.known.annotation) - 2))]
```

**(b) Selection of 90 cells**
```r
# Create a map between MCA and ARCHS4
map.df <- data.frame(mca.tissue = c("Embryonic-Mesenchyme", "Embryonic-Stem-Cell", "Trophoblast-Stem-Cell", "Fetal_Brain",
                                   "Neonatal-Calvaria","Fetal_Intestine", "Fetal-Liver", "Fetal_Lung", "Fetal_Stomache",
                                   "Neonatal-Heart", "Neonatal-Muscle",
                                   "Neonatal-Rib", "Neonatal-Skin",  "NeonatalPancreas"),
                     corresponding = c("frxn_embryo", "frxn_embryo", "frxn_embryo", "frxn_brain","frxn_brain",
                                       "frxn_small.intestine", "frxn_liver", 
                                       "frxn_lung", "frxn_stomach",  "frxn_heart", "frxn_muscle", "frxn_muscle", 
                                       "frxn_skin", "frxn_pancreas"),
                     stringsAsFactors = F)

# Identify top 90 cells for each tissue
tm.tissue <- unique(map.df$tm.tissue)
cell.selector <- c()
n.sample <- 90
for (i in 1:length(tm.tissue)) {
  curr.tissue <- tm.tissue[i]
  cell.names <- rownames(mca.meta)[which(mca.meta$tissue == curr.tissue)]
  curr.qp.subset <- full.qp.mtx.known.annotation.qp.score.only[cell.names, ]
  curr.map <- map.df$corresponding[which(map.df$tm.tissue == curr.tissue)]
  if (length(curr.map) <= 1){
    curr.qp.subset.sub <- data.frame(score = curr.qp.subset[,curr.map], cell.name = cell.names, stringsAsFactors = F)
  } else {
    curr.qp.subset.sub <- data.frame(score = rowSums(curr.qp.subset[,curr.map]), cell.name = cell.names, stringsAsFactors = F)
  }
  curr.qp.subset.sub.sort <- curr.qp.subset.sub[order(-curr.qp.subset.sub$score), ]
  cells.to.incl <- curr.qp.subset.sub.sort$cell.name[1:n.sample]
  
  cell.selector <- c(cell.selector, cells.to.incl)
}
saveRDS(full.qp.mtx.known.annotation.qp.score.only[cell.selector, ], "./MCA_embryonic_background.RDS")
```

*Note: This constructed QP background can be saved and reused and does not need to be reconstructed every time.*

### Identification of tissue correlate in the reference to the sample single-cell dataset

To find the correlated tissue in the reference to the sample single-cell dataset, we use a correlation based method. In brief, we calculate Pearson's correlation of the QP scores in a pairwise manner between each cell in the sample and each cell in the reference. Recall the assumption that cells in the sample that share similar combination of QP scores to those in MCA are marked to relate to the corresponding tissue in the MCA. If there is a significant percentage of reference cells of a tissue (over 70%) mapped to a cell, we record the tissue label. Then the frequency of each tissue label is calculated. Tissues with a frequency at least 0.5% (for cell number > 10,000) or at least 100 cells will be selected for further analysis. 

**1. Load the QP background matrix**
```r
background.qp.fpath <- system.file("extdata", "MCA Adult Background.Rds", package = "Capybara")
background.mtx <- readRDS(background.qp.fpath)
```

**2. Load the QP scores of the sample**

```r
## Load QP results
qp.rslt <- read.csv("./baron_bulk_classification_qp_scale.csv", row.names = 1, header = T, stringsAsFactors = F)

## Reshape the data
qp.rslt.sub <- qp.rslt[,c(1:(ncol(qp.rslt) - 2))]
```

**3. Correlation calculation**

*Note: we use WGCNA to calculate the correlation*

```r
mtx.test <- t(qp.rslt.sub[, colnames(background.mtx)])
ref.test <- t(background.mtx)

# Pearson's Correlation Calculation
corr.mtx <- WGCNA::cor(ref.test, mtx.test)
```

**4. Binarization based on correlation**

We perform binarization based on the correlation estimates.

```r
# Setup a correlation cutoff to the 90th quantile of the correlation matrix
correlation.cutoff <- quantile(corr.mtx, 0.90)

# Binarization based on the correlation
new.corr.bin <- corr.mtx
new.corr.bin[which(new.corr.bin >= correlation.cutoff)] <- 1
new.corr.bin[which(new.corr.bin < correlation.cutoff)] <- 0
new.corr.bin <- as.data.frame(new.corr.bin)
```

**5. Counting the tissues and select the final tissue types**

Count the frequency of occurrence of each tissue in the tissue list.

```r
# Count
count.in.cat <- c()
unique.cat <- unique(unlist(lapply(strsplit(rownames(new.corr.bin), "_"), function(x) x[1])))
for (uc in unique.cat) {
  curr.subset <- new.corr.bin[which(startsWith(rownames(new.corr.bin), uc)), c(1:1886)]
  count.in.cat[uc] <- sum(colSums(curr.subset) >= nrow(curr.subset) * 0.7)
}

count.in.cat <- as.data.frame(count.in.cat)
count.in.cat$perc <- round(count.in.cat$count.in.cat *100/sum(count.in.cat$count.in.cat), digits = 3)

# Check frequency
final.cell.types.adult <- rownames(count.in.cat)[which(count.in.cat$count.in.cat > 100)]
```

Below is a composition example for this pancreatic dataset, where we identify 3 major tissues, including stomach, pancreas, and small instestine.
<p align="center">
    <img src="/examples/bulk class mca pancreatic.png" height="800" width="400">
</p>


## Step 2: Generation of a High-Resolution Custom Reference, and Continuous Identity Measurement

After tissue-level classification, relevant cell types are selected from cell atlas and built as a single cell reference dataset. As an alternative, users could also use their own single-cell reference dataset to benchmark their samples.

### Systematic construction of a high-resolution reference

To alleviate the effect of technical variations, we construct pseudo-bulk references for each reference cell type. By default, 90 cells of each cell type would be used to build the reference. The construct.high.res.reference function returns a list containing expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference.

**Get the counts of the cell types involved in the tissues selected**

```r
pancreatic.all.meta <- mca.meta[which(mca.meta$cell.bc.tissue %in% final.cell.types.adult), ]

mca.counts.all.involved <- NULL
tissues.to.read <- unique(pancreatic.all.meta$tissue)
general.path <- "~/Box/Morris Lab/Classifier Analysis/Reference datasets/MCA/MCA_Counts/"
for (i in 1:length(tissues.to.read)) {
  curr.t <- tissues.to.read[i]
  curr.path.to.read <- paste0(general.path, curr.t, "/count.csv")
  curr.count <- read.csv(curr.path.to.read, header = T, row.names = 1, stringsAsFactors = F)
  if (is.null(mca.counts.all.involved)) {
    mca.counts.all.involved <- curr.count
  } else {
    mca.counts.all.involved <- cbind(mca.counts.all.involved, curr.count)
  }
}

## Meta data filtering
pancreatic.all.meta$cell.type <- gsub("Dendrtic cell", "Dendritic cell", pancreatic.all.meta$cell.type)
pancreatic.all.meta$cell.type.1 <- gsub("\\([^)]*\\)", "", pancreatic.all.meta$cell.type)
pancreatic.all.meta$cell.type.alone <- unlist(lapply(strsplit(pancreatic.all.meta$cell.type.1, "_"), function(x) x[1]))

## Filter out cell types with less than 30 cells
cell.type.alone.freq <- as.data.frame(table(pancreatic.all.meta$cell.type.alone))
cell.type.over.30 <- cell.type.alone.freq$Var1[which(cell.type.alone.freq$Freq >= 30)]
pancreatic.sub.meta <- pancreatic.all.meta[which(pancreatic.all.meta$cell.type.alone %in% as.character(cell.type.over.30)),]
coldata.df <- pancreatic.sub.meta
```

**Construction**

``` r
# Construction of a high-resolution reference
ref.list <- construct.high.res.reference(mca.counts.all.involved, coldata.df = coldata.df, criteria = "cell.type.alone")
# Get expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference
ref.df <- ref.construction(ref.list[[1]], ref.list[[2]], "cell.type")
```

### Application of quadratic programming on the self-established reference with the sample

``` r
# Measure cell identity in the reference dataset as a background 
single.round.QP.analysis(ref.df, ref.list[[1]], n.cores = 4, save.to.path = "./", save.to.filename = "01_MCA_Based_scClassifier_reference_mix90_normalize_select", unix.par = TRUE)

# Measure cell identity in the query dataset 
single.round.QP.analysis(ref.df, baron.expr, n.cores = 4, save.to.path = "./", save.to.filename = "02_MCA_Based_scClassifier_reference_mix90_test_normalize_select", unix.par = TRUE)
```

## Step 3: Discrete Cell Type Classification and Multiple Identity Scoring

### Empirical p-value calculation
With the constructed single-cell reference, we apply QP to both the sample and reference single-cell datasets to generate continuous measurements of cell identity. The result of this step includes two lists of p-value matrices: one for the reference and the other for the sample. For each cell, each column of the p-value matrix denotes a cell type, while each row describes each round of 50 (default).
``` r
# Read in background and testing identity scores

background.mtx <- read.csv("./01_MCA_Based_scClassifier_reference_mix90_normalize_select_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./02_MCA_Based_scClassifier_reference_mix90_test_normalize_select_scale.csv", header = T, row.names = 1, stringsAsFactors = F)

col.sub <- ncol(background.mtx) - 2

# Conduct reference randomization to get empirical p-value matrix
ref.perc.list <- percentage.calc(background.mtx[,c(1:col.sub)], background.mtx[,c(1:col.sub)])

# Conduct test randomization to get empirical p-value matrix
perc.list <- percentage.calc(as.matrix(mtx.test[,c(1:col.sub)]), as.matrix(background.mtx[,c(1:col.sub)]))
```

### Binarization with Mann-Whitney
A randomized test is performed using the background distributions as null to compute the occurrence probability or empirical p-values of each identity score. This test shapes the likelihood identity score occurrence as a continuous distribution, in which the cell type with the lowest likelihood rank is the classified identity. Capybara is also able to identify cells that harbor multiple identities, potentially representing cells transitioning between defined cell identities. To capture multiple cell identities, we use a Mann-Whitney (MW) test to compare the occurrence probabilities of the cell type with the lowest likelihood rank to that of other cell types, following the order from the second-lowest to the highest rank-sum. From this test, we calculate a p-value to determine whether two identities are equally likely to represent the identity of a specific cell. We stop our comparison when we identify the first cell type that is significantly (p-value < 0.05) less likely to represent one of the cell identities. A binarized matrix will be returned with each row representing a query cell and each column representing a possible cell type. 1 means inferred cell type in the matrix. 

``` r
# Binarization of inference results
bin.count <- binarization.mann.whitney(mtx = mtx.test[,c(1:col.sub)], ref.perc.ls = ref.perc.list, ref.meta = ref.list[[2]], perc.ls = perc.list)
```

### Classification

Finally, we return a classification table of each query cell and its inferred cell type. Cells with multiple inferred identities are marked as "Multi_ID". Cells with no significant inferred identity are marked as "unassigned".

``` r
classification <- binary.to.classification(bin.count[,c(1:col.sub)])
rownames(classification) <- classification$barcode
```

### Check the Classification Result

We check the classification results by comparing the labels that are shared between the reference and manual annotation of Baron et al., 2016. Further, we visualize the agreement using ggplot2.

```r
classification$actual <- baron.meta[rownames(classification), "cell.type"]

table.freq <- table(classification$actual, classification$call)
table.freq.perc <- apply(table.freq, 1, function(x) round(x * 100/sum(x), digits = 3))

rownames(table.freq.perc)[16] <- "beta"

table.freq.sub <- as.data.frame(table.freq.perc[c("B.cell", "beta", "Ductal.cell", "Endothelial.cell",
                                                  "Macrophage", "T.cell", "Dendritic.cell", "Stromal.cell", 
                                                  "Multi_ID", "Endocrine.cell"),
                                                c("B_cell", "beta", "ductal", "endothelial",
                                                  "macrophage", "T_cell", "immune_other", "activated_stellate", 
                                                  "alpha", "delta", "gamma")])
table.freq.sub$Capybara.Call <- rownames(table.freq.sub)
table.freq.melt <- melt(table.freq.sub)

table.freq.melt$Capybara.Call <- factor(table.freq.melt$Capybara.Call,
                                        levels = c("B.cell", "beta", "Ductal.cell", "Endothelial.cell",
                                                   "Macrophage", "T.cell", "Dendritic.cell", "Stromal.cell", 
                                                   "Multi_ID", "Endocrine.cell"),
                                        ordered = T)
table.freq.melt$variable <- factor(table.freq.melt$variable,
                                        levels = c("B_cell", "beta", "ductal", "endothelial",
                                                   "macrophage", "T_cell", "immune_other", "activated_stellate", 
                                                   "alpha", "delta", "gamma"),
                                        ordered = T)

ggplot(table.freq.melt, aes(x = Capybara.Call, y = variable, size=ifelse(value==0, NA,  value))) +
  geom_point(aes(colour = variable)) +
  scale_size_area(name = "Percentage", max_size=12) +
  scale_color_viridis_d(option = "A", begin = 0.15, end = 0.85) +
  ggtitle("Mouse Pancreatic Dataset (Baron et al., 2016)") +
  guides(fill = FALSE, color = FALSE) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12, face = "bold.italic", angle = 90),
        axis.text.y = element_text(size = 12, face = "bold.italic"),
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        title = element_text(face = "bold.italic", size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1))
```

Below is a dot plot example for this pancreatic dataset to show agreement.
<p align="center">
    <img src="/examples/pancreatic dot plot.png" height="450" width="400">
</p>

## Analysis of Cells with Multiple Identities

A unique aspect of Capybara is the classificaiton of cells with multiple identities, which are key to characterize cell fate transitions in a continuous process. Cells with multiple identities label transition harbors, while the discrete cell identities that connect these cells mark potential pivotal states/hallmarks during the continuous processes. In Capybara, we further develop a 'transition metric', a transition score, to measure the flux through the mixed cell identities. It is worth noting that the intention of this score is not to measure potential of each identity but to measure the dynamics going through each discrete state. For details, please refer to the paper. Here, we use an example of cardiomyocyte reprogramming (Stone et al., 2019) to demonstrate the preprocessing of data, classification, analysis of cells with multiple identities and calculation of transition scores.

### 1. Download the data

The dataset for the cardiomyocyte reprogramming can be found here under GEO: GSE131328. This dataset contains 6 timepoints of this reprogramming process, Day -1, 1, 2, 3, 7, and 14, where Day -1 marks the day of transduction of three transcription factors and Day 14 cells were sorted using a-MHC reporter (Stone et al., 2019). The data can be downloaded in terminal as well as in R.

```
wget https://www.ncbi.nlm.nih.gov/geo/download/acc=GSE133452&format=file&file=GSE133452%5Fm1%5F1%5F2%5F3%5F7%5F14P%5Fpaper%2Ecsv%2Egz
```

or

```r
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133452&format=file&file=GSE133452%5Fm1%5F1%5F2%5F3%5F7%5F14P%5Fpaper%2Ecsv%2Egz", "./cardiomyocyte_reprogramming_m1_14p.csv.gz")
unzip("./cardiomyocyte_reprogramming_m1_14p.csv.gz", overwrite = FALSE, exdir = ".")
```

### 2. Preprocessing of the data with Seurat

In this step, we preprocess the data with Seurat to filter the data and obtain a UMAP embedding of the data. For details of Seurat processing, please refer to the instructions or vignettes here - https://satijalab.org/seurat/vignettes.html.

```r
# Read in the file path for all features and genes
feature.file.path <- system.file("extdata", "features.tsv", package = "Capybara")

# Load the data
stone.et.al <- read.csv("./cardiomyocyte_reprogramming_m1_14p.csv", row.names = 1, header = T, stringsAsFactors = F)
feature.df <- read.table(feature.file.path, header = F, row.names = 1, stringsAsFactors = F)

# Map the gene names fr
gene.name.subset <- feature.df[intersect(stone.et.al$X, rownames(feature.df)), ]
stone.et.al.subset <- stone.et.al[which(stone.et.al$X %in% rownames(feature.df)), ]
stone.et.al.subset$gene.name <- gene.name.subset[stone.et.al.subset$X, "V2"]
stone.et.al.subset <- stone.et.al.subset[-which(duplicated(stone.et.al.subset$gene.name)), ]
rnm <- stone.et.al.subset$gene.name
stone.et.al.final <- stone.et.al.subset[, -c(1,ncol(stone.et.al.subset))]
rownames(stone.et.al.final) <- rnm

# Create Seurat object
sc.data.stone <- CreateSeuratObject(counts = stone.et.al.final, project = "cardiac.reprog", min.cells = 3, min.features = 200)

# Calculate mitochondria content
sc.data.stone[["percent.mt"]] <- PercentageFeatureSet(sc.data.stone, pattern = "mt-")

# Visualize QC metrics as a violin plot and scatter plot
VlnPlot(sc.data.stone, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(sc.data.stone, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc.data.stone, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Filter the dataset based on number of features
sc.data.stone <- subset(sc.data.stone, subset = nFeature_RNA > 200 & nFeature_RNA < 5500)

# Log normalize the data
sc.data.stone <- NormalizeData(sc.data.stone, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable gene identification
sc.data.stone <- FindVariableFeatures(sc.data.stone, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(sc.data.stone)
sc.data.stone <- ScaleData(sc.data.stone, features = all.genes)

# PCA
sc.data.stone <- RunPCA(sc.data.stone, features = VariableFeatures(object = sc.data.stone))

# JackStraw procedure and Elbow plot to select number of PCs
sc.data.stone <- JackStraw(sc.data.stone, num.replicate = 100)
sc.data.stone <- ScoreJackStraw(sc.data.stone, dims = 1:20)

JackStrawPlot(sc.data.stone, dims = 1:20)
ElbowPlot(sc.data.stone)

# Identify neighbors and clusters
sc.data.stone <- FindNeighbors(sc.data.stone, dims = 1:18)
sc.data.stone <- FindClusters(sc.data.stone, resolution = 0.8)

# UMAP embedding
sc.data.stone <- RunUMAP(sc.data.stone, dims = 1:18)
```

### 3. Classification

### Step 1. Tissue Classification

Here, we perform the same classification pipeline as described above in the first section, where we obtained four major tissues: neonatal skin, neonatal heart, fetal stomach, and fetal lung. 

**Load the bulk data**

```r
# File path
bulk.raw.path <- system.file("extdata", "Bulk Reference Raw.Rds", package = "Capybara")
bulk.rpkm.path <- system.file("extdata", "Bulk Reference RPKM.Rds", package = "Capybara")
# Read the matrices
bulk.raw <- readRDS(bulk.raw.path)
bulk.rpkm <- readRDS(bulk.rpkm.path)
```

**Application of Quadratic Programming using Bulk**

```r
single.round.QP.analysis(bulk.raw, stone.et.al, scale.bulk.sc = "scale", unix.par = TRUE, 
                         force.eq = 1, n.cores = 4, save.to.path = "./", 
                         save.to.filename = "stone_bulk_classification_qp")
```

**Correlation Analysist**

```r
## Load QP results
qp.rslt <- read.csv("./stone_bulk_classification_qp_scale.csv", row.names = 1, header = T, stringsAsFactors = F)

## Reshape the data
qp.rslt.sub <- qp.rslt[,c(1:(ncol(qp.rslt) - 2))]

## Background matrix
background.qp.fpath <- system.file("extdata", "MCA Embryonic Background.Rds", package = "Capybara")
background.mca <- readRDS(background.qp.fpath)
background.mtx <- background.mca[[2]]

## Correlation Analysis
mtx.test <- t(qp.rslt.sub[, colnames(background.mtx)])
ref.test <- t(background.mtx)

## Pearson's Correlation Calculation
corr.mtx <- WGCNA::cor(ref.test, mtx.test)

## Setup a correlation cutoff to the 90th quantile of the correlation matrix
correlation.cutoff <- quantile(corr.mtx, 0.90)

## Binarization based on the correlation
new.corr.bin <- corr.mtx
new.corr.bin[which(new.corr.bin >= correlation.cutoff)] <- 1
new.corr.bin[which(new.corr.bin < correlation.cutoff)] <- 0
new.corr.bin <- as.data.frame(new.corr.bin)
```

**Mapping to Tissues in Mouse Cell Atlas (MCA)**

```r
# Count
count.in.cat <- c()
unique.cat <- unique(unlist(lapply(strsplit(rownames(new.corr.bin), "_"), function(x) x[1])))
for (uc in unique.cat) {
  curr.subset <- new.corr.bin[which(startsWith(rownames(new.corr.bin), uc)), c(1:30729)]
  count.in.cat[uc] <- sum(colSums(curr.subset) >= nrow(curr.subset) * 0.80)
}

count.in.cat <- as.data.frame(count.in.cat)
count.in.cat$perc <- round(count.in.cat$count.in.cat *100/sum(count.in.cat$count.in.cat), digits = 3)

final.cell.types.fetal <- rownames(count.in.cat)[which(count.in.cat$count.in.cat > 100)]
```

Below is the composition for this cardiac reprogramming dataset, where we identify 4 major tissues.
<p align="center">
    <img src="/examples/cardiac_bulk_v2.png" height="800" width="400">
</p>

### Step 2. Construction of Reference at High-Resolution and Continuous Identity Measurements

**Get the counts of cell types in the selected tissues from MCA**

```r
# Background cells
mca <- read.csv("~/Box/Morris Lab/Classifier Analysis/Reference datasets/MCA/MCA_CellAssignments.csv",
                row.names = 1, header = T, stringsAsFactors = F)
mca.meta <- data.frame(row.names = mca$Cell.name, 
                       tissue = mca$Tissue,
                       cell.bc.tissue = unlist(lapply(strsplit(mca$Cell.name, "_"), function(x) x[1])),
                       cell.type = mca$Annotation,
                       stringsAsFactors = F)

cardiac.rp.all.meta <- mca.meta[which(mca.meta$cell.bc.tissue %in% final.cell.types.fetal), ]

mca.counts.all.involved <- NULL
tissues.to.read <- unique(cardiac.rp.all.meta$tissue)
general.path <- "~/Box/Morris Lab/Classifier Analysis/Reference datasets/MCA/MCA_Counts/"
for (i in 1:length(tissues.to.read)) {
  curr.t <- tissues.to.read[i]
  curr.path.to.read <- paste0(general.path, curr.t, "/count.csv")
  curr.count <- read.csv(curr.path.to.read, header = T, row.names = 1, stringsAsFactors = F)
  if (is.null(mca.counts.all.involved)) {
    mca.counts.all.involved <- curr.count
  } else {
    mca.counts.all.involved <- cbind(mca.counts.all.involved, curr.count)
  }
}

## meta data cleaning
cardiac.rp.all.meta$cell.type.1 <- gsub("\\([^)]*\\)", "", cardiac.rp.all.meta$cell.type)
cardiac.rp.all.meta$cell.type.alone <- unlist(lapply(strsplit(cardiac.rp.all.meta$cell.type.1, "_"), function(x) x[1]))

cardiac.rp.all.meta$cell.type.1 <- tolower(cardiac.rp.all.meta$cell.type.1)
coldata.df <- cardiac.rp.all.meta
```

**Construction**
```r
# Construction of a high-resolution reference
ref.list <- construct.high.res.reference(mca.counts.all.involved, coldata.df = coldata.df, criteria = "cell.type.1")
# Get expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference
ref.df <- ref.construction(ref.list[[1]], ref.list[[2]], "cell.type")
```

**Application of Quadratic Programming**
```r
single.round.QP.analysis(ref.df, ref.list[[1]], n.cores = 4, save.to.path = "./", save.to.filename = "stone_et_al_reference_MCA")
single.round.QP.analysis(ref.df, stone.et.al, n.cores = 4, save.to.path = "./", save.to.filename = "stone_et_al_test_MCA")
```

### Step 3. Discrete Cell Type Classification and Multiple Identity scoring


*Note: this will be continuously updating*
