###############################################
# Step 2: ARCHS4 Bulk Cleaning
###############################################
library(viridis)
library(pheatmap)
library(ggplot2)
library(reshape2)

### function to get the most similar n samples
get.most.connected <- function(mtx, n.sample) {
  corr.mtx <- WGCNA::cor(mtx)
  corr.mtx.upper <- corr.mtx * upper.tri(corr.mtx)
  corr.mtx.melt <- melt(corr.mtx.upper)
  corr.mtx.melt.pos <- corr.mtx.melt[which(corr.mtx.melt$value > 0), ]
  corr.mtx.melt.pos.sort <- corr.mtx.melt.pos[order(-corr.mtx.melt.pos$value), ]
  corr.mtx.melt.pos.sort$Var1 <- as.character(corr.mtx.melt.pos.sort$Var1)
  corr.mtx.melt.pos.sort$Var2 <- as.character(corr.mtx.melt.pos.sort$Var2)
  
  sample.list <- c()
  count.line <- 1
  while(length(sample.list) < n.sample) {
    sample.list <- unique(c(sample.list, unique(c(corr.mtx.melt.pos.sort$Var1[count.line], corr.mtx.melt.pos.sort$Var2[count.line]))))
    count.line <- count.line + 1
  }
  
  return(sample.list[1:90])
}

### Load the raw data/rpkm data from the step 1
raw.counts <- readRDS("03_raw_count_Tissue_bulk.Rds")
rpkm.count <- readRDS("03_rpkm_count_Tissue_bulk.Rds")
accession.full <- readRDS("01_full_accession_information_data_frame.Rds")
rownames(accession.full) <- accession.full$accession

### Subset the accession 
accession.sub <- accession.full[colnames(raw.counts), ]
accession.sub$source <- tolower(accession.sub$source)
accession.freq <- as.data.frame(table(accession.sub$source, accession.sub$series_id))

accession.sub$merge <- paste0(accession.sub$source, "_", accession.sub$series_id)
accession.sub.freq <- as.data.frame(table(accession.sub$merge))
accession.sub.freq$Var1 <- as.character(accession.sub.freq$Var1)
accession.sub.freq$tissue <- unlist(lapply(strsplit(accession.sub.freq$Var1, "_"), function(x) x[1]))

### Select for the most associated 90 samples across different GEO RNA-seq data
tissue.uniq <- unique(accession.sub$source)
tissue.uniq <- tissue.uniq[which(!is.na(tissue.uniq))]
series.to.keep <- c()
sample.accession.to.keep <- c()
new.accession <- data.frame()
new.sample.accession <- data.frame()
sample.count <- data.frame()
sample.rpkm <- data.frame()
n.sample <- 90
for (tis in tissue.uniq) {
  curr.freq.sub <- accession.sub.freq[which(accession.sub.freq$tissue == tis), ]
  curr.freq.sub.sort <- curr.freq.sub[order(-curr.freq.sub$Freq), ]
  curr.series.to.keep <- strsplit(curr.freq.sub.sort$Var1[1], "_")[[1]][2]
  curr.geo.to.keep <- accession.sub[which(accession.sub$series_id == curr.series.to.keep & accession.sub$source == tis), ]
  curr.geo.to.keep$tissue <- tis
  print(tis)
  print(nrow(curr.geo.to.keep))
  series.to.keep <- c(series.to.keep, curr.series.to.keep)
  #  sample.accession.to.keep <- c(sample.accession.to.keep, curr.geo.to.keep)
  
  if (nrow(curr.geo.to.keep) >= n.sample) {
    geo.sample <- get.most.connected(log2(rpkm.count[,curr.geo.to.keep$accession] + 1), n.sample = n.sample)
  } else {
    geo.sample <- sample(curr.geo.to.keep$accession, n.sample, replace = T)
  }
  
  curr.count.sample <- raw.counts[, geo.sample]
  curr.rpkm.sample <- rpkm.count[, geo.sample]
  
  colnames(curr.count.sample) <- paste0(tis, "_", seq(1,n.sample))
  colnames(curr.rpkm.sample) <- paste0(tis, "_", seq(1,n.sample))
  
  if (nrow(new.accession) <= 0) {
    new.accession <- curr.geo.to.keep
    new.sample.accession <- curr.geo.to.keep[geo.sample,]
    sample.count <- curr.count.sample
    sample.rpkm <- curr.rpkm.sample
  } else {
    new.accession <- rbind(new.accession, curr.geo.to.keep)
    new.sample.accession <- rbind(new.sample.accession, curr.geo.to.keep[geo.sample, ])
    sample.count <- cbind(sample.count, curr.count.sample)
    sample.rpkm <- cbind(sample.rpkm, curr.rpkm.sample)
  }
} 

accession.new.freq <- as.data.frame(table(new.accession$source))

### Look at the correlation between samples
corr.sample.rpkm <- WGCNA::cor(sample.rpkm)
corr.sample.rpkm.log <- WGCNA::cor(log2(sample.rpkm + 1))

my.col <- data.frame(row.names = rownames(corr.sample.rpkm), 
                     tissue = unlist(lapply(strsplit(rownames(corr.sample.rpkm), "_"), function(x) x[1])), 
                     stringsAsFactors = F)

### Heatmap plot to see if the samples selected are well correlated as well as well distinguished from other tissues
pheatmap::pheatmap(corr.sample.rpkm, color = viridis(20, option = "A"), show_colnames = F, show_rownames = F, cluster_cols = F, cluster_rows = F, clustering_method = "ward.D2", 
                   annotation_col = my.col, cellheight = 0.2, cellwidth = 0.2, file = "~/Desktop/test.pdf")
pheatmap::pheatmap(corr.sample.rpkm.log, color = viridis(20, option = "A"), show_colnames = F, show_rownames = F, cluster_cols = F, cluster_rows = F, clustering_method = "ward.D2", 
                   annotation_col = my.col, cellheight = 0.2, cellwidth = 0.2, file = "~/Desktop/test_log.pdf")

saveRDS(sample.count, "~/Box/Morris Lab/Classifier Analysis/ARCHS4 Reference/05_90_sampled_raw_counts.Rds")
saveRDS(sample.rpkm, "~/Box/Morris Lab/Classifier Analysis/ARCHS4 Reference/05_90_sampled_raw_rpkm.Rds")

### Compute the final averaged bulk dataset
final.rpkm.tissue <- data.frame()
final.raw.count.tissue <- data.frame()
uniq.tissue <- unique(my.col$tissue)

for (tis in uniq.tissue) {
  curr.rpkm.tis <- sample.rpkm[, rownames(my.col)[which(my.col$tissue == tis)]]
  curr.rpkm.mean <- as.data.frame(rowMeans(curr.rpkm.tis))
  
  curr.count.tis <- sample.count[, rownames(my.col)[which(my.col$tissue == tis)]]
  curr.count.mean <- as.data.frame(rowMeans(curr.count.tis))
  
  colnames(curr.rpkm.mean) <- tis
  colnames(curr.count.mean) <- tis
  
  if (ncol(final.rpkm.tissue) <= 0) {
    final.rpkm.tissue <- curr.rpkm.mean
    final.raw.count.tissue <- curr.count.mean
  } else {
    final.rpkm.tissue <- cbind(final.rpkm.tissue, curr.rpkm.mean)
    final.raw.count.tissue <- cbind(final.raw.count.tissue, curr.count.mean)
  }
}

### Check if the averaged bulk dataset has tissues that are distinguishable from each other
rpkm.tissue.cor <- WGCNA::cor(log2(final.rpkm.tissue + 1))
pheatmap::pheatmap(rpkm.tissue.cor, color = viridis(20, option = "A"), show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, clustering_method = "ward.D2", 
                   cellheight = 15, cellwidth = 15, file = "~/Desktop/test_log.pdf")

saveRDS(final.rpkm.tissue, "~/Box/Morris Lab/Classifier Analysis/ARCHS4 Reference/05_90_sample_final_rpkm_tissue_collapsed.Rds")
saveRDS(final.raw.count.tissue, "~/Box/Morris Lab/Classifier Analysis/ARCHS4 Reference/05_90_sample_final_raw_count_tissue_collapsed.Rds")

