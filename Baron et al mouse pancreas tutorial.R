library(Capybara)
baron.fns <- c("~/Box/Morris Lab/Classifier Analysis/Reference datasets/Baron et al human pancreatic dataset/Counts/GSM2230761_mouse1_umifm_counts.csv",
               "~/Box/Morris Lab/Classifier Analysis/Reference datasets/Baron et al human pancreatic dataset/Counts/GSM2230762_mouse2_umifm_counts.csv")

baron.expr <- data.frame()
baron.meta <- data.frame()
count <- 1
for (f in baron.fns) {
  curr.df <- read.csv(f, row.names = 1, stringsAsFactors = F, check.names = F)
  curr.bc <- paste0(rownames(curr.df), "_", curr.df$barcode, "_Sample_", count)
  rownames(curr.df) <- curr.bc
  if (nrow(baron.meta) <= 0) {
    baron.meta <- data.frame(row.names = curr.bc, cell.type = curr.df[curr.bc, "assigned_cluster"], barcode = curr.df[curr.bc, "barcode"], 
                             stringsAsFactors = F, check.names = F)
  } else {
    curr.meta <- data.frame(row.names = curr.bc, cell.type = curr.df[curr.bc, "assigned_cluster"], barcode = curr.df[curr.bc, "barcode"],
                            stringsAsFactors = F, check.names = F)
    baron.meta <- rbind(baron.meta, curr.meta)
  }
  curr.expr <- curr.df[, c(3:ncol(curr.df))]
  if (nrow(baron.expr) <= 0) {
    baron.expr <- curr.expr
  } else  {
    baron.expr <- rbind(baron.expr, curr.expr)
  }
  count <- count +1
}

count.mtx <- t(baron.expr)
write.csv(count.mtx, "~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/baron_et_al_pancreatic.csv", quote = F)
write.csv(baron.meta, "~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/baron_et_al_pancreatic_meta.csv", quote = F)

### Tissue Classification
# Read the bulk
bulk.rpkm.archs4 <- readRDS("~/Box/Morris Lab/Classifier Analysis/ARCHS4 Reference/05_90_sample_final_rpkm_tissue_collapsed.Rds")
bulk.count.archs4 <- readRDS("~/Box/Morris Lab/Classifier Analysis/ARCHS4 Reference/05_90_sample_final_raw_count_tissue_collapsed.Rds")

norm.bulk <- normalize.dt(bulk.count.archs4)
norm.sc.mtx <- normalize.dt(count.mtx)

# Calculate the scaling ratio
scale.ratio.sc <- calc.scale.ratio(bulk.count.archs4, count.mtx)
scale.norm.sc <- norm.bulk/scale.ratio.sc

# Calculate gene intersections
norm.ls.sc <- gene.intersect.sub(norm.bulk, norm.sc.mtx)
scale.ls.sc <- gene.intersect.sub(scale.norm.sc, norm.sc.mtx)

# Finish the log-normalization
log.norm.ls.sc <- lapply(norm.ls.sc, log1p)
log.scale.ls.sc <- lapply(scale.ls.sc, log1p)

qp.rslt <- sc.quad.prog.run(as.matrix(log.scale.ls.sc[[1]]), 
                            single.cell.transcriptome = log.scale.ls.sc[[2]], 
                            unix.parallel = TRUE, parallel.cores = 4, force.eq = 1)
write.csv(qp.rslt, "~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/Pancreatic Baron et al/bulk_classification_qp_whole_gene_set_rpkm.csv", quote = F, row.names = F)

# Read QP scores
qp.rslt <- read.csv("~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/Pancreatic Baron et al/bulk_classification_qp_whole_gene_set.csv", row.names = 1, header = T, stringsAsFactors = F)
qp.paga.rslt.sub <- qp.rslt[,c(1:(ncol(qp.rslt) - 2))]
qp.paga.rslt.stk <- stack(qp.paga.rslt.sub)
qp.paga.rslt.stk$pred <- unlist(lapply(strsplit(as.character(qp.paga.rslt.stk$ind), "_"), function(x) x[2]))

ggplot(qp.paga.rslt.stk, aes(x = values, y = pred, color = pred)) +
  geom_jitter() +
  labs(x = "Identity Score", y = "Prediction") +
  ggtitle("Direct Programming of mESC") +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(12, "Paired"), 
                                RColorBrewer::brewer.pal(8, "Set2"), 
                                RColorBrewer::brewer.pal(8, "Dark2"),
                                RColorBrewer::brewer.pal(4, "Pastel1"))) +
  theme(legend.position = "None",
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        title = element_text(face = "bold.italic", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_blank())

# Read MCA
background.mca <- readRDS("~/Box/Morris Lab/Classifier Analysis/ARCHS4 Reference/MCA ARCHS4/06_qp_top_90_each_category_backgroung.Rds")
background.mtx <- background.mca
mtx.test <- qp.paga.rslt.sub[,colnames(background.mtx)]

### cor test
ref.test <- t(background.mtx)
mtx.test.cor <- t(qp.paga.rslt.sub)
corr.mtx <- WGCNA::cor(ref.test, mtx.test.cor)

correlation.cutoff <- quantile(corr.mtx, 0.90)

new.corr.bin <- corr.mtx
new.corr.bin[which(new.corr.bin >= correlation.cutoff)] <- 1
new.corr.bin[which(new.corr.bin < correlation.cutoff)] <- 0
new.corr.bin <- as.data.frame(new.corr.bin)

new.corr.bin$cell.bc.ref <- rownames(new.corr.bin)
new.corr.bin.melt <- reshape2::melt(new.corr.bin)
new.corr.bin.melt.sub <- new.corr.bin.melt[which(new.corr.bin.melt$value > 0),]
new.corr.bin.melt.sub$cell.type <- unlist(lapply(strsplit(new.corr.bin.melt.sub$cell.bc.ref, "_"), function(x) x[1]))

count.in.cat <- c()
unique.cat <- unique(unlist(lapply(strsplit(rownames(new.corr.bin), "_"), function(x) x[1])))
for (uc in unique.cat) {
  curr.subset <- new.corr.bin[which(startsWith(rownames(new.corr.bin), uc)), c(1:1886)]
  count.in.cat[uc] <- sum(colSums(curr.subset) >= nrow(curr.subset) * 0.7)
}

count.in.cat <- as.data.frame(count.in.cat)
count.in.cat$perc <- round(count.in.cat$count.in.cat *100/sum(count.in.cat$count.in.cat), digits = 3)

final.cell.types.adult <- rownames(count.in.cat)[which(count.in.cat$count.in.cat > 100)]

comp.raw <- count.in.cat
comp.raw <- comp.raw[order(-comp.raw$perc), ]
comp.raw$Var1 <- rownames(comp.raw)
comp.raw$Var1 <- factor(comp.raw$Var1, comp.raw$Var1, ordered = T)

comp.raw$label <- "Mouse Pancreatic Cells"

ggplot(comp.raw, aes(x = comp.raw$label, y = comp.raw$perc, fill = comp.raw$Var1, label = comp.raw$Var1)) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5), fontface = "bold", aes(size = comp.raw$perc)) +
  scale_fill_manual(
    name = "Mapped MCA Tissue",
    values = c(RColorBrewer::brewer.pal(12, "Paired"),
               RColorBrewer::brewer.pal(8, "Set2"),
               RColorBrewer::brewer.pal(8, "Set3"))
  ) +
  labs(y = "Percentage of Cells") +
  ggtitle("Bulk Classification to MCA Tissues") +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold.italic", size = 12),
        title = element_text(face = "bold.italic", size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black", size =1))

### MCA construction
# Background cells
mca <- read.csv("~/Box/Morris Lab/Classifier Analysis/Reference datasets/MCA/MCA_CellAssignments.csv",
                row.names = 1, header = T, stringsAsFactors = F)
mca.meta <- data.frame(row.names = mca$Cell.name, 
                       tissue = mca$Tissue,
                       cell.bc.tissue = unlist(lapply(strsplit(mca$Cell.name, "_"), function(x) x[1])),
                       cell.type = mca$Annotation,
                       stringsAsFactors = F)

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

pancreatic.all.meta$cell.type <- gsub("Dendrtic cell", "Dendritic cell", pancreatic.all.meta$cell.type)
pancreatic.all.meta$cell.type.1 <- gsub("\\([^)]*\\)", "", pancreatic.all.meta$cell.type)
pancreatic.all.meta$cell.type.alone <- unlist(lapply(strsplit(pancreatic.all.meta$cell.type.1, "_"), function(x) x[1]))

cell.type.tissue.freq <- table(pancreatic.all.meta$cell.type.alone, pancreatic.all.meta$tissue)
cell.type.tissue.freq.bin <- cell.type.tissue.freq
cell.type.tissue.freq.bin[which(cell.type.tissue.freq.bin > 0)] <- 1
cell.type.expand <- as.data.frame(sort(rowSums(cell.type.tissue.freq.bin))/ncol(cell.type.tissue.freq.bin))

### High resolution construction
cell.type.alone.freq <- as.data.frame(table(pancreatic.all.meta$cell.type.alone))
cell.type.over.30 <- cell.type.alone.freq$Var1[which(cell.type.alone.freq$Freq >= 30)]
pancreatic.sub.meta <- pancreatic.all.meta[which(pancreatic.all.meta$cell.type.alone %in% as.character(cell.type.over.30)),]

cell.num.for.ref <- 90 
pancreatic.sub.meta$cell.type.1 <- tolower(pancreatic.sub.meta$cell.type.1)
coldata.df <- pancreatic.sub.meta
ct.freq.raw <- as.data.frame(table(coldata.df$cell.type.alone), stringsAsFactors = F)
ct.freq <- ct.freq.raw[which(ct.freq.raw$Freq >= 30), ]

# ref.sc, ref.meta, ref.df
ref.list <- construct.high.res.reference(mca.counts.all.involved, coldata.df = coldata.df, criteria = "cell.type.alone")
ref.df <- ref.construction(ref.list[[1]], ref.list[[2]], "cell.type")

# Run QP
single.round.QP.analysis(ref.df, ref.list[[1]], n.cores = 4, save.to.path = "~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/Pancreatic Baron et al/", save.to.filename = "01_MCA_Based_scClassifier_reference_mix90_normalize_select", unix.par = TRUE)
single.round.QP.analysis(ref.df, count.mtx, n.cores = 4, save.to.path = "~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/Pancreatic Baron et al/", save.to.filename = "02_MCA_Based_scClassifier_reference_mix90_test_normalize_select", unix.par = TRUE)

background.mtx <- read.csv("~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/Pancreatic Baron et al/01_MCA_Based_scClassifier_reference_mix90_normalize_select_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/Pancreatic Baron et al/02_MCA_Based_scClassifier_reference_mix90_test_normalize_select_scale.csv", header = T, row.names = 1, stringsAsFactors = F)

# reference Permutation
col.sub <- ncol(background.mtx) - 2
ref.perc.list <- percentage.calc(background.mtx[,c(1:col.sub)], background.mtx[,c(1:col.sub)])

# Test Permutation
perc.list <- percentage.calc(as.matrix(mtx.test[,c(1:col.sub)]), as.matrix(background.mtx[,c(1:col.sub)]))

bin.count <- binarization.mann.whitney(mtx = mtx.test[,c(1:col.sub)], ref.perc.ls = ref.perc.list, ref.meta = ref.list[[2]], perc.ls = perc.list)
classification <- binary.to.classification(bin.count[,c(1:col.sub)])

rownames(classification) <- classification$barcode
classification$actual <- baron.meta[rownames(classification), "cell.type"]

table.freq <- table(classification$actual, classification$call)
table.freq.perc <- apply(table.freq, 1, function(x) round(x * 100/sum(x), digits = 3))

rownames(table.freq.perc)[16] <- "beta"

table.freq.sub <- as.data.frame(table.freq.perc[c("B.cell", "beta", "Ductal.cell", "Endothelial.cell",
                                                "Macrophage", "T.cell", "Dendritic.cell", "Stromal.cell", "Endocrine.cell"), c(1,2,5,8,3,4,6,7,9,10,13)])
table.freq.sub$Capybara.Call <- rownames(table.freq.sub)
table.freq.melt <- melt(table.freq.sub)

table.freq.melt$Capybara.Call <- factor(table.freq.melt$Capybara.Call,
                                        levels = c("B.cell", "beta", "Ductal.cell", "Endothelial.cell",
                                                   "Macrophage", "T.cell", "Dendritic.cell", "Stromal.cell", "Multi_ID", "Endocrine.cell"),
                                        ordered = T)
table.freq.melt$variable <- factor(table.freq.melt$variable,
                                        levels = c("B_cell", "beta", "ductal", "endothelial",
                                                   "macrophage", "T_cell", "immune_other", "activated_stellate", "alpha", "delta", "gamma"),
                                        ordered = T)

pdf("~/Desktop/Morris Lab/Paper/Manuscript/scClassifier/Pancreatic Baron et al/dot plot.pdf", width = 8, height = 9, paper = "special")
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
dev.off()

