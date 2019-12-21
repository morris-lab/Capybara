#' Top Variable Gene Identificaiton
#'
#' This function identifies top n variable genes using method described in KeyGenes algorithm. In brief, this algorithm use cross validation based on LASSO regression. For detailed explanation, please refer to the paper listed in the note.
#' @param input.dir The path to where the dataset is saved. Dataset should be saved as a tab-delimited table with row = gene, column = cell.
#' @param output.dir The path to where the output should be saved. The output will be saved as a tab-delimited table with NO row names or column names.
#' @param top.number.count The number of top variable genes to extract. Default to 500 genes
#' @keywords top variable gene extraction
#' @note Code reference from Roost et. al., KeyGenes, a Tool to Probe Tissue Differentiation Using a Human Fetal Transcriptional Atlas, Cell Stem Cell Reports, 2015
#' @export
#' @examples
#' top.genes("~/Desktop/sample_data.txt", "~/Desktop/sample_gene_list_output.txt", top.number.count = 1000)
top.genes <- function(input.dir, output.dir, top.number.count = 500) {
  dataset <- input.dir
  top.dir <- output.dir

  DS <- read.table(dataset, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
  X <- voom(as.matrix(DS), normalize.method="none")$E
  vars <- apply(X, 1, var)
  top <- sort.list(-vars)[1:top.number.count]
  X <- X[top,]

  topl <- as.data.frame(rownames(X))
  colnames(topl) <- c(paste("Top", top.number.count, sep = "", collapse = ""))
  write.table(topl, top.dir, sep="\t", quote=F, row.names = FALSE, col.names = TRUE)
  return()
}

#' Normalization of Single-Cell RNA-Seq Data
#'
#' This function normalizes single-cell RNA-seq data using its raw counts. The normalization is performed to remove variation due to difference in read depth and coverage.
#' @param dt.st The dataset to normalize
#' @keywords normalization
#' @export
#' @examples
#' normalize.dt(single.cell.mtx)
normalize.dt <- function(dt.st) {
  # Calculate column sums
  csums <- colSums(dt.st)
  # Calculate averages of the sums
  cavg <- mean(csums)
  # Calculate normalized bulk
  norm.dt.st <- dt.st
  for (i in 1:length(csums)) {norm.dt.st[,i] <- (norm.dt.st[,i]/csums[i]) * cavg}
  return(norm.dt.st)
}

#' Scale Ratio Calculation
#'
#' This function calculate the scale ratio between the single-cell transcriptome and reference dataset to minimize the Lagrangian multiplier. In the other word, reduce the restriction.
#' @param bulk The reference dataset
#' @param sc The single-cell dataset
#' @keywords Scale ratio calculation
#' @export
#' @examples
#' calc.scale.ratio(reference.dt, single.cell.mtx)
calc.scale.ratio <- function(bulk, sc) {
  # Calculate column sums
  csums.bulk <- colSums(bulk)
  csums.sc <- colSums(sc)
  # Calculate averages of the sums
  cavg.bulk <- mean(csums.bulk)
  cavg.sc <- mean(csums.sc)
  # Calculate scaling ratio for each tissue
  scale.ratio.sc <- cavg.bulk/cavg.sc
  return(scale.ratio.sc)
}

#' Find intersection genes
#'
#' This function identifies the intersection genes between reference and single-cell data
#' @param bulk The reference dataset
#' @param sc The single-cell dataset
#' @export
#' @examples
#' gene.intersect.sub(reference.dt, single.cell.mtx)
gene.intersect.sub <- function(bulk, sc) {
  # Find genes with all zero-expression within sc/bulk
  rsums.bulk <- rowSums(bulk)
  rsums.sc <- rowSums(sc)
  # Identify blank genes
  bulk.blank.genes <- which(rsums.bulk == 0)
  sc.blank.genes <- which(rsums.sc == 0)
  # Remove those genes
  norm.bulk.sub <- bulk
  norm.sc.sub <- sc
  if (length(bulk.blank.genes) > 0) {norm.bulk.sub <- bulk[-(bulk.blank.genes),]}
  if (length(sc.blank.genes) > 0) {norm.sc.sub <- sc[-(sc.blank.genes),]}
  # Find intersection genes
  genes.inter <- intersect(rownames(norm.bulk.sub), rownames(norm.sc.sub))
  # Find the intersected sub dataset of bulk (normalized bulk)
  norm.bulk.sub.sub <- norm.bulk.sub[genes.inter, ]
  norm.sc.sub.sub <- norm.sc.sub[genes.inter, ]
  return(list(norm.bulk.sub.sub, norm.sc.sub.sub))
}

#' Reference Construction
#'
#' This function constructs reference from single-cell resolution reference data to be used for quadratic programming calculation
#' @param sc The single-cell resolution dataset
#' @param sc.aux The auxiliary data frame that annotate the single-cell resolutiond dataset
#' @param criteria The column name to use for construction of the reference
#' @export
#' @examples
#' ref.construction(single.ref.mtx, single.aux.df, "cell.type")
ref.construction <- function(sc, sc.aux, criteria) {
  criteria.col <- sc.aux[, criteria]
  ref.df <- data.frame()
  uniq.crit <- unique(criteria.col)

  for (i in 1:length(uniq.crit)) {
    curr.crit <- uniq.crit[i]
    sc.bc <- colnames(sc)[which(sc.aux[, criteria] == curr.crit)]
    curr.sc.sub <- sc[, sc.bc]
    curr.df <- as.data.frame(rowSums(curr.sc.sub))
    if (ncol(ref.df) <= 0) {
      ref.df <- curr.df
    } else {
      ref.df <- cbind(ref.df, curr.df)
    }
    colnames(ref.df)[i] <- paste0(criteria, "_", curr.crit)
  }
  return(ref.df)
}



