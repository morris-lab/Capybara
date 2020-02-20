#' Systematic construction of a high-resolution reference
#'
#' Create a pseudo-bulk reference by sampling 90-cells from each cell type to maintain cellular resolution while increasing transcriptional resolution
#' @param ref.mtx The single-cell reference dataset
#' @param coldata.df The metadata (cell type information) for cells in the high-resolution reference
#' @param cell.num.for.ref The number of cell numbers used to build the reference for each cell type. The default is cell.num.for.ref = 90
#' @keywords pseudo-bulk reference
#' @note
#' @export
#' @examples
#'
construct.high.res.reference <- function(ref.mtx, coldata.df, criteria,cell.num.for.ref = 90) {

  ct.freq <- as.data.frame(table(coldata.df[,criteria]), stringsAsFactors = F)
  ref.meta <- data.frame()
  ref.sc <- data.frame()
  mca.counts.all.involved.sub<-ref.mtx
  for (j in 1:nrow(ct.freq)) {
    curr.ct.at.test <- as.character(ct.freq[j,1])
    curr.cell.involve <- rownames(coldata.df)[which(coldata.df[,criteria] == curr.ct.at.test)]
    if (ct.freq$Freq[j] >= cell.num.for.ref) {
      sample.ref.cell <- c(get.most.connected(log1p(normalize.dt(mca.counts.all.involved.sub[,curr.cell.involve])), cell.num.for.ref/2),
                           get.least.connected(log1p(normalize.dt(mca.counts.all.involved.sub[,curr.cell.involve])), cell.num.for.ref/2))
    } else {
      sample.ref.cell <- sample(curr.cell.involve, size = cell.num.for.ref, replace = T)
    }
    new.bc <- paste0("Cell_", seq(((j-1) * cell.num.for.ref + 1), j * cell.num.for.ref))
    curr.meta <- data.frame(row.names = new.bc, cell.type = curr.ct.at.test, cell.bc = sample.ref.cell, stringsAsFactors = F)
    curr.sc <- mca.counts.all.involved.sub[, sample.ref.cell]
    colnames(curr.sc) <- new.bc
    if (ncol(ref.sc) <= 0) {
      ref.meta <- curr.meta
      ref.sc <- curr.sc
    } else {
      ref.meta <- rbind(ref.meta, curr.meta)
      ref.sc <- cbind(ref.sc, curr.sc)
    }
  }
  ref.df <- ref.construction(ref.sc, ref.meta, "cell.type")
  return(list(ref.sc, ref.meta, ref.df))
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

