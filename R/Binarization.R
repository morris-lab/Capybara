#' Binarization and Identity Calling from Identity Score
#'
#' This function calls single identity or multiple identities of query cells from empirical p-values. Inferred cell types are marked by 1 in a binarized matrix.
#' @param mtx The matrix of identity scores of query cells. The row number is the total number of query cells and the column number is the number of total number of possible cell types
#' @param ref.perc.ls Emprical p-values for reference cells
#' @param ref.meta The celltype meta information for reference cells
#' @param perc.ls Emprical p-values for query cells
#' @param bulk If the reference data type is bulk RNA-seq. The default is bulk = FALSE
#' @param map.df bulk mapping. The default is bulk = FALSE
#' @keywords binarization, identity calling
#' @note
#' @export
#' @examples
#'
binarization.mann.whitney <- function(mtx, ref.perc.ls, ref.meta, perc.ls, bulk = FALSE, map.df = NULL) {
  p_vals <- matrix(nrow = nrow(mtx), ncol = (ncol(mtx) - 1))
  type_rank_in_order_mtx <- c()
  which_type <- c()

  all.indx <- seq(1,length(perc.ls))
  # ref.indx <- which(unlist(lapply(perc.ls, function(x) x[[1]])) %in% rownames(ref.meta))
  # test.indx <- setdiff(all.indx, ref.indx)

  benchmark <- list()
  benchmark.final <- c()

  for (i in 1:length(ref.perc.ls)) {
    curr.cell <- names(ref.perc.ls[[i]])
    perc_mtx <- as.data.frame(ref.perc.ls[[i]][[1]])
    if (curr.cell %in% rownames(ref.meta)) {
      curr.ct <- ref.meta[curr.cell, "cell.type"]
      if (endsWith(curr.ct, " ")) curr.ct <- substr(curr.ct, start = 1, stop = (nchar(curr.ct) - 1))
      curr.cell.cell.type <- gsub(" ", ".", curr.ct)
      curr.cell.cell.type <- gsub("−", ".", curr.cell.cell.type)
      curr.cell.cell.type <- gsub("-", ".", curr.cell.cell.type)
      curr.cell.cell.type <- gsub("/", ".", curr.cell.cell.type)
      curr.cell.cell.type <- gsub("&", ".", curr.cell.cell.type)
      curr.cell.cell.type <- gsub("\\(", ".", curr.cell.cell.type)
      curr.cell.cell.type <- gsub("\\)", ".", curr.cell.cell.type)
      if (!bulk) {
        curr.bm <- paste0("frxn_cell.type_", curr.cell.cell.type)
      } else {
        if (is.null(map.df)) {
          stop("No bulk mapping found!")
        }
        curr.bm <- map.df[which(map.df$tissue == curr.ct), "corresponding"]
      }
      curr.celltype.perm <- perc_mtx[, curr.bm]
      if (length(benchmark) <= 0) {
        benchmark[[curr.bm]] <- curr.celltype.perm
      } else {
        benchmark[[curr.bm]] <- c(benchmark[[curr.bm]], curr.celltype.perm)
      }
    }
  }

  for (i in 1:length(benchmark)) {
    curr.name <- names(benchmark)[i]
    curr.benchmark <- benchmark[[i]]
    curr.benchmark <- curr.benchmark[which(curr.benchmark < 1)]
    curr.outliers <- boxplot.stats(curr.benchmark)$out
    higher.end.outliers <- curr.outliers[which(curr.outliers > getmode(curr.benchmark))]
    if (length(higher.end.outliers) > 0) benchmark[[i]] <- curr.benchmark[which(curr.benchmark <= min(higher.end.outliers))]
    benchmark.final[curr.name] <- max(benchmark[[i]][which(benchmark[[i]] < 1)])
  }
  # benchmark.final <- rep(quantile(benchmark.final, 0.9), length(benchmark.final))
  # names(benchmark.final) <- names(benchmark)

  for (i in all.indx) {
    perc_mtx<-as.matrix((perc.ls[[i]][[1]]))
    vec<-as.vector(perc_mtx)
    rank<-rank(vec,ties.method = "average")
    type_rank<-as.matrix(t(colSums(matrix(rank, nrow=50))))
    colnames(type_rank)<-colnames(perc_mtx)
    type_rank_in_order<-colnames(type_rank)[order(type_rank,decreasing = F)]
    # type_rank_in_order[which(endsWith(type_rank_in_order, "."))] <- substr(type_rank_in_order[which(endsWith(type_rank_in_order, "."))], 1, (nchar(type_rank_in_order[which(endsWith(type_rank_in_order, "."))]) - 1))
    type_rank_in_order_mtx<-rbind(type_rank_in_order_mtx,type_rank_in_order)

    if(mean(perc_mtx[,type_rank_in_order[1]]) >= (benchmark.final[type_rank_in_order[1]] + 0.01)) {
      which_type <- c(which_type, 0)
    } else {
      for (j in 1:(length(colnames(type_rank))-1)) {

        a<-type_rank_in_order[1]
        b<-type_rank_in_order[j+1]

        test<-wilcox.test(x=perc_mtx[,a],y=perc_mtx[,b],alternative = "less")
        p_vals[i,j]<-test$p.value
        if(test$p.value<.05){
          which_type<-c(which_type,j)
          break
        } else {
          if (j == (ncol(type_rank_in_order_mtx) - 1)) {
            which_type <- c(which_type, j+1)
          }
        }
      }
    }
  }

  type_rank_in_order_mtx <- as.matrix(type_rank_in_order_mtx)
  binary.mtx<-matrix(0,nrow = nrow(mtx), ncol = ncol(mtx))
  colnames(binary.mtx)<-colnames(mtx)
  cellnames<-c()
  for (i in 1:length(which_type)) {
    cellnames <- c(cellnames, names(perc.ls[[i]]))
    if (which_type[i] > 0) {
      for (j in 1:which_type[i]) {
        binary.mtx[i,type_rank_in_order_mtx[i,j]] <- 1
      }
    }
  }
  rownames(binary.mtx)<-cellnames

  return(binary.mtx)
}

