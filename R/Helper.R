#' Get the Most Connected Cells
#'
#' This function returns index of cells that are the most connected within a cell type
#' @param mtx The normalized reference count matrix
#' @param n.sample The number of reference cells to be included within each pseudo-bulk
get.most.connected <- function(mtx, n.sample) {
  corr.mtx <- WGCNA::cor(mtx)
  corr.mtx.upper <- corr.mtx * upper.tri(corr.mtx)
  corr.mtx.melt <- melt(corr.mtx.upper)
  corr.mtx.melt.pos <- corr.mtx.melt[which(corr.mtx.melt$value > 0), ]
  corr.mtx.melt.pos.sort <- corr.mtx.melt.pos[order(-corr.mtx.melt.pos$value), ]
  corr.mtx.melt.pos.sort$X1 <- as.character(corr.mtx.melt.pos.sort$X1)
  corr.mtx.melt.pos.sort$X2 <- as.character(corr.mtx.melt.pos.sort$X2)

  sample.list <- c()
  count.line <- 1
  while(length(sample.list) < n.sample) {
    sample.list <- unique(c(sample.list, unique(c(corr.mtx.melt.pos.sort$X1[count.line], corr.mtx.melt.pos.sort$X2[count.line]))))
    count.line <- count.line + 1
  }

  return(sample.list[1:n.sample])
}


#' Get the Least Connected Cells
#'
#' This function returns index of cells that are the least connected within a cell type
#' @param mtx The normalized reference count matrix
#' @param n.sample The number of reference cells to be included within each pseudo-bulk
get.lest.connected <- function(mtx, n.sample) {
  corr.mtx <- WGCNA::cor(mtx)
  corr.mtx.upper <- corr.mtx * upper.tri(corr.mtx)
  corr.mtx.melt <- melt(corr.mtx.upper)
  corr.mtx.melt.pos <- corr.mtx.melt[which(corr.mtx.melt$value > 0), ]
  corr.mtx.melt.pos.sort <- corr.mtx.melt.pos[order(corr.mtx.melt.pos$value), ]
  corr.mtx.melt.pos.sort$X1 <- as.character(corr.mtx.melt.pos.sort$X1)
  corr.mtx.melt.pos.sort$X2 <- as.character(corr.mtx.melt.pos.sort$X2)

  sample.list <- c()
  count.line <- 1
  while(length(sample.list) < n.sample) {
    sample.list <- unique(c(sample.list, unique(c(corr.mtx.melt.pos.sort$X1[count.line], corr.mtx.melt.pos.sort$X2[count.line]))))
    count.line <- count.line + 1
  }

  return(sample.list[1:n.sample])
}

#' Get the Medium Connected Cells
#'
#' This function returns index of cells that are the most connected within a cell type
#' @param mtx The normalized reference count matrix
#' @param n.sample The number of reference cells to be included within each pseudo-bulk
get.mid.connected <- function(mtx, n.sample) {
  corr.mtx <- WGCNA::cor(mtx)
  corr.mtx.upper <- corr.mtx * upper.tri(corr.mtx)
  corr.mtx.melt <- melt(corr.mtx.upper)
  corr.mtx.melt.pos <- corr.mtx.melt[which(corr.mtx.melt$value > 0), ]
  corr.mtx.melt.pos.sort <- corr.mtx.melt.pos[order(corr.mtx.melt.pos$value), ]
  corr.mtx.melt.pos.sort$X1 <- as.character(corr.mtx.melt.pos.sort$X1)
  corr.mtx.melt.pos.sort$X2 <- as.character(corr.mtx.melt.pos.sort$X2)

  corr.mtx.melt.pos.sort.most <- corr.mtx.melt.pos[order(-corr.mtx.melt.pos$value), ]
  corr.mtx.melt.pos.sort.most$X1 <- as.character(corr.mtx.melt.pos.sort.most$X1)
  corr.mtx.melt.pos.sort.most$X2 <- as.character(corr.mtx.melt.pos.sort.most$X2)

  sample.list <- c()
  count.line <- 1
  while(length(sample.list) < n.sample/2) {
    sample.list <- unique(c(sample.list, unique(c(corr.mtx.melt.pos.sort$X1[count.line], corr.mtx.melt.pos.sort$X2[count.line]))))
    count.line <- count.line + 1
  }
  count.line <- 1
  while(length(sample.list) < n.sample) {
    sample.list <- unique(c(sample.list, unique(c(corr.mtx.melt.pos.sort.most$X1[count.line], corr.mtx.melt.pos.sort.most$X2[count.line]))))
    count.line <- count.line + 1
  }

  return(sample.list[1:n.sample])
}


# Create the function.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



