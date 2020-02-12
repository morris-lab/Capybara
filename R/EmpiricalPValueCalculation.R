#' Empirical p-value of Identity Score Calculation
#'
#' This function resamples from a sample dataset and returns an empirical p-value of identity score
#' @param dens.x The background dataset that we compare the identity score of our sample with
#' @param curr.val The identity score of the sample dataset to be evaluated its significance 
#' @param n The number of times of resampling. The default is n=1000
#' @keywords resampling, empirical p-value
#' @note 
#' @export
#' @examples
#' 
sample.func <- function(dens.x, curr.val, prob, n = 1000) {
  samp <- sample(1:length(dens.x), n, replace = T, prob = prob)
  samp.val <- dens.x[samp]
  perc <- sum(samp.val > curr.val)/n
  return(perc)
}

#' Emprical p-values for query cells
#'
#' This function returns a list of emprical p-value matrices. Every matrix in the list contains emprical p-values of all possilbe celltypes for a single query cells under test.times resampled backgrounds.
#' @param mtx The matrix of identity scores of query cells. The row number is the total number of query cells and the column number is the number of total number of possible cell types
#' @param bkgd.mtx The matrix of identity score of reference cells. The row number is the total number of reference cells and the column number is the number of total number of possible cell types
#' @keywords resampling, empirical p-value, all query cells
#' @note 
#' @export
#' @examples
#' 
percentage.calc <- function(mtx, bkgd.mtx) {
  cnms <- colnames(mtx)
  cell.names <- rownames(mtx)
  mtx.cell.num <- nrow(mtx)
  mtx.ct.num <- ncol(mtx)
  test.times <- 50
  
  perc.ls <- pbmclapply(as.list(seq_len(mtx.cell.num)), perc.calc.aux, cell.names = cell.names, mtx.ct.num = mtx.ct.num, 
                        background.mtx = bkgd.mtx, mtx = mtx, cnms = cnms, mc.cores = 4)
  
  return(perc.ls)
}

#' Emprical p-values for a query cell
#'
#' This function returns emprical p-values of all possilbe celltypes for a query cells under test.times resampled backgrounds.
#' @param l The index of the sample cell of query 
#' @param cell.names A vector of sample cells' names
#' @param mtx.ct.num The number of all possible cell types
#' @param background.mtx The matrix of identity score of reference cells. The row number is the total number of reference cells and the column number is the number of total number of possible cell types
#' @param mtx The matrix of identity scores of query cells. The row number is the total number of query cells and the column number is the number of total number of possible cell types
#' @param cnms The names of all possible cell types
#' @param test.times The number of times resampling the background. The default is test.times=50
#' @keywords resampling, empirical p-value, single query
#' @note 
#' @export
#' @examples
#' 
perc.calc.aux <- function(l, cell.names, mtx.ct.num, background.mtx, mtx, cnms, test.times = 50) {
  curr.cell <- cell.names[l]
  perc.vec <- lapply(seq(1, mtx.ct.num),
                     function(x) {
                       dens <- density(background.mtx[,x], n = nrow(mtx))
                       curr.num <- mtx[l, x]
                       if (curr.num <= .Machine$double.eps * 2) {
                         pc.v <- rep(1, test.times)
                       } else {
                         pc.v <- apply(as.data.frame(seq_len(test.times)), 1, function(x) sample.func(dens.x = dens$x, curr.val = curr.num, prob = dens$y, n = 1000))
                       }
                       return(list(cnms[x], pc.v))
                     })
  comb.perc <- rbindlist(perc.vec)
  comb.perc$label <- rep(seq(1,test.times), length(unique(comb.perc$V1)))
  
  comb.df <- dcast(comb.perc, label~V1, value.var = "V2")
  comb.df <- comb.df[, -c(1)]
  curr.list <- list()
  curr.list[[curr.cell]] <- comb.df
  return(curr.list)
}




