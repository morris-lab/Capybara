#' Generate classification results from the binarized result matrix
#'
#' This function genrate classification table from the binarized result matrix
#' @param bin.count.rslt The binarized classification result matrix
#' @keywords calssification, result
#' @note 
#' @export
#' @examples
#' 
binary.to.classification <- function(bin.count.rslt) {
  class.rslt <- data.frame()
  bin.count.rsums <- rowSums(bin.count.rslt)
  for (i in 1:nrow(bin.count.rslt)) {
    curr.cell.bc <- rownames(bin.count.rslt)[i]
    if (bin.count.rsums[i] > 1) {
      curr.df <- data.frame(barcode = curr.cell.bc, call = "Multi_ID", stringsAsFactors = F)
    } else {
      if (bin.count.rsums[i] == 0) {
        curr.df <- data.frame(barcode = curr.cell.bc, call = "Unassigned", stringsAsFactors = F)
      } else {
        idx <- which(bin.count.rslt[i,] == 1)
        identity <- unlist(lapply(strsplit(colnames(bin.count.rslt), "frxn_cell.type_"), function(x) x[length(x)]))[idx]
        curr.df <- data.frame(barcode = curr.cell.bc, call = identity[length(identity)], stringsAsFactors = F)
      }
    }
    if (nrow(class.rslt) <= 0) {
      class.rslt <- curr.df
    } else {
      class.rslt <- rbind(class.rslt, curr.df)
    }
  }
  return(class.rslt)
}
