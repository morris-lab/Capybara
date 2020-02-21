#' Multi-ID Score-Based Filter
#'
#' This function filters the multiple identity cells based on their QP scores, where we assume that a low QP score (less than 10E-3) are not a true identity to consider
#' @param binary.counts The binary count matrix, which is the output from binarization with Mann Whitney. 
#' @param classification The classification result, which is the output from binary to classification.
#' @param qp.matrix The matrix that contains QP scores calculated for the sample cells
#' @param qp.threshold The threshold to cut off for the QP scores in the multiple identity listed cells
#' @return A list contain 2 elements, the first is the curated and filtered multiple identity data frame and the second is the new classification data frame.
#' @keywords Multiple identities
#' @export
#' @examples
#' multi.id.curate.qp(multi.id.meta)
multi.id.curate.qp <- function(binary.counts, classification, qp.matrix, qp.threshold = 10^-3) {
  
  bin.count.std.multi <- as.data.frame(binary.counts[classification$barcode[which(classification$call == "Multi_ID")], ])
  bin.count.std.multi$cell.bc <- rownames(bin.count.std.multi)
  
  bin.count.std.multi.melt <- reshape2::melt(bin.count.std.multi)
  bin.count.std.multi.melt <- bin.count.std.multi.melt[which(bin.count.std.multi.melt$value > 0), ]
  
  classification.new <- classification
  
  qp.column <- c()
  for (i in 1:nrow(bin.count.std.multi.melt)) {
    curr.cell <- as.character(bin.count.std.multi.melt$cell.bc[i])
    curr.cell.type <- as.character(bin.count.std.multi.melt$variable[i])
    qp.column <- c(qp.column, qp.matrix[curr.cell, curr.cell.type])
  }
  bin.count.std.multi.melt$qp.score <- qp.column
  bin.count.std.multi.melt.remain <- bin.count.std.multi.melt[which(bin.count.std.multi.melt$qp.score >= qp.threshold), ]
  
  cell.occurrence <- as.data.frame(table(bin.count.std.multi.melt.remain$cell.bc))
  single.id.cell <- as.character(cell.occurrence$Var1[which(cell.occurrence$Freq == 1)])
    
  actual.multi <- bin.count.std.multi.melt.remain
  if (length(single.id.cell) > 0) {
    single.id.cell.bin <- bin.count.std.multi.melt.remain[which(bin.count.std.multi.melt.remain$cell.bc %in% single.id.cell), ]
    rownames(single.id.cell.bin) <- single.id.cell.bin$cell.bc
    classification.new[rownames(single.id.cell.bin), "call"] <- unlist(lapply(strsplit(as.character(single.id.cell.bin$variable), "frxn_cell.type_"), function(x) x[2]))
    actual.multi <- bin.count.std.multi.melt.remain[-which(bin.count.std.multi.melt.remain$cell.bc %in% single.id.cell), ]
  }
  
  return(list(actual.multi, classification.new))
}

