#' Transition Score Calculation
#'
#' This function calculates the transition scores for each cell state that is connected by cells with multiple identities.
#' @param multi.id.meta A data frame that contains cells with multiple identities. Column 1 - cell barcode, Column 2 - cell type call, Column 3 - counts, Column 4 - corresponding QP scores.
#' @return A data frame that contains the calculated transition scores for each identity
#' @keywords Transition metric
#' @export
#' @examples
#' transition.score(multi.id.meta)
transition.score <- function(multi.id.meta) {
  
  unique.ct <- unique(as.character(multi.id.meta[,2]))
  ct.entropy <- data.frame()
  for (i in 1:length(unique.ct)) {
    curr.cell <- unique.ct[i]
    actual.sub <- multi.id.meta[which(as.character(multi.id.meta[,2]) == curr.cell), ]
    curr.entro <- sum(-actual.sub$qp.score*log(actual.sub$qp.score, base = max(2, nrow(actual.sub))))
    curr.df <- data.frame(cell.ct = curr.cell, entropy = curr.entro, stringsAsFactors = F)
    if (nrow(ct.entropy) <= 0) {
      ct.entropy <- curr.df
    } else {
      ct.entropy <- rbind(ct.entropy, curr.df)
    }
  }
  
  rownames(ct.entropy) <- unlist(lapply(strsplit(ct.entropy$cell.ct, "frxn_cell.type_"), function(x) x[2]))
  
  return(ct.entropy)
}

