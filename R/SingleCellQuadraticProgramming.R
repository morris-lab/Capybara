#' Single-Cell Quadratic Programming Calculation
#'
#' This function runs quadratic programming to identify the probability of the cells in single-cell RNA-seq belonging to cell types in the reference transcriptome
#' @param bulk.transcriptome The reference transcriptome taht contains the transcriptome of each potential cell type
#' @param single.cell.transcriptome the transcriptome profile of cells from single-cell RNA-sequencing
#' @param force.eq either 0 or 1. Setting to 0 assumes the 1st constraint as inequality. Setting to 1 assumes equality. Default to 0
#' @param unix.parallel boolean value, either TRUE or FALSE. If using unix/linux based systems, this command can be set to TRUE to parallelize use parallel package. Default to FALSE
#' @param windows.parallel boolean value, either TRUE or FALSE. If using Windows based systems, this command can be set to TRUE to parallelize use snow package. Default to FALSE
#' @param parallel.cores the number of cores to use for parallel processes. If no parallelization selected, no parallelization will be implemented. Only 1 core will be used
#' @keywords quadratic programming, scRNA-seq
#' @note Code reference from Treutlein et. al., Dissecting direct reprogramming from fibroblast to neuron using single-cell RNA-seq
#' @export
#' @examples
#' sc.quad.prog.run(ref.transcriptome, sc.transcriptome, force.eq = 1)
single.round.QP.analysis <- function(ref, sc.data, scale.bulk.sc = "scale", n.cores = 1, save.to.path, save.to.filename, norm = T, logging = T) {
  # Normalized the bulk and single-cell data
  if (norm) {
    norm.bulk <- normalize.dt(ref)
    norm.sc.mtx <- normalize.dt(sc.data)
  } else {
    norm.bulk <- as.matrix(ref)
    norm.sc.mtx <- as.matrix(sc.data)
  }


  # Calculate the scaling ratio
  scale.ratio.sc <- calc.scale.ratio(ref, sc.data)
  scale.norm.sc <- norm.bulk/scale.ratio.sc

  # Calculate gene intersections
  norm.ls.sc <- gene.intersect.sub(norm.bulk, norm.sc.mtx)
  scale.ls.sc <- gene.intersect.sub(scale.norm.sc, norm.sc.mtx)

  # Finish the log-normalization
  if (logging) {
    log.norm.ls.sc <- lapply(norm.ls.sc, log1p)
    log.scale.ls.sc <- lapply(scale.ls.sc, log1p)
  } else {
    log.norm.ls.sc <- norm.ls.sc
    log.scale.ls.sc <- scale.ls.sc
  }
  # Run QP and save to file
  if (scale.bulk.sc == "scale") {
    qp.rslt <- sc.quad.prog.run(as.matrix(log.scale.ls.sc[[1]]),
                                single.cell.transcriptome = log.scale.ls.sc[[2]],
                                unix.parallel = TRUE, parallel.cores = n.cores, force.eq = 1)
    write.csv(qp.rslt, paste0(save.to.path, save.to.filename, "_scale.csv"), quote = F, row.names = F)
  } else {
    if (scale.bulk.sc == "non-scale") {
      qp.rslt <- sc.quad.prog.run(as.matrix(log.norm.ls.sc[[1]]),
                                  single.cell.transcriptome = log.norm.ls.sc[[2]],
                                  unix.parallel = TRUE, parallel.cores = n.cores, force.eq = 1)
      write.csv(qp.rslt, paste0(save.to.path, save.to.filename, "_non_scale.csv"), quote = F, row.names = F)
    } else {
      qp.rslt.scl <- sc.quad.prog.run(as.matrix(log.scale.ls.sc[[1]]),
                                      single.cell.transcriptome = log.scale.ls.sc[[2]],
                                      unix.parallel = TRUE, parallel.cores = n.cores, force.eq = 1)
      qp.rslt.non.scl <- sc.quad.prog.run(as.matrix(log.norm.ls.sc[[1]]),
                                          single.cell.transcriptome = log.norm.ls.sc[[2]],
                                          unix.parallel = TRUE, parallel.cores = n.cores, force.eq = 1)
      write.csv(qp.rslt.scl, paste0(save.to.path, save.to.filename, "_scale.csv"), quote = F, row.names = F)
      write.csv(qp.rslt.non.scl, paste0(save.to.path, save.to.filename, "_non_scale.csv"), quote = F, row.names = F)
    }
  }
}

#' Single-Cell Quadratic Programming Calculation
#'
#' This function runs quadratic programming to identify the probability of the cells in single-cell RNA-seq belonging to cell types in the reference transcriptome
#' @param bulk.transcriptome The reference transcriptome taht contains the transcriptome of each potential cell type
#' @param single.cell.transcriptome the transcriptome profile of cells from single-cell RNA-sequencing
#' @param force.eq either 0 or 1. Setting to 0 assumes the 1st constraint as inequality. Setting to 1 assumes equality. Default to 0
#' @param unix.parallel boolean value, either TRUE or FALSE. If using unix/linux based systems, this command can be set to TRUE to parallelize use parallel package. Default to FALSE
#' @param windows.parallel boolean value, either TRUE or FALSE. If using Windows based systems, this command can be set to TRUE to parallelize use snow package. Default to FALSE
#' @param parallel.cores the number of cores to use for parallel processes. If no parallelization selected, no parallelization will be implemented. Only 1 core will be used
#' @keywords quadratic programming, scRNA-seq
#' @note Code reference from Treutlein et. al., Dissecting direct reprogramming from fibroblast to neuron using single-cell RNA-seq
#' @export
#' @examples
#' sc.quad.prog.run(ref.transcriptome, sc.transcriptome, force.eq = 1)
sc.quad.prog.run <- function(bulk.transcriptome, single.cell.transcriptome, force.eq = 0,
                             unix.parallel = FALSE, windows.parallel = FALSE, parallel.cores = 4) {
  # Check if the sizes of the matrices are comfortable with each other
  if (nrow(bulk.transcriptome) != nrow(single.cell.transcriptome)) {
    print("The number of genes included in the bulk expression is not the same as the single cell transcriptome")
    return()
  }

  # Initialize a data frame to hold the fraction identity matrix
  identity.matx <- data.frame()
  # Initialize the cell types given
  given.cell.typs <- colnames(bulk.transcriptome)

  # If using unix/linux based system, can use this process to parallel the the program
  # Set unix.parallel parameter to be TRUE, but set windows.parallel to be FALSE
  if (unix.parallel) {
    identy.mx.ls <- mclapply(as.list(seq_len(ncol(single.cell.transcriptome))), quad.prog.calc,
                             bulk.transcriptome = bulk.transcriptome,
                             single.cell.transcriptome = single.cell.transcriptome,
                             force.eq = force.eq, mc.cores = parallel.cores)
    identity.ls <- lapply(identy.mx.ls, function(identity) c(identity[[2]], identity[[1]]$solution,
                                                             identity[[1]]$Lagrangian[1], identity[[3]]))

    identity.matx <- matrix(unlist(identity.ls), byrow = TRUE,
                            nrow = length(identity.ls), ncol = length(identity.ls[[1]]))
  } else {
    # If using Windows based system, can use this process to parallel the the program
    # Set windows.parallel parameter to be TRUE, but set unix.parallel to be FALSE
    # NOTE: This might be slower than applying parallelization under Unix/Linux based system
    if (windows.parallel) {
      cl<-makeCluster(parallel.cores)
      identy.mx.ls <- clusterApply(cl, seq_len(ncol(single.cell.transcriptome)), quad.prog.calc,
                                   bulk.transcriptome = bulk.transcriptome,
                                   single.cell.transcriptome = single.cell.transcriptome,
                                   force.eq = force.eq)
      stopCluster(cl)
      identity.ls <- lapply(identy.mx.ls, function(identity) c(identity[[2]], identity[[1]]$solution,
                                                               identity[[1]]$Lagrangian[1], identity[[3]]))

      identity.matx <- matrix(unlist(identity.ls), byrow = TRUE,
                              nrow = length(identity.ls), ncol = length(identity.ls[[1]]))

      # If not using parallel for this case, just run the following process in a loop for serial
      # All parallel parameters option should be set to FALSE
    } else {
      for (i in 1:ncol(single.cell.transcriptome)){
        identity<-c()
        quad.rslt <- quad.prog.calc(i, bulk.transcriptome, single.cell.transcriptome, force.eq)
        QP <- quad.rslt[[1]]
        Error <- quad.rslt[[3]]

        identity <- c(colnames(single.cell.transcriptome)[i], QP$solution, QP$Lagrangian[1],Error)

        if (nrow(identity.matx) == 0) {
          identity.matx <- as.data.frame(t(as.matrix(identity)))
        } else {
          identity.matx <- rbind(identity.matx, as.data.frame(t(as.matrix(identity))))
        }
      }
    }
  }

  # Set the column names of the returning data frame
  col.frx.names <- paste("frxn_", given.cell.typs, sep = "")
  colnames(identity.matx)<-c("cell_name", col.frx.names,"Lagrangian","Error")
  return(identity.matx)
}

quad.prog.calc <- function(col.num, bulk.transcriptome, single.cell.transcriptome, force.eq) {
  library(quadprog)
  Y <- as.matrix(single.cell.transcriptome[, col.num])
  Rinv <- solve(chol(t(bulk.transcriptome) %*% bulk.transcriptome));
  C <- cbind(rep(-1, ncol(bulk.transcriptome)), diag(ncol(bulk.transcriptome)))
  n <- nrow(Rinv)
  q <- ncol(bulk.transcriptome)
  b <- c(-1, rep(0, q))
  d <- t(Y) %*% bulk.transcriptome

  QP<-solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, bvec = b, Amat = C, meq = force.eq)
  Error<-sum(abs(Y- bulk.transcriptome %*% QP$solution))

  return(list(QP, colnames(single.cell.transcriptome)[col.num], Error))
}

