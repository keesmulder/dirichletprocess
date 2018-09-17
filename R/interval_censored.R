
IntervalCensoredDraw <- function(dpobj) UseMethod("IntervalCensoredDraw")

IntervalCensoredDraw.vonMises <- function(dpobj) {

  y <- dpobj$data
  y_imp <- y[, 1]

  clusterLabels <- dpobj$clusterLabels
  clusterParams <- dpobj$clusterParameters

  is_censored <- which(y[, 1] != y[,2])

  cens_cl    <- clusterLabels[is_censored]
  cens_clust <- unique(cens_cl)


  for (clust_i in cens_clust) {
    y_this_clust <- which(clusterLabels == clust_i)
    y_imp[y_this_clust] <- rvmc(length(y_this_clust),
                          clusterParams[[1]][, , clust_i, drop = FALSE],
                          clusterParams[[2]][, , clust_i, drop = FALSE]
                          )
  }

  return(matrix(y_imp))
}
