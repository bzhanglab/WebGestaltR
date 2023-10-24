#' @title kMedoid
#' @description kMedoid clustering
#' @param idsInSet a list of sets of ids
#' @param score a vector of scores for each set
#' @param maxK maximum number of clusters
#' @importFrom cluster pam
kMedoid <- function(idsInSet, score, maxK = 10) {
  # first find out the union of sets, sorted
  all.genes <- sort(unique(unlist(idsInSet)))
  overlap.mat <- sapply(idsInSet, function(x) {
    as.integer(all.genes %in% x)
  })

  num <- length(idsInSet)
  if (num <= maxK) {
    maxK <- num - 1
  }
  sim.mat <- matrix(1, num, num)
  colnames(sim.mat) <- colnames(overlap.mat)

  if (num == 1) {
    return(list(sim.mat = sim.mat, ip.vec = c(1)))
  }

  for (i in 1:(num - 1)) {
    for (j in (i + 1):num) {
      x <- sum(bitwOr(overlap.mat[, i], overlap.mat[, j]))
      if (x == 0) { # if there is no overlap, set the similarity to -Inf
        sim.mat[i, j] <- -Inf
        sim.mat[j, i] <- -Inf
      } else {
        jaccardIndex <- sum(bitwAnd(overlap.mat[, i], overlap.mat[, j])) / x
        sim.mat[i, j] <- jaccardIndex
        sim.mat[j, i] <- jaccardIndex
      }
    }
  }


  if (max(sim.mat) == min(sim.mat)) {
    # this will generate error, so randomy add some noise to off diagonal elements
    mat.siz <- dim(sim.mat)[1]
    rand.m <- matrix(rnorm(mat.siz * mat.siz, 0, 0.01), mat.siz)
    # make it symmetric
    rand.m[lower.tri(rand.m)] <- t(rand.m)[lower.tri(rand.m)]
    sim.mat <- sim.mat + rand.m
    # make diagonal all 1
    diag(sim.mat) <- 1
  }

  # compute the k-medoid clustering
  kmRes <- pam(sim.mat, maxK, diss = TRUE, variant = "faster", medoids = 1:maxK) # TODO: Make parameter for number of clusters. Currently set to 5.

  # sort clusters to make exemplar the first member
  clusters <- vector(mode = "list", length(kmRes$id.med))
  if (length(kmRes$id.med) == 0) {
    return(NULL)
  }
  cluster_info <- kmRes$silinfo[["widths"]][, 1]
  for (i in seq_along(kmRes$id.med)) {
    clusters[[i]] <- names(cluster_info[cluster_info == i])
  }
  # print(kmRes$medoids)
  return(list(clusters = clusters, representatives = kmRes$medoids))
}
