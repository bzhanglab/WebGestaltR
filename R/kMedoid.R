kMedoid <- function(idsInSet, score){
    # first find out the union of sets, sorted
    all.genes <- sort(unique(unlist(idsInSet)))
    overlap.mat <- sapply(idsInSet, function(x) {as.integer(all.genes %in% x)})

    num <- length(idsInSet)
    sim.mat <- matrix(1, num, num)
    colnames(sim.mat) <- colnames(overlap.mat)

    if (num == 1) {
        return(list(sim.mat=sim.mat, ip.vec=c(1)))
    }

    for (i in 1:(num-1)) {
        for (j in (i+1):num) {
            jaccardIndex <- sum(bitwAnd(overlap.mat[, i], overlap.mat[, j])) / sum(bitwOr(overlap.mat[, i], overlap.mat[, j]))
            sim.mat[i, j] <- jaccardIndex
            sim.mat[j, i] <- jaccardIndex
        }
    }
    # if there is no overlap, set the similarity to -Inf
    for (i in 1:(num-1)) {
        for (j in (i+1):num) {
            if (sum(bitwOr(overlap.mat[, i], overlap.mat[, j])) == 0) {
                sim.mat[i, j] <- -Inf
                sim.mat[j, i] <- -Inf
            }
        }
    }
    # compute the k-medoid clustering
    kmRes <- pam(sim.mat, 5, diss=TRUE) # TODO: Make parameter for number of clusters. Currently set to 5.
    
    #sort clusters to make exemplar the first member
    clusters <- vector(mode="list", length(kmRes$medoids))
    if(length(kmRes$medoids) == 0){
        return(NULL)
    }
    for (i in 1:length(clusters)) {
        clusters[[i]] <- kmRes$clustering[[i]][order(kmRes$clustering[[i]] == i, decreasing=TRUE)]
    }
    # print(kmRes$medoids)
    return(list(clusters=sapply(clusters, names), representatives=kmRes$medoids))
}