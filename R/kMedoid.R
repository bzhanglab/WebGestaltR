kMedoid <- function(idsInSet, score, maxK = 10){
    # first find out the union of sets, sorted
    all.genes <- sort(unique(unlist(idsInSet)))
    overlap.mat <- sapply(idsInSet, function(x) {as.integer(all.genes %in% x)})

    num <- length(idsInSet)
    if (num <= maxK) {
        maxK <- num - 1
    }
    sim.mat <- matrix(1, num, num)
    colnames(sim.mat) <- colnames(overlap.mat)

    if (num == 1) {
        return(list(sim.mat=sim.mat, ip.vec=c(1)))
    }

    for (i in 1:(num-1)) {
        for (j in (i+1):num) {
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
		rand.m <- matrix(rnorm(mat.siz*mat.siz,0,0.01),mat.siz)
		# make it symmetric
		rand.m[lower.tri(rand.m)] = t(rand.m)[lower.tri(rand.m)]
		sim.mat <- sim.mat + rand.m
		# make diagonal all 1
		diag(sim.mat) <- 1
	}

    # compute the k-medoid clustering
    kmRes <- pam(sim.mat, maxK, diss=TRUE, variant = "faster") # TODO: Make parameter for number of clusters. Currently set to 5.
    
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
