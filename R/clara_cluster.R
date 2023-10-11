clara_cluster <- function(idsInSet, score){
    all.genes <- sort(unique(unlist(idsInSet)))
    overlap.mat <- sapply(idsInSet, function(x) {as.integer(all.genes %in% x)})
    kmRes <- clara(overlap.mat, 5, metric="manhattan", stand=FALSE, samples=1000, pamLike=TRUE)

    #sort clusters to make exemplar the first member
    clusters <- vector(mode="list", length(kmRes$medoids))
    # print(kmRes$clusinfo)
    if(length(kmRes$medoids) == 0){
        return(NULL)
    }
    for (i in 1:length(clusters)) {
        clusters[[i]] <- kmRes$clustering[[i]][order(kmRes$clustering[[i]] == i, decreasing=TRUE)]
    }
    # print(kmRes$medoids)
    return(list(clusters=sapply(clusters, names), representatives=kmRes$medoids))
}
