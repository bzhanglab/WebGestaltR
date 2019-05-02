#' Affinity Propagation
#'
#' Use affinity propagation to cluster similar gene sets to reduce redundancy in report.
#'
#' @param idsInSet A list of set names and their member IDs.
#' @param score A vector of addible scores with the same length used to assign input preference;
#'  higher score has larger weight, i.e. -logP.
#'
#' @return A list of \code{clusters} and \code{representatives} for each cluster.
#' \describe{
#'  \item{clusters}{A list of character vectors of set IDs in each cluster.}
#'  \item{representatives}{A character vector of representatives for each cluster.}
#' }
#'
#' @export
#' @importFrom apcluster apcluster
#' @importFrom stats rnorm
#' @author Zhiao Shi, Yuxing Liao
affinityPropagation <- function(idsInSet, score) {
	cat("Begin affinity propagation...\n")
	# compute the similiarity and input preference vector
	ret <- jaccardSim(idsInSet, score)

	sim.mat <- ret$sim.mat
	ip.vec <- ret$ip.vec

	apRes <- apcluster(sim.mat,p=ip.vec)
	#sort clusters to make exemplar the first member
	clusters <- vector(mode="list", length(apRes@clusters))
	if(length(apRes@clusters) == 0){
		return(NULL)
	}
	for (i in 1:length(apRes@clusters)) {
		exemplar <- apRes@exemplars[[i]]
		clusters[[i]] <- apRes@clusters[[i]][order(apRes@clusters[[i]] == exemplar, decreasing=TRUE)]
	}
	cat("End affinity propagation...\n")
	return(list(clusters=sapply(clusters, names), representatives=names(apRes@exemplars)))
}

#' Jaccard Similarity
#'
#' Calculate Jaccard Similarity.
#'
#' @inheritParams affinityPropagation
#'
#' @return A list of similarity matrix \code{sim.mat} and input preference vector \code{ip.vec}.
#'
#' @importFrom stats median
#' @author Zhiao Shi, Yuxing Liao
jaccardSim <- function(idsInSet, score){
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
	sim.mat[sim.mat == 0] <- -Inf
	# check sim.mat to see if it is identical for each pair
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
	# set the input preference (IP) for each gene set
	# give higher IP to gene set with larger -logP (remove sign)
	# IP <- maxScore for gene set with largest -logP value
	# IP <- minScore for gene set with smallest -logP value
	# other gene sets will have linearly interpolated IP value
	max.sig <- max(score)
	min.sig <- min(score)

	minScore <- 0
	tmp.sim.mat <- sim.mat
	tmp.sim.mat[!is.finite(tmp.sim.mat)] <- NA
	# get the median excluding -Inf
	maxScore <- median(tmp.sim.mat, na.rm=TRUE)
	if (abs(max.sig - min.sig) < .Machine$double.eps^0.5) {
		ip.vec <- NA
	} else{
		ip.vec <- minScore + (maxScore-minScore) * (score-min.sig) / (max.sig-min.sig)
	}
	return(list(sim.mat=sim.mat, ip.vec=ip.vec))
}
