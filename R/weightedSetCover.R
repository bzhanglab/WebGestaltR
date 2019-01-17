#' Weighted Set Cover
#'
#' Size constrained weighted set cover problem to find top N sets while
#' maximizing the coverage of all elements.
#'
#' @param idsInSet A list of set names and their member IDs.
#' @param costs A vector of the same length to add weights for penalty, i.e. 1/-logP.
#' @param topN The number of sets (or less when it completes early) to return.
#' @param nThreads The number of processes to use. In Windows, it fallbacks to 1.
#'
#' @return A list of \code{topSets} and \code{coverage}.
#' \describe{
#'  \item{topSets}{A list of set IDs.}
#'  \item{coverage}{The percentage of IDs covered in the top sets.}
#' }
#'
#' @importFrom dplyr bind_rows %>%
#' @importFrom parallel mclapply
#' @export
#' @author Zhiao Shi, Yuxing Liao
#'
weightedSetCover <- function(idsInSet, costs, topN, nThreads=4) {
  cat("Begin weighted set cover...\n")
  names(costs) <- names(idsInSet)
  if (.Platform$OS.type == "windows") {
    nThreads = 1
  }
  multiplier <- 10
  # we only start with the top (multiplier * topN) most
  # significant sets to reduce computational cost
  max_num_set <- multiplier * topN
  if (length(idsInSet) > max_num_set) {
    index <- order(costs)
    costs <- costs[index][1:max_num_set]
    idsInSet <- idsInSet[index][1:max_num_set]
  }

  s.hat <- 1.0
  # get all unique genes in all enriched sets
  all.genes <- unique(unlist(idsInSet))
  remain <- s.hat * length(all.genes)

  # final results, contains a list of gene set names
  cur.res <- c()
  # current candidates with marginal gain and size
  all.set.names <- names(idsInSet)
  mc_results <- mclapply(all.set.names, function(cur_name, cur_res, idsInSet, costs) {
      cur_gain <- marginalGain(cur_name, cur_res, idsInSet, costs)
      cur_size <- length(idsInSet[[cur_name]])
      return(data.frame(geneset.name=cur_name, gain=cur_gain, size=cur_size, stringsAsFactors=FALSE))
    }, cur_res=cur.res, idsInSet=idsInSet, costs=costs, mc.cores=nThreads)
  candidates <- mc_results %>% bind_rows()
  topN <- min(topN, nrow(candidates))
  for (i in seq(topN)) {
    # if there is no candidates, return
    if (nrow(candidates) == 0) {
      covered.genes <- unique(unlist(idsInSet[cur.res]))
      s.hat <- length(covered.genes) / length(all.genes)
      cat("No more candidates, ending weighted set cover\n")
      return(list(topSets=cur.res, coverage=s.hat))
    }
    # find the set with maximum marginal gain
    # tie breaker: for two sets with sname marginal gain, pick the one with
    # larger size
    candidates <- candidates[order(-candidates$gain, -candidates$size), ]
    # update remain
    remain <- remain - length(marginalBenefit(candidates[1, "geneset.name"], cur.res, idsInSet))
    cur.res <- c(cur.res, candidates[1,"geneset.name"])
    if (remain == 0) {
      covered.genes <- unique(unlist(idsInSet[cur.res]))
      s.hat <- length(covered.genes) / length(all.genes)
      cat("Remain is 0, ending weighted set cover\n")
      # full coverage solution
      return(list(topSets=cur.res, coverage=s.hat))
    }
    # update candidates
    # first remove the one just been selected
    candidates <- candidates[-1, ]
    # recalculate gain, remove rows with gain == 0
    mc_results <- mclapply(seq(nrow(candidates)), function(row, candidates, cur_res, idsInSet, costs){
         cur_name <- candidates[row, "geneset.name"]
         cur_gain <- marginalGain(cur_name, cur_res, idsInSet, costs)
         if(cur_gain != 0) {
           candidates[candidates$geneset.name == cur_name, "gain"] <- cur_gain
           tmp_candidate <- candidates[candidates$geneset.name == cur_name,]
           return(tmp_candidate)
         }
      }, candidates=candidates, cur_res=cur.res, idsInSet=idsInSet, costs=costs, mc.cores=nThreads)

    new_candidates <- mc_results %>% bind_rows()
    candidates <- new_candidates
  }
  # not fully covered, compute the current coverage and return
  covered.genes <- unique(unlist(idsInSet[cur.res]))
  s.hat <- length(covered.genes) / length(all.genes)
  cat("End weighted set cover...\n")
  return(list(topSets=cur.res, coverage=s.hat))
}


# return a list of genes from all.genes that has not been
# covered so far
# cur.set.name: name of the candidate gene set
# cur.res: vector of names of gene sets in current result
marginalBenefit <- function(cur.set.name, cur.res, idsInSet) {
  all.genes <- unique(unlist(idsInSet))
  cur.genes <- idsInSet[[cur.set.name]]
  if(length(cur.res) == 0){
    not.covered.genes <- cur.genes
  } else{
    covered.genes <- unique(unlist(idsInSet[cur.res]))
    not.covered.genes <- setdiff(cur.genes, covered.genes)
  }
  return(not.covered.genes)
}


marginalGain <- function(cur.set.name, cur.res, idsInSet, costs) {
  cur.cost <- costs[cur.set.name]
  cur.mben <- marginalBenefit(cur.set.name, cur.res, idsInSet)
  return(length(cur.mben) / cur.cost)
}
