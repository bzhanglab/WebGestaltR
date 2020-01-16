#' Wrapper for batch WebGestaltR runs
#'
#' @param interestGeneFolder Run WebGestaltR for gene list files in the folder.
#' @param isParallel If jobs are run parallelly in the batch.
#'
#' @return The WebGestaltRBatch function returns a list of enrichment results.
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#'
#' @export
#' @aliases WebGestaltR_batch
#' @rdname WebGestaltR
#'
WebGestaltRBatch <- function(interestGeneFolder=NULL, enrichMethod="ORA", isParallel=FALSE, nThreads=3, ...) {
	args <- list(...)
	if(enrichMethod=="ORA" || enrichMethod=="NTA"){
		interestGeneFiles <- list.files(interestGeneFolder,pattern="\\.txt",full.names=TRUE)
	}

	if(enrichMethod=="GSEA"){
		interestGeneFiles <- list.files(interestGeneFolder,pattern="\\.rnk",full.names=TRUE)
	}

	projectNames <- unlist(lapply(strsplit(basename(interestGeneFiles),split=".",fixed=TRUE),function(e){return(paste(e[-length(e)],collapse="."))}))

	resultList <- list()

	if(isParallel==TRUE){
		cl <- makeCluster(nThreads)
		registerDoParallel(cl)
		resultList <- foreach(i=1:length(interestGeneFiles), .packages="WebGestaltR") %dopar% {
			args$interestGeneFile <- interestGeneFiles[i]
			args$projectName <- projectNames[i]
			args$enrichMethod <- enrichMethod
			sig <- do.call(WebGestaltR, args)
			re <- list(filename=interestGeneFiles[i], enrichResult=sig)
			return(re)
		}
		stopCluster(cl)
	}else{
		for(i in c(1:length(interestGeneFiles))){
			cat("Process file: ",interestGeneFiles[i],"\n",sep="")
			args$interestGeneFile <- interestGeneFiles[i]
			args$projectName <- projectNames[i]
			args$enrichMethod <- enrichMethod
			sig <- do.call(WebGestaltR, args)
			re <- list(filename=interestGeneFiles[i], enrichResult=sig)
			resultList[[i]] <- re
		}
	}
	return(resultList)
}

#' @export
WebGestaltR_batch <- function(is.output=TRUE, ...) {
	warning("Function WebGestaltR_batch is deprecated and changed to WebGestaltRBatch!\n")
	return(WebGestaltRBatch(isOutput=is.output, ...))
}
