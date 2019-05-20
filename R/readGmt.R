#' Read GMT File
#'
#' @inheritParams WebGestaltR
#' @param gmtFile The file path or URL of the GMT file.
#'
#' @return A data frame with columns of "geneSet", "description", "gene".
#'
#' @importFrom httr content
#' @importFrom tools file_ext
#' @export
#'
readGmt <- function(gmtFile, cache=NULL) {
#####Change a gmt file to a three column matrix (gene set name, gene set description and genes)#######
	if (startsWith(gmtFile, "http://") || startsWith(gmtFile, "https://")) {
		response <- cacheUrl(gmtFile, cache)
		if (response$status_code == 200) {
			data <- unlist(strsplit(content(response), "\n", fixed=TRUE))
		} else {
			stop(webRequestError(response))
		}
	} else {
		if (file_ext(gmtFile) != "gmt") {
			stop(gmtFormatError("empty"))
		}
		# remove BOM in some windows files
		data <- gsub("\xEF\xBB\xBF", "", readLines(gmtFile, skipNul=TRUE), useBytes=TRUE)
	}
	data <- strsplit(data, "\t", useBytes=TRUE)
	data <- lapply(data,.toList)
	data <- do.call("rbind",data)

	if (is.null(data)) {
		stop(gmtFormatError("incorrect"))
	} else {
		data <- as.data.frame(data, stringsAsFactors=FALSE)
		colnames(data) <- c("geneSet", "description", "gene")
		return(data)
	}
}
readGMT <- readGmt


.toList <- function(data) {
	if (length(data)>2) {
		data <- data[!is.na(data) && !is.null(data)]
		# replace % in gene set names to _, because png treats % in filename specially
		data1 <- cbind(rep(gsub('%', '_', data[1], fixed=TRUE), length(data)-2), rep(data[2], length(data)-2), data[c(-1,-2)])
		return(data1)
	} else {
		return(NULL)
	}
}

#' Prepare Input Matrix for GSEA
#'
#' @param rank A 2 column Data Frame of gene and score
#' @param gmt 3 column Data Frame of geneSet, description, and gene
#'
#' @return A matrix used for input to \code{swGsea}.
#'
#' @importFrom dplyr filter select distinct inner_join %>%
#' @export
#'
prepareInputMatrixGsea <- function(rank, gmt) {
	genes <- rank$gene
	gmt <- gmt %>% filter(.data$gene %in% genes)
	geneSets <- (gmt %>% select(.data$geneSet, .data$description) %>% distinct())$geneSet
	# 0 or 1 matrix indicating gene and gene set relationship
	rel <- fillInputDataFrame(gmt, genes, geneSets)
	# R implementation
	# rel <- matrix(0, nrow=length(genes), ncol=length(geneSets), dimnames=list(genes, geneSets))
	#
	# for (i in 1:nrow(gmt)) {
	# 	rel[gmt[i, "gene"], gmt[i, "geneSet"]] <- 1
	# }
	# rel <- as.data.frame(rel)
	# rel$gene <- genes
	return(inner_join(rank, rel, by="gene"))
}

readGMT <- function(...) {
	warning("Function readGMT is deprecated and changed to readGmt!\n")
	return(readGmt(...))
}
