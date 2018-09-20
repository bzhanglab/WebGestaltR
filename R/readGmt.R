readGmt <- function(gmtFile){
#####Change a gmt file to a three column matrix (gene set name, gene set description and genes)#######
	if (startsWith(gmtFile, "http")) {
		response <- GET(gmtFile)
		if (response$status_code == 200) {
			data <- unlist(strsplit(content(response), "\n", fixed=TRUE))
		} else {
			return(webRequestError(response))
		}
	} else {
		if(file_extension(gmtFile) != "gmt"){
			return(gmtFormatError("empty"))
		}
		data <- readLines(gmtFile)
	}
	data <- strsplit(data,"\t")
	data <- lapply(data,.toList)
	data <- do.call("rbind",data)

	if(is.null(data)){
		return(gmtFormatError("incorrect"))
	}else{
		data <- as.data.frame(data, stringsAsFactors=FALSE)
		colnames(data) <- c("geneSet", "description", "gene")
		return(data)
	}
}
readGMT <- readGmt


.toList <- function(data){
	if(length(data)>2){
		data <- data[!is.na(data)]
		data1 <- cbind(rep(data[1],length(data)-2),rep(data[2],length(data)-2),data[c(-1,-2)])
		return(data1)
	}else{
		return(NULL)
	}
}

prepareInputMatrixGsea <- function(rank, gmt) {
	# rank is 2 column Data Frame of gene and score
	# gmt is 3 column Data Frame of geneSet, geneSetLink, and gene
	geneSets <- (gmt %>% select(geneSet, description) %>% distinct())$geneSet
	genes <- rank$gene
	# 0 or 1 matrix indicating gene and gene set relationship
	relDf <- as.data.frame(matrix(0, nrow=length(genes), ncol=length(geneSets), dimnames=list(genes, geneSets)))
	for (i in 1:nrow(gmt)) {
		if (gmt[i, "gene"] %in% genes) {
			relDf[gmt[i, "gene"], gmt[i, "geneSet"]] <- 1
		}
	}
	relDf$gene <- genes
	return(inner_join(rank, relDf, by="gene"))
}

readGMT <- function(...) {
	cat("WARNING: Function readGMT is deprecated and changed to readGmt!\n")
	return(readGmt(...))
}
