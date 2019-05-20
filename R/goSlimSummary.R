#' GO Slim Summary
#'
#' Outputs a brief summary of input genes based on GO Slim data.
#'
#' @inheritParams WebGestaltR
#' @param geneList A list of input genes.
#' @param outputFile Output file name.
#' @param isOutput Boolean if a plot is save to \code{outputFile}.
#' @param outputType File format of the plot: \code{pdf}, \code{bmp} or \code{png}.
#'
#' @return A list of the summary result.
#' @export
#' @aliases GOSlimSummary
#'
goSlimSummary <- function(organism="hsapiens", geneList, outputFile, outputType="pdf", isOutput=TRUE, cache=NULL, hostName="http://www.webgestalt.org") {
	organisms <- listOrganism(hostName=hostName, cache=cache)
	if (!organism %in% organisms) {
		stop("ERROR: ", organism, " can not be supported.")
	}

	if (!is.vector(geneList)) {
		stop("ERROR: Please upload a list of entrez gene ids.")
	} else {
		geneList <- as.character(geneList)
	}

	outputTypeList <- c("pdf","png","bmp")
	if (!outputType %in% outputTypeList) {
	stop("ERROR: The output Type ", outputType ," is invalid. Please select one from pdf, png, and bmp.")
	}

	bpRes <- .processData(organism, hostName, geneList, "Biological_Process")
	bpGoCnts <- bpRes$goTermCounts
	bpUnclassified <- bpRes$dataUnclassified
	if (isOutput) {
		write.table(bpGoCnts, paste0(outputFile, "_bp.txt"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	}

	ccRes <- .processData(organism, hostName, geneList, "Cellular_Component")
	ccGoCnts <- ccRes$goTermCounts
	ccUnclassified <- ccRes$dataUnclassified
	if (isOutput) {
		write.table(ccGoCnts, paste0(outputFile, "_cc.txt"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	}

	mfRes <- .processData(organism, hostName, geneList, "Molecular_Function")
	mfGoCnts <- mfRes$goTermCounts
	mfUnclassified <- mfRes$dataUnclassified
	if (isOutput) {
		write.table(mfGoCnts, paste0(outputFile, "_mf.txt"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	}

	if (outputType=="pdf") {
		pdf(paste0(outputFile, ".pdf"), height=8, width=16)
	} else if (outputType=="png") {
		png(paste0(outputFile, ".png"), width=2800, height=1300, res=200)
	} else if (outputType=="bmp") {
		bmp(paste0(outputFile, ".bmp"), width=2800, height=1300, res=200)
	}

	layout(matrix(1:3,1,3))

	.plotData(geneList, bpGoCnts, bpUnclassified, "Biological Process", "red")
	.plotData(geneList, ccGoCnts, ccUnclassified, "Cellular Component", "blue")
	.plotData(geneList, mfGoCnts, mfUnclassified, "Molecular Function", "green")

	dev.off()
	return(NULL)
}

#' @importFrom httr POST content
#' @importFrom dplyr select distinct inner_join arrange filter %>%
#' @importFrom readr read_tsv
.processData <- function(organism, hostName, geneList, ontology) {
	names <- c("entrezgene", "accession", "name")
	if (startsWith(hostName, "file://")) {
		goPath <- removeFileProtocol(file.path(hostName, "goslim", paste0(organism, "_GOSlim_", ontology, ".table")))
		goSlimData <- read_tsv(goPath, col_names=names, col_types="ccc", quote="")
		goSlimData <- filter(goSlimData, .data$entrezgene %in% geneList)
	} else {
		goUrl <- file.path(hostName, "api", "goslim")
		response <- POST(goUrl, body=list(organism=organism, ontology=ontology, entrezgenes=geneList), encode="json")
		if (response$status_code != 200) {
			stop(webRequestError(response))
		}
		goSlimData <- content(response)[["goslim"]]
		goSlimData <- data.frame(matrix(unlist(lapply(goSlimData, FUN=function(x) { x[names] })), nrow=length(goSlimData), byrow=TRUE), stringsAsFactors=FALSE)
		colnames(goSlimData) <- names
	}


	if (nrow(goSlimData) == 0) {
		stop("ERROR: The ID type of the uploaded list is not entrez gene id. Please use IDMapping function to map other id types to entrez gene id.")
	}

	dataUnclassified <- setdiff(geneList, unique(goSlimData[["entrezgene"]]))

	goTermCount <- tapply(goSlimData[["entrezgene"]], goSlimData[["accession"]], length)
	goTermCount <- data.frame(accession=names(goTermCount), geneNum=unname(goTermCount), stringsAsFactors=FALSE)
	uniqueGoCounts <- goSlimData %>% select(.data$accession, .data$name) %>% distinct() %>%
		inner_join(goTermCount, by="accession") %>%
		arrange(desc(.data$geneNum))

	re <- list(goTermCounts=uniqueGoCounts,dataUnclassified=dataUnclassified)
	return(re)
}

.plotData <- function(geneList,goCounts,dataUnclassified,ontology,color){
	par(mar=c(20,5,2,2))
	c <- c(length(geneList), goCounts$geneNum, length(dataUnclassified))
	names(c) <- c("all", goCounts$name, "unclassified")
	maxC <- max(c)
	xx <- barplot(c,main=paste("Bar chart of ",ontology," categories",sep=""),col=color,xaxt="n",xlab="",ylim=c(0,1.2*maxC),cex.axis=1.5,las=2,cex.main=1.7)
	text(x=xx+0.3,y=c+maxC*0.05,labels=c,pos=3,cex=1.5,col="black",srt=90)
	axis(side=1,at=xx,labels=names(c),las=2,cex.axis=1.5)
}


#' @export
GOSlimSummary <- function(...) {
	warning("Function GOSlimSummary is deprecated and changed to goSlimSummary!")
	return(goSlimSummary(...))
}
