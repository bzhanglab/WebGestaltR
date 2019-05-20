#' Wrapper for mapping ids to "genesymbol"
#'
#' @inheritParams idMapping
#'
#' @importFrom tools file_ext
#'
#' @export
#' @rdname idMapping
idToSymbol <- function(organism="hsapiens", dataType="list", inputGeneFile=NULL,
                       inputGene=NULL, sourceIdType="ensembl_gene_id", collapseMethod="mean",
                       mappingOutput=FALSE, outputFileName=NULL, cache=NULL,
                       hostName="http://www.webgestalt.org/") {
	# various error checking
	errorTest <- .hostNameError(hostName)
	if(!is.null(errorTest)){
		return(errorTest)
	}
	errorTest <- .organismError(organism, hostName)
	if(!is.null(errorTest)){
		return(errorTest)
	}
	if(!(dataType %in% c("list", "rnk", "matrix"))){
		error <- paste("ERROR: Data type ",dataType," is not supported by idToSymbol.
						Please select from 'list', 'rnk' and 'matrix'.",sep="")
		cat(error)
		return(error)
	}
	errorTest <- .collapseMethodError(collapseMethod)
	if(!is.null(errorTest)){
		return(errorTest)
	}
	errorTest <- idTypeError(idType=sourceIdType, organism=organism, hostName=hostName, cache=cache)
	if(!is.null(errorTest)){
		return(errorTest)
	}
	# actual mapping
	if(dataType == "list" || dataType == "rnk"){
		if(mappingOutput && is.null(outputFileName)){
			outputFileName <- paste0("wgr_", dataType)
		}
		geneMap <- idMappingGene(organism=organism, dataType=dataType, inputGeneFile=inputGeneFile, inputGene=inputGene, sourceIdType=sourceIdType, targetIdType="genesymbol",  collapseMethod=collapseMethod, mappingOutput=mappingOutput, outputFileName=outputFileName, hostName=hostName)
		return(geneMap)
	} else if(dataType == "matrix"){
		if(mappingOutput && is.null(outputFileName)){
			outputFileName <- paste0("wgr_", dataType, "_converted.txt")
		}
		# cct or cbt
		if(!is.null(inputGeneFile)){
			inputMat <- .testMatrixFormat(inputGeneFile, collapseMethod)
		} else if(!is.null(inputGene)){
			inputMat <- .testMatrixFormat(inputGene, collapseMethod)
		}
		inputId <- rownames(inputMat)
		geneMap <- idMappingGene(organism=organism, dataType="list", inputGeneFile=NULL, inputGene=inputId, sourceIdType=sourceIdType, targetIdType="genesymbol",  mappingOutput=FALSE, hostName=hostName)
		idMap <- geneMap$mapped[,c(1,2)]
		id <- as.vector(idMap[,2])
		inputMat <- inputMat[idMap[,1],]
		inputMat <- mergeDuplicate(id,inputMat,collapseMethod)
		re <- list(data=inputMat,idMap=idMap)
		# write to output file
		if(mappingOutput){
			mtrx <- re$data
			mtrx <- cbind(rownames(mtrx),mtrx)
			colnames(mtrx)[1] <- "GeneSymbol"
			write.table(mtrx,file=outputFileName,row.names=F,col.names=T,sep="\t",quote=F)
		}
		return(re)
	}
}

.testMatrixFormat <- function(inputMat, collapseMethod="maxSD"){
	if(class(inputMat)=="character"){
		if(file_ext(inputMat)!="cct" && file_ext(inputMat)!="cbt"){
			stop("The extension of the input file should be 'cct' or 'cbt'. \n")
		}else{
			inputMat <- read.table(inputMat,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
			geneId <- inputMat[,1]
			if(length(geneId)!=length(unique(geneId))){
				cat("The input data contain duplicate Id. The function will use ",collapseMethod," to collapse duplicate Id in each sample!\n",sep="")
				inputMat <- mergeDuplicate(geneId,inputMat[,2:ncol(inputMat)] ,collapseMethod)
			}else{
				inputMat <- inputMat[,c(2:ncol(inputMat))]
				rownames(inputMat) <- geneId
			}
		}
	}else{
		if(class(inputMat) != "matrix" && class(inputMat) != "data.frame"){
			stop("The type of input data should be a matrix or data.frame object. Other types of data are not supported by this package.!\n")
		}else{
			x <- apply(inputMat,2,function(e) return(class(e)=="numeric" || class(e)=="integer"))
			y <- all(x==TRUE)
			if(y==FALSE){
				stop("The input matrix or data.frame object should only contain numeric or integer values.\n")
			}
		}
	}
	if(ncol(inputMat)<6){
		stop("The data should contain at least six samples!\n")
	}
	if(length(which(inputMat %in% c("Inf","-Inf")))>0){
		stop("The input data contain Inf which may be gernated by some wrong operation, such as log(0) or 1/0. Please re-process the data and remove the Inf\n")
	}
	return(inputMat)
}
