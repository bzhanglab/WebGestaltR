parameterErrorMessage <- function(hostName="http://www.webgestalt.org/", cache=NULL, ...){
	errorTest <- NULL
	argList <- list(...)

	argNames <- names(argList)
	# special cases, assuming minNum, maxNum appear together
	if("minNum" %in% argNames && "maxNum" %in% argNames){
		errorTest <- .minMaxNumError(minNum=argList$minNum,maxNum=argList$maxNum)
		if(!is.null(errorTest)){
			return(errorTest)
		}
		argList$minNum <- NULL
		argList$maxNum <- NULL
	}

	if("organism" %in% argNames){
		errorTest <- .organismError(organism=argList$organism, hostName=hostName, cache=cache)
		if(!is.null(errorTest)){
			return(errorTest)
		}
		argList$organism <- NULL
	}

	# individual checking function has name with format:
	# "." + parameterName + "Error"
	for(i in seq_along(argList)){
		errorTest <- do.call(paste0(".", names(argList[i]),"Error"), argList[i])
		if(!is.null(errorTest)){
			return(errorTest)
		}
	}
	return(errorTest)
}

.enrichMethodError <- function(enrichMethod){ #####Input method error
	existingMethods <- c("ORA", "GSEA", "NTA")

	if(!(enrichMethod %in% existingMethods)){
		error <- paste0("ERROR: The enrichment method '",enrichMethod, "' is not supported.")
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.organismError <- function(organism, hostName, cache) { ####Input organism error
	organisms <- listOrganism(hostName=hostName, cache=cache)
	organisms <- c(organisms, "others")
	if(!(organism %in% organisms)){
		error <- paste0("ERROR: The organism '",organism,"' is not supported.")
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.collapseMethodError <- function(collapseMethod){ ###Input collapseMethod error
	collapseMethodList <- c("mean","median","min","max")
	if(!(collapseMethod %in% collapseMethodList)){
		error <- "ERROR: WebGesalt only supports the following collapse methods: mean, median, min and max."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.minMaxNumError <- function(minNum,maxNum){  #####Input the minimum and maximum number of category error
	if(!is.wholenumber(minNum) || minNum<1){
		error <- "ERROR: The minimum number of genes annotated to the category should be an positive integer."
		cat(error)
		return(error)
	}else{
		if(!is.wholenumber(maxNum) || maxNum<1){
			error <- "ERROR: The maximum number of genes annotated to the category should be an positive integer."
			cat(error)
			return(error)
		}else{
			if(minNum>=maxNum){
				error <- "ERROR: The minimum number of genes annotated to the category should be less than the maximum number."
				cat(error)
				return(error)
			}else{
				return(NULL)
			}
		}
	}
}

.fdrMethodError <- function(fdrMethod){   ##Input fdrMethod error
	fdrMethodList <- c("holm","hochberg","hommel","bonferroni","BH","BY")

	if(!(fdrMethod %in% fdrMethodList)){
		error <- "ERROR: WebGesalt only supports the following FDR methods: holm, hochberg, hommel, bonferroni, BH and BY."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.sigMethodError <- function(sigMethod){   ##Input sig method error
	sigMethodList <- c("fdr","top")
	if(!(sigMethod %in% sigMethodList)){
		error <- "ERROR: WebGesalt only supports two methods to identify the enriched categories: 'fdr' method identifies the categories with FDR less than 'fdrThr' and 'top' method ranks all categories based on the FDR and selects the 'topThr' categories."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.fdrThrError <- function(fdrThr){   ##Input fdr threshold error
	if(!is.numeric(fdrThr) || fdrThr<0 || fdrThr>1){
		error <- "ERROR: The FDR threshold should be a numberic from 0 to 1."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.topThrError <- function(topThr){   ##Input top threshold error
	if(!is.wholenumber(topThr) || topThr<1){
		error <- "ERROR: The number of top categories should be a positive integer."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.reportNumError <- function(reportNum){   ##Input the number of visualized enriched category error
	if (!is.wholenumber(reportNum) || reportNum<0) {
		error <- "ERROR: The number of enriched categories shown in the final report should be a positive integer."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.perNumError <- function(perNum){   ##Input the number of permutation error
	if(!is.wholenumber(perNum) || perNum<0 || perNum>10000){
		error <- "ERROR: The number of permutation for GSEA method should be a positive integer and less than 10000."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.isOutputError <- function(isOutput){  ##Input isOutput error
	if(!is.logical(isOutput)){
		error <- "ERROR: isOutput should be an R logical object (TRUE or FALSE)."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.outputDirectoryError <- function(outputDirectory){   ##Input directory error
	if(!dir.exists(outputDirectory)){
		error <- paste("ERROR: The output directory ",outputDirectory," does not exist. please change another directory or create the directory.",sep="")
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.dagColorError <- function(dagColor){   ##Input dag color type error
	dagColorList <- c("binary","continuous")
	if(!(dagColor %in% dagColorList)){
		error <- "ERROR: dagColor should be binary or continuous."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.hostNameError <- function(hostName){ ####Input hostName error
		if(grepl("^file://", hostName, perl=TRUE)){
			return(NULL)
		}
		archives <- listArchiveUrl()
		if(!(hostName %in% archives[,2])){
			error <- paste("ERROR: The host name ",hostName," is incorrect. Please use listArchiveURL function to find the correct host name.",sep="")
			cat(error)
			return(error)
		}else{
			return(NULL)
		}
}

.dataTypeError <- function(dataType){  ######Input data type error for ID Mapping function
	dataTypeA <- c("list","rnk","gmt")
	if(!(dataType %in% dataTypeA)){
		error <- paste("ERROR: Data type ",dataType," can not be supported. Please select from 'list', 'rnk' and 'gmt'.",sep="")
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

.mappingOutputError <- function(mappingOutput){  ##Input mappingOutput error for ID Mapping function
	if(!is.logical(mappingOutput)){
		error <- "ERROR: mappingOutput should be an R logical object (TRUE or FALSE)."
		cat(error)
		return(error)
	}else{
		return(NULL)
	}
}

##############Other functions for the parameter error test############
testNull <- function(parameter){
	if(!is.null(parameter)){
		if(length(parameter)==1 && length(which(parameter=="NULL"))==1){
			parameter <- NULL
		}
	}
	return(parameter)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
	re <- tryCatch(abs(x - round(x)) < tol,warning=function(e){return(FALSE)},error=function(e){return(FALSE)})
	return(re)
}


##################Other errors###############

descriptionFileError <- function(type){
	if(type=="format"){
		error <- "ERROR: The description file for the functional categetories should have a 'des' extension."
	}

	if(type=="columnNum"){
		error <- "ERROR: The description file should contain two columns: the first column is the ID of the gene sets that should be the same with the uploaded enrichment gene sets, the second column is the description of the gene sets."
	}

	if(type=="overlap"){
		error <- "ERROR: The ID types of the uploaded functional database file and description file are different. Please check the uploaded files."
	}

	cat(error)
	return(error)
}

gmtFormatError <- function(type){
	if(type=="empty"){
		error <- "ERROR: Please upload a file with extension 'gmt'."
	}
	if(type=="incorrect"){
		error <- "ERROR: Invalid GMT file format. Please check the format of the GMT file from http://www.webgestalt.org/WebGestalt_2017_Manual.pdf."
	}
	cat(error)
	return(error)
}


enrichDatabaseError <- function(type,enrichDatabase="",organism=""){
	if(type=="unsupported"){
		error <- paste("ERROR: ",enrichDatabase," can not be supported for organism ",organism,".",sep="")
	}

	if(type=="empty"){
		error <- "ERROR: Please select the functional database or select 'others' to upload the functional database."
	}

	if(type=="others"){
		error <- "ERROR: Please upload a 'gmt' file as the functional database for 'others' organism."
	}

	cat(error)
	return(error)
}


idTypeError <- function(idType, organism, hostName, cache) {
		idTypes <- listIdType(organism=organism, hostName=hostName, cache=cache)
		if(!(idType %in% idTypes)){
				error <- paste("ERROR: The ID type ",idType," can not be supported for organism ",organism,".",sep="")
				cat(error)
				return(error)
		}else{
				return(NULL)
		}
}

targetIdTypeError <- function(idType, organism, hostName, cache) {
		idTypes <- listIdType(organism=organism, hostName=hostName, cache=cache)
		if(!(idType %in% idTypes)){
				error <- paste("ERROR: The target ID type ",idType," can not be supported for organism ",organism,".",sep="")
				cat(error)
				return(error)
		}else{
				return(NULL)
		}
}

stardardDiffError <- function(standardSource,standardTarget){
	if(standardSource!=standardTarget){
		error <- paste(standardSource," based ID type can not map to ",standardTarget," based ID type.",sep="")
		cat(error)
		return(error)
	}else{
		return(NULL)
	}

}

idMappingError <- function(type,idType="",topF=""){
	if(type=="unmapped"){
		error <- paste("ERROR: The ID type of the uploaded list is not consistent with the input ID type ",idType,". Examples of the input ID type: ",topF,".",sep="")
	}
	if(type=="empty"){
		error <- "ERROR: No IDs are mapped. Please check your input."
	}

	cat(error)
	return(error)
}

interestGeneError <- function(type){
	if(type=="empty"){
		error <- "ERROR: Please upload an interesting ID list."
	}
	if(type=="emptyType"){
		error <- "ERROR: Please select the ID type of the interesting ID list."
	}
	if(type=="unmatch"){
		error <- "The ID type of the input ID list can not match to the ID type of the functional database."
	}

	if(type=="unannotated"){
		error <- "ERROR: All IDs in the uploaded list are not annotated to any category of the functional database."
	}

	if(type=="onlyOne"){
		error <- "ERROR: All IDs in the uploaded list can only annotate to one category, which can not perform the enrichment analysis."
	}

	cat(error)
	return(error)
}

referenceGeneError <- function(type){
	if(type=="empty"){
		error <- "ERROR: Please upload a reference list or select an existing reference set."
	}
	if(type=="emptyType"){
		error <- "ERROR: Please select the ID type of the reference list."
	}
	if(type=="unmatch"){
		error <- "The ID type of the input reference list can not match to the ID type of the functional database."
	}
	if(type=="unmatch_set"){
		error <- "The ID type of the selected reference set can not match to the ID type of the functional database."
	}
	if(type=="existingRef"){
		error <- "ERROR: Please select an existing reference set."
	}
	if(type=="interestEmpty"){
		error <- "ERROR: All IDs in the interesting list are not included in the overlapping IDs between IDs in the reference list and IDs annotated to the functional database."
	}

	cat(error)
	return(error)
}

webRequestError <- function(response) {
	error <- paste(response["status_code"], content(response)["message"], "at", response["url"], "\n")
	cat(error)
	return(error)
}

webApiError <- function(response) {
	error <- response[["message"]]
	cat(paste0(error, "\n"))
	return(error)

}

.hasError <- function(obj){
	if(is.character(obj) && length(obj)==1 && length(grep("ERROR:",obj))>0) {
		return(TRUE)
	}
	return(FALSE)
}
