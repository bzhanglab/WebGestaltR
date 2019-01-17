#' @importFrom stats sd IQR
mergeDuplicate <- function(id,data,collapseMode="maxSD"){
	data <- as.data.frame(data, stringsAsFactors=FALSE)
	if(!is.vector(id)){
		stop("The input id should be a vector object!\n")
	}

	if(!is.matrix(data) && !is.data.frame(data)){
		stop("The input data matrix should be a matrix or data.frame object!\n")
	}

	if(length(which(collapseMode %in% c("mean","median","maxSD","maxIQR","max","min")))==0){
		stop("The input 'collapseMode' is not valide! Please select an option from 'mean', 'median', 'maxSD',  'maxIQR', 'max' and 'min' (use mean, median, max standard derivation, max interquartile range, maximum and minimum of duplicate genes for each sample)!\n")
	}

	if(ncol(data)==1){
		if(collapseMode=="maxSD" || collapseMode=="maxIQR"){
			collapseMode <- "mean"
		}
	}

	mergeR <- c()

	if(ncol(data)==1){
		data <- data.frame(id=id,data=data[,1],stringsAsFactors=F)
		data <- tapply(data[,2],data[,1],collapseMode,na.rm=TRUE)
		return(data)

	}else{
		if(collapseMode=="mean" || collapseMode=="median" || collapseMode=="max" || collapseMode=="min"){
			mergeR <- lapply(split(c(1:nrow(data)),id),function(u){return(apply(data[u,],2,collapseMode,na.rm=TRUE))})
			data <- do.call(rbind,mergeR)
			return(data)
		}

		if(collapseMode=="maxSD" || collapseMode=="maxIQR"){
			if(collapseMode=="maxSD"){
				mergeR <- apply(data,1,sd,na.rm=TRUE)
			}

			if(collapseMode=="maxIQR"){
				mergeR <- apply(data,1,IQR,na.rm=TRUE)
			}

			mergeR <- data.frame(id=id,score=mergeR,pos=c(1:length(mergeR)),stringsAsFactors=F)
			mergeR <- mergeR[order(mergeR[,1],-mergeR[,2]),]
			finalPos <- tapply(mergeR[,3],mergeR[,1],.returnMaxPos)
			data <- data[finalPos,]
			id <- id[finalPos]

			rownames(data) <- id
			return(data)
		}
	}
}

.returnMaxPos <- function(pos){
	return(pos[1])
}
