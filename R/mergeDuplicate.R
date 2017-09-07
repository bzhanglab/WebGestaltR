mergeDuplicate <- function(id,data,collapse_mode="maxSD"){
    
    data <- as.data.frame(data)
    if(!is.vector(id)){
        stop("The input id should be a vector object!\n")
    }
    
    if(!is.matrix(data) && !is.data.frame(data)){
        stop("The input data matrix should be a matrix or data.frame object!\n")
    }
    
    
    if(length(which(collapse_mode %in% c("mean","median","maxSD","maxIQR","max","min")))==0){
        stop("The input 'collapse_mode' is not valide! Please select an option from 'mean', 'median', 'maxSD',  'maxIQR', 'max' and 'min' (use mean, median, max standard derivation, max interquartile range, maximum and minimum of duplicate genes for each sample)!\n")
    }
    
    
    if(ncol(data)==1){
    	if(collapse_mode=="maxSD" || collapse_mode=="maxIQR"){
    		collapse_mode <- "mean"
    	}
    }
    
    mergeR <- c()
    
    if(ncol(data)==1){
    	data <- data.frame(id=id,data=data[,1],stringsAsFactors=F)
    	data <- tapply(data[,2],data[,1],collapse_mode,na.rm=TRUE)
    	return(data)
    	
   	}else{
    
    	if(collapse_mode=="mean" || collapse_mode=="median" || collapse_mode=="max" || collapse_mode=="min"){
    	
    		mergeR <- lapply(split(c(1:nrow(data)),id),function(u){return(apply(data[u,],2,collapse_mode,na.rm=TRUE))})
        
        data <- do.call(rbind,mergeR)
        return(data)
    	}
    
    
   	 if(collapse_mode=="maxSD" || collapse_mode=="maxIQR"){
        if(collapse_mode=="maxSD"){
            mergeR <- apply(data,1,sd,na.rm=TRUE)
        }
        
        if(collapse_mode=="maxIQR"){
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