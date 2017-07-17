formatCheck <- function(dataType="list",inputGeneFile=NULL,inputGene=NULL){

		dataTypeA <- c("list","rnk")
		if(length(dataType)>1 || !is.character(dataType) || length(which(dataTypeA==dataType))==0){
			error <- "ERROR: dataType parameter can only be one of 'list' or 'rnk'."
      cat(error)
      return(error)
		}
	
		if(dataType=="list"){
			if(!is.null(inputGeneFile)){
				if(file_extension(inputGeneFile)!="txt"){
					error <- "ERROR: For the user ID list, please upload a 'txt' file with only one column."
          cat(error)
          return(error)
				}else{
					inputGene <- tryCatch(fread(input=inputGeneFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE,showProgress=FALSE),error=function(e){return("ERROR: The format of the uploaded gene list is incorrect, the file name contains the special characters or the character encoding in the file is not UTF-8. Please check the file format, file name or character encoding.")})
					
          if(.hasError(inputGene)){
						cat(inputGene)
						return(inputGene)
					}
					
					
					if(ncol(inputGene)!=1){
						error <- "ERROR: For the user ID list, please upload a 'txt' file with only one column."
          	cat(error)
          	return(error)
					}else{
						inputGene <- as.character(inputGene[,1])
						return(inputGene)
					}
				}
			}else{
				if(!is.null(inputGene)){
					if(!is.vector(inputGene)){
						error <- "ERROR: For the user ID list, please upload an R vector object."
          	cat(error)
          	return(error)
					}else{
						inputGene <- as.character(inputGene)
						return(inputGene)
					}
				}else{
					error <- "ERROR: Please upload a file or an R object for the ID list."
         	cat(error)
          return(error)
				}
			}
		}
		
		if(dataType=="rnk"){
			if(!is.null(inputGeneFile)){
				if(file_extension(inputGeneFile)!="rnk"){
					error <- "ERROR: For the ranked list, please upload a 'rnk' file with two columns (ids and scores)."
          cat(error)
          return(error)
				}else{
					inputGene <- tryCatch(fread(input=inputGeneFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE,showProgress=FALSE),error=function(e){return("ERROR: The format of the uploaded ranked list is incorrect, the file name contains the special characters or the character encoding in the file is not UTF-8. Please check the file format, file name or character encoding.")})
					
          if(.hasError(inputGene)){
						cat(inputGene)
						return(inputGene)
					}
					
					
					if(ncol(inputGene)!=2){
						error <- "ERROR: For the ranked list, please upload a 'rnk' file with two columns (ids and scores)."
          	cat(error)
          	return(error)
					}else{
						if(!is.numeric(inputGene[,2]) && !is.integer(inputGene[,2])){
							error <- "ERROR: The second column of the ranked list should be the numeric scores."
          		cat(error)
          		return(error)
						}else{
							#########GSEA do not allow the second column contains NA. Thus, we should remove NA first############
							inputGene <- inputGene[!is.na(inputGene[,2]),]
							inputGene[,1] <- as.character(inputGene[,1])
							return(inputGene)
						}
					}
				}
			}else{
				if(!is.null(inputGene)){
					if(!is.data.frame(inputGene)){
						error <- "ERROR: For the ranked list, please upload an R data.frame object."
          	cat(error)
          	return(error)
					}else{
						if(!is.numeric(inputGene[,2]) && !is.integer(inputGene[,2])){
							error <- "ERROR: The second column of the ranked list should be the numeric scores."
          		cat(error)
          		return(error)
						}else{
							#########GSEA do not allow the second column contains NA. Thus, we should remove NA first############
							inputGene <- inputGene[!is.na(inputGene[,2]),]
							
							inputGene[,1] <- as.character(inputGene[,1])
							return(inputGene)
						}
					}
				}else{
					error <- "ERROR: Please upload a file or an R object for the ranked list."
         	cat(error)
          return(error)
				}
			}
		}
}
