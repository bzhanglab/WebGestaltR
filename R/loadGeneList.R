loadInterestGene <- function(organism="hsapiens",dataType="list",inputGeneFile=NULL,inputGene=NULL,geneType="entrezgene",collapseMethod="mean",methodType="R",hostName="http://www.webgestalt.org/",geneSet){
		
    if(is.null(inputGeneFile) && is.null(inputGene)){
    	return(interestGeneError(type="empty"))
    }else{
    	if(organism!="others"){
    		if(is.null(geneType)){
    			return(interestGeneError(type="emptyType"))
    		}else{
    			mapRe <- .uploadGene_existingOrganism(organism=organism,dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene,geneType=geneType,collapseMethod=collapseMethod,geneSet=geneSet,methodType=methodType,hostName=hostName)
    			if(is.character(mapRe) && length(mapRe)==1 && length(grep("ERROR:",mapRe))>0){
    				return(mapRe)
    			}
    		}
    	}else{
    		mapRe <- .uploadGene_Others(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene,geneSet=geneSet)
    		if(is.character(mapRe) && length(mapRe)==1 && length(grep("ERROR:",mapRe))>0){
    			return(mapRe)
    		}
    	}
    }
    
    #if organism is not others, the function will return a mapping result with mapped and unmapped list
    #if organism is others, the function will return a matrix with gene list
    return(mapRe)
}

loadReferneceGene <- function(organism="hsapiens",referenceGeneFile=NULL,referenceGene=NULL,referenceGeneType="entrezgene",referenceSet=NULL,collapseMethod="mean",methodType="R",hostName="http://www.webgestalt.org/",geneSet,interestGene_List){
    
    referenceGene_List <- NULL
    referenceGeneMap <- NULL
    	 
   	if(is.null(referenceGeneFile) && is.null(referenceGene) && is.null(referenceSet)){
    	 	return(referenceGeneError(type="empty"))
    }else{
    	if(organism!="others"){
  			if(!is.null(referenceGeneFile) || !is.null(referenceGene)){
  				if(is.null(referenceGeneType)){
  					return(referenceGeneError(type="emptyType"))
  				}else{
    				mapRe <- .uploadGene_existingOrganism(organism=organism,dataType="list",inputGeneFile=referenceGeneFile,inputGene=referenceGene,geneType=referenceGeneType,collapseMethod=collapseMethod,geneSet=geneSet,methodType=methodType,hostName=hostName)
    				if(is.character(mapRe) && length(mapRe)==1 && length(grep("ERROR:",mapRe))>0){
    					return(mapRe)
    				}
    				gene_standardId <- identifyStandardId(hostName=hostName,idtype=referenceGeneType,organism=organism,type="interest")
    				referenceGene_List <- mapRe$mapped[,gene_standardId]
    			}			
    		}else{ ###referenceGeneFile and referenceGene are both NULL. But referenceSet is not NULL
    			refS <- listReferenceSet(organism=organism,hostName=hostName)
    			if(length(which(refS==referenceSet))==0){
    				return(referenceGeneError(type="existingRef"))
    			}
    			
    			ref_standardId <- identifyStandardId(hostName=hostName,idtype=referenceSet,organism=organism,type="reference")
    			referenceGene_List <- fread(input=file.path(hostName,"data","reference",paste(organism,"_",referenceSet,"_",ref_standardId,".table",sep="")),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE)
      		referenceGene_List <- as.character(unique(referenceGene_List[,1]))
      	}
    	}else{ ##For other organisms
    		if(!is.null(referenceGeneFile) || !is.null(referenceGene)){    	
    				referenceGene_List <- .uploadGene_Others(dataType="list",inputGeneFile=referenceGeneFile,inputGene=referenceGene,geneSet=geneSet)
    				if(is.character(referenceGene_List) && length(referenceGene_List)==1 && length(grep("ERROR:",referenceGene_List))>0){
    					return(referenceGene_List)
    				}
    		}else{
    				return(referenceGeneError(type="empty"))
    		}
    	}
    }
    	 
    ##compare interest gene list and reference gene list
    if(length(intersect(interestGene_List,intersect(referenceGene_List,geneSet[,3])))==0){
    		return(referenceGeneError(type="interestEmpty"))
    }
    return(referenceGene_List)
}


.uploadGene_existingOrganism <- function(organism,dataType,inputGeneFile,inputGene,geneType,collapseMethod,geneSet,methodType,hostName){
					
					geneMap <- IDMapping(organism=organism,dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene,sourceIdType=geneType,targetIdType=NULL,collapseMethod=collapseMethod,mappingOutput=FALSE,methodType=methodType,hostName=hostName)
					
					if(is.character(geneMap) && length(geneMap)==1 && length(grep("ERROR:",geneMap))>0){
						return(geneMap)
					}
					
					#gene_standardId <- identifyStandardId(hostName=hostName,idtype=geneType,organism=organism,type="interest")  ##identifyStandardId in IDMapping_component.R 
    			#if(gene_standardId!=databaseStandardId){  ###the standardId of the input genes should be the same with the standardarId of the functional database
    			#	return(interestGeneError(type="unmatch"))
    			#}
					
					geneMap_MappedList <- geneMap$mapped
					standardId <- geneMap$standardId
					
					gene_List <- as.character(unique(geneMap_MappedList[,standardId]))
    			ov <- intersect(gene_List,geneSet[,3])
					
					if(length(ov)==0){
    				return(interestGeneError(type="unannotated"))
    			}
					
					###Because if all genes are annotated to only one category, GSEA will return the error, we need to avoid this error by reporting the error in the R#
    			gL <- geneSet[geneSet[,3] %in% gene_List,,drop=FALSE]
    			gL <- tapply(gL[,3],gL[,1],length)
    			if(length(gL)==1){
    				return(interestGeneError(type="onlyOne"))
    			}
    			return(geneMap)
}


.uploadGene_Others <- function(dataType,inputGeneFile,inputGene,geneSet){
		
		gene_List <- formatCheck(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
		
		if(is.character(gene_List) && length(gene_List)==1 && length(grep("ERROR:",gene_List))>0){
			return(gene_List)
		}
		gene_List <- as.data.frame(gene_List,stringsAsFactors=F)
		ov <- intersect(gene_List[,1],geneSet[,3])
    if(length(ov)==0){
    		return(interestGeneError(type="unannotated"))
    }
    
    ###Because if all genes are annotated to only one category, GSEA will return the error, we need to avoid this error by reporting the error in the R#
    gL <- geneSet[geneSet[,3] %in% gene_List[,1],,drop=FALSE]
    gL <- tapply(gL[,3],gL[,1],length)
    if(length(gL)==1){
    		return(interestGeneError(type="onlyOne"))
    }
    ####gene_list is a matrix with one or two columns#####
    return(gene_List)
}	

