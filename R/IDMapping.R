IDMapping <- function(organism="hsapiens",dataType="list",inputGeneFile=NULL,inputGene=NULL,sourceIdType,targetIdType=NULL,collapseMethod="mean",mappingOutput=FALSE, outputFileName="",methodType="R",hostName="http://www.webgestalt.org/"){
    
     #############Check general parameters########
    errorTest <- parameterErrorMessage(organism=organism,dataType=dataType,collapseMethod=collapseMethod,methodType=methodType,hostName=hostName,mappingOutput=mappingOutput)
		
		if(!is.null(errorTest)){
			return(errorTest)
		}
    
     ############Check source id type#########
     errorTest <- IDTypeERROR(idType=sourceIdType,organism=organism,hostName=hostName)
     if(!is.null(errorTest)){
     	return(errorTest)
     }
     
    ##########Identify the standardID for the input ID type###########
     standardSource <- identifyStandardId(hostName=hostName,idtype=sourceIdType,organism=organism,type="interest")
     
     ############Check target id type######### 
     if(!is.null(targetIdType)){
     		errorTest <- targetIDTypeERROR(idType=targetIdType,organism=organism,hostName=hostName)
     		if(!is.null(errorTest)){
     			return(errorTest)
     		}else{
     			standardTarget <- identifyStandardId(hostName=hostName,idtype=targetIdType,organism=organism,type="interest")
     			errorTest <- stardardDiffError(standardSource=standardSource,standardTarget=standardTarget)
     			if(!is.null(errorTest)){
     					return(errorTest)
     			}
     		}
     }else{
     	targetIdType <- standardSource
     }
     
     ##########gene level ID Mapping##########
     if(standardSource=="entrezgene"){
     	idmap <- IDMapping_gene(organism=organism,dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene,sourceIdType=sourceIdType,standardID=standardSource,targetIdType=targetIdType, collapseMethod=collapseMethod,mappingOutput=mappingOutput,outputFileName=outputFileName,methodType=methodType,hostName=hostName)
     }
    
     if(standardSource=="phosphositeSeq"){
     	idmap <- IDMapping_phosphosite(organism=organism,dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene,sourceIdType=sourceIdType,standardID=standardSource,collapseMethod=collapseMethod,mappingOutput=mappingOutput,outputFileName=outputFileName,methodType=methodType,hostName=hostName)
     }
     
     if(.hasError(idmap)){
    			return(idmap)
     }else{
     	 idmap$standardId <- standardSource
		 	 return(idmap)
     }
}
