loadGeneSet <- function(organism="hsapiens",enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,collapseMethod="mean",methodType="R",hostName="http://www.webgestalt.org/"){
	geneSet <- NULL    ##gene sets
	geneSetDes <- NULL ##gene set description file
	geneSetDAG <- NULL ##gene set DAG file
	geneSetNet <- NULL ##gene set network file
	standardId <- NULL

	if(organism!="others"){
		geneSets <- listGeneSet(organism=organism,hostName=hostName)
		if(!(enrichDatabase %in% geneSets[,1])){  ###if the upload gene set name is not in the database
			if(!is.null(enrichDatabase)){
				if(enrichDatabase=="others"){    ##if the user upload their own data
					#cat("Because 'enrichDatabase' is 'others', user can upload their own gene sets using GMT file and WebGestaltR will transform ids in the gene sets to entrez ids based on the parameter 'enrichDatabaseType'!\n")
					#######Read GMT File and transform id##########
					if(!is.null(enrichDatabaseFile)){
						geneSet <- IDMapping(organism=organism,dataType="gmt",inputGeneFile=enrichDatabaseFile,sourceIdType=enrichDatabaseType,targetIdType=NULL,collapseMethod=collapseMethod,mappingOutput=FALSE,methodType=methodType,hostName=hostName)
						if(.hasError(geneSet)){
							return(geneSet)
						}
						standardId <- geneSet$standardid
						geneSet <- geneSet$mapped
						geneSet <- unique(geneSet[,c("geneset","link",standardId)])
						if(!is.null(enrichDatabaseDescriptionFile)){     ##upload description file
							geneSetDes <- .loadEnrichDatabaseDescriptionFile(geneSet,enrichDatabaseDescriptionFile)
							if(.hasError(geneSetDes)){
								return(geneSetDes)
							}
						}
					}else{  #no enrichDatabaseFile for 'others' enrichDatabase
						return(gmtFormatError("empty"))
					}
				}else{ ##input a wrong enrichDatabase name
					return(enrichDataBaseError(type="unsupported",enrichDatabase=enrichDatabase,organism=organism))
				}
			}else{  #enrichDatabase is NULL
				return(enrichDataBaseError(type="empty"))
			}
		}else{  #input a correct enrichDatabase
			standardId <- geneSets[geneSets[,1]==enrichDatabase,3]   ###get the ID type of the enriched database, such as entrezgene or phosphositeSeq

			#########Read GMT file from the existing database###########
			geneSet <- readGMT(file.path(hostName,"data","geneset",paste(organism,"_",enrichDatabase,"_",standardId,".gmt",sep="")))

			if(.hasError(geneSet)){
				return(geneSet)
			}
			#########Read the description file#############
			geneSetDesFile <- file.path(hostName,"data","geneset",paste(organism,"_",enrichDatabase,"_",standardId,".des",sep=""))
			geneSetDes <- tryCatch(fread(input=geneSetDesFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE),warning=function(e){return(NULL)},error=function(e){return(NULL)})  #####read the des file. If no des file, return NULL. For the des file, First column is the category id and the second is the description

			###########Try to load the DAG file#################
			geneSetDAGFile <- file.path(hostName,"data","geneset",paste(organism,"_",enrichDatabase,"_",standardId,".dag",sep=""))
			geneSetDAG <- tryCatch(fread(input=geneSetDAGFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE),warning=function(e){return(NULL)},error=function(e){return(NULL)})  #####read the dag file. If no dag file, return NULL. For the dag file, First column is the parent term and the second is the child term

			###########Try to load the network file if the gene sets are generated from the network##########
			geneSetNetFile <- file.path(hostName,"data","geneset",paste(organism,"_",enrichDatabase,"_",standardId,".net",sep=""))
			geneSetNet <- tryCatch(fread(input=geneSetNetFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE),warning=function(e){return(NULL)},error=function(e){return(NULL)})    ####Read the net file. If no net file, return NULL. For the net file, First column is the gene set name, the second is the network related to the gene set and the third column is the related functions.
		}
	}else{
		#########Read GMT file for other orgnisms from user files###########
		if(!is.null(enrichDatabaseFile)){
			geneSet <- readGMT(enrichDatabaseFile)
		if(.hasError(geneSet)){
			return(geneSet)
		}
			if(!is.null(enrichDatabaseDescriptionFile)){     ##upload description file
				geneSetDes <- .loadEnrichDatabaseDescriptionFile(geneSet,enrichDatabaseDescriptionFile)
				if(.hasError(geneSetDes)){
					return(geneSetDes)
				}
			}
		}else{  ##enrichDatabaseFile is NULL
			return(enrichDataBaseError(type="others"))
		}
	}

	re <- list(geneSet=geneSet,geneSetDes=geneSetDes,geneSetDAG=geneSetDAG,geneSetNet=geneSetNet,standardId=standardId)
	return(re)
}

.loadEnrichDatabaseDescriptionFile <- function(geneSet,enrichDatabaseDescriptionFile){
	if(file_extension(enrichDatabaseDescriptionFile)!="des"){
		return(descriptionFileError("format"))
	}else{
		geneSetDes <- fread(input=enrichDatabaseDescriptionFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
		if(ncol(geneSetDes)!=2){
			return(descriptionFileError("columnNum"))
		}else{
			if(length(intersect(unique(geneSet[,1]),geneSetDes[,1]))<0.6*length(unique(geneSet[,1]))){
				return(descriptionFileError("overlap"))
			}else{
				colnames(geneSetDes) <- c("geneset","description")
				return(geneSetDes)
			}
		}
	}
}
