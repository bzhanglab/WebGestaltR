#' @importFrom dplyr select distinct filter %>%
#' @importFrom httr modify_url
loadGeneSet <- function(organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, collapseMethod="mean", hostName="http://www.webgestalt.org/") {
	geneSet <- NULL    ##gene sets
	geneSetDes <- NULL ##gene set description file
	geneSetDag <- NULL ##gene set DAG file
	geneSetNet <- NULL ##gene set network file
	standardId <- NULL

	if(organism!="others"){
		geneSetInfo <- listGeneSet(organism=organism,hostName=hostName)
		if (!(enrichDatabase %in% geneSetInfo$name)) {  ###if the upload gene set name is not in the database
			if(!is.null(enrichDatabase)){
				if(enrichDatabase=="others"){    ##if the user upload their own data
					#cat("Because 'enrichDatabase' is 'others', user can upload their own gene sets using GMT file and WebGestaltR will transform ids in the gene sets to entrez ids based on the parameter 'enrichDatabaseType'!\n")
					#######Read GMT File and transform id##########
					if(!is.null(enrichDatabaseFile)){
						geneSet <- idMapping(organism=organism, dataType="gmt", inputGeneFile=enrichDatabaseFile, sourceIdType=enrichDatabaseType, targetIdType=NULL, collapseMethod=collapseMethod, mappingOutput=FALSE, hostName=hostName)
						if(.hasError(geneSet)){
							return(geneSet)
						}
						standardId <- geneSet$standardId
						geneSet <- geneSet$mapped %>% select(.data$geneSet, .data$description, gene=.data$entrezgene) %>% distinct()
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
					return(enrichDatabaseError(type="unsupported",enrichDatabase=enrichDatabase,organism=organism))
				}
			}else{  #enrichDatabase is NULL
				return(enrichDatabaseError(type="empty"))
			}
		}else{  #input a correct enrichDatabase
			standardId <- filter(geneSetInfo, .data$name==enrichDatabase)[[1, "idType"]]  # get the ID type of the enriched database, such as entrezgene or phosphsiteSeq

			#########Read GMT file from the existing database###########
			gmtUrl <- modify_url(file.path(hostName,"api","geneset"), query=list(organism=organism, database=enrichDatabase, standardId=standardId, fileType="gmt"))
			geneSet <- readGmt(gmtUrl)

			if(.hasError(geneSet)){
				return(geneSet)
			}
			#########Read the description file#############
			geneSetDes <- .loadGeneSetData(hostName, organism, enrichDatabase, standardId, "des")

			###########Try to load the DAG file#################
			geneSetDag <- .loadGeneSetData(hostName, organism, enrichDatabase, standardId, "dag")

			###########Try to load the network file if the gene sets are generated from the network##########
			geneSetNet <- .loadGeneSetData(hostName, organism, enrichDatabase, standardId, "net")
		}
	}else{
		#########Read GMT file for other orgnisms from user files###########
		if(!is.null(enrichDatabaseFile)){
			geneSet <- readGmt(enrichDatabaseFile)
			if(.hasError(geneSet)){
				return(geneSet)
			}
			if(!is.null(enrichDatabaseDescriptionFile)){     ##upload description file
				geneSetDes <- .loadEnrichDatabaseDescriptionFile(geneSet, enrichDatabaseDescriptionFile)
				if(.hasError(geneSetDes)){
					return(geneSetDes)
				}
			}
		}else{  ##enrichDatabaseFile is NULL
			return(enrichDatabaseError(type="others"))
		}
	}

	re <- list(geneSet=geneSet,geneSetDes=geneSetDes,geneSetDag=geneSetDag,geneSetNet=geneSetNet,standardId=standardId)
	return(re)
}

#' @importFrom httr GET content
#' @importFrom readr read_tsv
.loadGeneSetData <- function(hostName, organism, database, standardId, fileType) {
	# read gene set files from API or returns NULL
	geneSetUrl <- file.path(hostName,"api","geneset")
	response <- GET(geneSetUrl, query=list(organism=organism, database=database, standardId=standardId, fileType=fileType))
	if (response$status_code == 200) {
		geneSetData <- read_tsv(content(response), col_names=FALSE, col_types="cc")
	} else {
		geneSetData <- NULL
	}
	return(geneSetData)
}

#' @importFrom readr read_tsv
#' @importFrom tools file_ext
.loadEnrichDatabaseDescriptionFile <- function(geneSet, enrichDatabaseDescriptionFile){
	if(file_ext(enrichDatabaseDescriptionFile)!="des"){
		return(descriptionFileError("format"))
	}else{
		geneSetDes <- read_tsv(enrichDatabaseDescriptionFile, col_names=c("geneSet", "description"), col_types="cc")
		if(ncol(geneSetDes)!=2){
			return(descriptionFileError("columnNum"))
		}else{
			if(length(intersect(unique(geneSet$geneSet), geneSetDes$geneSet)) < 0.6 * length(unique(geneSet[,1]))){
				return(descriptionFileError("overlap"))
			}else{
				return(geneSetDes)
			}
		}
	}
}
