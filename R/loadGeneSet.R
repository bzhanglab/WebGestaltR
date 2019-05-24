#' Load gene set data
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of \code{geneSet}, \code{geneSetDes}, \code{geneSetDag}, \code{geneSetNet}, \code{standardId}.
#' \describe{
#'  \item{geneSet}{Gene set: A data frame with columns of "geneSet", "description", "genes"}
#'  \item{geneSetDes}{Description: A data frame with columns of two columns of gene set ID and description}
#'  \item{geneSetDag}{DAG: A edge list data frame of two columns of parent and child. Or a list of data frames if multilple databases are given.}
#'  \item{geneSetNet}{Network: A edge list data frame of two columns connecting nodes. Or a list of data frames if multilple databases are given.}
#'  \item{standardId}{The standard ID of the gene set}
#' }
#'
#' @importFrom dplyr select distinct filter %>%
#' @importFrom httr modify_url
#' @export
loadGeneSet <- function(organism="hsapiens", enrichDatabase=NULL, enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, cache=NULL, hostName="http://www.webgestalt.org/") {
	# TODO: multiple custom database ID types?
	geneSet <- NULL    ##gene sets
	geneSetDes <- NULL ##gene set description file
	geneSetDag <- list() ##gene set DAG file
	geneSetNet <- list() ##gene set network file
	standardId <- NULL

	if (organism != "others" && !is.null(enrichDatabaseFile) && is.null(enrichDatabaseType)) {
		stop("The ID type should be given in enrichDatabaseType for custom GMT files, e.g. genesymbol.")
	}
	# necessary because loop with length used below. enrichDatabase will skip NULL
	if (!is.vector(enrichDatabaseFile)) {
		enrichDatabaseFile = list(enrichDatabaseFile)
	}
	if (!is.vector(enrichDatabaseDescriptionFile)) {
		enrichDatabaseDescriptionFile = list(enrichDatabaseDescriptionFile)
	}
	if (length(enrichDatabaseFile) != length(enrichDatabaseDescriptionFile)) {
		stop("The number of custom database and its description files should be equal. Use NULL for placeholder.")
	}
	if (organism != "others") {  # supported organism
		geneSetInfo <- listGeneSet(organism=organism, hostName=hostName, cache=cache)
		#  load build-in databases
		for (enrichDb in enrichDatabase) {
			if (is.null(enrichDb) || enrichDb == "others") { next }  # just for backward compatibility
			if (!enrichDb %in% geneSetInfo$name) {
				warning("Database ", enrichDb, " is not supported")
				next
			}
			# get the ID type of the enriched database, such as entrezgene or phosphsiteSeq
			thisStandardId <- filter(geneSetInfo, .data$name==enrichDb)[[1, "idType"]]
			if (!is.null(standardId) && standardId != thisStandardId) {
				stop("Databases have inconsistent ID types. Mixed gene annotation databases with phosphosite databases?")
			}
			standardId <- thisStandardId

			#########Read GMT file from the existing database###########
			if (startsWith(hostName, "file://")) {
				gmtPath <- removeFileProtocol(file.path(hostName, "geneset", paste0(paste(organism, enrichDb, standardId, sep="_"), ".gmt")))
				thisGeneSet <- readGmt(gmtPath)
				thisGeneSet$database <- enrichDb # add a column for database source
				geneSet <- rbind(geneSet, thisGeneSet)
			} else {
				gmtUrl <- modify_url(file.path(hostName, "api", "geneset"), query=list(organism=organism, database=enrichDb, standardId=standardId, fileType="gmt"))
				thisGeneSet <-  readGmt(gmtUrl, cache=cache)
				thisGeneSet$database <- enrichDb
				oriCnt <- nrow(thisGeneSet)
				thisGeneSet <- thisGeneSet %>% filter(!(.data$geneSet %in% unique(!!geneSet$geneSet)))
				if (nrow(thisGeneSet) < oriCnt) {
				  warning(paste("Duplicate gene set names in", enrichDb, "have been ignored."))
				}
				geneSet <- rbind(geneSet, thisGeneSet)
			}

			#########Read the description file#############
			#geneSetDes <- rbind(geneSetDes, .loadGeneSetData(hostName, organism, enrichDb, standardId, "des", cache))
			thisGeneSetDes <- .loadGeneSetData(hostName, organism, enrichDb, standardId, "des", cache)
			if (!is.null(thisGeneSetDes) && !is.null(geneSetDes)) {
			  thisGeneSetDes <- thisGeneSetDes %>% filter(!(.data$geneSet %in% unique(!!geneSetDes$geneSet)))
			}
			geneSetDes <- rbind(geneSetDes, thisGeneSetDes)

			###########Try to load the DAG file#################
			# assignment considering possible return of NULL
			# list[] <- NULL will delete the element
			geneSetDag[enrichDb] <- list(.loadGeneSetData(hostName, organism, enrichDb, standardId, "dag", cache))

			###########Try to load the network file if the gene sets are generated from the network##########
			geneSetNet[enrichDb] <- list(.loadGeneSetData(hostName, organism, enrichDb, standardId, "net", cache))
		}

		# load local database files
		for (i in 1:length(enrichDatabaseFile)) {
			enrichDbFile <- enrichDatabaseFile[[i]]
			if (is.null(enrichDbFile)) { next }

			thisGeneSet <- idMapping(organism=organism, dataType="gmt", inputGeneFile=enrichDbFile, sourceIdType=enrichDatabaseType, targetIdType=NULL, mappingOutput=FALSE, cache=cache, hostName=hostName)
			thisStandardId <- thisGeneSet$standardId  # should be just enrichDatabaseType here
			if (!is.null(standardId) && standardId != thisStandardId) {
				stop("Databases have inconsistent ID types. Mixed gene annotation databases with phosphosite databases?")
			}
			standardId <- thisStandardId

			thisGeneSet <- thisGeneSet$mapped %>% select(.data$geneSet, .data$description, gene=standardId) %>% distinct()

			# load local description files
			enrichDbDesFile <- enrichDatabaseDescriptionFile[[i]]
			if (!is.null(enrichDbDesFile)) {
				thisGeneSetDes <- .loadEnrichDatabaseDescriptionFile(thisGeneSet, enrichDbDesFile)
				if (!is.null(thisGeneSetDes) && !is.null(geneSetDes)) {
				  thisGeneSetDes <- thisGeneSetDes %>% filter(!(.data$geneSet %in% unique(!!geneSetDes$geneSet)))
				}
				geneSetDes <- rbind(geneSetDes, thisGeneSetDes)
			}

			fileName <- gsub(".gmt", "", basename(enrichDbFile), fixed=TRUE)
			thisGeneSet$database <- fileName
			oriCnt <- nrow(thisGeneSet)
			thisGeneSet <- thisGeneSet %>% filter(!(.data$geneSet %in% unique(!!geneSet$geneSet)))
			if (nrow(thisGeneSet) < oriCnt) {
			  warning(paste("Duplicate gene set names in", fileName, "have been ignored."))
			}
			geneSet <- rbind(geneSet, thisGeneSet)

			geneSetDag[fileName] <- list(NULL) # correct way to assign NULL to list element
			geneSetNet[fileName] <- list(NULL)
		}
	} else { # custom organisms
		for (i in 1:length(enrichDatabaseFile)) {
			enrichDbFile <- enrichDatabaseFile[[i]]
			if (is.null(enrichDbFile)) { next }
			thisGeneSet <- readGmt(enrichDbFile)

			enrichDbDesFile <- enrichDatabaseDescriptionFile[[i]]
			if (!is.null(enrichDbDesFile)) {
				thisGeneSetDes <- .loadEnrichDatabaseDescriptionFile(thisGeneSet, enrichDbDesFile)
				if (!is.null(thisGeneSetDes) && !is.null(geneSetDes)) {
				  thisGeneSetDes <- thisGeneSetDes %>% filter(!(.data$geneSet %in% unique(!!geneSetDes$geneSet)))
				}
				geneSetDes <- rbind(geneSetDes, thisGeneSetDes)
			}
			fileName <- gsub(".gmt", "", basename(enrichDbFile), fixed=TRUE)
			thisGeneSet$database <- fileName
			oriCnt <- nrow(thisGeneSet)
			thisGeneSet <- thisGeneSet %>% filter(!(.data$geneSet %in% unique(!!geneSet$geneSet)))
			if (nrow(thisGeneSet) < oriCnt) {
			  warning(paste("Duplicate gene set names in", enrichDb, "have been ignored."))
			}
			geneSet <- rbind(geneSet, thisGeneSet)

			geneSetDag[fileName] <- list(NULL)
			geneSetNet[fileName] <- list(NULL)
		}
	}
	if (is.null(geneSet)) { stop(enrichDatabaseError(type="empty")) }

	# unlist if just one database, for backward compatibility
	if (length(geneSetDag) == 1) { geneSetDag <- geneSetDag[[1]] }
	if (length(geneSetNet) == 1) { geneSetNet <- geneSetNet[[1]] }

	if (length(unique(geneSet$database)) == 1) {
		# remove database column for single source
		geneSet <- select(geneSet, -.data$database)
	}
	re <- list(geneSet=geneSet, geneSetDes=geneSetDes, geneSetDag=geneSetDag, geneSetNet=geneSetNet,standardId=standardId)
	return(re)
}

#' @importFrom httr content
#' @importFrom readr read_tsv
.loadGeneSetData <- function(hostName, organism, database, standardId, fileType, cache=NULL) {
	# read gene set files from API or returns NULL
	if (startsWith(hostName, "file://")) {
		geneSetPath <- removeFileProtocol(file.path(hostName, "geneset", paste(paste(organism, database, standardId, sep="_"), fileType, sep=".")))
		if (file.exists(geneSetPath)) {
			geneSetData <- read_tsv(geneSetPath, col_names=FALSE, col_types="cc")
		} else {
			geneSetData <- NULL
		}
	} else {
		geneSetUrl <- file.path(hostName,"api","geneset")
		response <- cacheUrl(geneSetUrl, cache=cache, query=list(organism=organism, database=database, standardId=standardId, fileType=fileType))
		if (response$status_code == 200) {
			geneSetData <- read_tsv(content(response), col_names=FALSE, col_types="cc")
		} else {
			geneSetData <- NULL
		}
	}
	if (!is.null(geneSetData) && fileType == "des") {
		colnames(geneSetData) <- c("geneSet", "description")
	}
	return(geneSetData)
}

#' @importFrom readr read_tsv
#' @importFrom tools file_ext
.loadEnrichDatabaseDescriptionFile <- function(geneSet, enrichDatabaseDescriptionFile) {
	if (file_ext(enrichDatabaseDescriptionFile) != "des") {
		warning(descriptionFileError("format"))
		return(NULL)
	} else {
		geneSetDes <- read_tsv(enrichDatabaseDescriptionFile, col_names=c("geneSet", "description"), col_types="cc")
		if (ncol(geneSetDes)!=2) {
			warning(descriptionFileError("columnNum"))
			return(NULL)
		} else {
			if (length(intersect(unique(geneSet$geneSet), geneSetDes$geneSet)) < 0.6 * length(unique(geneSet[,1]))) {
				warning(descriptionFileError("overlap"))
				return(NULL)
			} else {
				return(geneSetDes)
			}
		}
	}
}
