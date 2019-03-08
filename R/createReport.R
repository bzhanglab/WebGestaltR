#' createReport
#'
#' Generate HTML report for ORA and GSEA
#'
#' @importFrom rjson toJSON
#' @importFrom whisker whisker.render rowSplit
#'
#' @keywords internal
#'
createReport <- function(hostName, outputDirectory, organism="hsapiens", projectName, enrichMethod, geneSet, geneSetDes, geneSetDag, geneSetNet, interestingGeneMap, referenceGeneList, enrichedSig, geneTables, clusters, background, enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, reportNum=20, perNum=1000, dagColor="binary"){
	outputHtmlFile <- file.path(outputDirectory, paste0("Project_", projectName), paste0("Report_", projectName, ".html"))

	# if hostname starts with "file://", it is used as WebGestaltReporter
	if (startsWith(hostName, "file://")) {
		# change back hostName for web assets and browsers will cache it.
		hostName <- "http://www.webgestalt.org"
	}

	numAnnoRefUserId <- NULL
	dagJson <- list()
	allEnrichedSig <- enrichedSig
	repAdded <- FALSE
	if(organism!="others"){
		if (!is.null(enrichedSig) && reportNum < nrow(enrichedSig)) {
			if (enrichMethod == "ORA") {
				enrichedSig <- enrichedSig[1:reportNum, ]
			} else if (enrichMethod == "GSEA") {
				enrichedSig <- getTopGseaResults(enrichedSig, reportNum / 2)[[1]]
			}
			# Add representatives if they are not in top ReportNum. So could be more if ReportNum.is small and high redundancy in top
			numRes <- nrow(enrichedSig)
			enrichedSig <- keepRep(enrichedSig, allEnrichedSig, clusters$ap$representatives)
			enrichedSig <- keepRep(enrichedSig, allEnrichedSig, clusters$wsc$representatives)
			repAdded <- nrow(enrichedSig) > numRes
		}

		#####Summary Tab########
		bodyContent <- summaryDescription(projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase, enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, allEnrichedSig, reportNum, perNum, geneSet, repAdded, hostName)

		standardId <- interestingGeneMap$standardId
		if (enrichMethod == 'ORA') {
			interestGeneList <- unique(interestingGeneMap$mapped[[standardId]])
			numAnnoRefUserId <- length(intersect(interestGeneList, intersect(referenceGeneList, geneSet$gene)))
		}
		###########GOSlim summary#########################
		if(standardId=="entrezgene"){
			bodyContent <- paste(bodyContent, goSlimReport(projectName), sep='\n')
		}

		############Enrichment result##################
		if(!is.null(enrichedSig)){
			bodyContent <- paste(bodyContent, enrichResultSection(enrichMethod, geneSetDes, geneSetDag, geneSetNet, clusters), seq='\n')
			if (!is.null(geneSetDag)) {
				dagRes <- expandDag(enrichedSig$geneSet, geneSetDag)
				dagEdges <- dagRes$edges
				dagNodes <- getDagNodes(enrichedSig, dagRes$allNodes, geneSetDes, enrichMethod, dagColor)
				dagJson <- c(dagEdges, dagNodes)
			}
		}
	}else{
		###########Organism is others. No mapping information#############
		#############summary for the analysis###################
		if (enrichMethod == 'ORA') {
			numAnnoRefUserId <- length(intersect(interestingGeneMap, intersect(referenceGeneList, geneSet$gene)))
		}
		if (!is.null(enrichedSig) && reportNum < nrow(enrichedSig)) {
			if (enrichMethod == "ORA") {
				enrichedSig <- enrichedSig[1:reportNum, ]
			} else if (enrichMethod == "GSEA") {
				enrichedSig <- getTopGseaResults(enrichedSig, reportNum / 2)[[1]]
			}
			# Add representatives if they are not in top ReportNum. So could be more if ReportNum.is small and high redundancy in top
			numRes <- nrow(enrichedSig)
			enrichedSig <- keepRep(enrichedSig, allEnrichedSig, clusters$ap$representatives)
			enrichedSig <- keepRep(enrichedSig, allEnrichedSig, clusters$wsc$representatives)
			repAdded <- nrow(enrichedSig) > numRes
		}

		bodyContent <- summaryDescription(projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase, enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, allEnrichedSig, reportNum, perNum, repAdded, geneSet)

		##############Enrich Result################
		if(!is.null(enrichedSig)){
			bodyContent <- paste(bodyContent, enrichResultSection(enrichMethod, geneSetDes, geneSetDag, geneSetNet, clusters), seq='\n')
		}
		standardId <- NULL
	}
	if (is.null(enrichedSig)) {
		enrichedSig <- data.frame()
	}
	if (is.null(background)) {
		background <- data.frame()
	}
	version <- packageVersion("WebGestaltR")
	# use major and minor version numbers. If API changes, version should be bumped
	# patch number should not matter
	version <- paste(version[1, 1], version[1, 2], sep=".")
	hasGeneSetDag = !is.null(geneSetDag)
	hasCytoscape <- hasGeneSetDag || !is.null(geneSetNet) # DAG or network needs cytoscape

	header <- readLines(system.file("templates/header.mustache", package="WebGestaltR"))
	footer <- readLines(system.file("templates/footer.mustache", package="WebGestaltR"))
	template <- readLines(system.file("templates/template.mustache", package="WebGestaltR"))
	data <- list(hostName=hostName, geneSetNet=geneSetNet, bodyContent=bodyContent,
				organism=organism, enrichDatabase=enrichDatabase,
				sigJson=toJSON(unname(rowSplit(enrichedSig))), insigJson=toJSON(unname(rowSplit(background))),
				dagJson=toJSON(unname(dagJson)), hasGeneSetDag=hasGeneSetDag, version=version,
				clusterJson=toJSON(clusters), hasCytoscape=hasCytoscape,
				geneTableJson=toJSON(geneTables), standardId=standardId, numAnnoRefUserId=numAnnoRefUserId,
				methodIsGsea=enrichMethod=="GSEA", hasGeneSetDes=!is.null(geneSetDes)
				)
	cat(whisker.render(template, data, partials=list(header=header, footer=footer)), file=outputHtmlFile)
}
