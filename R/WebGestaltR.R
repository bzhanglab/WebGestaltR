#' Comprehensive R function for the enrichment analysis
#'
#' Main function for enrichment analysis
#'
#' WebGestaltR function can perform three enrichment analyses:
#' ORA (Over-Representation Analysis) and GSEA (Gene Set Enrichment Analysis).and
#' NTA (Network Topology Analysis). Based on the user-uploaded gene list or gene list
#' with scores, WebGestaltR function will first map the gene list to the entrez gene
#' ids and then summarize the gene list based on the GO (Gene Ontology) Slim. After
#' performing the enrichment analysis, WebGestaltR function also returns a user-friendly
#' HTML report containing GO Slim summary and the enrichment analysis result. If functional
#' categories have DAG (directed acyclic graph) structure or genes in the functional
#' categories have network structure, those relationship can also be visualized in the
#' report.
#'
#' @param enrichMethod Enrichment methods: \code{ORA}, \code{GSEA} or \code{NTA}.
#' @param organism Currently, WebGestaltR supports 12 organisms. Users can use the function
#'   \code{listOrganism} to check available organisms. Users can also input \code{others} to
#'   perform the enrichment analysis for other organisms not supported by WebGestaltR. For
#'   other organisms, users need to provide the functional categories, interesting list and
#'   reference list (for ORA method). Because WebGestaltR does not perform the ID mapping for
#'   the other organisms, the above data should have the same ID type.
#' @param enrichDatabase The functional categories for the enrichment analysis. Users can use
#'   the function \code{listGeneSet} to check the available functional databases for the
#'   selected organism. Users can also input \code{others} to provide a custom functional
#'   databases not supported by WebGestaltR for the selected organism.
#' @param enrichDatabaseFile If users set \code{organism} as \code{others} or set
#'   \code{enrichDatabase} as \code{others}, users need to provide a GMT file as the functional
#'   category for enrichment analysis. The extension of the file should be \code{gmt} and the
#'   first column of the file is the category ID, the second one is the external link for the
#'   category. Genes annotated to the category are from the third column. All columns are
#'   separated by tabs.
#' @param enrichDatabaseType If users set \code{enrichDatabase} as \code{others}, WebGestaltR
#'   will also perform ID mapping for the supplied GMT file. Thus, users need to set the ID
#'   type of the genes in the \code{enrichDatabaseFile}. If users set \code{organism} as
#'   \code{others}, users do not need to set this ID type because WebGestaltR will not perform
#'   ID mapping for other organisms. The supported ID types of WebGestaltR for the selected
#'   organism can be found by the function \code{listIdType}.
#' @param enrichDatabaseDescriptionFile Users can also provide a description file for the custom
#'   \code{enrichDatabaseFile}. The extension of the description file should be \code{des}. The
#'   description file contains two columns: the first column is the category ID that should be
#'   exactly the same as the category ID in the custom \code{enrichDatabaseFile} and the second
#'   column is the description of the category. All columns are separated by tabs.
#' @param interestGeneFile If \code{enrichMethod} is \code{ORA} or \code{NTA}, the extension of
#'   the \code{interestGeneFile} should be \code{txt} and the file can only contain one column:
#'   the interesting gene list. If \code{enrichMethod} is \code{GSEA}, the extension of the
#'   \code{interestGeneFile} should be \code{rnk} and the file should contain two columns
#'   separated by tab: the gene list and the corresponding scores.
#' @param interestGene Users can also use an R object as the input. If \code{enrichMethod} is
#'   \code{ORA} or \code{NTA}, \code{interestGene} should be an R \code{vector} object
#'   containing the interesting gene list. If \code{enrichMethod} is \code{GSEA},
#'   \code{interestGene} should be an R \code{data.frame} object containing two columns: the
#'   gene list and the corresponding scores.
#' @param interestGeneType The ID type of the interesting gene list. The supported ID types of
#'   WebGestaltR for the selected organism can be found by the function \code{listIdType}. If
#'   the \code{organism} is \code{others}, users do not need to set this parameter.
#' @param collapseMethod The method to collapse duplicate IDs with scores. \code{mean},
#'   \code{median}, \code{min} and \code{max} represent the mean, median, minimum and maximum
#'   of scores for the duplicate IDs.
#' @param referenceGeneFile For the ORA method, the users need to upload the reference gene
#'   list. The extension of the \code{referenceGeneFile} should be \code{txt} and the file can
#'   only contain one column: the reference gene list.
#' @param referenceGene For the ORA method, users can also use an R object as the reference
#'   gene list. \code{referenceGene} should be an R \code{vector} object containing the
#'   reference gene list.
#' @param referenceGeneType The ID type of the reference gene list. The supported ID types
#'   of WebGestaltR for the selected organism can be found by the function \code{listIdType}.
#'   If the \code{organism} is \code{others}, users do not need to set this parameter.
#' @param referenceSet Users can directly select the reference set from existing platforms in
#'   WebGestaltR and do not need to provide the reference set through \code{referenceGeneFile}.
#'   All existing platforms supported in WebGestaltR can be found by the function
#'   \code{listReferenceSet}. If \code{referenceGeneFile} and \code{refereneceGene} are
#'   \code{NULL}, WebGestaltR will use the \code{referenceSet} as the reference gene set.
#'   Otherwise, WebGestaltR will use the user supplied reference set for enrichment analysis.
#' @param minNum WebGestaltR will exclude the categories with the number of annotated genes
#'   less than \code{minNum} for enrichment analysis. The default is \code{10}.
#' @param maxNum WebGestaltR will exclude the categories with the number of annotated genes
#'   larger than \code{maxNum} for enrichment analysis. The default is \code{500}.
#' @param sigMethod Two methods of significance are available in WebGestaltR: \code{fdr} and
#'   \code{top}. \code{fdr} means the enriched categories are identified based on the FDR and
#'   \code{top} means all categories are ranked based on FDR and then select top categories
#'   as the enriched categories. The default is \code{fdr}.
#' @param fdrMethod For the ORA method, WebGestaltR supports five FDR methods: \code{holm},
#'   \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH} and \code{BY}. The default
#'   is \code{BH}.
#' @param fdrThr The significant threshold for the \code{fdr} method. The default is \code{0.05}.
#' @param topThr The threshold for the \code{top} method. The default is \code{10}.
#' @param reportNum The number of enriched categories visualized in the final report. The default
#'   is \code{20}. A larger \code{reportNum} may be slow to render in the report.
#' @param perNum The number of permutations for the GSEA method. The default is \code{1000}.
#' @param isOutput If \code{isOutput} is TRUE, WebGestaltR will create a folder named by
#'   the \code{projectName} and save the results in the folder. Otherwise, WebGestaltR will
#'   only return an R \code{data.frame} object containing the enrichment results. If
#'   hundreds of gene list need to be analyzed simultaneously, it is better to set
#'   \code{isOutput} to \code{FALSE}. The default is \code{TRUE}.
#' @param outputDirectory The output directory for the results.
#' @param projectName The name of the project. If \code{projectName} is \code{NULL},
#'   WebGestaltR will use time stamp as the project name.
#' @param dagColor If \code{dagColor} is \code{binary}, the significant terms in the DAG
#'   structure will be colored by steel blue for ORA method or steel blue (positive related)
#'   and dark orange (negative related) for GSEA method. If \code{dagColor} is \code{continous},
#'   the significant terms in the DAG structure will be colored by the color gradient based on
#'   corresponding FDRs.
#' @param setCoverNum The number of expected gene sets after set cover to reduce redundancy.
#'   It could get fewer sets if the coverage reaches 100\%. The default is \code{10}.
#' @param networkConstructionMethod Netowrk construction method for NTA. Either
#'   \code{Network_Retrieval_Prioritization} or \code{Network_Expansion}. Network Retrieval &
#'   Prioritization first uses random walk analysis to calculate random walk probabilities
#'   for the input seeds, then identifies the relationships among the seeds in the selected
#'   network and returns a retrieval sub-network. The seeds with the top random walk
#'   probabilities are highlighted in the sub-network. Network Expansion first uses random
#'   walk analysis to rank all genes in the selected network based on their network
#'   proximity to the input seeds and then return an expanded sub-network in which nodes
#'   are the input seeds and their top ranking neighbors and edges represent their
#'   relationships.
#' @param neighborNum The number of neighbors to include in NTA Network Expansion method.
#' @param highlightType The type of nodes to highlight in the NTA Network Expansion method,
#'   either \code{Seeds} or \code{Neighbors}.
#' @param highlightSeedNum The number of top input seeds to highlight in NTA Network Retrieval
#'   & Prioritizaiton method.
#' @param nThreads The number of cores to use for GSEA and set cover, and in batch function.
#' @param hostName The server URL for accessing data. Mostly for development purposes.
#' @param ... In batch function, passes parameters to WebGestaltR function.
#'   Also handles backward compatibility for some parameters in old versions.
#'
#' @return The WebGestaltR function returns a data frame containing the enrichment analysis
#'   result and also outputs an user-friendly HTML report if \code{isOutput} is \code{TRUE}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ####### ORA example #########
#' geneFile <- system.file("extdata", "interestingGenes.txt", package="WebGestaltR")
#' refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
#' outputDirectory <- getwd()
#' enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
#'   enrichDatabase="pathway_KEGG", interestGeneFile=geneFile,
#'   interestGeneType="genesymbol", referenceGeneFile=refFile,
#'   referenceGeneType="genesymbol", isOutput=TRUE,
#'   outputDirectory=outputDirectory, projectName=NULL)
#'
#' ####### GSEA example #########
#' rankFile <- system.file("extdata", "GeneRankList.rnk", package="WebGestaltR")
#' outputDirectory <- getwd()
#' enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
#'   enrichDatabase="pathway_KEGG", interestGeneFile=rankFile,
#'   interestGeneType="genesymbol", sigMethod="top", topThr=10, minNum=5,
#'   outputDirectory=outputDirectory)
#'
#' ####### NTA example #########
#' enrichResult <- WebGestaltR(enrichMethod="NTA", organism="hsapiens",
#'   enrichDatabase="network_PPI_BIOGRID", interestGeneFile=geneFile,
#'   interestGeneType="genesymbol", sigMethod="top", topThr=10,
#'   outputDirectory=getwd(), highlightSeedNum=10,
#'   networkConstructionMethod="Network_Retrieval_Prioritization")
#' }
#'
WebGestaltR <- function(enrichMethod="ORA", organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, sigMethod="fdr", fdrMethod="BH", fdrThr=0.05, topThr=10, reportNum=20, perNum=1000, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, dagColor="continuous", setCoverNum=10, networkConstructionMethod=NULL, neighborNum=10, highlightType="Seeds", highlightSeedNum=10, nThreads=1, hostName="http://www.webgestalt.org/",  ...) {
	extraArgs <- list(...)
	if ('keepGSEAFolder' %in% names(extraArgs) | 'keepGseaFolder' %in% names(extraArgs)) {
		cat("WARNING: Parameter keepGSEAFolder is obsolete.\n")
	}
	if ('is.output' %in% names(extraArgs)) {
		isOutput <- extraArgs$is.output
		cat("WARNING: Parameter is.output is deprecated and changed to isOutput!\n")
		warning("Column names of the result data frame are modified.")
	}
	if ('methodType' %in% names(extraArgs)) {
		cat("WARNING: Parameter methodType is obsolete.\n")
	}
	if ('lNum' %in% names(extraArgs)) {
		warning("Parameter lNum is obsolete.\n")
	}
	if ('dNum' %in% names(extraArgs)) {
		cat("WARNING: Parameter dNum is deprecated and changed to reportNum.\n")
		reportNum <- extraArgs$dNum
	}

	## TODO: add para test for NTA
	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, perNum=perNum, isOutput=isOutput, outputDirectory=outputDirectory, dagColor=dagColor, hostName=hostName)
	if(!is.null(errorTest)){
		return(errorTest)
	}

	if(is.null(projectName)){
		projectName <- as.character(as.integer(Sys.time()))
	}
	projectName <- sanitizeFileName(projectName) # use for GOSlim summary file name, convert punct to _
	if (enrichMethod == "ORA") {
		enrichR <- WebGestaltROra(organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,  interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, referenceGeneFile=referenceGeneFile, referenceGene=referenceGene, referenceGeneType=referenceGeneType, referenceSet=referenceSet, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, setCoverNum=setCoverNum, isOutput=isOutput, outputDirectory=outputDirectory, projectName=projectName, dagColor=dagColor, nThreads=nThreads, hostName=hostName)
	} else if (enrichMethod == "GSEA") {
		enrichR <- WebGestaltRGsea(organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,  interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, setCoverNum=setCoverNum, perNum=perNum, isOutput=isOutput, outputDirectory=outputDirectory, projectName=projectName, dagColor=dagColor, nThreads=nThreads, hostName=hostName)
	} else if (enrichMethod == "NTA") {
		enrichR <- WebGestaltRNta(organism=organism, network=enrichDatabase, method=networkConstructionMethod, neighborNum=neighborNum, highlightSeedNum=highlightSeedNum, inputSeed=interestGene, inputSeedFile=interestGeneFile, interestGeneType=interestGeneType, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, outputDirectory=outputDirectory, projectName=projectName, highlightType=highlightType, hostName=hostName)
	}

	return(enrichR)
}
