#' @title WebGestaltRMultiOmics
#' @description Perform multi-omics analysis using WebGestaltR
#' @param enrichMethod Enrichment methods: \code{ORA}or \code{GSEA}.
#' @param organism Currently, WebGestaltR supports 12 organisms. Users can use the function
#'   \code{listOrganism} to check available organisms. Users can also input \code{others} to
#'   perform the enrichment analysis for other organisms not supported by WebGestaltR. For
#'   other organisms, users need to provide the functional categories, interesting list and
#'   reference list (for ORA method). Because WebGestaltR does not perform the ID mapping for
#'   the other organisms, the above data should have the same ID type.
#' @param enrichDatabase The functional categories for the enrichment analysis. Users can use
#'   the function \code{listGeneSet} to check the available functional databases for the
#'   selected organism. Multiple databases in a vector are supported for ORA and GSEA.
#' @param enrichDatabaseFile Users can provide one or more GMT files as the functional
#'   category for enrichment analysis. The extension of the file should be \code{gmt} and the
#'   first column of the file is the category ID, the second one is the external link for the
#'   category. Genes annotated to the category are from the third column. All columns are
#'   separated by tabs. The GMT files will be combined with \code{enrichDatabase}.
#' @param enrichDatabaseType The ID type of the genes in the \code{enrichDatabaseFile}.
#'   If users set \code{organism} as \code{others}, users do not need to set this ID type because
#'   WebGestaltR will not perform ID mapping for other organisms. The supported ID types of
#'   WebGestaltR for the selected organism can be found by the function \code{listIdType}.
#' @param enrichDatabaseDescriptionFile Users can also provide description files for the custom
#'   \code{enrichDatabaseFile}. The extension of the description file should be \code{des}. The
#'   description file contains two columns: the first column is the category ID that should be
#'   exactly the same as the category ID in the custom \code{enrichDatabaseFile} and the second
#'   column is the description of the category. All columns are separated by tabs.
#' @param analyteListFiles If \code{enrichMethod} is \code{ORA}, the extension of
#'   the \code{analyteListFiles} should be \code{txt} and each file can only contain one column:
#'   the interesting analyte list. If \code{enrichMethod} is \code{GSEA}, the extension of the
#'   \code{analyteListFiles} should be \code{rnk} and the files should contain two columns
#'   separated by tab: the analyte list and the corresponding scores.
#' @param analyteLists Users can also use an R object as the input. If \code{enrichMethod} is
#'   \code{ORA}, \code{analyte} should be an R \code{vector} contain multiple R \code{vector} object
#'   containing the interesting analyte lists. If \code{enrichMethod} is \code{GSEA},
#'   \code{analyteLists} should be an \code{vector} of R \code{data.frame} objects containing two columns: the
#'   gene list and the corresponding scores.
#' @param analyteLists \code{vector} of the ID type of the corresponding interesting analyte list. The supported ID types of
#'   WebGestaltR for the selected organism can be found by the function \code{listIdType}. If
#'   the \code{organism} is \code{others}, users do not need to set this parameter. The length of \code{analyteLists} should be
#'   the same as the length of \code{analyteListFiles} or \code{analyteLists}.
#' @param collapseMethod The method to collapse duplicate IDs with scores. \code{mean},
#'   \code{median}, \code{min} and \code{max} represent the mean, median, minimum and maximum
#'   of scores for the duplicate IDs.
#' @param referenceListFiles For the ORA method, the users need to upload the reference gene
#'   list. The extension of the \code{referenceListFile} should be \code{txt} and the file can
#'   only contain one column: the reference gene list.
#' @param referenceLists For the ORA method, users can also use an R object as the reference
#'   gene list. \code{referenceLists} should be an R \code{vector} object containing the
#'   reference gene list.
#' @param referenceGeneType The ID type of the reference gene list. The supported ID types
#'   of WebGestaltR for the selected organism can be found by the function \code{listIdType}.
#'   If the \code{organism} is \code{others}, users do not need to set this parameter.
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
#' @param gseaP The exponential scaling factor of the phenotype score. The default is \code{1}.
#'   When p=0, ES reduces to standard K-S statistics (See original paper for more details).
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
#' @param saveRawGseaResult Whether the raw result from GSEA is saved as a RDS file, which can be
#'   used for plotting. Defaults to \code{FALSE}. The list includes
#'   \describe{
#'     \item{Enrichment_Results}{A data frame of GSEA results with statistics}
#'     \item{Running_Sums}{A matrix of running sum of scores for each gene set}
#'     \item{Items_in_Set}{A list with ranks of genes for each gene set}
#'  }
#' @param gseaPlotFormat The graphic format of GSEA enrichment plots. Either \code{svg},
#'   \code{png}, or \code{c("png", "svg")} (default).
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
#' @param cache A directory to save data cache for reuse. Defaults to \code{NULL} and disabled.
#' @param hostName The server URL for accessing data. Mostly for development purposes.
#' @export
WebGestaltRMultiOmics <- function(analyteLists = NULL, analyteListFiles = NULL, analyteTypes = NULL, enrichMethod = "ORA", organism = "hsapiens",
                                  enrichDatabase = NULL, enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL,
                                  collapseMethod = "mean", minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05,
                                  topThr = 10, reportNum = 20, setCoverNum = 10, perNum = 1000, p = 1, isOutput = TRUE, outputDirectory = getwd(),
                                  projectName = NULL, dagColor = "binary", saveRawGseaResult = FALSE, plotFormat = "png", nThreads = 1, cache = NULL,
                                  hostName = "https://www.webgestalt.org/", useWeightedSetCover = TRUE, useAffinityPropagation = FALSE,
                                  usekMedoid = FALSE, kMedoid_k = 10, isMetaAnalysis = TRUE, mergeMethod = "mean", normalizationMethod = "rank",
                                  referenceLists = NULL, referenceListFiles = NULL, referenceTypes = NULL) {
  VALID_MERGE_METHODS <- c("mean", "max")
  VALID_NORM_METHODS <- c("rank", "median", "mean")
  VALID_ENRICH_METHODS <- c("ORA", "GSEA")

  # Null Check
  analyteLists <- testNull(analyteLists)
  analyteListFiles <- testNull(analyteListFiles)
  analyteTypes <- testNull(analyteTypes)
  referenceLists <- testNull(referenceLists)
  referenceTypes <- testNull(referenceTypes)
  referenceListFiles <- testNull(referenceListFiles)
  enrichMethod <- testNull(enrichMethod)
  organism <- testNull(organism)
  enrichDatabase <- testNull(enrichDatabase)
  enrichDatabaseFile <- testNull(enrichDatabaseFile)
  enrichDatabaseType <- testNull(enrichDatabaseType)
  enrichDatabaseDescriptionFile <- testNull(enrichDatabaseDescriptionFile)
  error_msg <- parameterErrorMessage(
    enrichMethod = enrichMethod, organism = organism, collapseMethod = collapseMethod, minNum = minNum, maxNum = maxNum,
    fdrMethod = fdrMethod, sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, reportNum = reportNum, isOutput = isOutput,
    outputDirectory = outputDirectory, dagColor = dagColor, hostName = hostName, cache = cache
  )
  if (!is.null(error_msg)) {
    stop(error_msg)
  }

  # Verify parameters
  mergeMethod <- tolower(mergeMethod)
  normalizationMethod <- tolower(normalizationMethod)
  enrichMethod <- toupper(enrichMethod)
  if (enrichMethod == "ORA" && !isMetaAnalysis) {
    stop("ORA only supports meta-analysis. isMetaAnalysis must be set to TRUE")
  }
  if (!(mergeMethod %in% VALID_MERGE_METHODS)) {
    stop(paste0(mergeMethod, " is not a valid merge method.\nValid options are: ", paste(VALID_MERGE_METHODS, collapse = ", ")))
  }
  if (!(normalizationMethod %in% VALID_NORM_METHODS)) {
    stop(paste0(normalizationMethod, " is not a valid normalization method.\nValid options are: ", paste(VALID_NORM_METHODS, collapse = ", ")))
  }
  if (length(analyteLists) != length(analyteTypes) && length(analyteListFiles) != length(analyteTypes)) {
    stop("analyte lists and analyteTypes must be the same length.")
  }
  if (!(enrichMethod %in% VALID_ENRICH_METHODS)) {
    stop(paste0(enrichMethod, " is not a valid enrichment method for multiomics.\nValid options are: ", paste(VALID_ENRICH_METHODS, collapse = ", ")))
  }
  if (length(enrichDatabase) > 1 || length(enrichDatabaseFile) > 1) {
    stop("Only one enrichDatabase or enrichDatabaseFile can be specified for multiomics.")
  }

  if (enrichMethod == "ORA") {
    cat("Performing multi-omics ORA\nLoading the functional categories...\n")
    databases <- c()
    if (!is.null(enrichDatabase)) { # Need to get correct name for metabolite databases
      if (length(unique(analyteTypes)) == 1) {
        databases <- enrichDatabase
      } else {
        for (i in seq_along(analyteTypes)) {
          databases <- c(databases, get_gmt_file(hostName, analyteTypes[i], enrichDatabase[i], organism, cache))
        }
      }
    } else {
      databases <- NULL
    }
    all_sets <- loadGeneSet(
      organism = organism, enrichDatabase = databases, enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType,
      enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile, cache = cache, hostName = hostName
    )
    if (length(all_sets) > 1) {
      geneSet <- all_sets[[1]]$geneSet
      geneSetDes <- all_sets[[1]]$geneSetDes
      geneSetNet <- all_sets[[1]]$geneSetNet
      for (i in 2:length(all_sets)) {
        geneSet <- rbind(geneSet, all_sets[[i]]$geneSet)
        geneSetDes <- rbind(geneSetDes, all_sets[[i]]$geneSetDes)
        geneSetDag <- rbind(geneSetDag, all_sets[[i]]$geneSetDag)
        geneSetNet <- rbind(geneSetNet, all_sets[[i]]$geneSetNet)
        databaseStandardId <- "multiomics"
      }
    } else {
      geneSet <- all_sets$geneSet
      geneSetDag <- all_sets$geneSetDag
      geneSetNet <- all_sets$geneSetNet
      databaseStandardId <- all_sets$standardId
    }

    rm(all_sets)

    cat("Loading the ID lists...\n")
    interest_lists <- list()
    if (is.null(analyteLists)) {
      for (i in seq_along(analyteListFiles)) {
        interestingGeneMap <- loadInterestGene(
          organism = organism, dataType = "list", inputGeneFile = analyteListFiles[i], inputGene = NULL,
          geneType = analyteTypes[i], collapseMethod = collapseMethod, cache = cache,
          hostName = hostName, geneSet = geneSet
        )
        if (organism == "others") {
          interestGeneList <- unique(interestingGeneMap)
          interest_lists[[i]] <- interestGeneList
        } else {
          interestStandardId <- interestingGeneMap$standardId
          interestGeneList <- unique(interestingGeneMap$mapped[[interestStandardId]])
          interest_lists[[i]] <- interestGeneList
        }
      }
    } else {
      for (i in seq_along(analyteLists)) {
        interestingGeneMap <- loadInterestGene(
          organism = organism, dataType = "list", inputGeneFile = NULL, inputGene = analyteLists[i],
          geneType = analyteTypes[i], collapseMethod = collapseMethod, cache = cache,
          hostName = hostName, geneSet = geneSet
        )
        if (organism == "others") {
          interestGeneList <- unique(interestingGeneMap)
          interest_lists[[i]] <- interestGeneList
        } else {
          interestStandardId <- interestingGeneMap$standardId
          interestGeneList <- unique(interestingGeneMap$mapped[[interestStandardId]])
          interest_lists[[i]] <- interestGeneList
        }
      }
    }

    # Load Gene Sets
    cat("Loading the reference lists...\n")

    referenceGeneList <- loadReferenceGene(organism = organism, referenceGeneFile = referenceGeneFile, referenceGene = referenceGene, referenceGeneType = referenceGeneType, referenceSet = referenceSet, collapseMethod = collapseMethod, hostName = hostName, geneSet = geneSet, interestGeneList = interestGeneList, cache = cache)

    ## Meta-analysis
  }
}


load_combined_gmt <- function(enrichDatabase, enrichDatabaseFile, enrichDatabaseDescriptionFile, enrichDatabaseType, analyteLists, analyteListFiles, analyteTypes, organism, cache, hostName) {
  databases <- c()
  if (!is.null(enrichDatabase)) { # Need to get correct name for metabolite databases
    if (length(unique(analyteTypes)) == 1) {
      databases <- enrichDatabase
    } else {
      types_processed <- c()
      for (i in seq_along(analyteTypes)) {
        if (analyteTypes[i] %in% types_processed) {
          next
        }
        databases <- c(databases, get_gmt_file(hostName, analyteTypes[i], enrichDatabase[i], organism, cache))
        types_processed <- c(types_processed, analyteTypes[i])
      }
      databases <- unique(databases)
    }
  } else {
    databases <- NULL
  }
  all_sets <- loadGeneSet(
    organism = organism, enrichDatabase = databases, enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType,
    enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile, cache = cache, hostName = hostName
  )
  return(all_sets)
}
