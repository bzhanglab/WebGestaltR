get_gmt_file <- function(hostName, id_type, db_name, organism, cache) {
  std_id <- identifyStandardId(hostName, id_type, organism, "interest", cache)
  if (is.null(std_id)) {
    stop(idTypeError(id_type, organism, hostName, cache))
  }
  geneSetInfo <- listGeneSet(organism = organism, hostName = hostName, cache = cache)
  if (std_id == "rampc") {
    if (grepl("pathway", db_name, fixed = TRUE) && organism == "hsapiens") {
      if (grepl("RampDB_Metabolomics", db_name, fixed = TRUE)) {
        return(geneSetInfo[geneSetInfo$idType == "rampc" & geneSetInfo$name == db_name, 1][[1]])
      }
      db_name <- gsub("pathway_TAGMETABOLITE", "", db_name)
      db_name <- gsub("pathway_", "", db_name)
      return(geneSetInfo[geneSetInfo$idType == "rampc" & geneSetInfo$name == paste0("pathway_TAGMETABOLITE", db_name), 1][[1]])
    } else if (organism != "hsapiens") {
      stop("Only the organism 'hsapiens' is supported for metabolites.")
    } else {
      stop(paste0("The database ", db_name, " is not supported for metabolites."))
    }
  } else {
    return(geneSetInfo[geneSetInfo$idType == std_id & geneSetInfo$name == db_name, 1][[1]])
  }
}
