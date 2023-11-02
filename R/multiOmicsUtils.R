get_gmt_file <- function(hostName, id_type, db_name, organism, cache) {
  std_id <- identifyStandardId(hostName, id_type, organism, "interest", cache)
  if (is.null(std_id)) {
    stop(idTypeError(id_type, organism, hostName, cache))
  }
  geneSetInfo <- listGeneSet(organism = organism, hostName = hostName, cache = cache)
  if (std_id == "rampc") {
    db_name <- gsub("pathway_TAGMETABOLITE", "", db_name)
    db_name <- gsub("pathway_", "", db_name)
    return(geneSetInfo[geneSetInfo$idType == "rampc" & geneSetInfo$name == paste0("pathway_TAGMETABOLITE", db_name), 1][[1]])
  } else {
    return(geneSetInfo[geneSetInfo$idType == std_id & geneSetInfo$name == db_name, 1][[1]])
  }
}
