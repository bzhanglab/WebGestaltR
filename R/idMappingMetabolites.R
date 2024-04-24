#' @importFrom httr POST content
#' @importFrom dplyr right_join select left_join %>%
idMappingMetabolites <- function(organism = "hsapiens", dataType = "list", inputGeneFile = NULL, inputGene = NULL, sourceIdType, targetIdType, standardId, collapseMethod = "mean", mappingOutput = FALSE, outputFileName = "", hostName = "https://www.webgestalt.org/") {
  ########### Check input data type###############
  inputGene <- idMappingInput(dataType = dataType, inputGeneFile = inputGeneFile, inputGene = inputGene)

  ########## ID Mapping Specify to phosphosite level###############
  if (dataType == "list") {
    inputGeneL <- unique(inputGene)
  }

  if (dataType == "rnk") {
    ###### Collapse the gene ids with multiple scores##########
    x <- tapply(inputGene$score, inputGene$gene, collapseMethod)
    inputGene <- data.frame(gene = names(x), score = as.numeric(x), stringsAsFactors = FALSE)
    inputGeneL <- inputGene$gene
    colnames(inputGene) <- c(sourceIdType, "score")
  }

  if (dataType == "gmt") {
    colnames(inputGene) <- c("geneSet", "link", sourceIdType)
    inputGeneL <- unique(inputGene$gene)
  }
  mappedInputGene <- NULL
  if (startsWith(hostName, "file://")) {
    sourceMap <- read_tsv(
      removeFileProtocol(file.path(hostName, "xref", paste(organism, sourceIdType, paste0(standardId, ".table"), sep = "_"))),
      col_names = c(standardId, "userId"), col_types = "cc", quote = ""
    ) %>% filter(.data$userId %in% inputGeneL)
    if (targetIdType == standardId || targetIdType == sourceIdType) {
      mappedInputGene <- sourceMap
    } else {
      targetMap <- read_tsv(
        removeFileProtocol(file.path(hostName, "xref", paste(organism, targetIdType, paste0(standardId, ".table"), sep = "_"))),
        col_names = c(standardId, targetIdType), col_types = "cc", quote = ""
      )
      mappedInputGene <- inner_join(sourceMap, targetMap, by = c(standardId))
    }
    if (nrow(mappedInputGene) == 0) {
      return(idMappingError("empty"))
    }
    mappedInputGene <- select(mappedInputGene, .data$userId, targetIdType)
    unmappedIds <- setdiff(inputGeneL, mappedInputGene$userId)
  } else {
    response <- POST(file.path(hostName, "api", "idmapping"),
      encode = "json",
      body = list(
        organism = organism, sourceType = sourceIdType,
        targetType = targetIdType, ids = inputGeneL, standardId = standardId
      )
    )
    if (response$status_code != 200) {
      stop(webRequestError(response))
    }
    mapRes <- content(response)
    if (mapRes$status == 1) {
      stop(webApiError(mapRes))
    }

    mappedIds <- mapRes$mapped

    unmappedIds <- unlist(mapRes$unmapped)

    if (length(mappedIds) == 0) {
      stop(idMappingError("empty"))
    }
    names <- c("sourceId", "targetId")
    mappedInputGene <- data.frame(matrix(unlist(lapply(replace_null(mappedIds), FUN = function(x) {
      x[names]
    })), nrow = length(mappedIds), byrow = TRUE), stringsAsFactors = FALSE)

    colnames(mappedInputGene) <- c("userId", "rampc")
    unmappedIds <- append(unmappedIds, mappedInputGene[duplicated(mappedInputGene$rampc), ]$userId)
    mappedInputGene <- mappedInputGene[!duplicated(mappedInputGene$rampc), ]
    response <- POST(file.path(hostName, "api", "idmapping"),
      encode = "json",
      body = list(
        organism = organism, sourceType = "rampc",
        targetType = "metabolite_name", ids = mappedInputGene$rampc, standardId = standardId
      )
    )
    mapRes <- content(response)
    if (mapRes$status == 1) {
      stop(webApiError(mapRes))
    }
    newMapped <- mapRes$mapped
    meta_names <- data.frame(matrix(unlist(lapply(replace_null(newMapped), FUN = function(x) {
      x[names]
    })), nrow = length(newMapped), byrow = TRUE), stringsAsFactors = FALSE)
    colnames(meta_names) <- c("id", "meta")
  }
  mappedInputGene$geneSymbol <- mappedInputGene$userId
  mappedInputGene$geneName <- meta_names$meta

  # Link

  old_id_type <- sourceIdType
  sourceIdType <- tolower(sourceIdType)
  if (sourceIdType == "hmdb") {
    mappedInputGene$gLink <- paste0("https://www.hmdb.ca/metabolites/", replace_prefix(mappedInputGene$userId, "hmdb:"))
  } else if (sourceIdType == "chebi") {
    mappedInputGene$gLink <- paste0("https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:", replace_prefix(mappedInputGene$userId, "chebi:"))
  } else if (sourceIdType == "cas") {
    mappedInputGene$gLink <- paste0("https://commonchemistry.cas.org/detail?cas_rn=", replace_prefix(mappedInputGene$userId, "cas:"))
  } else if (sourceIdType == "chemspider") {
    mappedInputGene$gLink <- paste0("https://www.chemspider.com/Chemical-Structure.", replace_prefix(mappedInputGene$userId, "chemspider:"), ".html")
  } else if (sourceIdType == "kegg") {
    mappedInputGene$gLink <- paste0("https://www.genome.jp/entry/", replace_prefix(mappedInputGene$userId, "kegg:"))
  } else if (sourceIdType == "kegg_glycan") {
    mappedInputGene$gLink <- paste0("https://www.genome.jp/entry/", replace_prefix(mappedInputGene$userId, "kegg_glycan:"))
  } else if (sourceIdType == "lipidbank") {
    cuts <- replace_prefix(mappedInputGene$userId, "lipidbank:")
    mappedInputGene$gLink <- paste0("https://lipidbank.jp/", sapply(cuts, function(x) {
      return(substring(x, 1, 3))
    }))
  } else if (sourceIdType == "lipidmaps") {
    mappedInputGene$gLink <- paste0("https://www.lipidmaps.org/databases/lmsd/", replace_prefix(mappedInputGene$userId, "lipidmaps:"))
  } else if (sourceIdType == "plantfa") {
    mappedInputGene$gLink <- paste0("https://plantfadb.org/fatty_acids/", replace_prefix(mappedInputGene$userId, "plantfa:"))
  } else if (sourceIdType == "pubchem") {
    mappedInputGene$gLink <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", replace_prefix(mappedInputGene$userId, "pubchem:"))
  } else if (sourceIdType == "swisslipids") {
    mappedInputGene$gLink <- paste0("http://www.swisslipids.org/#/entity/", replace_prefix(mappedInputGene$userId, "swisslipids:"), "/")
  } else if (sourceIdType == "wikidata") {
    mappedInputGene$gLink <- paste0("https://www.wikidata.org/wiki/", replace_prefix(mappedInputGene$userId, "wikidata:"))
  } else if (sourceIdType == "rampc") {
    mappedInputGene$gLink <- replicate(length(mappedInputGene$userId), "https://rampdb.nih.gov") # No links for rampc. Links to home page
  } else {
    mappedInputGene$gLink <- paste0("URL NOT FOUND FOR TYPE ", sourceIdType)
  }

  if (dataType == "list") {
    inputGene <- select(mappedInputGene, .data$userId, .data$geneSymbol, .data$geneName, targetIdType, .data$gLink)
  }

  if (dataType == "rnk") {
    old_ids <- inputGene[[sourceIdType]]
    inputGene <- data.frame("score" = inputGene$score, "userId" = add_prefix(old_ids, old_id_type))
    inputGene <- inner_join(mappedInputGene, inputGene, by = "userId")
  }

  if (dataType == "gmt") {
    old_ids <- inputGene[[sourceIdType]]
    inputGene <- data.frame("score" = inputGene$score, "userId" = add_prefix(old_ids, old_id_type))
    inputGene <- mappedInputGene %>%
      inner_join(inputGene, by = "userId") %>%
      select(.data$geneSet, .data$link, .data$userId, .data$geneSymbol, .data$geneName, targetIdType, .data$gLink)
  }
  # inputGene <- mappedInputGene
  ############# Output#######################
  if (mappingOutput) {
    idMappingOutput(outputFileName, inputGene, unmappedIds, dataType, old_id_type, targetIdType = targetIdType)
  }
  r <- list(mapped = inputGene, unmapped = unmappedIds)
  return(r)
}

add_prefix <- function(x, sourceIdType) {
  no_prefixes <- c("rampc")
  if (sourceIdType %in% no_prefixes) {
    return(x)
  }
  uppers <- c("LIPIDMAPS", "CAS")
  if (toupper(sourceIdType) %in% uppers) {
    return(unlist(sapply(x, function(y) {
      if (grepl(":", y, fixed = TRUE)) {
        return(y)
      }
      return(paste0(toupper(sourceIdType), ":", toupper(y)))
    })))
  } else {
    return(unlist(sapply(x, function(y) {
      if (grepl(":", y, fixed = TRUE) && sourceIdType != "swisslipids") {
        return(y)
      }
      return(paste0(sourceIdType, ":", y))
    })))
  }
}

replace_prefix <- function(x, prefix) {
  return(sapply(x, function(y) {
    return(sub(prefix, "", y, ignore.case = TRUE))
  }))
}

.combineG <- function(e) {
  e <- e[-length(e)]
  e <- paste(e, collapse = "_")
  return(e)
}
