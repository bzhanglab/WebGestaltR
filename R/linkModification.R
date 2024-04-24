#' Modify the link to highlight the genes in the pathways
#'
#' Currently, we only have wikipathway and kegg pathways that need to modify the link
#'
#' @keywords internal
linkModification <- function(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap, hostName = "https://www.webgestalt.org/") {
    if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE) && (interestingGeneMap$standardId == "rampc")) {
        link <- keggMetaboliteLinkModification(enrichPathwayLink, geneList, interestingGeneMap, hostName)
        return(link)
    } else if (grepl("www.wikipathways.org", enrichPathwayLink, fixed = TRUE) && (interestingGeneMap$standardId == "rampc")) {
        link <- wikiMetaboliteLinkModification(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap, hostName)
        return(link)
    } else if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
        link <- keggLinkModification(enrichPathwayLink, geneList)
        return(link)
    } else if (grepl("www.wikipathways.org", enrichPathwayLink, fixed = TRUE)) {
        link <- wikiLinkModification(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap)
        return(link)
    }
    return(enrichPathwayLink)
}

keggLinkModification <- function(enrichPathwayLink, geneList) {
    geneList <- gsub(";", "+", geneList)
    enrichPathwayLink <- paste(enrichPathwayLink, "+", geneList, sep = "")
    return(enrichPathwayLink)
}

keggMetaboliteLinkModification <- function(enrichPathwayLink, geneList, interestingGeneMap, hostName) {
    geneList <- simple_mapping(unlist(strsplit(geneList, ";")), "hsapiens", "rampc", "kegg", "rampc", hostName)
    geneList <- sapply(geneList, function(x) x <- gsub("kegg:", "", x, ignore.case = TRUE))
    geneList <- paste(geneList, collapse = "+")
    enrichPathwayLink <- paste(enrichPathwayLink, "+", geneList, sep = "")
    return(enrichPathwayLink)
}

wikiMetaboliteLinkModification <- function(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap, hostName) {
    geneMap <- interestingGeneMap$mapped
    hmdbGeneList <- simple_mapping(unlist(strsplit(geneList, ";")), "hsapiens", "rampc", "hmdb", "rampc", hostName, no_dups = TRUE)
    # hmdbGeneList <- sapply(hmdbGeneList, function(x) x <- gsub("hmdb:", "HMDB_", x, ignore.case = TRUE))
    geneMap <- filter(geneMap, .data$rampc %in% geneList)
    geneList <- unlist(strsplit(geneList, ";"))
    if (grepl("PathwayWidget", enrichPathwayLink, fixed = FALSE)) {
        enrichPathwayLink <- gsub("www.wikipathways.org/wpi/PathwayWidget.php?id=", "pathway-viewer.toolforge.org/embed/", enrichPathwayLink, fixed = TRUE)
        enrichPathwayLink <- paste0(enrichPathwayLink, "?")
    }
    if (enrichMethod == "ORA") {
        enrichPathwayLink <- paste0(enrichPathwayLink, colorPos, "=", paste(hmdbGeneList, collapse = ",", sep = ""))
    } else if (enrichMethod == "GSEA") {
        scores <- filter(interestingGeneMap$mapped, .data$rampc %in% geneList)[["score"]]
        if (length(unlist(scores)) == 0) {
            return(enrichPathwayLink)
        }
        maxScore <- max(scores)
        minScore <- min(scores)
        tmp <- getPaletteForGsea(maxScore, minScore)
        palette <- tmp[[1]]
        breaks <- tmp[[2]]
        colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
        colorStr <- paste(gsub("#", "%23", colors, fixed = TRUE), collapse = ",")
        for (i in seq_along(hmdbGeneList)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, URLencode(colors[i]), "=", hmdbGeneList[i], "&")
        }
        enrichPathwayLink <- paste0(enrichPathwayLink, "&colors=", colorStr)
    }
    return(enrichPathwayLink)
}

wikiLinkModification <- function(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap) {
    geneMap <- interestingGeneMap$mapped
    # print(geneMap)
    geneList <- unlist(strsplit(geneList, ";"))
    geneMap <- filter(geneMap, .data$entrezgene %in% geneList)
    enrichPathwayLink <- paste0(
        enrichPathwayLink,
        paste0(sapply(geneMap$geneSymbol, function(x) paste0("&label[]=", x)), collapse = "")
        # not many pathway have entrezgene xref. Using both also seem to interfere with coloring
        # paste0(sapply(geneMap$entrezgene, function(x) paste0("&xref[]=", x, ",Entrez Gene")), collapse="")
    )
    if (enrichMethod == "ORA") {
        enrichPathwayLink <- paste0(enrichPathwayLink, "&colors=", colorPos)
    } else if (enrichMethod == "GSEA") {
        scores <- filter(interestingGeneMap$mapped, .data$entrezgene %in% geneList)[["score"]]
        if (length(unlist(scores)) == 0) {
            return(enrichPathwayLink)
        }
        maxScore <- max(scores)
        minScore <- min(scores)
        tmp <- getPaletteForGsea(maxScore, minScore)
        palette <- tmp[[1]]
        breaks <- tmp[[2]]
        colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
        colorStr <- paste(gsub("#", "%23", colors, fixed = TRUE), collapse = ",")
        enrichPathwayLink <- paste0(enrichPathwayLink, "&colors=", colorStr)
    }
    return(enrichPathwayLink)
}


simple_mapping <- function(id_list, organism, source_id, target_id, standard_id, hostName, no_dups = FALSE) {
    if (source_id == target_id) {
        return(id_list)
    }
    response <- POST(file.path(hostName, "api", "idmapping"),
        encode = "json",
        body = list(
            organism = organism, sourceType = source_id,
            targetType = target_id, ids = id_list, standardId = standard_id
        )
    )
    if (response$status_code != 200) {
        stop(webRequestError(response))
    }
    mapRes <- content(response)
    mappedIds <- mapRes$mapped
    names <- c("sourceId", "targetId")
    if (is.null(mappedIds)) {
        return(c(""))
    }
    tryCatch(
        {
            mappedInputGene <- data.frame(matrix(unlist(lapply(replace_null(mappedIds), FUN = function(x) {
                x[names]
            })), nrow = length(mappedIds), byrow = TRUE), stringsAsFactors = FALSE)

            colnames(mappedInputGene) <- c("sourceId", "targetId")
            if (no_dups) {
                mappedInputGene <- mappedInputGene[!duplicated(mappedInputGene$sourceId), ]
                mappedInputGene <- mappedInputGene[!duplicated(mappedInputGene$targetId), ]
            }
            return(mappedInputGene$targetId)
        },
        error = function(e) {
            warning("No mapping result found. May be caused by empty sets chosen by significance method.")
            return(c(""))
        }
    )
}
