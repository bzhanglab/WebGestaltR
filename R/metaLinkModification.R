#' Modify the link to highlight the genes in the pathways
#'
#' Currently, we only have wikipathway and kegg pathways that need to modify the link
#'
#' @keywords internal
metaLinkModification <- function(enrichMethod, enrichPathwayLink, geneList_list, interestingGeneMap_list, hostName = "https://www.webgestalt.org/") {
    list_standard_ids <- sapply(interestingGeneMap_list, function(x) x$standardId)
    sets <- unique(list_standard_ids)
    all_sets <- list()
    for (i in seq_along(sets)) {
        all_sets[[i]] <- list()
        for (j in seq_along(list_standard_ids)) {
            if (sets[[i]] == list_standard_ids[[j]]) {
                all_sets[[i]] <- c(all_sets[[i]], j)
            }
        }
    }
    all_sets <- lapply(all_sets, unlist)
    is_modified <- FALSE
    modified_links <- c("www.kegg.jp", "www.wikipathways.org")
    is_modified <- any(sapply(modified_links, function(x) any(grepl(x, enrichPathwayLink, fixed = TRUE))))
    if (grepl("PathwayWidget", enrichPathwayLink, fixed = FALSE)) {
        print("here")
        enrichPathwayLink <- gsub("www.wikipathways.org/wpi/PathwayWidget.php?id=", "pathway-viewer.toolforge.org/embed/", enrichPathwayLink, fixed=TRUE) 
    }
    if (is_modified) {
        enrichPathwayLink <- paste0(enrichPathwayLink, "?")
        for (i in seq_along(all_sets)) {
            if (sets[[i]] == "rampc") { # metabolite
                all_genes <- unique(unlist(lapply(all_sets[[i]], function(x) interestingGeneMap_list[[x]]$mapped$rampc)))
                mapped_genes <- simple_mapping(all_genes, "hsapiens", "rampc", "hmdb", "rampc", hostName, no_dups = TRUE)
                mapping_table <- data.frame(mapped_genes, all_genes, stringsAsFactors = FALSE)
                all_genes <- mapped_genes
                for (j in seq_along(all_sets[[i]])) {
                    print(head(geneList_list[[all_sets[[i]][[j]]]]))
                    hmdb_geneList <- mapping_table[mapping_table$all_genes %in% geneList, ]$mapped_genes
                    interestingGeneMap <- interestingGeneMap_list[[all_sets[[i]][[j]]]]
                    if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                        set_addition <- meta_keggMetaboliteLinkModification(enrichPathwayLink, hmdb_geneList, interestingGeneMap, hostName)
                    } else if (grepl("www.wikipathways.org", enrichPathwayLink, fixed = TRUE)) {
                        set_addition <- meta_wikiMetaboliteLinkModification(enrichMethod, hmdb_geneList, all_genes, interestingGeneMap, hostName, j)
                    }
                    enrichPathwayLink <- paste0(enrichPathwayLink, set_addition, "&")
                }
            } else { # entrezgene

                all_genes <- unique(unlist(lapply(all_sets[[i]], function(x) interestingGeneMap_list[[x]]$mapped$standardId)))
                # mapped_genes <- simple_mapping(all_genes, "hsapiens", "rampc", "hmdb", "rampc", hostName, no_dups = TRUE)
                # mapping_table <- data.frame(mapped_genes, all_genes, stringsAsFactors = FALSE)
                # all_genes <- mapped_genes
                for (j in seq_along(all_sets[[i]])) {
                    geneList <- geneList_list[[all_sets[[i]][[j]]]]
                    print(geneList)
                    # hmdb_geneList <- mapping_table[mapping_table$all_genes %in% geneList, ]$mapped_genes
                    interestingGeneMap <- interestingGeneMap_list[[all_sets[[i]][[j]]]]
                    if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                        set_addition <- meta_keggLinkModification(enrichPathwayLink, geneList, interestingGeneMap, hostName)
                    } else if (grepl("www.wikipathways.org", enrichPathwayLink, fixed = TRUE)) {
                        set_addition <- meta_wikiLinkModification(enrichMethod, geneList, all_genes, interestingGeneMap, hostName, j)
                    }
                    enrichPathwayLink <- paste0(enrichPathwayLink, set_addition, "&")
                }
            }
        }
    }
    return(enrichPathwayLink)
}

meta_keggMetaboliteLinkModification <- function(enrichPathwayLink, geneList, interestingGeneMap, hostName) {
    geneList <- simple_mapping(unlist(strsplit(geneList, ";")), "hsapiens", "rampc", "kegg", "rampc", hostName)
    geneList <- sapply(geneList, function(x) x <- gsub("kegg:", "", x, ignore.case = TRUE))
    geneList <- paste(geneList, collapse = "+")
    enrichPathwayLink <- paste(enrichPathwayLink, "+", geneList, sep = "")
    return(enrichPathwayLink)
}

hmdb_map <- function(geneList, interestingGeneMap, hostName) {
    geneMap <- interestingGeneMap$mapped
    hmdbGeneList <- simple_mapping(unlist(strsplit(geneList, ";")), "hsapiens", "rampc", "hmdb", "rampc", hostName, no_dups = TRUE)
    hmdbGeneList <- sapply(hmdbGeneList, function(x) x <- gsub("hmdb:", "", x, ignore.case = TRUE))
    geneMap <- filter(geneMap, .data$rampc %in% geneList)
    return(geneMap)
}

meta_wikiMetaboliteLinkModification <- function(enrichMethod, geneList, all_genes, interestingGeneMap, hostName, color_index) {
    enrichPathwayLink <- ""
    found <- sapply(geneList, function(x) x <- gsub("hmdb:", "", x, ignore.case = TRUE))
    not_found <- all_genes[!(geneList %in% all_genes)]
    not_found <- sapply(not_found, function(x) x <- gsub("hmdb:", "", x, ignore.case = TRUE))
    if (enrichMethod == "ORA") {
        ora_color <- get_ora_colors(color_index)
        enrichPathwayLink <- paste0(ora_color, "=")
        for (i in seq_along(found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, "HMDB_", found[[i]], ",")
        }
        enrichPathwayLink <- substr(enrichPathwayLink, 1, nchar(enrichPathwayLink) - 1)
        ora_white <- get_white(color_index)
        enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
        for (i in seq_along(not_found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, "HMDB_", not_found[[i]], ",")
        }
    } else if (enrichMethod == "GSEA") {
        scores <- filter(interestingGeneMap$mapped, .data$rampc %in% geneList)[["score"]]
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

meta_wikiLinkModification <- function(enrichMethod, geneList, all_genes, interestingGeneMap, hostName, color_index) {
    enrichPathwayLink <- ""
    found <- geneList
    not_found <- all_genes[!(geneList %in% all_genes)]
    if (enrichMethod == "ORA") {
        ora_color <- get_ora_colors(color_index)
        enrichPathwayLink <- paste0(ora_color, "=")
        for (i in seq_along(found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, found[[i]], ",")
        }
        enrichPathwayLink <- substr(enrichPathwayLink, 1, nchar(enrichPathwayLink) - 1)
        ora_white <- get_white(color_index)
        enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
        for (i in seq_along(not_found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, not_found[[i]], ",")
        }
    } else if (enrichMethod == "GSEA") {
        scores <- filter(interestingGeneMap$mapped, .data$rampc %in% geneList)[["score"]]
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

get_ora_colors <- function(set_number) {
    colors <- c("%23ffadad", "%23ffd6a5", "%23fdffb6", "%23caffbf", "%2396e3ff", "%23b5baff", "%23ffb3ff", "%23ffdfc4", "%23fffffc")
    if (set_number > length(colors)) {
        colors <- c(colors, rep(colors, ceiling(set_number / length(colors)) - 1))
    }
    return(colors[1:set_number])
}

get_white <- function(offset) {
    blue_offset <- ceiling(offset / 3)
    green_offset <- floor(offset / 3) + (offset %% 3 == 2)
    red_offset <- floor(offset / 3)
    red_hex <- sprintf("%02x", 255 - red_offset)
    blue_hex <- sprintf("%02x", 255 - blue_offset)
    green_hex <- sprintf("%02x", 255 - green_offset)
    return(paste0("%23", red_hex, green_hex, blue_hex))
}
