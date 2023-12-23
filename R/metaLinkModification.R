#' Modify the link to highlight the genes in the pathways
#'
#' Currently, we only have wikipathway and kegg pathways that need to modify the link
#'
#' @keywords internal
metaLinkModification <- function(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap_list, hostName = "https://www.webgestalt.org/", geneSet) {
    if (enrichMethod == "GSEA") {
        return(enrichPathwayLink)
    }
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
        enrichPathwayLink <- gsub("www.wikipathways.org/wpi/PathwayWidget.php?id=", "pathway-viewer.toolforge.org/embed/", enrichPathwayLink, fixed=TRUE)
        enrichPathwayLink <- paste0(enrichPathwayLink, "?")
    } else if (grepl("kegg.jp", enrichPathwayLink, fixed = TRUE)) {
            enrichPathwayLink <- gsub("show_pathway?", "show_pathway?map=", enrichPathwayLink, fixed=TRUE)
    }
    if (is_modified) {
        for (i in seq_along(all_sets)) {
            if (sets[[i]] == "rampc") { # metabolite
                kegg_colordbs <- list()
                all_genes <- unique(unlist(lapply(all_sets[[i]], function(x) interestingGeneMap_list[[x]]$mapped$rampc)))
                all_genes <- all_genes[all_genes %in% geneList]
                if (is.null(all_genes) || length(all_genes) == 0) {
                    next
                }
                for (j in seq_along(all_sets[[i]])) {
                    interestingGeneMap <- interestingGeneMap_list[[all_sets[[i]][[j]]]]
                    rampc_geneList <- interestingGeneMap$mapped$rampc[interestingGeneMap$mapped$rampc %in% all_genes]
                    interestingGeneMap <- interestingGeneMap_list[[all_sets[[i]][[j]]]]
                    if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                        mapped_genes <- simple_mapping(all_genes, "hsapiens", "rampc", "kegg", "rampc", hostName, no_dups = TRUE)
                        if (is.null(mapped_genes)) {
                            next
                        }
                        mapping_table <- data.frame(mapped_genes, all_genes, stringsAsFactors = FALSE)
                        kegg_geneList <- sapply(rampc_geneList, function(x) mapping_table[mapping_table$mapped_genes == x, "mapped_genes"])
                        set_addition <- meta_keggMetaboliteLinkModification(enrichMethod, kegg_geneList, rampc_geneList, all_genes, interestingGeneMap, hostName, j)
                        kegg_colordbs[[length(kegg_colordbs) + 1]] <- set_addition
                    } else if (grepl("toolforge.org", enrichPathwayLink, fixed = TRUE)) {
                        mapped_genes <- simple_mapping(all_genes, "hsapiens", "rampc", "hmdb", "rampc", hostName, no_dups = TRUE)
                        if (is.null(mapped_genes)) {
                            next
                        }
                        mapping_table <- data.frame(mapped_genes, all_genes, stringsAsFactors = FALSE)
                        hmdb_geneList <- sapply(rampc_geneList, function(x) mapping_table[mapping_table$mapped_genes == x, "mapped_genes"])
                        all_genes <- mapped_genes
                        set_addition <- meta_wikiMetaboliteLinkModification(enrichMethod, hmdb_geneList, rampc_geneList, all_genes, interestingGeneMap, hostName, j)
                        enrichPathwayLink <- paste0(enrichPathwayLink, set_addition, "&")
                    }
                }
                if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                    if (i == 1) {
                        enrichPathwayLink <- paste0(enrichPathwayLink, "&multi_query=")
                    }
                    for (j in seq_along(all_genes)) {
                        gene <- all_genes[[j]]
                        color_string <- ""
                        for (k in seq_along(kegg_colordbs)) {
                            color_db <- kegg_colordbs[[k]]
                            if (gene %in% color_db$analyte) {
                                color <- unique(unlist(color_db$color[color_db$analyte == gene]))
                                color <- color[1]
                                color <- gsub("#", "%23", color, fixed = TRUE)
                                color_string <- paste0(color_string, "+", color)
                            }
                        }
                        split_string <- "%0A"
                        if (j == 1) {
                            split_string <- ""
                        }
                        if (color_string != "") {
                            enrichPathwayLink <- paste0(enrichPathwayLink, gene, "+", color_string, split_string)
                        }
                    }
                }
            } else { # entrezgene
                kegg_colordbs <- list()
                all_genes <- unique(unlist(lapply(all_sets[[i]], function(x) interestingGeneMap_list[[x]]$mapped$entrezgene)))
                all_genes <- all_genes[all_genes %in% geneList]
                if (is.null(all_genes)) {
                    next;
                }
                for (j in seq_along(all_sets[[i]])) {
                    interestingGeneMap <- interestingGeneMap_list[[all_sets[[i]][[j]]]]
                    entrez_geneList <- interestingGeneMap$mapped$entrezgene[interestingGeneMap$mapped$entrezgene %in% all_genes]
                    if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                        set_addition <- meta_keggLinkModification(enrichMethod, entrez_geneList, all_genes, interestingGeneMap, hostName, j)
                         kegg_colordbs[[length(kegg_colordbs) + 1]] <- set_addition
                    } else if (grepl("toolforge.org", enrichPathwayLink, fixed = TRUE)) {
                        set_addition <- meta_wikiLinkModification(enrichMethod, entrez_geneList, all_genes, interestingGeneMap, hostName, j)
                        enrichPathwayLink <- paste0(enrichPathwayLink, set_addition, "&")
                    }
                    
                }
                if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                    if (i == 1) {
                        enrichPathwayLink <- paste0(enrichPathwayLink, "&multi_query=")
                    }
                    for (j in seq_along(all_genes)) {
                        gene <- all_genes[[j]]
                        color_string <- ""
                        for (k in seq_along(kegg_colordbs)) {
                            color_db <- kegg_colordbs[[k]]
                            if (gene %in% color_db$analyte) {
                                 color <- unique(unlist(color_db$color[color_db$analyte == gene]))
                                color <- color[1]
                                color <- gsub("#", "%23", color, fixed = TRUE)
                                color_string <- paste0(color_string, "+", color)
                            }
                        }
                        split_string <- "%0A"
                        if (j == 1) {
                            split_string <- ""
                        }
                        if (color_string != "") {
                            enrichPathwayLink <- paste0(enrichPathwayLink, gene, "+", color_string, split_string)
                        }
                    }
                }
            }
        }
    }
    return(enrichPathwayLink)
}


meta_keggMetaboliteLinkModification <- function(enrichMethod, geneList, rampc_geneList, all_genes, interestingGeneMap, hostName, color_index) {
     enrichPathwayLink <- ""
    found <- sapply(geneList, function(x) x <- gsub("kegg:", "", x, ignore.case = TRUE))
    not_found <- all_genes[!(all_genes %in% geneList)]
    not_found <- sapply(not_found, function(x) x <- gsub("kegg:", "", x, ignore.case = TRUE))
    color_db <- list(color = c(), analyte = c())
    if (enrichMethod == "ORA") {
        ora_color <- get_ora_colors(color_index)
        enrichPathwayLink <- paste0(ora_color, "=")
        for (i in seq_along(found)) {
            color_db$analyte[[i]] <- found[[i]]
            color_db$color[[i]] <- ora_color
        }
        ora_white <- get_white(color_index)
        for (i in seq_along(not_found)) {
            color_db$analyte[[i]] <- not_found[[i]]
            color_db$color[[i]] <- ora_white
        }
        color_db$analyte <- unlist(color_db$analyte)
        color_db$color <- unlist(color_db$color)
    } else if (enrichMethod == "GSEA") {
        scores <- filter(interestingGeneMap$mapped, .data$rampc %in% rampc_geneList)[["score"]]
        maxScore <- max(scores)
        minScore <- min(scores)
        tmp <- getPaletteForGsea(maxScore, minScore)
        palette <- tmp[[1]]
        palette <- shift_colors(palette, color_index)
        breaks <- tmp[[2]]
        colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
        for (i in seq_along(colors)) {
            colors[[i]] <- gsub("#", "%23", colors[[i]], fixed = TRUE)
        }
        unique_colors = unique(colors)
        for (i in seq_along(unique_colors)) {
            color <- unique_colors[[i]]
            all_colored_genes <- found[[colors == color]]
            all_colored_genes <- paste(all_colored_genes, collapse = ",")
            enrichPathwayLink <- paste0(enrichPathwayLink, "&", color, "=", all_colored_genes)
        }
        ora_white <- get_white(color_index)
        enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
        for (i in seq_along(not_found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, not_found[[i]], ",")
        }
    }
    return(as.data.frame(color_db))
}

meta_keggLinkModification <- function(enrichMethod, geneList, all_genes, interestingGeneMap, hostName, color_index) {
    enrichPathwayLink <- ""
    found <- geneList
    not_found <- all_genes[!(all_genes %in% geneList)]
    color_db <- list(color = c(), analyte = c())
    if (enrichMethod == "ORA") {
        ora_color <- get_ora_colors(color_index, reverse = TRUE)
        enrichPathwayLink <- paste0(ora_color, "=")
        for (i in seq_along(found)) {
            color_db$analyte[[i]] <- found[[i]]
            color_db$color[[i]] <- ora_color
        }
        ora_white <- get_white(color_index)
        for (i in seq_along(not_found)) {
            color_db$analyte[[i + length(found)]] <- not_found[[i]]
            color_db$color[[i + length(found)]] <- ora_white
        }
        color_db$analyte <- unlist(color_db$analyte)
        color_db$color <- unlist(color_db$color)
    } 
    # else if (enrichMethod == "GSEA") {
    #     # scores <- filter(interestingGeneMap$mapped, .data$entrezgene %in% rampc_geneList)[["score"]]
    #     maxScore <- max(scores)
    #     minScore <- min(scores)
    #     tmp <- getPaletteForGsea(maxScore, minScore)
    #     palette <- tmp[[1]]
    #     palette <- shift_colors(palette, color_index)
    #     breaks <- tmp[[2]]
    #     colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
    #     for (i in seq_along(colors)) {
    #         colors[[i]] <- gsub("#", "%23", colors[[i]], fixed = TRUE)
    #     }
    #     unique_colors = unique(colors)
    #     for (i in seq_along(unique_colors)) {
    #         color <- unique_colors[[i]]
    #         all_colored_genes <- found[[colors == color]]
    #         all_colored_genes <- paste(all_colored_genes, collapse = ",")
    #         enrichPathwayLink <- paste0(enrichPathwayLink, "&", color, "=", all_colored_genes)
    #     }
    #     ora_white <- get_white(color_index)
    #     enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
    #     for (i in seq_along(not_found)) {
    #         enrichPathwayLink <- paste0(enrichPathwayLink, not_found[[i]], ",")
    #     }
    # }
    return(color_db)
}

hmdb_map <- function(geneList, interestingGeneMap, hostName) {
    geneMap <- interestingGeneMap$mapped
    hmdbGeneList <- simple_mapping(unlist(strsplit(geneList, ";")), "hsapiens", "rampc", "hmdb", "rampc", hostName, no_dups = TRUE)
    hmdbGeneList <- sapply(hmdbGeneList, function(x) x <- gsub("hmdb:", "", x, ignore.case = TRUE))
    geneMap <- filter(geneMap, .data$rampc %in% geneList)
    return(geneMap)
}

meta_wikiMetaboliteLinkModification <- function(enrichMethod, geneList, rampc_geneList, all_genes, interestingGeneMap, hostName, color_index) {
    enrichPathwayLink <- ""
    found <- sapply(geneList, function(x) x <- gsub("hmdb:", "HMDB_", x, ignore.case = TRUE))
    not_found <- all_genes[!(all_genes %in% geneList)]
    not_found <- sapply(not_found, function(x) x <- gsub("hmdb:", "HMDB_", x, ignore.case = TRUE))
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
        scores <- filter(interestingGeneMap$mapped, .data$rampc %in% rampc_geneList)[["score"]]
        maxScore <- max(scores)
        minScore <- min(scores)
        tmp <- getPaletteForGsea(maxScore, minScore)
        palette <- tmp[[1]]
        palette <- shift_colors(palette, color_index)
        breaks <- tmp[[2]]
        colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
        for (i in seq_along(colors)) {
            colors[[i]] <- gsub("#", "%23", colors[[i]], fixed = TRUE)
        }
        unique_colors = unique(colors)
        for (i in seq_along(unique_colors)) {
            color <- unique_colors[[i]]
            all_colored_genes <- found[[colors == color]]
            all_colored_genes <- paste(all_colored_genes, collapse = ",")
            enrichPathwayLink <- paste0(enrichPathwayLink, "&", color, "=", all_colored_genes)
        }
        ora_white <- get_white(color_index)
        enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
        for (i in seq_along(not_found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, not_found[[i]], ",")
        }
    }
    
    return(enrichPathwayLink)
}

meta_wikiLinkModification <- function(enrichMethod, geneList, all_genes, interestingGeneMap, hostName, color_index) {
    enrichPathwayLink <- ""
    found <- geneList
    not_found <- all_genes[!(all_genes %in% geneList)]
    if (enrichMethod == "ORA") {
        ora_color <- get_ora_colors(color_index, reverse = TRUE)
        enrichPathwayLink <- paste0(ora_color, "=")
        for (i in seq_along(found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, "Entrez Gene_",found[[i]], ",")
        }
        enrichPathwayLink <- substr(enrichPathwayLink, 1, nchar(enrichPathwayLink) - 1)
        ora_white <- get_white(color_index)
        enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
        for (i in seq_along(not_found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, "Entrez Gene_", not_found[[i]], ",")
        }
    } else if (enrichMethod == "GSEA") {
        scores <- filter(interestingGeneMap$mapped, .data$entrezgene %in% geneList)[["score"]]
        maxScore <- max(scores)
        minScore <- min(scores)
        tmp <- getPaletteForGsea(maxScore, minScore)
        palette <- tmp[[1]]
        palette <- shift_colors(palette, color_index)
        breaks <- tmp[[2]]
        colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
        for (i in seq_along(colors)) {
            colors[[i]] <- gsub("#", "%23", colors[[i]], fixed = TRUE)
        }
        unique_colors = unique(colors)
        for (i in seq_along(unique_colors)) {
            color <- unique_colors[[i]]
            all_colored_genes <- found[[colors == color]]
            all_colored_genes <- paste(all_colored_genes, collapse = ",")
            enrichPathwayLink <- paste0(enrichPathwayLink, "&", color, "=", all_colored_genes)
        }
        ora_white <- get_white(color_index)
        enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
        for (i in seq_along(not_found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, not_found[[i]], ",")
        }
    }
    return(enrichPathwayLink)
}

get_ora_colors <- function(set_number, reverse = FALSE) {
    colors <- c("%23ffadad", "%23ffd6a5", "%23fdffb6", "%23caffbf", "%2396e3ff", "%23b5baff", "%23ffb3ff", "%23ffdfc4")
    if (reverse) {
        colors <- rev(colors)
    }
    if (set_number > length(colors)) {
        colors <- c(colors, rep(colors, ceiling(set_number / length(colors)) - 1))
    }
    return(colors[set_number])
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

shift_colors <- function(colors, shift) {
    colors <- gsub("#", "0x", colors, fixed = TRUE)
    colors <- paste0("#", sapply(colors, function(x) sprintf("%06x", strtoi(x) + shift)))
    return(colors)
}