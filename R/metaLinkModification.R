#' Modify the link to highlight the genes in the pathways
#'
#' Currently, we only have wikipathway and kegg pathways that need to modify the link
#'
#' @keywords internal
metaLinkModification <- function(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap_list, hostName = "https://www.webgestalt.org/", geneSet) {
    geneList <- unlist(unique(unlist(geneList)))
    print(geneList)
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
    original_link <- enrichPathwayLink
    modified_links <- c("www.kegg.jp", "www.wikipathways.org")
    is_modified <- any(sapply(modified_links, function(x) any(grepl(x, enrichPathwayLink, fixed = TRUE))))
    pathway_id <- geneSet
    if (grepl("PathwayWidget", enrichPathwayLink, fixed = FALSE)) {
        enrichPathwayLink <- gsub("www.wikipathways.org/wpi/PathwayWidget.php?id=", "pathway-viewer.toolforge.org/embed/", enrichPathwayLink, fixed = TRUE)
        enrichPathwayLink <- paste0(enrichPathwayLink, "?")
    } else if (grepl("kegg.jp", enrichPathwayLink, fixed = TRUE)) {
        enrichPathwayLink <- gsub("show_pathway?", "show_pathway?map=", enrichPathwayLink, fixed = TRUE)
        reg_pat <- "(?<=\\?map=)([A-z]{3,5})(?=\\d)"
        enrichPathwayLink <- gsub(reg_pat, "map", enrichPathwayLink, perl = TRUE)
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
                if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                    mapping_table <- full_simple_mapping(all_genes, "hsapiens", "rampc", "kegg", "rampc", hostName, no_dups = TRUE)
                } else if (grepl("toolforge.org", enrichPathwayLink, fixed = TRUE)) {
                    mapping_table <- full_simple_mapping(all_genes, "hsapiens", "rampc", "hmdb", "rampc", hostName, no_dups = TRUE)
                }
                if (is.null(mapping_table)) {
                    next
                }
                colnames(mapping_table) <- c("all_genes", "mapped_genes")
                for (j in seq_along(all_sets[[i]])) {
                    interestingGeneMap <- interestingGeneMap_list[[all_sets[[i]][[j]]]]
                    rampc_geneList <- interestingGeneMap$mapped$rampc[interestingGeneMap$mapped$rampc %in% mapping_table$all_genes]
                    if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                        if (is.null(mapping_table$mapped_genes)) {
                            next
                        }
                        kegg_geneList <- sapply(rampc_geneList, function(x) mapping_table[mapping_table$all_genes == x, "mapped_genes"])
                        set_addition <- meta_keggMetaboliteLinkModification(enrichMethod, kegg_geneList, rampc_geneList, all_genes, interestingGeneMap, hostName, j, mapping_table)
                        kegg_colordbs[[length(kegg_colordbs) + 1]] <- set_addition
                    } else if (grepl("toolforge.org", enrichPathwayLink, fixed = TRUE)) {
                        hmdb_geneList <- sapply(rampc_geneList, function(x) mapping_table[mapping_table$all_genes == x, "mapped_genes"])
                        hmdb_all_genes <- mapping_table$mapped_genes
                        set_addition <- meta_wikiMetaboliteLinkModification(enrichMethod, hmdb_geneList, rampc_geneList, hmdb_all_genes, interestingGeneMap, hostName, j)
                        enrichPathwayLink <- paste0(enrichPathwayLink, set_addition, "&")
                    }
                }
                if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                    if (!grepl("&multi_query=", enrichPathwayLink, fixed = TRUE)) {
                        enrichPathwayLink <- paste0(enrichPathwayLink, "&multi_query=")
                    }
                    for (j in seq_along(mapping_table$mapped_genes)) {
                        gene <- mapping_table$mapped_genes[[j]]
                        gene <- gsub("kegg:", "", gene, ignore.case = TRUE)
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
                        split_string <- "%0d%0a"
                        if (color_string != "") {
                            enrichPathwayLink <- paste0(enrichPathwayLink, gene, color_string, split_string)
                        }
                    }
                }
            } else { # entrezgene
                kegg_colordbs <- list()
                all_genes_symbol <- unique(unlist(lapply(all_sets[[i]], function(x) interestingGeneMap_list[[x]]$mapped$geneSymbol)))
                all_genes <- unique(unlist(lapply(all_sets[[i]], function(x) interestingGeneMap_list[[x]]$mapped$entrezgene)))
                mapping_table <- data.frame(all_genes, all_genes_symbol, stringsAsFactors = FALSE)
                all_genes <- all_genes[all_genes %in% geneList]
                if (is.null(all_genes) || length(all_genes) == 0) {
                    next
                }
                if (grepl("kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                    kegg_ontology <- full_simple_mapping(all_genes, "hsapiens", "entrezgene", "kegg_ontology", "entrezgene", hostName, no_dups = TRUE)
                }
                all_genes_symbol <- mapping_table[mapping_table$all_genes %in% all_genes, "all_genes_symbol"]
                if (is.null(all_genes)) {
                    next
                }
                for (j in seq_along(all_sets[[i]])) {
                    interestingGeneMap <- interestingGeneMap_list[[all_sets[[i]][[j]]]]
                    genesymbol_geneList <- interestingGeneMap$mapped$geneSymbol[interestingGeneMap$mapped$geneSymbol %in% all_genes_symbol]
                    if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                        set_addition <- meta_keggLinkModification(enrichMethod, genesymbol_geneList, all_genes_symbol, interestingGeneMap, hostName, j)
                        kegg_colordbs[[length(kegg_colordbs) + 1]] <- set_addition
                    } else if (grepl("toolforge.org", enrichPathwayLink, fixed = TRUE)) {
                        set_addition <- meta_wikiLinkModification(enrichMethod, genesymbol_geneList, all_genes_symbol, interestingGeneMap, hostName, j)
                        enrichPathwayLink <- paste0(enrichPathwayLink, set_addition, "&")
                    }
                }
                if (grepl("www.kegg.jp", enrichPathwayLink, fixed = TRUE)) {
                    if (!is.null(kegg_ontology) && nrow(kegg_ontology) > 0) {
                        if (!grepl("&multi_query=", enrichPathwayLink, fixed = TRUE)) {
                            enrichPathwayLink <- paste0(enrichPathwayLink, "&multi_query=")
                        }
                        all_displayed_genes <- kegg_ontology$sourceId
                        all_genes_kegg <- kegg_ontology$targetId
                        for (j in seq_along(all_displayed_genes)) {
                            kegg_gene <- all_genes_kegg[[j]]
                            gene_symbol <- mapping_table[mapping_table$all_genes == all_displayed_genes[[j]], "all_genes_symbol"]
                            color_string <- ""
                            for (k in seq_along(kegg_colordbs)) {
                                color_db <- kegg_colordbs[[k]]
                                if (gene_symbol %in% color_db$analyte) {
                                    color <- unique(unlist(color_db$color[color_db$analyte == gene_symbol]))
                                    color <- color[1]
                                    color <- gsub("#", "%23", color, fixed = TRUE)
                                    color_string <- paste0(color_string, "+", color)
                                }
                            }
                            split_string <- "%0d%0a"
                            if (color_string != "") {
                                enrichPathwayLink <- paste0(enrichPathwayLink, kegg_gene, color_string, split_string)
                            }
                        }
                    }
                }
            }
        }
    }
    if (nchar(enrichPathwayLink) > 2000) { # URL length limit
        enrichPathwayLink <- paste0(hostName, "/long_url.html?pathway_url=", URLencode(original_link), "&pathway_id=", URLencode(pathway_id))
    }
    if (is.na(enrichPathwayLink) || is.null(enrichPathwayLink) || enrichPathwayLink == "") {
        enrichPathwayLink <- original_link
    }
    return(enrichPathwayLink)
}


meta_keggMetaboliteLinkModification <- function(enrichMethod, kegg_geneList, rampc_geneList, all_genes, interestingGeneMap, hostName, color_index, mapping_table) {
    found <- sapply(kegg_geneList, function(x) x <- gsub("kegg:", "", x, ignore.case = TRUE))
    rampc_not_found <- all_genes[!(all_genes %in% rampc_geneList)]
    not_found <- list()
    for (i in seq_along(rampc_not_found)) {
        not_found[[i]] <- mapping_table[mapping_table$all_genes == rampc_not_found[[i]], "mapped_genes"]
    }
    not_found <- sapply(not_found, function(x) x <- gsub("kegg:", "", x, ignore.case = TRUE))
    color_db <- list(color = c(), analyte = c())
    if (enrichMethod == "ORA") {
        ora_color <- get_ora_colors(color_index)
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
        scores <- list()
        for (i in seq_along(found)) {
            found_rampc <- rampc_geneList[[i]]
            scores[[i]] <- interestingGeneMap$mapped[interestingGeneMap$mapped$rampc == found_rampc, "score"]
        }
        scores <- unlist(scores)
        if (length(unlist(scores)) != 0) {
            maxScore <- max(unlist(scores))
            minScore <- min(unlist(scores))
            tmp <- getPaletteForGsea(maxScore, minScore)
            palette <- tmp[[1]]
            palette <- shift_colors(palette, color_index)
            breaks <- tmp[[2]]
            colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
            for (i in seq_along(colors)) {
                colors[[i]] <- gsub("#", "%23", colors[[i]], fixed = TRUE)
            }
            for (i in seq_along(found)) {
                color <- colors[[i]]
                color_db$analyte[[i]] <- found[[i]]
                color_db$color[[i]] <- color
            }

            ora_white <- get_white(color_index)
            for (i in seq_along(not_found)) {
                color_db$analyte[[i + length(color_db$analyte)]] <- not_found[[i]]
                color_db$color[[i + length(color_db$color)]] <- ora_white
            }
        }
        color_db$analyte <- unlist(color_db$analyte)
        color_db$color <- unlist(color_db$color)
    }
    return(color_db)
}

meta_keggLinkModification <- function(enrichMethod, geneList, all_genes, interestingGeneMap, hostName, color_index) {
    found <- geneList
    not_found <- all_genes[!(all_genes %in% geneList)]
    color_db <- list(color = c(), analyte = c())
    if (enrichMethod == "ORA") {
        ora_color <- get_ora_colors(color_index, reverse = TRUE)
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
    } else if (enrichMethod == "GSEA") {
        scores <- list()
        for (i in seq_along(found)) {
            scores[[i]] <- interestingGeneMap$mapped[interestingGeneMap$mapped$geneSymbol == found[[i]], "score"]
        }
        scores <- unlist(scores)
        if (length(unlist(scores)) != 0) {
            maxScore <- max(unlist(scores))
            minScore <- min(unlist(scores))
            tmp <- getPaletteForGsea(maxScore, minScore)
            palette <- tmp[[1]]
            palette <- shift_colors(palette, color_index)
            breaks <- tmp[[2]]
            colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
            for (i in seq_along(colors)) {
                colors[[i]] <- gsub("#", "%23", colors[[i]], fixed = TRUE)
            }
            for (i in seq_along(found)) {
                color <- colors[[i]]
                color_db$analyte[[i]] <- found[[i]]
                color_db$color[[i]] <- color
            }

            ora_white <- get_white(color_index)
            for (i in seq_along(not_found)) {
                color_db$analyte[[i + length(found)]] <- not_found[[i]]
                color_db$color[[i + length(found)]] <- ora_white
            }
            color_db$analyte <- unlist(color_db$analyte)
            color_db$color <- unlist(color_db$color)
        }
    }
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
        scores <- list()
        for (i in seq_along(found)) {
            scores[[i]] <- interestingGeneMap$mapped[interestingGeneMap$mapped$rampc == rampc_geneList[[i]], "score"]
        }
        scores <- unlist(scores)
        if (length(unlist(scores)) != 0) {
            maxScore <- max(scores)
            minScore <- min(scores)
            tmp <- getPaletteForGsea(maxScore, minScore, 100)
            palette <- tmp[[1]]
            palette <- shift_colors(palette, color_index)
            breaks <- tmp[[2]]
            colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
            for (i in seq_along(colors)) {
                colors[[i]] <- gsub("#", "%23", colors[[i]], fixed = TRUE)
            }
            unique_colors <- unique(colors)
            for (i in seq_along(unique_colors)) {
                color <- unique_colors[[i]]
                all_colored_genes <- c()
                for (j in seq_along(colors)) {
                    if (colors[[j]] == color) {
                        all_colored_genes <- c(all_colored_genes, geneList[[j]])
                    }
                }
                all_colored_genes <- paste(all_colored_genes, collapse = ",")
                enrichPathwayLink <- paste0(enrichPathwayLink, "&", color, "=", all_colored_genes)
            }
        }
        ora_white <- get_white(color_index)

        if (length(not_found) != 0) {
            enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
            for (i in seq_along(not_found)) {
                enrichPathwayLink <- paste0(enrichPathwayLink, not_found[[i]], ",")
            }
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
            enrichPathwayLink <- paste0(enrichPathwayLink, found[[i]], ",")
        }
        enrichPathwayLink <- substr(enrichPathwayLink, 1, nchar(enrichPathwayLink) - 1)
        ora_white <- get_white(color_index)
        enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
        for (i in seq_along(not_found)) {
            enrichPathwayLink <- paste0(enrichPathwayLink, not_found[[i]], ",")
        }
    } else if (enrichMethod == "GSEA") {
        scores <- list()
        for (i in seq_along(found)) {
            scores[[i]] <- interestingGeneMap$mapped[interestingGeneMap$mapped$geneSymbol == found[[i]], "score"]
        }
        scores <- unlist(scores)
        if (length(unlist(scores)) != 0) {
            maxScore <- max(scores)
            minScore <- min(scores)
            tmp <- getPaletteForGsea(maxScore, minScore, 100)
            palette <- tmp[[1]]
            palette <- shift_colors(palette, color_index)
            breaks <- tmp[[2]]
            colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
            for (i in seq_along(colors)) {
                colors[[i]] <- gsub("#", "%23", colors[[i]], fixed = TRUE)
            }
            unique_colors <- unique(colors)
            for (i in seq_along(unique_colors)) {
                color <- unique_colors[[i]]
                all_colored_genes <- c()
                for (j in seq_along(colors)) {
                    if (colors[[j]] == color) {
                        all_colored_genes <- c(all_colored_genes, geneList[[j]])
                    }
                }
                all_colored_genes <- paste(all_colored_genes, collapse = ",")
                enrichPathwayLink <- paste0(enrichPathwayLink, "&", color, "=", all_colored_genes)
            }
        }
        ora_white <- get_white(color_index)
        if (length(not_found) != 0) {
            enrichPathwayLink <- paste0(enrichPathwayLink, "&", ora_white, "=")
            for (i in seq_along(not_found)) {
                enrichPathwayLink <- paste0(enrichPathwayLink, not_found[[i]], ",")
            }
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
    # Add shift to blue part of hex code for each color in colors

    colors <- gsub("#", "", colors, fixed = TRUE) # Remove '#' from color codes
    colors <- strsplit(colors, "") # Split each color code into individual characters
    colors <- lapply(colors, function(x) {
        # Convert each color code to RGB
        r <- strtoi(paste0(x[1], x[2]), 16)
        g <- strtoi(paste0(x[3], x[4]), 16)
        b <- strtoi(paste0(x[5], x[6]), 16)
        # Add shift to blue part
        b <- b + shift
        # Ensure b is within 0-255
        b <- ifelse(b > 255, 255, ifelse(b < 0, 0, b))
        # Convert RGB back to hex
        return(paste0("#", sprintf("%02x", r), sprintf("%02x", g), sprintf("%02x", b)))
    })
    # Unlist to get back a vector of color codes
    colors <- unlist(colors)
    return(colors)
}


#     # colors <- gsub("#", "0x", colors, fixed = TRUE)
#     # colors <- paste0("#", sapply(colors, function(x) sprintf("%06x", strtoi(x) + shift)))
#     # return(colors)
# }

full_simple_mapping <- function(id_list, organism, source_id, target_id, standard_id, hostName, no_dups = FALSE) {
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
        return(NULL)
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
            return(mappedInputGene)
        },
        error = function(e) {
            warning("No mapping result found. May be caused by empty sets chosen by significance method.")
            return(NULL)
        }
    )
}
