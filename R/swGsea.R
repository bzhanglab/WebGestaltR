#' Site Weighted Gene Set Enrichment Analysis
#'
#' Performs site weighted gene set enrichment analysis or standard GSEA when
#' likelihood/weight columns in \code{input_df} are 1 or 0, \code{p=1},
#' \code{q=1} and \code{thresh_type="val"}.
#'
#' The formula for weighting is as follows
#' \deqn{\frac{s_{j}^{q}|r_{j}|^{p}}{\sum s^{q}|r|^{p}}}{(s_j^q * abs(r_j)^p) / (\sum s^q * abs(r)^p)}
#' Where r is log ratio score, s is likelihood score, j is the index of the gene.
#'
#' @param input_df A data frame in which first column is name of item of interest
#' (gene, protein, phosphosite, etc.), the second is the correlation of that item
#' of interest with the phenotype (typically log ratio of expression for phenotype
#' vs. normal), and the remaining columns are the scores for the likelihood that
#' the item belongs in each set (one column per set).
#' @param thresh_type The type of \code{thresh}. Use 'percentile' to include all
#' scores over that percentile given in \code{thresh} (i.e., 0.9 would be all items
#' in 90th percentile, or top 10 percent); 'list' to include a list of set lists
#' where the set lists are in the same order as the corresponding set columns in
#' the \code{input_df}; 'val' to apply a single threshold value to all sets; or
#' 'values' to use a vector of unique cutoffs for each set (needs to be in the
#' same order as the sets are specified in the columns of \code{input_df}")
#' @param thresh Depends on \code{thresh_type}. A list of lists of the items in
#' each set (with same names as colnames of the scores); a numeric vector of
#' threshold scores for each set (in the same order as the colnames of the scores
#' in the input_df), or a single percentile value between 0 and 1 (i.e., if
#' \code{thresh}=0.9, the 90th percentile of the score or the highest scoring 10%
#' of of the items are included in the set for each scoring regimen) (\code{thresh}
#' ="all" is not supported at this time, as it doesn't result in a Kolgorov-Smirnoff
#' statistic; this may be worked in as an alternate scoring method later on).
#' @param thresh_action Either "include", "exclude (default)", or "adjust";
#' this specifies how to treat each set if it doesn't contain a minimum number of
#' items or contains all of the items; this option cannot be used with predefined
#' lists of items in sets (if the number of items in a given set doesn't meet
#' requirements, that set will be skipped).
#' @param min_set_size,max_set_size The minimum/maximum number of items each
#' set needs for the analysis to proceed.
#' @param max_score,min_score A optional numeric vector of minimum/maximum boundaries
#' to clip scores for each set.
#' @param psuedocount Psuedocount (pc) is used for rescaling set scores:
#' \code{(score - min_score + pc)/(max_score - min_score +pc)}; this is needed to
#' prevent division by 0 if \code{max_score==min_score} (in this case, all scores
#' for items in set will be 1, which is equivalent to standard GSEA); it also allows
#' users to adjust weights for scores that are close to the minimum for the scores in
#' the set (unless min_score==max_score): as psuedocount value approaches 0, scaled
#' minimum scores also approach 0; as psuedocount approaches infinity, scaled minimum
#' scores approach the scaled maximum scores (which equal 1); this value must be
#' larger than 0.
#' @param perms The number of permutations.
#' @param p The exponential scaling factor of the phenotype score (second column in
#' \code{input_df}).
#' @param q The exponential scaling factor of the likelihood score (weights).
#' @param nThreads The number of threads to use in calculating permutaions.
#' @param rng_seed Random seed.
#' @param fork A boolean. Whether pass "fork" to \code{type} parameter of
#' \code{makeCluster} on Unix-like machines.
#'
#' @return A list of \code{Enrichment_Results}, \code{Items_in_Set} and \code{Running_Sums}.
#' \describe{
#'  \item{Enrichment_Results}{A data frame with row names of gene set and columns of
#'  "ES", "NES", "p_val", "fdr".}
#'  \item{Items_in_Set}{A list of one-column data frames. Describes genes and their
#'  ranks in each set.}
#'  \item{Running_Sums}{Running sum scores along genes sorted by ranked scores,
#'  with gene sets as columns.}
#' }
#'
#' @importFrom stats quantile
#' @importFrom dplyr arrange desc
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
#' @export
#' @author Eric Jaehnig
#'
swGsea <- function(input_df, thresh_type = "percentile", thresh = 0.9, thresh_action = "exclude", min_set_size = 10, max_set_size = 500, max_score = "max", min_score = "min", psuedocount = 0.001, perms = 1000, p = 1, q = 1, nThreads = 1, rng_seed = 1, fork = FALSE) {
    # check input parameters
    if (thresh_type != "percentile" & thresh_type != "list" & thresh_type != "val" & thresh_type != "values") {
        stop("invalid thresh_type specified; needs to be set to 'percentile' to include all scores over that percentile (i.e., 0.9 would be all items in 90th percentile, or top 10 percent) 'list' to include a list of set lists where the set lists are in the same order as the corresponding set columns in the input_df, 'val' to apply a single threshold value to all sets, or 'values' to use a vector of unique cutoffs for each set (needs to be in the same order as the sets are specified in the columns of input_df")
    }
    if (thresh_action != "exclude" & thresh_action != "include" & thresh_action != "adjust") {
        stop("invalid thresh_action specified; needs to be set to 'exclude' to skip set if it contains no items after applying score threshold (or contains all items), or 'include' to include values for the set at the end of the results (ES and NES automatically set to 0 and pval to 1) or 'adjust' to adjust threshold to add at least min_set_size items below thresh to set (or remove all items equal to the minimum set score value from the set)")
    }
    if (min_set_size < 3) {
        stop("please set 'min_set_size' to 3 or greater (default is 5)")
    }
    if (max_score != "max") {
        if (length(max_score) != (ncol(input_df) - 2) || (!is.numeric(max_score))) {
            stop("max_score needs to be set to max or contain a numeric vector of maximum scores for each set")
        }
    }
    if (min_score != "min") {
        if (length(min_score) == (ncol(input_df) - 2) || (!is.numeric(min_score))) {
            stop("min_score needs to be set to min or contain a numeric vector of minimum scores for each set")
        }
    }
    if (!psuedocount > 0) {
        stop("psuedocount must be greater than 0")
    }
    pc <- psuedocount

    # re-order by log ratios (rank order from highest to lowest)
    expt <- colnames(input_df)[2]
    enr_test <- colnames(input_df)[3:ncol(input_df)]
    colnames(input_df)[c(1, 2)] <- c("item", "expression_val")
    input_df <- arrange(input_df, desc(.data$expression_val))

    # get and check size of set items; build in-set matrix of 1's for items in set and 0's for items not in set
    inset_mat <- matrix(0, nrow = length(input_df$item), ncol = length(enr_test))
    dimnames(inset_mat) <- list(input_df$item, enr_test)
    # if list of items provided for each set, check to make sure each item is in dataset and set inset_mat to 1 if it is and 0 if not
    if (thresh_type == "list" & is.list(thresh)) {
        thresh_action <- "exclude"
        for (a in 1:length(thresh)) {
            skip <- ""
            if (length(setdiff(thresh[[a]], input_df$item)) > 0) {
                skip <- setdiff(thresh[[a]], input_df$item)
                warning(paste0(paste(skip, collapse = ", "), " are not items in the input dataframe and will be skipped.\n"))
            }
            for (b in seq_along(thresh[[a]])) {
                if (!thresh[[a]][b] %in% skip) {
                    inset_mat[thresh[[a]][b], a] <- 1
                }
            }
        }
        # if numeric threshold provided, set inset_mat to 1 for items that meet threshold or 0 for those that don't
    } else if (is.numeric(thresh)) {
        if (thresh_type == "values" && length(thresh) == length(enr_test)) {
            for (a in 1:length(thresh)) {
                items_in_set <- input_df$item[input_df[, enr_test[a]] >= thresh[a]]
                inset_mat[items_in_set, enr_test[a]] <- 1
            }
        } else if (thresh_type == "percentile" && thresh > 0 & thresh < 1) {
            thresh1 <- vector(mode = "numeric", length = length(enr_test))
            for (a in 1:length(enr_test)) {
                thresh1[a] <- quantile(input_df[, enr_test[a]], probs = thresh)
                items_in_set <- input_df$item[input_df[, enr_test[a]] >= thresh1[a]]
                inset_mat[items_in_set, enr_test[a]] <- 1
            }
            thresh <- thresh1
        } else if (thresh_type == "val" & length(thresh) == 1) {
            for (a in 1:ncol(inset_mat)) {
                items_in_set <- input_df$item[input_df[, enr_test[a]] >= thresh]
                inset_mat[items_in_set, enr_test[a]] <- 1
            }
            thresh <- rep(thresh, times = ncol(inset_mat))
        }
    } else {
        stop("improper threshold specified")
    }

    # check to make sure minimum number of items present in each set; if not, skip or adjust according to thresh_action
    skipped_sets <- character(0)
    for (c in 1:ncol(inset_mat)) {
        if (sum(inset_mat[, c]) < min_set_size) {
            too_small <- colnames(inset_mat)[c]
            warning(paste0(expt, " does not contain minimum number of items in set for ", too_small, "\n"))
            if (thresh_action == "exclude" | thresh_action == "include") {
                skipped_sets <- c(skipped_sets, too_small)
            } else {
                # lower threshold to first value that would include minimum number of items in set
                lower_scores <- sort(input_df[input_df[too_small] < thresh, too_small], decreasing = T)
                new_thresh <- lower_scores[min_set_size - sum(inset_mat[, c])]
                # set all inset_mat values that meet the lower threshold to 1
                inset_mat[inset_mat[, c] >= new_thresh, c] <- 1
                # check to make sure the set meets size requirements after applying new threshold, if not exclude set
                if ((sum(inset_mat[, c]) < min_set_size) | (sum(inset_mat[, c]) == length(input_df$item))) {
                    skipped_sets <- c(skipped_sets, too_small)
                    warning(paste0("cannot adjust threshold for ", too_small, " to meet size requirements for analysis in ", expt), "; skipping set\n")
                }
            }
        }
    }

    # now check to make sure each set contains fewer items than the max_set_size threshold (default=500); adjust score threshold to obtain valid set size or skip according to thresh_action
    for (c in 1:ncol(inset_mat)) {
        check_col <- colnames(inset_mat)[c]
        if (sum(inset_mat[, c]) > max_set_size) {
            warning(paste0(expt, " has more than ", max_set_size, " items in set ", check_col, "\n"))
            if (thresh_action == "exclude" | thresh_action == "include") {
                skipped_sets <- c(skipped_sets, check_col)
            } else {
                # set items with minimum value to 0 to reduce set size (note: while loop not necessary here because max set size is all items in the dataset, but may be useful if we decide to impose a maximum set size later)
                while (sum(inset_mat[, c]) > max_set_size) {
                    inset_mat[input_df$item[input_df[, check_col] == min(input_df[inset_mat[, c] == 1, check_col])], check_col] <- 0
                }
                # make sure this adjustment didn't remove too many items
                if (sum(inset_mat[, c]) < min_set_size) {
                    skipped_sets <- c(skipped_sets, check_col)
                    warning(paste0("cannot adjust threshold for ", check_col, " to meet size requirements for analysis in ", expt), "; skipping set\n")
                }
            }
        }
    }

    # remove skipped columns
    # print(c("sets that don't contain proper number of items:", skipped_sets))
    inset_mat <- inset_mat[, !(colnames(inset_mat) %in% skipped_sets)]

    if (is.null(ncol(inset_mat))) {
        stop("All gene sets are skipped! Please try to decrease the minimum set size.\n")
    }

    # generate list containing names of items in each set and ranks of those items
    items_in_set <- list()
    rust_parts <- list()
    for (it in 1:ncol(inset_mat)) {
      items_in_set[[colnames(inset_mat)[it]]] <- data.frame(which(inset_mat[, it] == 1), stringsAsFactors = F)
      rust_parts[[it]] <- rownames(items_in_set[[colnames(inset_mat)[it]]])
      colnames(items_in_set[[colnames(inset_mat)[it]]]) <- "rank"
    }
    rust_analytes <- input_df[, 1]
    rust_ranks <- input_df[, 2]
    rust_sets <- colnames(inset_mat)
    rust_result <- gsea_rust(min_set_size, max_set_size, perms, rust_sets, rust_parts, rust_analytes, rust_ranks)
    output_df <- data.frame(fdr = rust_result$fdr, p_val = rust_result$p_val, ES = rust_result$ES, NES = rust_result$NES, leading_edge = rust_result$leading_edge)
    rownames(output_df) <- rust_result$gene_sets
    running_sum <- do.call('cbind', rust_result$running_sum)
    dimnames(running_sum) <- list(rust_analytes, names(rust_result$running_sum))
    if ((thresh_action == "include") & (length(skipped_sets) > 0)) {
        new_row <- data.frame(matrix(0, nrow = 1, ncol = 4), stringsAsFactors = F)
        colnames(new_row) <- colnames(output_df)
        new_row$p_val <- 1
        new_row$fdr <- 1
        for (i in 1:length(skipped_sets)) {
            rownames(new_row) <- skipped_sets[i]
            output_df <- rbind(output_df, new_row)
        }
    }
    output_df$fdr[output_df$fdr > 1] <- 1
    return(list(Enrichment_Results = output_df, Running_Sums = running_sum, Items_in_Set = items_in_set))
}


#' Prepare input for standard GSEA
#'
#' A helper to read files for performing standard GSEA.
#'
#' @param rankFile Path of the rnk file
#' @param gmtFile Path of the GMT file
#'
#' @return a data frame to be used in \code{swGsea}
#'
#' @importFrom readr read_tsv
#'
#' @export
#'
prepareGseaInput <- function(rankFile, gmtFile) {
    rank <- read_tsv(rankFile, col_names = c("gene", "score"), col_types = "cd")
    gmt <- readGmt(gmtFile)
    inputDf <- prepareInputMatrixGsea(rank, gmt)
    return(inputDf)
}
