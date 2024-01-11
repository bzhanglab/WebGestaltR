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
#' @param input_df_list A data frame in which first column is name of item of interest
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
#' @param fdrMethod For the ORA method, WebGestaltR supports five FDR methods: \code{holm},
#'   \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH} and \code{BY}. The default
#'   is \code{BH}.
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
#' @importFrom poolr stouffer
#' @importFrom doRNG %dorng%
#' @export
#' @author John Elizarraras
#'
multiswGsea <- function(input_df_list, thresh_type = "percentile", thresh = 0.9, thresh_action = "exclude", min_set_size = 10,
                        max_set_size = 500, max_score = "max", min_score = "min", psuedocount = 0.001, perms = 1000, p = 1,
                        q = 1, nThreads = 1, rng_seed = 1, fork = FALSE, fdrMethod = "BH") {
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


    output_df_list <- list()
    gseaRes_list <- list()
    gseaRes_list[[1]] <- NULL
    # output_df_list[[1]] <- NULL # Set null for now. Will be filled later.
    # running_sum_list[[1]] <- NULL
    # items_in_set_list[[1]] <- NULL

    for (i in seq_along(input_df_list)) {
        inputDf <- input_df_list[[i]]
        gseaRes <- swGsea(inputDf,
            thresh_type = "val", perms = perms,
            min_set_size = min_set_size, max_set_size = max_set_size, p = p,
            nThreads = nThreads, rng_seed = rng_seed
        )
        gseaRes_list[[i + 1]] <- gseaRes
        output_df_list[[i]] <- gseaRes$Enrichment_Results
        # running_sum_list[[i + 1]] <- gseaRes$Running_Sums
        # items_in_set_list[[i + 1]] <- gseaRes$Items_in_Set
    }

    all_gene_sets <- unique(unlist(lapply(output_df_list, rownames)))
    meta_ps <- list()
    biggest_p <- 1 - .Machine$double.eps
    meta_items_in_sets <- list()
    for (i in seq_along(all_gene_sets)) {
        gene_set <- all_gene_sets[[i]]
        p_vals <- c()
        for (j in seq_along(output_df_list)) {
            if (gene_set %in% rownames(output_df_list[[j]])) {
                list_p <- output_df_list[[j]][gene_set, "p_val"]
                if (list_p == 0.0) {
                    list_p <- .Machine$double.eps
                } else if (list_p > biggest_p) {
                    list_p <- biggest_p
                }
                p_vals <- append(p_vals, list_p)
                relevant_items <- gseaRes_list[[j + 1]]$Items_in_Set[[gene_set]]
                if (length(meta_items_in_sets) < i) {
                    meta_items_in_sets[[i]] <- relevant_items
                } else {
                    meta_items_in_sets[[i]] <- rbind(meta_items_in_sets[[i]], relevant_items)
                }
            }
        }
        if (length(p_vals) < 2) {
            meta_ps[[i]] <- p_vals[1]
        } else {
            meta_ps[[i]] <- stouffer(p_vals)$p[1]
        }
    }

    meta_fdrs <- sapply(meta_ps, function(x) {
        p.adjust(x, method = fdrMethod)
    })

    meta_output_df <- data.frame(
        fdr = unlist(meta_fdrs), p_val = unlist(meta_ps), ES = numeric(length(all_gene_sets)), NES = numeric(length(all_gene_sets)),
        leading_edge = numeric(length(all_gene_sets)), stringsAsFactors = FALSE
    )
    rownames(meta_output_df) <- all_gene_sets
    gseaRes_list[[1]] <- list(Enrichment_Results = meta_output_df, Running_Sums = numeric(nrow(meta_output_df)), Items_in_Set = meta_items_in_sets)

    return(gseaRes_list)
}