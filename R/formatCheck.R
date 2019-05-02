#' Check Format and Read Data
#'
#' @inheritParams idMapping
#'
#' @return A list of data frame
#'
#' @importFrom tools file_ext
#' @importFrom readr read_tsv stop_for_problems
#'
#' @export
#'
formatCheck <- function(dataType="list", inputGeneFile=NULL, inputGene=NULL) {
	dataTypeA <- c("list", "rnk")
	if (length(dataType)>1 || !is.character(dataType) || length(which(dataTypeA==dataType))==0) {
		stop("ERROR: DataType parameter can only be one of 'list' or 'rnk'.")
	}

	if (dataType=="list") {
		if (!is.null(inputGeneFile)) {
			if (file_ext(inputGeneFile)!="txt") {
				stop("ERROR: For the user ID list, please upload a 'txt' file with only one column.")
			} else {
				inputGene <- read_tsv(inputGeneFile, col_names=FALSE, col_types="c")
				error <- tryCatch(stop_for_problems(inputGene),
					error=function(e) stop("ERROR: The format of the uploaded gene list is incorrect. Please make sure the file only contains one column of gene IDs.")
				)
				return(inputGene[[1]])
			}
		} else {
			if (!is.null(inputGene)) {
				if (!is.vector(inputGene)) {
					stop("ERROR: For the user ID list, please upload an R vector object.")
				} else {
					inputGene <- unname(sapply(inputGene, as.character)) # as.character will convert length-1 vector to primitive
					return(inputGene)
				}
			} else {
				stop("ERROR: Please upload a file or an R object for the ID list.")
			}
		}
	}

	if (dataType=="rnk") {
		if (!is.null(inputGeneFile)) {
			if (file_ext(inputGeneFile)!="rnk") {
				stop("ERROR: For the ranked list, please upload a 'rnk' file with two columns (ids and scores).")
			} else {
				inputGene <- read_tsv(inputGeneFile, col_names=c("gene", "score"), col_types="cd")
				# readr just gives warning
				if (all(is.na(inputGene$gene)) || all(is.na(inputGene$score))) {
					stop("ERROR: For the ranked list, please upload a 'rnk' file with two columns (ids and scores).")
				}

				error <- tryCatch(stop_for_problems(inputGene),
					error=function(e) stop("ERROR: The format of the uploaded ranked list is incorrect. Please make sure the file contains two column of gene IDs and scores.")
				 )

				#########GSEA do not allow the second column contains NA. Thus, we should remove NA first############
				return(filter(inputGene, !is.na(.data$score)))
			}
		} else {
			if (!is.null(inputGene)) {
				if (!is.data.frame(inputGene)) {
					stop("ERROR: For the ranked list, please upload an R data.frame object.")
				} else {
					if (!is.numeric(inputGene[[2]]) && !is.integer(inputGene[[2]])) {
						stop("ERROR: The second column of the ranked list should be the numeric scores.")
					} else {
						#########GSEA do not allow the second column contains NA. Thus, we should remove NA first############
						inputGene <- inputGene[!is.na(inputGene[[2]]), ]
						inputGene[[1]] <- as.character(inputGene[[1]])
						colnames(inputGene) <- c("gene", "score")
						return(inputGene)
					}
				}
			} else {
				stop("ERROR: Please upload a file or an R object for the ranked list.")
			}
		}
	}
}
