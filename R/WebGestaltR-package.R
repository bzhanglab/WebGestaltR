#' WebGestaltR: The R interface for enrichment analysis with WebGestalt.
#'
#' @docType package
#' @name WebGestaltR
#' @import methods
#' @import grDevices
#' @import graphics
#' @import utils
#' @importFrom Rcpp sourceCpp
#' @useDynLib WebGestaltR
#'
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

.onAttach <- function(lib, pkg) {
	packageStartupMessage("******************************************\n")
	packageStartupMessage("*                                        *\n")
	packageStartupMessage("*          Welcome to WebGestaltR !      *\n")
	packageStartupMessage("*                                        *\n")
	packageStartupMessage("******************************************\n")
}
