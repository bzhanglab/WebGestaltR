#' WebGestaltR: The R interface for enrichment analysis with WebGestalt.
#'
#' @docType package
#' @name WebGestaltR
#' @import methods
#' @import grDevices
#' @import graphics
#' @import utils
#' @importFrom stats p.adjust phyper
#' @importFrom Rcpp sourceCpp
#' @useDynLib WebGestaltR
#'
NULL

.onAttach <- function(lib, pkg) {
	packageStartupMessage("******************************************\n")
	packageStartupMessage("*                                        *\n")
	packageStartupMessage("*          Welcome to WebGestaltR !      *\n")
	packageStartupMessage("*                                        *\n")
	packageStartupMessage("******************************************\n")
}
