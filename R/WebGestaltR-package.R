#' WebGestaltR: The R interface for enrichment analysis with WebGestalt.
#'
#' @docType package
#' @name WebGestaltR
#' @import methods
#' @import grDevices
#' @import graphics
#' @import utils
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @useDynLib WebGestaltR
#'
NULL


.onAttach <- function(lib, pkg) {
  packageStartupMessage("******************************************
")
  packageStartupMessage("*                                        *
")
  packageStartupMessage("*       Welcome to WebGestaltR-rust!     *
")
  packageStartupMessage("*                                        *
")
  packageStartupMessage("*                                        *
")
  packageStartupMessage("******************************************
")
}
