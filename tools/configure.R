
# MIT License
# Copyright (c) 2020 Hiroaki Yutani

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# The following fields in DESCRIPTION can be used for customizing the behavior.
#
# Config/<package name>/MSRV (optional):
#   Minimum Supported Rust version (e.g. 1.41.0). If this is specified, errors
#   when the installed cargo is newer than the requirement.

SYSINFO_OS      <- tolower(Sys.info()[["sysname"]])
SYSINFO_MACHINE <- Sys.info()[["machine"]]
HAS_32BIT_R     <- dir.exists(file.path(R.home(), "bin", "i386"))
USE_UCRT        <- identical(R.version$crt, "ucrt")


# Utilities ---------------------------------------------------------------

#' Read a field of the package's DESCRIPTION file
#'
#' The field should have the prefix
#'
#' @param field
#'   Name of a field without prefix (e.g. `"check_cargo"`).
#' @param prefix
#'   Prefix of the field (e.g. `"Config/rextendr/`).
#' @param optional
#'   If `TRUE`, return `NA` when there's no field. Otherwise raise an error.
#'
get_desc_field <- function(field, prefix = DESC_FIELD_PREFIX, optional = TRUE) {
  field <- paste0(prefix, field)
  if (length(field) != 1) {
    stop("Field must be length one of character vector")
  }

  # `read.dcf()` always succeeds even when the field is missing.
  # Detect the failure by checking NA
  x <- read.dcf("DESCRIPTION", fields = field)[[1]]

  if (isTRUE(is.na(x)) && !isTRUE(optional)) {
    stop("Failed to get the field ", field, " from DESCRIPTION")
  }

  x
}

# This is tricky; while DESC_FIELD_PREFIX is used in get_desc_field()'s default,
# this variable is defined by get_desc_field(). It's no problem as long as the
# default is not used before it exists!
DESC_FIELD_PREFIX <- paste0("Config/", get_desc_field("Package", prefix = ""), "/")


safe_system2 <- function(cmd, args) {
  result <- list(success = FALSE, output = "")

  output_tmp <- tempfile()
  on.exit(unlink(output_tmp, force = TRUE))

  suppressWarnings(ret <- system2(cmd, args, stdout = output_tmp))

  if (!identical(ret, 0L)) {
    return(result)
  }

  result$output  <- readLines(output_tmp)
  result$success <- TRUE
  result
}

# check_cargo -------------------------------------------------------------

#' Check if the cargo command exists and the version is above the requirements
#'
#' @return
#'   `TRUE` invisibly if no error was found.
check_cargo <- function() {
  ### Check if cargo command works without error ###

  cat("*** Checking if cargo is installed\n")

  cargo_cmd <- "cargo"
  cargo_args <- "version"

  res_version <- safe_system2(cargo_cmd, cargo_args)

  if (!isTRUE(res_version$success)) {
    stop(errorCondition("cargo command is not available", class = c("webgestaltr_error_cargo_check", "error")))
  }

  ### Check the version ###

  msrv <- get_desc_field("MSRV", optional = TRUE)

  if (isTRUE(!is.na(msrv))) {
    cat("*** Checking if cargo is newer than the required version\n")

    version <- res_version$output

    ptn <- "cargo\\s+(\\d+\\.\\d+\\.\\d+)"
    m <- regmatches(version, regexec(ptn, version))[[1]]

    if (length(m) != 2) {
      stop(errorCondition("cargo version returned unexpected result", class = c("webgestaltr_error_cargo_check", "error")))
    }

    if (package_version(m[2]) < package_version(msrv)) {
      msg <- sprintf("The installed version of cargo (%s) is older than the requirement (%s)", m[2], msrv)
      stop(errorCondition(msg, class = c("webgestaltr_error_cargo_check", "error")))
    }
  }

  ### Check the targets ###

  if (identical(SYSINFO_OS, "windows")) {
    cat("*** Checking if the required Rust target is installed\n")

    targets <- safe_system2("rustup", c("target", "list", "--installed"))

    # rustup might not exist if Rust is installed directly via the .msi installer
    # in that case, just ignore and pray that the compilation will succeed.
    if (!isTRUE(targets$success)) {
      return(invisible(TRUE))
    }

    if (!isTRUE("x86_64-pc-windows-gnu" %in% targets$output)) {
      msg <- "The required target x86_64-pc-windows-gnu is not installed"
      stop(errorCondition(msg, class = c("webgestaltr_error_cargo_check", "error")))
    }
  }

  invisible(TRUE)
}

# MAIN --------------------------------------------------------------------

### Check cargo toolchain ###

cargo_check_result <- tryCatch(
  check_cargo(),
  # Defer errors if it's raised by functions here
  webgestaltr_error_cargo_check = function(e) e$message
)

# If cargo is confirmed fine, exit here. But, even if the cargo is not available
# or too old, it's not the end of the world. There might be a pre-compiled
# binary available for the platform.
if (isTRUE(cargo_check_result)) {
  cat("*** cargo is ok\n")
  quit("no", status = 0)
}

cat(sprintf("
-------------- ERROR: CONFIGURATION FAILED --------------------

[cargo check result]
%s

Please refer to <https://www.rust-lang.org/tools/install> to install Rust.

---------------------------------------------------------------

", cargo_check_result))
quit("no", status = 2)