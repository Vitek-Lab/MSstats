## I am assuming but am not certain that when R byte-compiles this directory, it
## will read this first and create its environment with this material first,
## since I prefixed it with '01'.  I should probably test this and see if it is
## true.  In any event, maybe it does not matter as the important stuff is in a
## .onLoad?

## The following was taken from ggplot2's ggplot2.r
## I presume it is a blanket importer cue for roxygen2 to add
## import statements to the NAMESPACE file so that when ggplot2 is used
## it will ensure that these libraries are available.
## I checked the roxygen documentation and it appears that
## imports are saved as the exclusive set, as a result repeating these
## at each function declaration serves to make explicit what each function
## requires while not (I think) adding excessive cruft to the NAMESPACE


#' MSstats: a suite of tools to make MS analyses easier.
#'
#' This provides a series of helpers for working with MS/MS data, DDA/DIA/etc.
#'
#' @docType package
#' @name MSstats
#' @import Biobase
#' @import data.table
#' @import doSNOW
#' @import dplyr
#' @import foreach
#' @import ggplot2
#' @import ggrepel
#' @import gplots
#' @import graphics
#' @import grDevices
#' @import grid
#' @import gtable
#' @import knitr
#' @import lme4
#' @import limma
#' @import logging
#' @import magrittr
#' @import marray
#' @import methods
#' @import minpack.lm
#' @import MSnbase
#' @import preprocessCore
#' @import Rcpp
#' @import reshape2
#' @import scales
#' @import snow
#' @import stats
#' @import survival
#' @import utils
NULL

#' Pipe operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' data.table's funky column assignment operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name :=
#' @rdname column_assignment
#' @keywords internal
#' @export
#' @importFrom data.table :=
NULL

#' dopar
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %dopar%
#' @rdname dopar
#' @keywords internal
#' @export
#' @importFrom foreach %dopar%
NULL

#' Set up a log4j-ish log.
#'
#' There are a few places in this code which do various logging tasks.
#' There exist a few packages in R which handle this type of work reasonably well,
#' I arbitrarily chose 'logging'.
#' @name .onLoad
#' @rdname onLoad
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  log_start <- try(suppressWarnings(file.remove("msstats.log")), silent=TRUE)
  log_start <- logging::logReset()
  log_start <- logging::basicConfig()
  log_start <- logging::addHandler(logging::writeToFile,
                                   file="msstats.log",
                                   level="WARN")
  log_start <- logging::addHandler(logging::writeToFile,
                                   file="msstats.log",
                                   level="INFO")
  log_start <- logging::addHandler(logging::writeToFile,
                                   file="msstats.log",
                                   level="ERROR")
  logging::writeToFile("Writing sessionInfo, this should be at the end I think.", log_start)
  logging::writeToFile(as.matrix(sessionInfo(), header=TRUE, sep="\t"), log_start)
}
