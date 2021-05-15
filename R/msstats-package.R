#' @import data.table
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
#' @useDynLib MSstats, .registration=TRUE
NULL


#' Set default logging object when package is loaded
#' @param ... ignored
#' @return none, sets options called MSstatsLog and MSstatsMsg
#' @keywords internal
.onLoad = function(...) {
    MSstatsConvert::MSstatsLogsSettings()
}
