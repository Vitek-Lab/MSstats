#' @import data.table
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
#' @useDynLib MSstatsdev, .registration=TRUE
NULL

#' log4r appender used not to write messages 
#' 
#' A convenience function written to save time on checking if messages should
#' be printed or logs should be written to a file.
#' 
#' @param level log level 
#' @param ... messages - ignored
#' @return NULL invisibly
#' @keywords internal
.nullAppender = function(level, ...) {
    invisible(NULL)
}


#' Set default logging object when package is loaded
#' @param ... ignored
#' @importFrom log4r file_appender console_appender
#' @return none, sets options called MSstatsLog and MSstatsMsg
#' @keywords internal
.onLoad = function(...) {
    logs = getOption("MSstatsLog")
    msgs = getOption("MSstatsMsg")
    time_now = Sys.time()
    path = paste0("./MSstats_log_", gsub("[ :\\-]", "_", time_now), ".log")
    
    if (is.null(logs)) {
        ms_logs = file_appender(path)
        options(MSstatsLog = ms_logs)
    }
    if (is.null(msgs)) {
        ms_messages = console_appender()
        options(MSstatsMsg = ms_messages)
    }
}
