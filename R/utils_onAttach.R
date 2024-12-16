#' Load core packages when loading MSstats
#' @noRd
#' @keywords internal
.onAttach = function(libname, pkgname) {
    # Load core packages dynamically
    to_load <- c("MSstatsConvert")
    
    # Attach them quietly
    if (length(to_load) > 0) {
        lapply(to_load, function(pkg) {
            suppressPackageStartupMessages(library(pkg, character.only = TRUE))
        })
    }
}