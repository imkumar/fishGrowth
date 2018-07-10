
#' @import TMB
#' @useDynLib fishGrowth

.onUnload <- function(lib) {
    library.dynam.unload("fishGrowth", lib)
}
