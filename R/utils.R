### ============================================================================
### Some low-level utilities
### Functions in this file are not exported
###

#' @importFrom S4Vectors isSingleString wmsg
#' @importFrom tools file_path_as_absolute
.normarg_filepath <- function (filepath, what) {
    if (!isSingleString(filepath)) {
        stop(wmsg(what, " must be a single string specifying the path ",
            "to the matter file where the dataset is located"))
    }
    file_path_as_absolute(filepath)
}
