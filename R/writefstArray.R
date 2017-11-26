### =========================================================================
### writeHDF5Array()
### -------------------------------------------------------------------------
###

# TODO: Do this properly with realization sink (see writeHDF5Array)
#' @importFrom fst write.fst
#' @export
writefstArray <- function(x, file, compress = 0) {
  dim_x <- dim(x)
  if (length(dim_x) != 2L) {
    stop("'x' must be 2-dimensional")
  }
  x <- as.data.frame(x)
  write.fst(x, path = file, compress = compress)
  fstArray(file)
}
