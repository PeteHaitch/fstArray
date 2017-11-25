### =========================================================================
### writeHDF5Array()
### -------------------------------------------------------------------------
###

# TODO: Do this properly with realization sink (see writeHDF5Array)
#' @importFrom fst write.fst
#' @export
writefstArray <- function(x, file, compress = 0) {
  x <- as.data.frame(x)
  write.fst(x, path = file, compress = compress)
  fstArray(file)
}
