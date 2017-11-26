### ============================================================================
### writefstArray()
### ----------------------------------------------------------------------------
###

#' @include utils.R

### ----------------------------------------------------------------------------
### fstRealizationSink objects
###
### The fstRealizationSink class is a concrete RealizationSink subclass that
### implements an fst realization sink.
###

#' @export
setClass("fstRealizationSink",
         contains = "RealizationSink",
         representation(
           dim = "integer",
           file = "character",
           compress = "integer"
         )
)

# TODO: dimnames,fstRealizationSink-method (could just be a slot of the class
#       like HDF5RealizationSink)

#' @importFrom S4Vectors new2
#' @export
fstRealizationSink <- function(dim, file = NULL, compress = NULL) {
  if (is.null(file)) {
    file <- getfstDumpFile(for.use = TRUE)
  } else {
    file <- normalize_dump_file(file)
  }
  if (is.null(compress)) {
    compress <- getfstDumpCompressionLevel()
  } else {
    compress <- normalize_compression_level(compress)
  }
  appendDatasetCreationTofstDumpLog(file, dim, compress)
  if (is.null(dimnames)) {
    dimnames <- vector("list", length(dim))
  } else {
    ## TODO: Write the rownames to the fst file?
  }
  new2("fstRealizationSink", dim = dim, file = file, compress = compress)
}

# TODO: This needs to append to the existing sink@file, otherwise only the last
#       block will be in the fst file after running writefstArray()
#' @importFrom DelayedArray makeNindexFromArrayViewport
#' @importFrom fst write.fst
#' @export
setMethod("write_block_to_sink", "fstRealizationSink",
          function(block, sink, viewport) {
            stopifnot(identical(dim(sink), refdim(viewport)),
                      identical(dim(block), dim(viewport)))
            index <- makeNindexFromArrayViewport(viewport,
                                                 expand.RangeNSBS = TRUE)
            # TODO: Is there a way to avoid the coercion to a data frame before
            #       writing to the fst file?
            block <- as.data.frame(block)
            write.fst(x = block, path = sink@file, compress = sink@compress)
          }
)

### ----------------------------------------------------------------------------
### Coercing a fstRealizationSink object.
###

# TODO: This coercion needs to propagate the rownames *thru* the fst file.
setAs("fstRealizationSink", "fstArraySeed",
      function(from) fstArraySeed(from@file)
)

# NOTE: This coercion currently drops the rownames but will naturally propagate
#       them when coercion from fstRealizationSink to fstArraySeed propagates
#       them. See TODO above.
setAs("fstRealizationSink", "fstArray",
      function(from) fstArray(as(from, "fstArraySeed"))
)

setAs("fstRealizationSink", "DelayedArray",
      function(from) {
        ans <- fstArray(as(from, "fstArraySeed"))
        ## Temporarily needed because coercion from fstRealizationSink to
        ## HDF5ArraySeed does not propagate the rownames at the moment. See
        ## TODO above.
        ## TODO: Remove line below when TODO above is addressed.
        rownames(ans) <- rownames(from)
        ans
      }
)

### ----------------------------------------------------------------------------
### writefstArray()
###

# NOTE: Write the array to the current dump if 'file' is not specified.
# Return a fstArray object pointing to the newly written fst dataset on disk.
# TODO: How to handle rownames?
#' @importFrom S4Vectors isTRUEorFALSE
#' @export
writefstArray <- function(x, file = NULL, compress = NULL, verbose = FALSE) {
  x_dim <- dim(x)
  if (length(x_dim) != 2L) {
    stop("'x' must be 2-dimensional")
  }
  if (!isTRUEorFALSE(verbose)) {
    stop("'verbose' must be TRUE or FALSE")
  }
  sink <- fstRealizationSink(dim = x_dim, file = file, compress = compress)
  if (verbose) {
    old_verbose <- DelayedArray:::set_verbose_block_processing(verbose)
    on.exit(DelayedArray:::set_verbose_block_processing(old_verbose))
  }
  write_array_to_sink(x, sink)
  as(sink, "fstArray")
}

### ----------------------------------------------------------------------------
### Coercion to fstArray
###
### The methods below write the object to disk. Note that coercion from
### fstRealizationSink to fstArray is already taken care of by the specific
### method above and doesn't write anything to disk. So coercing to fstArray
### in general writes the object to disk *except* when the object to coerce is
### a fstRealizationSink object.
###

# NOTE: writes to current dump
.as_fstArray <- function(from) writefstArray(from)

setAs("ANY", "fstArray", .as_fstArray)

# NOTE: Automatic coercion method from DelayedArray to fstArray silently
#       returns a broken object (unfortunately these dummy automatic coercion
#       methods don't bother to validate the object they return). So we
#       overwrite it.
setAs("DelayedArray", "fstArray", .as_fstArray)
setAs("DelayedMatrix", "fstMatrix", .as_fstArray)
