### ============================================================================
### fstArray objects
###

setOldClass("fst_table")
# TODO: Provide access to other args of fst::read.fst() as slots in class?
#' @rdname fstArray-class
#' @export
setClass("fstArraySeed",
         contains = "Array",
         representation(
           fst_table = "fst_table"
         )
)

### ----------------------------------------------------------------------------
### dim()
###

#' @importFrom fst metadata_fst
#' @export
setMethod("dim", "fstArraySeed", function(x) {
  # TODO: dim.fst_table returns a double because `NrOfRows` is stored as a
  #       double rather than an int; is this a a real issue (E.g., construct a
  #       fst with > .Machine$integer.max rows and see what breaks) or is it
  #       just cosmetic
  as.integer(dim(x@fst_table))
})

### ----------------------------------------------------------------------------
### dimnames()
###

# TODO: fst file has no concept of rownames; how to handle in fstArray?
#       1. Don't support rownames
#       2. Store rownames in _fstArraySeed_ (will require care when subsetting)
#' @importFrom fst metadata_fst
#' @export
setMethod("dimnames", "fstArraySeed", function(x) {
  dimnames(x@fst_table)
})

### ----------------------------------------------------------------------------
### path()
###

#' @importFrom BiocGenerics path
#' @export
setMethod("path", "fstArraySeed", function(object, ...) {
  .subset2(object@fst_table, "meta")$path
})

### ----------------------------------------------------------------------------
### extract_array()
###

.extract_array_from_fstArraySeed <- function(x, index) {
  # TODO: Avoid coercion to matrix which incurs a copy; see fstRead()
  as.matrix(fstRead(x@fst_table, index))
}

#' @importMethodsFrom DelayedArray extract_array
#' @export
setMethod("extract_array", "fstArraySeed", .extract_array_from_fstArraySeed)

### ----------------------------------------------------------------------------
### fstArraySeed() constructer
###

#' @importFrom fst fst
#' @importFrom S4Vectors isTRUEorFALSE new2
#' @rdname fstArray-class
#' @export
fstArraySeed <- function(filepath, old_format = FALSE) {
  if (!isTRUEorFALSE(old_format)) {
    stop("'old_format' must be TRUE or FALSE")
  }
  if (is(filepath, "fst_table")) {
    if (old_format != .subset2(filepath, "old_format")) {
      stop("The specified 'old_format' does not match that of the 'fst_table'")
    }
    fst_table <- filepath
  } else {
    filepath <- .normarg_filepath(filepath)
    fst_table <- fst(path = filepath, old_format = old_format)
  }
  # TODO: Check that all columns of fst_table are of same type (or coercible)?
  new2("fstArraySeed", fst_table = fst_table)
}

### ----------------------------------------------------------------------------
### fstArray and fstMatrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate fstArray and fstMatrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

#' fst datasets as array-like objects
#'
#' @description We provide 2 classes for representing an (on-disk) fst dataset
#' as an array-like object in R:
#'
#' - _fstArray_: A high-level class, _fstArray_ extends
#' [DelayedArray::DelayedArray-class]. All the operations available on
#' [DelayedArray::DelayedArray-class] objects work on _fstArray_ objects.
#' - _fstArraySeed_: A low-level class for pointing to a fst dataset. No
#' operation can be performed directly on an _fstArraySeed_ object. It needs to
#' be wrapped in a [DelayedArray::DelayedArray-class] or _fstArray_ object first.
#' An _fstArray_ object is just an _fstArraySeed_ wrapped in a
#' [DelayedArray::DelayedArray-class] object.
#'
#' @param filepath The path (as a single character string) to the fst file where
#' the dataset is located or an [fst::fst_table][fst::fst] instance.
#' @param old_format use `TRUE` to read fst files generated with a fst package
#' version lower than v0.8.0.
#'
#' @return An _fstArray_ object for `fstArray()`.
#'
#' An _fstArraySeed_ object for `fstArraySeed()`.
#'
#' @seealso
#' - [fst::write_fst()] for writing a data frame to disk as an 'fst' file.
#' - [DelayedArray::DelayedArray-class] objects.
#' - [DelayedArray::DelayedArray-utils] for common operations on [DelayedArray::DelayedArray-class] objects
#' - [base::array] objects in base R.
#'
#' @examples
#' # ----------------------------------------------------------------------
#' # CONSTRUCTION
#' # ----------------------------------------------------------------------
#'
#' # Simulate some data in a data frame and write to disk as an 'fst' file
#' library(fst)
#' x <- as.data.frame(replicate(10, runif(10000)))
#' fst_file <- tempfile(fileext = ".fst")
#' write_fst(x, fst_file)
#'
#' # Construct an fstArray from an fst file
#' fst_array <- fstArray(fst_file)
#' fst_array
#' # Construct an fstArray from an fst::fst_table
#' fst_table <- fst(fst_file)
#' class(fst_table)
#' fstArray(fst_table)
#'
#' # ----------------------------------------------------------------------
#' # dim/dimnames
#' # ----------------------------------------------------------------------
#'
#' dim(fst_array)
#'
#' dimnames(fst_array)
#' dimnames(fst_array) <- list(paste0("gene", seq_len(nrow(fst_array))),
#'                             paste0("S", seq_len(ncol(fst_array))))
#' fst_array
#'
#' # ----------------------------------------------------------------------
#' # SLICING (A.K.A. SUBSETTING)
#' # ----------------------------------------------------------------------
#'
#' fst_array2 <- drop(fst_array[31:40, c("S3", "S6")])
#' fst_array2
#'
#' dim(fst_array2)
#' as.array(fst_array2)
#' stopifnot(identical(dim(as.array(fst_array2)), dim(fst_array2)))
#' stopifnot(identical(dimnames(as.array(fst_array2)), dimnames(fst_array2)))
#' # ----------------------------------------------------------------------
#' # SummarizedExperiment OBJECTS WITH DELAYED ASSAYS
#' # ----------------------------------------------------------------------
#'
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(fst_array)
#' se
#'
#' stopifnot(validObject(se, complete = TRUE))
#'
# NOTE: These aliases are necessary because I want them in the .Rd but not
#         with a S4method{} in the 'Usage' section (roxygen2 is doing this)
#' @rawRd \alias{dim,fstArraySeed-method}
#' @rawRd \alias{dimnames,fstArraySeed-method}
#' @rawRd \alias{extract_array,fstArraySeed-method}
#' @rawRd \alias{DelayedArray,fstArraySeed-method}
#' @rawRd \alias{path,fstArraySeed-method}
#'
#' @export
setClass("fstArray", contains = "DelayedArray")

#' @rdname fstArray-class
#' @export
setClass("fstMatrix", contains = c("DelayedMatrix", "fstArray"))

setAs("fstArray", "fstMatrix", function(from) new("fstMatrix", from))

# NOTE: For internal use only.
setMethod("matrixClass", "fstArray", function(x) "fstMatrix")

#' @importFrom DelayedArray seed
#' @importFrom S4Vectors wmsg
.validate_fstArray <- function(x){
  if (!is(seed(x), "fstArraySeed")) {
    return(wmsg("'x@seed' must be a fstArraySeed object"))
  }
  if (!DelayedArray:::is_pristine(x)) {
    return(wmsg("'x' carries delayed operations"))
  }
  TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("fstArray", .validate_fstArray)

setAs("ANY", "fstMatrix",
      function(from) as(as(from, "fstArray"), "fstMatrix")
)

### ----------------------------------------------------------------------------
### Constructor
###

#' @export
setMethod("DelayedArray", "fstArraySeed", function(seed) {
  DelayedArray:::new_DelayedArray(seed, Class = "fstArray")
})

# NOTE: Works directly on a fstArraySeed object, in which case it must be
#       called with a single argument.
#' @importMethodsFrom DelayedArray DelayedArray
#' @rdname fstArray-class
#' @export
fstArray <- function(filepath, old_format = FALSE) {
  if (is(filepath, "fstArraySeed")) {
    seed <- filepath
  } else {
    seed <- fstArraySeed(filepath, old_format)
  }
  DelayedArray(seed)
}
