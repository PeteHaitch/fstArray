### ============================================================================
### fstArray objects
### ----------------------------------------------------------------------------

# TODO: Provide access to other args of fst::read.fst() as slots in class?
#' @rdname fstArray-class
#' @export
setClass("fstArraySeed",
         contains = "Array",
         slots = c(
           file = "character",
           old_format = "logical"
         )
)

### ----------------------------------------------------------------------------
### dim()
###

#' @export
setMethod("dim", "fstArraySeed", function(x) {
  metadata <- fst::metadata_fst(path = x@file, old_format = x@old_format)
  # TODO: `NrOfRows` is stored as a double rather than an int; is this an
  #       issue? E.g., construct a fst with > .Machine$integer.max rows and see
  #       what breaks
  c(as.integer(metadata$nrOfRows), length(metadata$columnBaseTypes))
})

### ----------------------------------------------------------------------------
### dimnames()
###

# TODO: fst file has no concept of rownames; how to handle in fstArray?
#       1. Don't support rownames
#       2. Store rownames in fstArraySeed (will require care when subsetting)
#' @importFrom fst metadata_fst
#' @export
setMethod("dimnames", "fstArraySeed", function(x) {
  metadata <- metadata_fst(path = x@file, old_format = x@old_format)
  list(NULL, metadata$columnNames)
})

### ----------------------------------------------------------------------------
### extract_array()
###

# TODO: `DelayedArray:::expand_Nindex_RangeNSBS()` may be useful, here.
# TODO: See some of the tricks used by `[`,DelayedArray-method when called with
#       a linear index. Briefly, it loads column-oriented blocks of the data,
#       subsets the relevant elements, and then joins these up in the correct
#       order at the end. This should be possible here to minimise the number
#       of disk I/O ops
#' @importFrom IRanges findOverlaps IRanges reduce
#' @importFrom S4Vectors .Call2 subjectHits wmsg
#' @importFrom BiocGenerics end start
#' @importFrom fst metadata_fst read_fst
# NOTE: rows must be a integer vector (although I think an fst file can have
#       > .Machine$integer.max rows)
# NOTE: cols must be a character vector
.extract_data_frame_from_fstArraySeed <-
  function(x, index, as.data.table = FALSE) {
    if (as.data.table) {
      # NOTE: Would need to be careful if supporting as.data.table = TRUE when
      #       multiple chunks of rows are read from the .fst file and need to
      #       be combined
      stop("Currently, 'as.data.table' must be FALSE")
    }
    # TODO: Sanity check rows and cols
    stopifnot(length(index) == 2L)
    rows <- index[[1L]]
    columns <- index[[2L]]
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L)) {
      metadata <- metadata_fst(path = x@file, old_format = x@old_format)
      # TODO: `types` adapted fst:::print.fst.metadata(); find a more robust
      #       way to infer column types
      types <- c("unknown", "character", "factor", "ordered factor",
                 "integer", "POSIXct", "difftime", "IDate", "ITime", "double",
                 "Date", "POSIXct", "difftime", "ITime", "logical", "integer64",
                 "nanotime", "raw")
      column_types <- metadata$columnTypes
      # TODO: This is perhaps too strict. E.g., if column types are 'integer'
      #       and 'double' then could promote result to 'double' and everything
      #       should work. But, for now, we'll keep
      if (any(column_types != column_types[1L])) {
        stop(wmsg(x@file, " contains multiple column types: ",
                  paste0(types[unique(column_types)], collapse = ", ")))
      }
      type <- types[column_types[[1L]]]
      ans <- switch(EXPR = type,
                    "integer" = integer(),
                    "double" = numeric(),
                    last = NULL)
      if (is.null(ans)) {
        stop(wmsg(x@file, " contains column type unsupported by fstArray: ",
                  type))
      }
      dim(ans) <- ans_dim
    } else {
      if (is.numeric(columns)) {
        # NOTE: fst::read.fst() requires column names rather than indices
        columns <- colnames(x)[columns]
      }
      if (is.null(rows)) {
        ans <- read_fst(path = x@file,
                        columns = columns,
                        from = 1,
                        to = NULL,
                        as.data.table = FALSE,
                        old_format = x@old_format)
      } else {
        rows_as_ranges <- IRanges(rows, width = 1L)
        reduced_ranges <- reduce(rows_as_ranges)
        froms <- start(reduced_ranges)
        tos <- end(reduced_ranges)
        ans_list <- mapply(function(from, to) {
          read_fst(path = x@file,
                   columns = columns,
                   from = from,
                   to = to,
                   as.data.table = FALSE,
                   old_format = x@old_format)
        }, from = froms, to = tos, SIMPLIFY = FALSE)
        ans <- do.call(rbind, ans_list)
        # NOTE: Need to rearrange if row indices are unsorted and/or contain
        #       duplicates
        if (is.unsorted(rows) || anyDuplicated(rows)) {
          ol <- findOverlaps(rows_as_ranges, sort(rows_as_ranges))
          ans <- ans[subjectHits(ol), , drop = FALSE]
        }
      }
    }
    ans
  }

#' @importMethodsFrom DelayedArray extract_array
#' @export
setMethod("extract_array", "fstArraySeed", function(x, index) {
  # TODO: Avoid coercion to matrix which incurs a copy; can I read from a
  #       fst file directly into a matrix of appropriate dimensions?
  as.matrix(.extract_data_frame_from_fstArraySeed(x = x,
                                                  index = index,
                                                  as.data.table = FALSE))
})

### ----------------------------------------------------------------------------
### fstArraySeed() constructer
###

#' @importFrom S4Vectors isSingleString isTRUEorFALSE new2 wmsg
#' @importFrom tools file_path_as_absolute
#' @rdname fstArray-class
#' @export
fstArraySeed <- function(file, old_format = FALSE) {
  if (!isSingleString(file)) {
    stop(wmsg("'file' must be a single string specifying the path to ",
              "the fst file where the dataset is located"))
  }
  if (!isTRUEorFALSE(old_format)) {
    stop("'old_format' must be TRUE or FALSE")
  }
  # TODO: Check that all columns of fst are of same type (or coercible)?
  file <- file_path_as_absolute(file)
  new2("fstArraySeed", file = file, old_format = old_format)
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
#' - fstArray: A high-level class, fstArray extends
#' [DelayedArray::DelayedArray-class]. All the operations available on
#' [DelayedArray::DelayedArray-class] objects work on fstArray objects.
#' - fstArraySeed: A low-level class for pointing to a fst dataset. No
#' operation can be performed directly on an fstArraySeed object. It needs to
#' be wrapped in a [DelayedArray::DelayedArray-class] or fstArray object first.
#' A fstArray object is just a fstArraySeed wrapped in a
#' [DelayedArray::DelayedArray-class] object.
#'
#' @param file The path (as a single character string) to the fst file where
#' the dataset is located.
#' @param old_format use `TRUE` to read fst files generated with a fst package
#' version lower than v0.8.0
#'
#' @return An fstArray object for `fstArray()`.
#'
#' An fstArraySeed object for `fstArraySeed()`.
#'
#' @seealso
#' - [DelayedArray::DelayedArray-class] objects.
#' - [DelayedArray::DelayedArray-utils] for common operations on [DelayedArray::DelayedArray-class] objects
#' - [base::array] objects in base R.
#' 
#'
# NOTE: These aliases are necessary because I want them in the .Rd but not 
#         with a S4method{} in the 'Usage' section (roxygen2 is doing this)
#' @rawRd \alias{dim,fstArraySeed-method}
#' @rawRd \alias{dimnames,fstArraySeed-method}
#' @rawRd \alias{extract_array,fstArraySeed-method}
#' @rawRd \alias{DelayedArray,fstArraySeed-method}
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
### as.data.frame
###

# TODO: Export?
#' @importFrom DelayedArray seed
setMethod("as.data.frame", "fstMatrix", function(x, ...) {
  .extract_data_frame_from_fstArraySeed(x = seed(x),
                                        index = unname(x@index),
                                        as.data.table = FALSE)
})

# TODO: Would this be useful?
#       This would introduce a dependency on data.table. If this were supported
#       it would nice to only have data.table in Suggests rather than Imports.
#       I don't know if it's possible to have data.table in Suggests and still
#       define an as.data.table,fstMatrix-method. Also, data.table may soon be
#       non-optional dependency of fst
#       (https://github.com/fstpackage/fst/blob/develop/NEWS.md). See
#       https://github.com/tidyverse/dbplyr/blob/6be777d8b23d588f19c98de52f4e58f16c2ef67e/R/zzz.R for example of conditional S3 registration and
#       https://stat.ethz.ch/pipermail/r-package-devel/2017q4/002159.html
#       for discussion of various options
# #' @importFrom data.table
# #' @importFrom DelayedArray seed
# #' @export
# setMethod("as.data.table", "fstMatrix", function(x, ...) {
#   .extract_data_frame_from_fstArraySeed(x = seed(x),
#                                         index = unname(x@index),
#                                         as.data.table = TRUE)
# })

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
fstArray <- function(file, old_format = FALSE) {
  if (is(file, "fstArraySeed")) {
    # TODO: Check old_format agrees with file@old_format?
    seed <- file
  } else {
    seed <- fstArraySeed(file, old_format)
  }
  DelayedArray(seed)
}
