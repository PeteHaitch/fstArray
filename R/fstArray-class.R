### ============================================================================
### fstArray objects
### ----------------------------------------------------------------------------

#' @include utils.R

#' @export
setClass("fstArraySeed",
         contains = "Array",
         slots = c(
           file = "character"
         )
)

### ----------------------------------------------------------------------------
### dim()
###

#' @importFrom fst fst.metadata
#' @export
setMethod("dim", "fstArraySeed", function(x) {
  metadata <- fst.metadata(x@file)
  # TODO: `NrOfRows` is stored as a double rather than an int; is this an
  #       issue? E.g., construct a fst with > .Machine$integer.max rows and see
  #       what breaks
  c(as.integer(metadata$NrOfRows), length(metadata$ColumnTypes))
})

### ----------------------------------------------------------------------------
### dimnames()
###

# TODO: fst file has no concept of rownames; how to handle in fstArray?
#       1. Don't support rownames
#       2. Store rownames in fstArraySeed (will require care when subsetting)
#' @importFrom fst fst.metadata
#' @export
setMethod("dimnames", "fstArraySeed", function(x) {
  metadata <- fst.metadata(x@file)
  list(NULL, metadata$ColumnNames)
})

### ----------------------------------------------------------------------------
### extract_array()
###

#' @importFrom fst fst.metadata read.fst
#' @importFrom IRanges findOverlaps IRanges reduce
#' @importFrom S4Vectors .Call2 subjectHits
#' @importFrom BiocGenerics end start
# NOTE: rows must be a integer vector (although I think an fst file can have
#       > .Machine$integer.max rows)
# NOTE: cols must be a character vector
.extract_data_frame_from_fstArraySeed <-
  function(x, index, as.data.table = FALSE) {
    # TODO: Sanity check rows and cols
    stopifnot(length(index) == 2L)
    rows <- index[[1L]]
    columns <- index[[2L]]
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L)) {
      metadata <- fst.metadata(x@file)
      # TODO: `types` adapted fst:::print.fst.metadata(); find a more robust
      #       way to infer column types
      types <- c("character", "integer", "numeric", "logical", "factor")
      ans <- new(types[metadata$ColumnTypes[1L]])
      dim(ans) <- ans_dim
    } else {
      if (is.numeric(columns)) {
        # NOTE: fst::read.fst() requires column names rather than indices
        columns <- colnames(x)[columns]
      }
      if (is.null(rows)) {
        ans <- read.fst(path = x@file,
                        columns = columns,
                        from = 1,
                        to = NULL,
                        as.data.table = as.data.table)
      } else {
        rows_as_ranges <- IRanges(rows, width = 1L)
        reduced_ranges <- reduce(rows_as_ranges)
        froms <- start(reduced_ranges)
        tos <- end(reduced_ranges)
        ans_list <- mapply(function(from, to) {
          read.fst(path = x@file,
                   columns = columns,
                   from = from,
                   to = to,
                   as.data.table = FALSE)
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

# TODO: Provide access to other args of fst::read.fst()?
#' @importFrom S4Vectors isSingleString new2 wmsg
#' @importFrom tools file_path_as_absolute
#' @export
fstArraySeed <- function(file) {
  if (!isSingleString(file)) {
    stop(wmsg("'file' must be a single string specifying the path to ",
              "the fst file where the dataset is located"))
  }
  # TODO: Check that all columns of fst are of same type (or coercible)?
  file <- file_path_as_absolute(file)
  new2("fstArraySeed", file = file)
}

### ----------------------------------------------------------------------------
### fstArray and fstMatrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate fstArray and fstMatrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

#' @export
setClass("fstArray", contains = "DelayedArray")

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
### as.matrix
###

# NOTE: A shortcut to coerce a fstMatrix to a matrix; this avoids
#       several checks in the as.array,DelayedArray-method, as well as some
#       overhead from S4 method dispatch
# NOTE: Defining as.matrix() rather than as.array() because fstArray currently
#       only support 2-dimensional arrays (matrices) and because as.array()
#       called on a data.frame will break whereas as.matrix() will work
#' @importFrom DelayedArray seed
#' @export
setMethod("as.matrix", "fstMatrix", function(x, ...) {
  # TODO: Avoid coercion to matrix which incurs a copy; can I read from a
  #       fst file directly into a matrix of appropriate dimensions?
  as.matrix(.extract_data_frame_from_fstArraySeed(x = seed(x),
                                                  index = unname(x@index),
                                                  as.data.table = FALSE))
})

#' @importFrom DelayedArray seed
#' @export
setMethod("as.data.frame", "fstMatrix", function(x, ...) {
  .extract_data_frame_from_fstArraySeed(x = seed(x),
                                        index = unname(x@index),
                                        as.data.table = FALSE)
})

# TODO: This would introduce a dependency on data.table, but it would nice to
#       only have this in Suggests rather than Imports. I don't know if it's
#       possible to have data.table in Suggests and still define an
#       as.data.table,fstMatrix-method
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
#' @export
fstArray <- function(file) {
  if (is(file, "fstArraySeed")) {
    seed <- file
  } else {
    seed <- fstArraySeed(file)
  }
  DelayedArray(seed)
}
