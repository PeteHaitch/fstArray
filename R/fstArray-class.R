### ============================================================================
### fstArray objects
### ----------------------------------------------------------------------------

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
#' @importFrom S4Vectors .Call2 wmsg
# NOTE: rows must be a integer vector (although I think an fst file can have
#       > .Machine$integer.max rows)
# NOTE: cols must be a character vector
.extract_array_from_fstArraySeed <- function(x, index) {
  # TODO: Sanity check rows and cols
  stopifnot(length(index) == 2L)
  rows <- index[[1L]]
  columns <- index[[2L]]
  ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
  if (any(ans_dim == 0L)) {
    metadata <- fst.metadata(x@file)
    # TODO: `types` copied from fst:::print.fst.metadata(); find a more robust
    #       way to infer column types
    types <- c("character", "integer", "real", "logical", "factor")
    ans <- new(types[metadata$ColumnTypes[1L]])
    dim(ans) <- ans_dim
  } else {
    if (is.numeric(columns)) {
      # NOTE: read.fst() requires column names rather than indices
      columns <- colnames(x)[columns]
    }
    if (is.null(rows)) {
      # TODO: Avoid coercion to matrix which incurs a copy; can I read from a
      #       fst file directly into a matrix of appropriate dimensions?
      ans <- as.matrix(read.fst(x@file, columns, from = 1, to = NULL))
    } else {
      # TODO: Is it possible to read arbitrary rows of fst file instead of
      #       iterating over [from, to]-ranges?
      # NOTE: This .Call2() is the relevant internal code called when running
      #       `reduce(IRanges(rows, width = 1L))`
      rows_as_ranges <- .Call2("Ranges_reduce",
                               as.integer(rows), # start
                               rep(1L, length(rows)), # width
                               FALSE, # drop.empty.ranges
                               1L, # min.gapwidth
                               FALSE, # with.revmap
                               FALSE, # with.inframe.attrib
                               PACKAGE = "IRanges")
      froms <- rows_as_ranges$start
      tos <- rows_as_ranges$start + rows_as_ranges$width - 1L
      ans_list <- mapply(function(from, to) {
        # TODO: Avoid coercion to matrix which incurs a copy; can I read from a
        #       fst file directly into a matrix of appropriate dimensions?
        as.matrix(read.fst(x@file, columns = columns, from = from, to = to))
      }, from = froms, to = tos, SIMPLIFY = FALSE)
      ans <- do.call(rbind, ans_list)
      # TODO: Allow for duplicated and unsorted row indices
      if (is.unsorted(rows) || anyDuplicated(rows)) {
        ans <- ans[rows, , drop = FALSE]
      }
    }
  }
  ans
}

#' @importMethodsFrom DelayedArray extract_array
#' @export
setMethod("extract_array", "fstArraySeed", .extract_array_from_fstArraySeed)

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

#' @importFrom S4Vectors wmsg
.validate_fstArray <- function(x){
  if (!is(x@seed, "fstArraySeed")) {
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
#' @export
fstArray <- function(file) {
  if (is(file, "fstArraySeed")) {
    seed <- file
  } else {
    seed <- fstArraySeed(file)
  }
  DelayedArray(seed)
}
