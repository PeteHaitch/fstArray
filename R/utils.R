### ============================================================================
### Some low-level utilities
### Functions in this file are not exported
###

#' @importFrom S4Vectors isSingleString wmsg
#' @importFrom tools file_path_as_absolute
.normarg_filepath <- function(filepath, what) {
  if (!isSingleString(filepath)) {
    stop(wmsg(what, " must be a single string specifying the path ",
              "to the fst file where the dataset is located"))
  }
  file_path_as_absolute(filepath)
}

# NOTE: fstRead() returns a data.frame
# TODO: A version of fstRead() that can read from a fst file directly into a
#       matrix of appropriate dimensions
fstRead <- function(fst_table, Nindex) {
  # NOTE: Need special case for when all elements of Nindex are NULL because
  #       of how subsetting is implemented for fst_table objects;
  #       `fst_table[NULL, ..., NULL]` hangs.
  #       NULL subscripts come through here when a fstArray has a 'missing'
  #       element of Nindex
  missing_idx <- S4Vectors:::sapply_isNULL(Nindex)
  if (all(missing_idx)) {
    return(fst_table[])
  }
  # TODO: Fix the NOTE below in fst package
  # NOTE: Can't subset a fst_table object by Nindex containing zero-length
  #       element of Nindex
  #       E.g. this works
  #         data.frame(x = 1:8, y = 8:1)[integer(0), integer(0)]
  #       but
  #         write_fst(data.frame(x = 1:8, y = 8:1), "a.fst")
  #         fst("a.fst")[integer(0), integer(0)]
  #       errors
  empty_idx <- lengths(Nindex) == 0
  if (isTRUE(all(empty_idx))) {
    val <- fst_table[1, ][integer(0), integer(0)]
    return(val)
  }
  DelayedArray:::subset_by_Nindex(fst_table, Nindex)
}
