### ============================================================================
### fst dump management
### ----------------------------------------------------------------------------
###

#' @include utils.R

### ----------------------------------------------------------------------------
### A global internal counters  for the dump files
###
### The 2 counters are safe to use in the context of parallel execution e.g.
###
###   library(BiocParallel)
###   bplapply(1:5, function(i) .get_dump_files_global_counter(increment=TRUE))
###

.get_dump_files_global_counter_filepath <- function() {
  file.path(tempdir(), "fstArray_dump_files_global_counter")
}

# NOTE: Called by .onLoad() hook (see zzz.R file).
init_fst_dump_files_global_counter <- function() {
  filepath <- .get_dump_files_global_counter_filepath()
  init_global_counter(filepath)
}

.get_dump_files_global_counter <- function(increment = FALSE) {
  filepath <- .get_dump_files_global_counter_filepath()
  get_global_counter(filepath, increment = increment)
}

### ----------------------------------------------------------------------------
### Normalization (with basic checking) of an fst file path or dataset name
###

# Return the *absolute path* to the dump file.
# Has the side effect of creating the file as an empty fst file if it does
# not exist yet.
#' @importFrom S4Vectors isSingleString wmsg
#' @importFrom tools file_path_as_absolute
normalize_dump_file <- function(file) {
  if (!isSingleString(file) || file == "")
    stop(wmsg("'file' must be a non-empty string specifying the path ",
              "to a new fst file"))
  if (!file.exists(file)) {
    fstCreateFile(file)
  }
  file_path_as_absolute(file)
}

### ----------------------------------------------------------------------------
### Very low-level stuff used in this file only
###

fstCreateFile <- function(file) {
  if (is.character(file)) {
    file <- normalizePath(file, mustWork = FALSE)
    if (file.exists(file)) {
      message("file '", file, "' already exists.")
    } else {
      file.create(file)
    }
  } else {
    stop("file has to be a valid filename.")
  }
}

.dump_settings_envir <- new.env(parent = emptyenv())

# NOTE: Create directory 'dir' if it doesn't exist yet.
#' @importFrom S4Vectors wmsg
#' @importFrom tools file_path_as_absolute
.set_dump_dir <- function(dir) {
  # NOTE: Even though file_path_as_absolute() will trim the trailing slashes,
  #       we need to do this early. Otherwise, checking for the existence of a
  #       file of the same name as the to-be-created directory will fail.
  if (nchar(dir) > 1L) {
    dir <- trim_trailing_slashes(dir)
  }
  if (!dir.exists(dir)) {
    if (file.exists(dir)) {
      stop(wmsg("\"", dir, "\" already exists and is a file, ",
                "not a directory"))
    }
    if (!suppressWarnings(dir.create(dir))) {
      stop("cannot create directory \"", dir, "\"")
    }
  }
  dir <- file_path_as_absolute(dir)
  assign("dir", dir, envir = .dump_settings_envir)
}

.set_dump_autofiles_mode <- function() {
  suppressWarnings(rm(list = "specfile", envir = .dump_settings_envir))
}

# Create file as an empty fst file if it doesn't exist yet.
.set_dump_specfile <- function(file) {
  file <- normalize_dump_file(file)
  assign("specfile", file, envir = .dump_settings_envir)
}

# Return the user-specified file of the dump or an error if the user didn't
# specify a file.
.get_dump_specfile <- function() {
  get("specfile", envir = .dump_settings_envir)
}

### ----------------------------------------------------------------------------
### get/setfstDumpDir()
###

# #' @export
getfstDumpDir <- function() {
  get("dir", envir = .dump_settings_envir)
}

# Create auto file as an empty fst file if it doesn't exist yet.
.get_dump_autofile <- function(increment = FALSE) {
  counter <- .get_dump_files_global_counter(increment = increment)
  file <- file.path(getfstDumpDir(), sprintf("auto%06d.fst", counter))
  if (!file.exists(file)) {
    fstCreateFile(file)
  }
  file
}

# Called by .onLoad() hook (see zzz.R file).
#' @importFrom S4Vectors isSingleString wmsg
# #' @export
setfstDumpDir <- function(dir) {
  if (missing(dir)) {
    dir <- file.path(tempdir(), "fstArray_dump")
  } else if (!isSingleString(dir) || dir == "") {
    stop(wmsg("'dir' must be a non-empty string specifying the path ",
              "to a new or existing directory"))
  }
  dir <- .set_dump_dir(dir)
  .set_dump_autofiles_mode()
  .get_dump_autofile()
  invisible(dir)
}

### ----------------------------------------------------------------------------
### get/setfstDumpFile()
###

# Set the current fst dump file. Create it as an empty fst file if it doesn't
# exist yet.
#' @importFrom S4Vectors isSingleString
# #' @export
setfstDumpFile <- function(file) {
  if (missing(file)) {
    .set_dump_autofiles_mode()
    file <- .get_dump_autofile()
  } else {
    if (!isSingleString(file) || file == "") {
      stop("'file' must be a non-empty string")
    }
    if (has_trailing_slash(file)) {
      setfstDumpDir(file)
      file <- .get_dump_autofile()
    } else {
      file <- .set_dump_specfile(file)
    }
  }
}

# Return the *absolute path* to the dump file.
#' @importFrom S4Vectors isTRUEorFALSE
# #' @export
getfstDumpFile <- function(for.use = FALSE) {
  if (!isTRUEorFALSE(for.use)) {
    stop("'for.use' must be TRUE or FALSE")
  }
  file <- try(.get_dump_specfile(), silent = TRUE)
  if (is(file, "try-error")) {
    file <- .get_dump_autofile(increment = for.use)
  }
  file
}

### ----------------------------------------------------------------------------
### set/getfstDumpCompressionLevel
###

#' @importFrom S4Vectors isSingleNumber
normalize_compression_level <- function(compress) {
  if (!isSingleNumber(compress)) {
    stop("'compress' must be a single number")
  }
  if (!is.integer(compress)) {
    compress <- as.integer(compress)
    if (compress < 0L || compress > 100L) {
      stop(wmsg("'compress' must be between 0 (no compression) ",
                "and 100 (highest and slowest compression)"))
    }
  }
  compress
}

# TODO: Update default compresion level to fst::write_fst default
# Called by .onLoad() hook (see zzz.R file).
# #' @export
setfstDumpCompressionLevel <- function(compress = 0L) {
  compress <- normalize_compression_level(compress)
  assign("compress", compress, envir = .dump_settings_envir)
}

# #' @export
getfstDumpCompressionLevel <- function() {
  get("compress", envir = .dump_settings_envir)
}

### ----------------------------------------------------------------------------
### Dump log
###

# Called by .onLoad() hook (see zzz.R file).
get_fst_dump_logfile <- function() {
  file.path(tempdir(), "fstArray_dump_log")
}

.get_dataset_creation_global_counter_filepath <- function() {
  file.path(tempdir(), "fstArray_dataset_creation_global_counter")
}

# Called by .onLoad() hook (see zzz.R file).
init_fst_dataset_creation_global_counter <- function() {
  filepath <- .get_dataset_creation_global_counter_filepath()
  init_global_counter(filepath)
}

.get_dataset_creation_global_counter <- function(increment = FALSE) {
  filepath <- .get_dataset_creation_global_counter_filepath()
  get_global_counter(filepath, increment = increment)
}

# TODO: Consider adding type to log
# NOTE: Use a lock mechanism so is safe to use in the context of parallel
#       execution.
appendDatasetCreationTofstDumpLog <- function(file, dim, compress) {
  logfile <- get_fst_dump_logfile()
  locked_path <- lock_file(logfile)
  on.exit(unlock_file(logfile))
  counter <- .get_dataset_creation_global_counter(increment = TRUE)
  dims <- paste0(dim, collapse = "x")
  cat(as.character(Sys.time()), counter, file, dims, compress,
      sep = "\t", file = locked_path, append = TRUE)
  cat("\n", file = locked_path, append = TRUE)
}


showfstDumpLog <- function() {
  COLNAMES <- c("time", "counter", "file", "dims", "compress")
  # NOTE: The number of lines in the log file is the current value of the
  #       dataset creation counter minus one.
  counter <- .get_dataset_creation_global_counter()
  if (counter == 1L) {
    dump_log <- data.frame(time = character(0),
                           counter = integer(0),
                           file = character(0),
                           dims = character(0),
                           compress = integer(0),
                           stringsAsFactors = FALSE)
  } else {
    logfile <- get_fst_dump_logfile()
    locked_path <- lock_file(logfile)
    on.exit(unlock_file(logfile))
    dump_log <- read.table(locked_path, sep = "\t", stringsAsFactors = FALSE)
    colnames(dump_log) <- COLNAMES
    fmt <- "[%s] #%d In file '%s': creation of dataset (%s, compress=%d)"
    message(paste0(sprintf(fmt,
                           dump_log$time, dump_log$counter,
                           dump_log$file,
                           dump_log$dims,
                           dump_log$level),
                   "\n"),
            appendLF = FALSE)
  }
  invisible(dump_log)
}
