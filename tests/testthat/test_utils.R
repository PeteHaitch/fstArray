### ============================================================================
### Test utility functions
###

context(".normarg_filepath()")

test_that("Normalises filepath", {
  path <- tempfile()
  file.create(path)
  path <- tools::file_path_as_absolute(path)
  funky_path <- file.path(
    dirname(path), "..", basename(dirname(path)), basename(path))
  expect_identical(.normarg_filepath(funky_path), path)
})

test_that("Errors on bad filepath", {
  # TODO: Test exact match of error message?
  expect_error(.normarg_filepath(c("file1.fst", "file2.fst"), "'filepath'"))
})
