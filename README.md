
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fstArray

[![Travis-CI Build
Status](https://travis-ci.org/PeteHaitch/fstArray.svg?branch=master)](https://travis-ci.org/PeteHaitch/fstArray)

**fstArray** provides an on-disk backend for a
[*DelayedArray*](http://bioconductor.org/packages/DelayedArray/) object
using the [**fst**](https://cran.r-project.org/package=fst) file format.
**fstArray** provides a convenient way to store matrix-like data in an
`.fst` file and operate on it using delayed operations and block
processing.

The `fst` file format is designed for serializing data frames.
Therefore, **fstArray** can only handle 2-dimensional arrays (matrices),
which it does by (ab)using the `fst` file format to store columns of the
matrix as columns of the serialized data frame. Consequently,
**fstArray** is best used for tasks where data are processed column-wise
(or by all columns at once) and, ideally, by processing contiguous rows.

[**HDF5Array**](http://bioconductor.org/packages/HDF5Array/) is the
obvious alternative to **fstArray**. The `HDF5` file format provides
arrays of arbitrary dimension and allows alternative chunking and
serialization strategies to better support arbitrary data access
patterns.

In initial testing, I find that **fstArray** is generally faster at
reading data from disk, surprisingly often even when reading data using
sub-optimal access patterns.

**Limitations**: Currently, there is no equivalent to
`HDF5Array::writeHDF5Array()`, i.e. `writefstArray()`. This will be
added in a future version of **fstArray** once the core **fst** library
supports [appending data directly on/to
disk](https://github.com/fstpackage/fst/issues/91). For now, an
*fstArray* can be created by wrapping a `fst` file created by
`fst::write_fst()`.

## Installation

You can install fstArray from github with:

``` r
# install.packages("devtools")
devtools::install_github("PeteHaitch/fstarray/")
```

## Example

``` r
library(fstArray)
library(fst)
library(HDF5Array)
library(microbenchmark)
```

Here’s a simple example comparing the time to read subsets of the data
back into memory using **fstArray** and **HDF5Array**. **fstArray** and
**HDF5Array** are compared without and with data compression\[1\]. To
make the results more comparable, we chunk the *HDF5Array* by column
when using compression, mimicing the storage scheme of the *fstArray*.

The data are a ‘long’ matrix with 1 million rows and 10 columns filled
with pseudorandom numbers between 0 and 1.

``` r
# Simulate a 'long' matrix with numeric data
nrow <- 1e6
ncol <- 10
x <- matrix(runif(nrow * ncol), 
            ncol = ncol, 
            dimnames = list(NULL, letters[seq_len(ncol)]))
```

  - **TODO**: Also demonstate on data with better compressibility (e.g.
    sequencing coverage-like data), and ‘wide’ data. So, 4 datasets:
    long runif, wide runif, long coverage, wide scRNA-seq
  - **TODO**: *HDF5Array* with default chunking

### Summary of results

**fstArray** is slightly-to-noticeably faster than **HDF5Array** when
reading contiguous rows of data into memory as an ordinary *matrix*,
with the performance gap increasing when data compression is used. This
is the best-case data access pattern for a *fstArray*.

As expected, **fstArray** is slower than **HDF5Array** when reading
random rows of uncompressed data into memory as an ordinary matrix.
Surprisingly, however, **fstArray** is faster than **HDF5Array** on this
test when data compression is used.

Finally, the worst-case data access pattern for **fstArray** is reading
random elements of data into memory as an ordinary vector. Surprisingly,
**fstArray** is 2x faster than **HDF5Array** on this benchmark.

### Loading data into memory

**TODO**: Revisit text after re-running benchmarking

We now compare the performance of loading the data into memory as an
ordinary matrix. We first create *fstArray* and *HDF5Array*
representations of the data to be used in this benchmarking:

``` r
# Construct fstArray and HDF5Array instances

# Compression levels
# NOTE: Amount of HDF5 compression (`level`) is on a scale of {0, 1, ..., 9}
hdf5_compression_level <- getHDF5DumpCompressionLevel()
hdf5_compression_level
#> [1] 6
# NOTE: Amount of fst compression (`compress`) is on a scale of {0, 1, ..., 100}
fst_compression_level <- as.integer(hdf5_compression_level / 9 * 101)

# NOTE: Each fst file can only contain one dataset whereas a HDF5 file can 
#       contain multiple datasets
fst_no_compression <- tempfile(fileext = ".fst")
fst::write_fst(x = as.data.frame(x),
               path = fst_no_compression,
               compress = 0)
fstarray_no_compression <- fstArray(fst_no_compression)
fst_with_compression <- tempfile(fileext = ".fst")
fst::write_fst(x = as.data.frame(x),
               path = fst_with_compression,
               compress = fst_compression_level)
fstarray_with_compression <- fstArray(fst_with_compression)
hdf5array_no_compression <- writeHDF5Array(x = x,
                                           chunk_dim = c(nrow(x), 1L),
                                           level = 0)
hdf5array_with_compression <- writeHDF5Array(x = x,
                                             chunk_dim = c(nrow(x), 1L),
                                             level = hdf5_compression_level)
```

All objects contain the same data\[2\]:

``` r
all.equal(as.matrix(fstarray_no_compression), 
          as.matrix(fstarray_with_compression))
#> [1] TRUE
all.equal(as.matrix(fstarray_no_compression), 
          as.matrix(hdf5array_no_compression),
          check.attributes = FALSE)
#> [1] TRUE
all.equal(as.matrix(fstarray_no_compression), 
          as.matrix(hdf5array_with_compression),
          check.attributes = FALSE)
#> [1] TRUE
```

Similar file sizes are achieved for *fstArray* and *HDF5Array*
instances:

| name                         | class       | compression | size (bytes) |
| :--------------------------- | :---------- | :---------- | -----------: |
| fstarray\_no\_compression    | *fstArray*  | FALSE       |     80000622 |
| fstarray\_with\_compression  | *fstArray*  | TRUE        |     70537581 |
| hdf5array\_no\_compression   | *HDF5Array* | FALSE       |     80000467 |
| hdf5array\_with\_compression | *HDF5Array* | TRUE        |     52043398 |

#### Loading all data into memory

It is faster to load an entire *fstArray* than an entire *HDF5Array*.
This is despite **fstArray** currently having to do a coercion from a
data frame to a matrix when loading data into memory\[^coercion\].

``` r
# Benchmark loading all data from disk as an ordinary matrix
# TODO: Make a plot of result
microbenchmark(
  fstarray_no_compression = as.matrix(fstarray_no_compression),
  fstarray_with_compression = as.matrix(fstarray_with_compression),
  hdf5array_no_compression = as.matrix(hdf5array_no_compression),
  hdf5array_with_compression = as.matrix(hdf5array_with_compression),
  times = 10)
#> Unit: milliseconds
#>                        expr      min       lq     mean   median       uq
#>     fstarray_no_compression 102.5663 127.5894 140.7778 135.0653 150.5365
#>   fstarray_with_compression 140.0467 140.6319 186.4727 157.8618 192.6200
#>    hdf5array_no_compression 293.3690 304.2665 370.3637 365.1818 424.6711
#>  hdf5array_with_compression 666.2197 726.2155 752.4901 753.8953 782.6305
#>       max neval cld
#>  198.1223    10 a  
#>  306.8345    10 a  
#>  494.7491    10  b 
#>  819.0964    10   c
```

#### Loading contiguous rows of data into memory

Now, loading 10,000 contiguous rows from the middle of the matrix.
Again, using an *fstArray* is a faster than an *HDF5Array*.

``` r
rows <- 500001:600000
# TODO: Make a plot of result (either ggplot2::autoplot() or what's done in fst
#       README)
microbenchmark(
  fstarray_no_compression = as.matrix(fstarray_no_compression[rows, ]),
  fstarray_with_compression = as.matrix(fstarray_with_compression[rows, ]),
  hdf5array_no_compression = as.matrix(hdf5array_no_compression[rows, ]),
  hdf5array_with_compression = as.matrix(hdf5array_with_compression[rows, ]),
  times = 10)
#> Unit: milliseconds
#>                        expr       min        lq      mean    median
#>     fstarray_no_compression  57.57056  58.95272  76.35915  69.73875
#>   fstarray_with_compression  76.18813  80.73812  88.81016  85.99229
#>    hdf5array_no_compression 106.34253 116.09887 138.27331 146.25949
#>  hdf5array_with_compression 470.16741 474.82125 495.50887 489.72496
#>         uq      max neval cld
#>   90.47825 123.1876    10 a  
#>   91.30427 118.1874    10 a  
#>  156.64767 162.8812    10  b 
#>  518.77466 531.4070    10   c
```

#### Loading random rows of data into memory

We load 10,000 randomly selected rows, the worst data access pattern for
an *fstArray*. Unsurprisingly, an *fstArray* is slower than *HDF5Array*,
so algorithms should be designed accordingly when targetting
**fstArray**.

``` r
rows <- sample(nrow(x), 10000)
microbenchmark(
  fstarray_no_compression = as.matrix(fstarray_no_compression[rows, ]),
  fstarray_with_compression = as.matrix(fstarray_with_compression[rows, ]),
  hdf5array_no_compression = as.matrix(hdf5array_no_compression[rows, ]),
  hdf5array_with_compression = as.matrix(hdf5array_with_compression[rows, ]),
  times = 10)
#> Unit: seconds
#>                        expr       min        lq     mean    median
#>     fstarray_no_compression  7.712612  7.808474  7.94896  7.963342
#>   fstarray_with_compression 11.067966 11.187916 11.37206 11.289247
#>    hdf5array_no_compression 13.737002 13.769707 13.80272 13.805314
#>  hdf5array_with_compression 14.156510 14.236498 14.29448 14.296227
#>         uq       max neval  cld
#>   7.989697  8.275694    10 a   
#>  11.532162 11.854357    10  b  
#>  13.831551 13.903240    10   c 
#>  14.324118 14.550491    10    d
```

#### Loading random elements of data into memory

Benchmark loading 100,000 random elements of data from disk as an
ordinary vector

**TODO**: How does `[,DelayedArray-method` work?

``` r
elements <- sample(length(x), 100000)
microbenchmark(
  fstarray_no_compression = fstarray_no_compression[elements],
  fstarray_with_compression = fstarray_with_compression[elements],
  hdf5array_no_compression = hdf5array_no_compression[elements],
  hdf5array_with_compression = hdf5array_with_compression[elements],
  times = 10)
#> Unit: seconds
#>                        expr      min       lq     mean   median       uq
#>     fstarray_no_compression 1.305587 1.317115 1.369593 1.358739 1.394148
#>   fstarray_with_compression 1.977945 1.997501 2.066805 2.047448 2.126624
#>    hdf5array_no_compression 2.397374 2.415154 2.492379 2.454718 2.594813
#>  hdf5array_with_compression 3.104193 3.190160 3.269943 3.221916 3.349427
#>       max neval  cld
#>  1.514849    10 a   
#>  2.244289    10  b  
#>  2.626823    10   c 
#>  3.527146    10    d
```

1.  The level of compression is the default used by **HDF5Array**,
    scaled for use with **fstArray**.

2.  Currently, a *HDF5Array* cannot store the dimnames of the array
    whereas an *fstArray* can stored column names but not rownames.
