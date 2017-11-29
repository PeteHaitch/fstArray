
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fstArray

[![Travis-CI Build Status](https://travis-ci.org/PeteHaitch/fstArray.svg?branch=master)](https://travis-ci.org/PeteHaitch/fstArray)

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
`fst::write.fst()`.

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
fst::write.fst(x = as.data.frame(x),
               path = fst_no_compression,
               compress = 0)
fstarray_no_compression <- fstArray(fst_no_compression)
fst_with_compression <- tempfile(fileext = ".fst")
fst::write.fst(x = as.data.frame(x),
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
| fstarray\_no\_compression    | *fstArray*  | FALSE       |     80000262 |
| fstarray\_with\_compression  | *fstArray*  | TRUE        |     51549064 |
| hdf5array\_no\_compression   | *HDF5Array* | FALSE       |     80000467 |
| hdf5array\_with\_compression | *HDF5Array* | TRUE        |     52044232 |

#### Loading all data into memory

It is faster to load an entire *fstArray* than an entire *HDF5Array*.
This is despite **fstArray** currently having to do a coercion from a
data frame to a matrix when loading data into memory\[^coercion\].

``` r
# Benchmark loading all data from disk as an ordinary matrix
microbenchmark(
  fstarray_no_compression = as.matrix(fstarray_no_compression),
  fstarray_with_compression = as.matrix(fstarray_with_compression),
  hdf5array_no_compression = as.matrix(hdf5array_no_compression),
  hdf5array_with_compression = as.matrix(hdf5array_with_compression),
  times = 10)
#> Unit: milliseconds
#>                        expr       min       lq      mean    median
#>     fstarray_no_compression  87.22516 125.3228  195.0390  227.4134
#>   fstarray_with_compression 194.48685 205.0601  237.7612  215.9094
#>    hdf5array_no_compression 474.02809 560.4660  606.5597  635.2825
#>  hdf5array_with_compression 925.15713 934.6240 1026.2005 1033.9560
#>         uq       max neval cld
#>   250.7414  271.6434    10 a  
#>   237.0523  334.9566    10 a  
#>   667.1981  670.9407    10  b 
#>  1052.3974 1155.8843    10   c
```

#### Loading contiguous rows of data into memory

Now, loading 10,000 contiguous rows from the middle of the matrix.
Again, using an *fstArray* is a faster than an *HDF5Array*.

``` r
rows <- 500001:600000
microbenchmark(
  fstarray_no_compression = as.matrix(fstarray_no_compression[rows, ]),
  fstarray_with_compression = as.matrix(fstarray_with_compression[rows, ]),
  hdf5array_no_compression = as.matrix(hdf5array_no_compression[rows, ]),
  hdf5array_with_compression = as.matrix(hdf5array_with_compression[rows, ]),
  times = 10)
#> Unit: milliseconds
#>                        expr       min        lq      mean    median
#>     fstarray_no_compression  71.61974  77.18624  87.64469  88.38161
#>   fstarray_with_compression  86.87007  87.23516 101.77223  98.74093
#>    hdf5array_no_compression 129.37389 138.05747 146.56411 142.37873
#>  hdf5array_with_compression 480.86832 492.38430 511.25762 511.38382
#>         uq      max neval cld
#>   95.26093 104.0776    10 a  
#>  107.86351 129.1088    10 a  
#>  156.56137 169.8482    10  b 
#>  523.76236 553.9024    10   c
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
#>                        expr       min        lq      mean    median
#>     fstarray_no_compression  8.002406  8.230121  8.341033  8.377256
#>   fstarray_with_compression 10.829110 11.145667 11.162961 11.171536
#>    hdf5array_no_compression 13.232171 13.264697 13.309880 13.300659
#>  hdf5array_with_compression 13.709490 13.717944 13.806906 13.784387
#>         uq       max neval  cld
#>   8.484604  8.606015    10 a   
#>  11.257842 11.344139    10  b  
#>  13.345654 13.404611    10   c 
#>  13.834393 14.095220    10    d
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
#>     fstarray_no_compression 1.471993 1.492619 1.501118 1.496295 1.504998
#>   fstarray_with_compression 1.551851 1.587602 1.616991 1.594200 1.636589
#>    hdf5array_no_compression 2.519001 2.551852 2.610002 2.580096 2.691903
#>  hdf5array_with_compression 3.251076 3.325241 3.432805 3.441513 3.534551
#>       max neval  cld
#>  1.541528    10 a   
#>  1.763489    10  b  
#>  2.728017    10   c 
#>  3.609503    10    d
```

1.  The level of compression is the default used by **HDF5Array**,
    scaled for use with **fstArray**.

2.  Currently, a *HDF5Array* cannot store the dimnames of the array
    whereas an *fstArray* can stored column names but not rownames.
