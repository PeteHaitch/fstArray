
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fstArray

**fstArray** provides an on-disk backend for a
[*DelayedArray*](http://bioconductor.org/packages/DelayedArray/) object
using the [**fst**](https://cran.r-project.org/package=fst) file format.
**fstArray** provides a convenient way to store matrix-like data in an
`.fst` file and operate on it using delayed operations and block
processing.

**fstArray** is currently an experimental package. The `fst` file format
is designed for serializing data frames. Therefore, **fstArray** can
only handle 2-dimensional arrays (matrices), which it does by (ab)using
the `fst` file format to store columns of the matrix as columns of the
serialized data frame. Consequently, **fstArray** is best used for tasks
where data are processed column-wise and, ideally, contiguously by row.

[**HDF5Array**](http://bioconductor.org/packages/HDF5Array/) is the
obvious alternative to **fstArray**. It provides arrays of arbitrary
dimension and allows alternative chunking and serialization strategies
to better support arbitrary data access patterns.

## Installation

You can install fstArray from github with:

``` r
# install.packages("devtools")
devtools::install_github("PeteHaitch/fstarray/")
```

## Example

``` r
library(fstArray)
library(HDF5Array)
library(microbenchmark)
```

Here’s a simple example comparing the time to write the data to disk and
read it back into memory. The data are a ‘long’ matrix with 10 million
rows and 10 columns filled with pseudorandom numbers between 0 and 1.

**fstArray** and **HDF5Array** are compared both without compression and
using their default compression levels. To make the results more
comparable, we chunk the *HDF5Array* by column when using compression,
mimicing the storage scheme of the *fstArray*.

### Writing data to disk

First, serialization. `writefstArray()` is many times faster than
`writeHDF5Array()`, both without or with compression. Remarkably, in
this case using compression in `writefstArray()` comes at no additional
time cost, and potentially even a speed-up. Also remarkable is that
`writefstArray()` must currently coerce the matrix, `x`, to a data frame
before writing to disk, and yet it is still many times faster than
`writeHDF5Array()`\[1\].

**TODO:** Demonstate on data with better compressibility?

``` r
# Simulate a 'long' matrix with numeric data
nrow <- 1e7
ncol <- 10
x <- matrix(runif(nrow * ncol), 
            ncol = ncol, 
            dimnames = list(NULL, letters[seq_len(ncol)]))

# Compression levels
# NOTE: Amount of HDF5 compression (`level`) is on a scale of {0, 1, ..., 9}
hdf5_compression <- 6L
# NOTE: Amount of fst compression (`compress`) is on a scale of {0, 1, ..., 100}
fst_compression <- as.integer(hdf5_compression / 9 * 101)

# Benchmark writing data to disk
tmp_fst_file <- tempfile(fileext = ".fst")
microbenchmark(
  fst_no_compression = writefstArray(x = x, 
                                     file = tmp_fst_file,
                                     compress = 0),
  fst_compression = writefstArray(x = x, 
                                  file = tmp_fst_file,
                                  compress = fst_compression),
  hdf5_no_compression = writeHDF5Array(x = x, 
                                       chunk_dim = c(nrow(x), 1L),
                                       level = 0),
  hdf5_compression = writeHDF5Array(x = x, 
                                    chunk_dim = c(nrow(x), 1L),
                                    level = hdf5_compression),
  times = 5)
#> Unit: seconds
#>                 expr         min          lq        mean      median
#>   fst_no_compression    3.083715    3.192469    3.365107    3.332872
#>      fst_compression    8.852152    9.133636   11.262680   10.815920
#>  hdf5_no_compression   46.534255   46.867626   47.340029   47.496249
#>     hdf5_compression 1105.155307 1107.806281 1112.367453 1112.129620
#>           uq         max neval  cld
#>     3.502163    3.714315     5 a   
#>    11.870003   15.641688     5  b  
#>    47.668909   48.133106     5   c 
#>  1116.475221 1120.270836     5    d
```

### Loading data into memory

We now compare the performance of loading the data into memory as an
ordinary matrix. We first create *fstArray* and *HDF5Array*
representations of the data to be used in this benchmarking:

``` r
# Construct fstArray and HDF5Array instances
# NOTE: Each fst file can only contain one dataset whereas a HDF5 file can 
#       contain multiple datasets
fst_file_no_compression <- tempfile(fileext = ".fst")
fst_array_no_compression <- writefstArray(
  x = x,
  file = fst_file_no_compression,
  compress = 0)
file.size(seed(fst_array_no_compression)@file)
#> [1] 800000262

fst_file_compression <- tempfile(fileext = ".fst")
fst_array_compression <- writefstArray(
  x = x,
  file = fst_file_compression,
  compress = fst_compression)
file.size(seed(fst_array_compression)@file)
#> [1] 514599615

hdf5_array_no_compression <- writeHDF5Array(x = x,
                                            chunk_dim = c(nrow(x), 1L),
                                            level = 0)
file.size(seed(hdf5_array_no_compression)@file)
#> [1] 800000467

hdf5_array_compression <- writeHDF5Array(x = x,
                                         chunk_dim = c(nrow(x), 1L),
                                         level = hdf5_compression)
file.size(seed(hdf5_array_compression)@file)
#> [1] 520384071

# All objects contain the same data (noting that HDF5Array cannot currently 
# store dimnames)
all.equal(as.matrix(fst_array_no_compression), 
          as.matrix(fst_array_compression))
#> [1] TRUE
all.equal(as.matrix(fst_array_no_compression), 
          as.matrix(hdf5_array_no_compression),
          check.attributes = FALSE)
#> [1] TRUE
all.equal(as.matrix(fst_array_no_compression), 
          as.matrix(hdf5_array_compression),
          check.attributes = FALSE)
#> [1] TRUE
```

It is faster to load an entire *fstArray* than an entire *HDF5Array*.
This is despite **fstArray** currently having to do a coercion from a
data frame to a matrix when loading data into memory\[2\].

``` r
# Benchmark loading all data from disk as an ordinary matrix
microbenchmark(
  fst_no_compression = as.matrix(fst_array_no_compression),
  fst_compression = as.matrix(fst_array_compression),
  hdf5_no_compress = as.matrix(hdf5_array_no_compression),
  hdf5_compression = as.matrix(hdf5_array_compression),
  times = 10)
#> Unit: seconds
#>                expr       min        lq      mean    median        uq
#>  fst_no_compression  3.694870  3.737811  4.130642  4.168446  4.258776
#>     fst_compression  4.784667  4.910638  5.107753  5.059689  5.314448
#>    hdf5_no_compress  7.391479  7.434107  7.610014  7.530068  7.670305
#>    hdf5_compression 11.088528 11.622033 11.692814 11.696720 11.917724
#>        max neval  cld
#>   5.020059    10 a   
#>   5.433516    10  b  
#>   8.287026    10   c 
#>  12.056967    10    d
```

Now, loading 10,000 contiguous rows from the middle of the matrix.
Again, using an *fstArray* is a faster than an *HDF5Array*.

``` r
rows <- 56001:66000
# Benchmark loading 10,000 contiguous rows of data from disk as an ordinary 
# matrix
microbenchmark(
  fst_no_compression = as.matrix(fst_array_no_compression[rows, ]),
  fst_compression = as.matrix(fst_array_compression[rows, ]),
  hdf5_no_compress = as.matrix(hdf5_array_no_compression[rows, ]),
  hdf5_compression = as.matrix(hdf5_array_compression[rows, ]),
  times = 10)
#> Unit: milliseconds
#>                expr        min         lq       mean     median         uq
#>  fst_no_compression   19.03678   19.92618   25.91896   28.11214   30.68091
#>     fst_compression   20.78093   21.78739   62.63242   27.86290   34.17939
#>    hdf5_no_compress   23.85341   24.97857   30.16047   28.86553   34.80519
#>    hdf5_compression 3942.72205 3953.07229 4009.50013 4004.99938 4035.87510
#>         max neval cld
#>    32.68903    10  a 
#>   379.83055    10  a 
#>    39.74725    10  a 
#>  4113.37204    10   b
```

Finally, we load 1,000 randomly selected rows, the worst data access
pattern for an *fstArray*. Unsurprisingly, an *fstArray* is slower than
*HDF5Array*, so algorithms should be designed accordingly when
targetting **fstArray**.

``` r
# Benchmark loading 1,000 random rows of data from disk as an ordinary 
# matrix
rows <- sample(nrow(x), 1000)
microbenchmark(
  fst_no_compression = as.matrix(fst_array_no_compression[rows, ]),
  fst_compression = as.matrix(fst_array_compression[rows, ]),
  hdf5_no_compress = as.matrix(hdf5_array_no_compression[rows, ]),
  hdf5_compression = as.matrix(hdf5_array_compression[rows, ]),
  times = 10)
#> Unit: milliseconds
#>                expr       min        lq      mean    median        uq
#>  fst_no_compression  870.2085  895.9988  909.6758  907.8969  923.3426
#>     fst_compression 1185.4158 1224.2673 1243.2411 1232.3954 1279.3689
#>    hdf5_no_compress  174.6490  179.0273  184.6093  181.4484  183.5894
#>    hdf5_compression 3953.6065 3968.0939 4011.7968 4004.4323 4050.7211
#>        max neval  cld
#>   948.5513    10  b  
#>  1291.7582    10   c 
#>   215.7316    10 a   
#>  4101.2778    10    d
```

1.  This coercion and its associated costs may be avoided in future
    versions of **fstArray**, bringing even better performance.

2.  This coercion and its associated costs may be avoided in future
    versions of **fstArray**, bringing even better performance.
