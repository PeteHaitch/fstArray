
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
read it back into memory. The data are a ‘long’ matrix with 1 million
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
nrow <- 1e6
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
  times = 4)
#> Unit: milliseconds
#>                 expr        min        lq       mean     median         uq
#>   fst_no_compression   332.5177   347.755   400.3892   380.8732   453.0234
#>      fst_compression   969.4770   981.049  1019.5806   993.4713  1058.1122
#>  hdf5_no_compression  4720.5070  4780.557  5148.6595  4918.7139  5516.7616
#>     hdf5_compression 23385.0729 24477.333 28632.0234 25846.8578 32786.7134
#>         max neval cld
#>    507.2928     4  a 
#>   1121.9026     4  a 
#>   6036.7034     4  a 
#>  39449.3051     4   b
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
#> [1] 80000262

fst_file_compression <- tempfile(fileext = ".fst")
fst_array_compression <- writefstArray(
  x = x,
  file = fst_file_compression,
  compress = fst_compression)
file.size(seed(fst_array_compression)@file)
#> [1] 51549510

hdf5_array_no_compression <- writeHDF5Array(x = x,
                                            chunk_dim = c(nrow(x), 1L),
                                            level = 0)
file.size(seed(hdf5_array_no_compression)@file)
#> [1] 80000467

hdf5_array_compression <- writeHDF5Array(x = x,
                                         chunk_dim = c(nrow(x), 1L),
                                         level = hdf5_compression)
file.size(seed(hdf5_array_compression)@file)
#> [1] 52044162

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
#> Unit: milliseconds
#>                expr       min        lq      mean    median        uq
#>  fst_no_compression  499.5862  508.9117  548.0469  528.3886  550.4041
#>     fst_compression  529.0804  615.3113  634.0544  650.5486  658.7690
#>    hdf5_no_compress  795.0818  833.3832  871.1172  864.8767  899.5042
#>    hdf5_compression 1237.5181 1317.3233 1366.2545 1332.0140 1395.8388
#>        max neval cld
#>   683.1246    10 a  
#>   680.7218    10 a  
#>   958.0919    10  b 
#>  1644.5266    10   c
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
#>                expr       min        lq      mean    median        uq
#>  fst_no_compression  11.54032  13.12189  28.68029  16.79507  23.16494
#>     fst_compression  14.22772  15.70217  18.48921  19.55733  20.86771
#>    hdf5_no_compress  18.76924  20.13804  23.26911  21.73978  23.96825
#>    hdf5_compression 449.42446 487.47269 507.26217 502.56514 529.22231
#>        max neval cld
#>  129.21273    10  a 
#>   22.84729    10  a 
#>   35.94599    10  a 
#>  568.90620    10   b
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
#>  fst_no_compression 1259.6833 1319.2925 1420.6746 1393.0589 1515.6856
#>     fst_compression 1710.8086 1777.4868 1923.7452 1839.3203 2094.8817
#>    hdf5_no_compress  131.3710  134.1296  145.8661  141.3843  150.4381
#>    hdf5_compression  582.0486  593.4916  633.6262  604.3689  666.6740
#>        max neval  cld
#>  1751.2337    10   c 
#>  2465.8563    10    d
#>   172.7642    10 a   
#>   733.6317    10  b
```

1.  This coercion and its associated costs may be avoided in future
    versions of **fstArray**, bringing even better performance.

2.  This coercion and its associated costs may be avoided in future
    versions of **fstArray**, bringing even better performance.
