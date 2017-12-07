- [ ] The 'proper' `writefstArray()` implementation will only work once **fst** supports appending directly on/to disk. Fortunately, [this is planned for `v0.9.0`](https://github.com/fstpackage/fst/issues/91)
- [ ] s/fstArray/fstMatrix/ ? After all, only 2-dimensional arrays are supported because an fst file is designed to support data frames. Hervé notes in the 'Implementing A DelayedArray Backend' that only the fstMatrix form is needed, but it seems a good idea to future-proof the package name by calling fstArray
- [ ] Avoid using copy of HDF5Array internal lock file stuff? Might use a IPCMutex (in BiocParallel or https://github.com/mtmorgan/IPCMutex)
- [ ] Can speed up `[,fstArray-method` using something like the below. This can probably be further optimised by digging into the guts of `[,DelayedArray-method` when called with a 'linear index'. In fact, this will probably also work for random row access

```r
rows <- sample(nrow(x), 10000)
elements <- DelayedArray:::to_linear_index(list(rows, NULL), dim(x))

system.time(a <- as.matrix(fstarray_no_compression[rows, ]))
system.time(b <- matrix(fstarray_no_compression[elements], ncol = ncol(x),
                        dimnames = list(NULL, colnames(x))))
all.equal(a, b)
```

- [ ] Reconcile number of cores used by fst with BiocParallel?
