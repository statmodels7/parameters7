# The Spectral Decomposition a Parameter's Quantities Are Read From

The eigenvalues and eigenvectors of \\M(\eta)\\, with the eigenvalues
the declared rank says are zero identified by position rather than by
size.

## Usage

``` r
param_spectrum(s, eta)
```

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A list with `values`, `vectors`, and `keep`, a logical vector marking
the `s@rank` directions that carry the matrix.

## Details

The rank is taken from the object and never re-derived here. Counting
eigenvalues above a relative tolerance is not scale invariant, so a
parameter whose components differ by many orders of magnitude would be
assigned a different rank at different \\\eta\\ – and a fitted model
with smoothing parameters that far apart is ordinary. The object settled
the question once, at construction, from the components.
