# Validate a Matrix Handed Back to a Parameter

Checks that `m` is a square symmetric numeric matrix of the parameter's
dimension.

## Usage

``` r
check_matrix(s, m, tol = 1e-08)
```

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- m:

  The matrix supplied by the caller.

- tol:

  The relative tolerance for the symmetry check.

## Value

`m`, symmetrized.
